"""
Module 10 — Machine Learning Scoring (XGBoost + Optuna + SHAP).
Snakemake script: called via script: directive.

PURPOSE
Build two classifiers that combine all analytical outputs into a unified
clinical scoring framework:

  Classifier A — 4-class staging
    Input : all samples (Ctrl / NEW / MS-Rituximab-Stable / MS-Rituximab-Progressive)
    Output: per-sample probability vector over the four classes
    Derived: MEAS score = P(Progressive) − P(Ctrl), a continuous disease
              severity metric that correlates with EDSS, NfL, and CDI.

  Classifier B — binary Rituximab response
    Input : treated patients only (Stable vs Progressive)
    Output: binary classification + ROC-AUC
    Purpose: identify epigenomic features that predict who will progress
              despite Rituximab treatment — a clinically critical question.

FEATURE BLOCKS (assembled per sample)
  1. Top RRE features from DESeq2 (significant in any contrast, ranked by
     variance and capped at max_feats/2 to avoid high-dimensional overfitting).
  2. Cell-type deconvolution z-scores from Module 4 (one per cell type).
  3. TF AUC scores from Module 8 (one per TF).
  4. Composite indices from Module 4: CDI, NII, IAI.
  5. Scaled B Cell Index (bci_scaled) from Module 3.

MODEL: XGBoost (gradient-boosted trees)
  Chosen because: handles mixed-type features, robust to correlated features,
  natively multi-class, compatible with SHAP for interpretability.

HYPERPARAMETER TUNING: Optuna (Bayesian optimisation with TPE sampler)
  Tunes max_depth, learning_rate, n_estimators, subsample, colsample_bytree,
  min_child_weight via the same CV as evaluation (no separate inner loop to
  avoid the prohibitive cost of nested CV on small clinical cohorts).

CROSS-VALIDATION
  LOOCV (Leave-One-Out) is used if any group has < 20 samples (default in
  small MS cohorts) because it maximises training data per fold.
  5-fold stratified CV is used for larger cohorts.

INTERPRETABILITY: SHAP TreeExplainer
  SHAP (SHapley Additive exPlanations) values quantify each feature's
  contribution to each prediction. The summary plot reveals which features
  drive the classifier's decisions and in which direction.
"""

import warnings
warnings.filterwarnings("ignore")

from pathlib import Path
import numpy as np
import pandas as pd
from sklearn.model_selection import StratifiedKFold, LeaveOneOut
from sklearn.metrics import (
    accuracy_score, roc_auc_score, confusion_matrix, ConfusionMatrixDisplay
)
from sklearn.preprocessing import LabelEncoder
import xgboost as xgb
import optuna
optuna.logging.set_verbosity(optuna.logging.WARNING)
import shap
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.backends.backend_pdf as pdf_backend
import seaborn as sns


# ── Step 1: Load all input data ───────────────────────────────────────────────
meta      = pd.read_csv(snakemake.input.meta,      sep="\t", comment="#")
composite = pd.read_csv(snakemake.input.composite,  sep="\t")   # CDI, NII, IAI per sample
deconv    = pd.read_csv(snakemake.input.deconv,     sep="\t")   # cell-type z-scores (long format)
tf_auc    = pd.read_csv(snakemake.input.tf_auc,     sep="\t", index_col=0)  # TF × samples matrix
bci       = pd.read_csv(snakemake.input.bci,        sep="\t")   # BCI scores from Module 3

# Pipeline configuration parameters passed from Snakemake
rng         = np.random.default_rng(snakemake.params.random_seed)
n_trials    = int(snakemake.params.n_optuna_trials)  # Optuna tuning iterations
cv_strategy = snakemake.params.cv_strategy           # "loocv" or "stratified_kfold"
cv_folds    = int(snakemake.params.cv_folds)         # used only for stratified kfold
max_feats   = int(snakemake.params.max_features)     # cap on total features

# Restrict to QC-passing samples (exclude QC_Mix and flagged samples)
meta   = meta[meta["qc_include"].astype(str).str.upper() == "TRUE"]
GROUPS = ["Ctrl", "NEW", "MS-Rituximab-Stable", "MS-Rituximab-Progressive"]

# ── Step 2: Assemble the feature matrix ──────────────────────────────────────

# ── Block 1: Top DESeq2 RRE features ─────────────────────────────────────────
# Load the full RRE count matrix and collect all region IDs that were
# significant in at least one contrast (from the _significant.tsv files).
# Then select the top max_feats/2 by variance — high-variance features carry
# the most discriminatory information for XGBoost.
rre_counts   = pd.read_csv(snakemake.input.rre_counts, sep="\t", index_col=0)
sig_features = set()
for f in snakemake.input.deseq2_full:
    try:
        df = pd.read_csv(f, sep="\t", index_col=0)
        sig_features.update(df.index.tolist())
    except Exception:
        pass
if sig_features:
    top_rre = rre_counts.loc[rre_counts.index.isin(sig_features)].T.copy()
    # Keep the highest-variance features up to the cap
    top_rre = top_rre[top_rre.std().nlargest(min(max_feats // 2, len(sig_features))).index]
else:
    # No significant features found (pipeline early-stage); use empty frame
    top_rre = pd.DataFrame(index=meta["sample_id"])

# ── Block 2: Cell-type z-scores ────────────────────────────────────────────────
# Convert from long format (sample_id, cell_type, zscore) to wide format
# (sample_id as index, one column per cell type). Prefix 'CT_' to avoid
# name collisions with RRE features.
ct_pivot = deconv.pivot_table(
    index="sample_id", columns="cell_type", values="zscore"
)
ct_pivot.columns = [f"CT_{c}" for c in ct_pivot.columns]

# ── Block 3: TF AUC scores ────────────────────────────────────────────────────
# tf_auc is currently TF × samples; transpose to samples × TF and prefix 'TF_'.
tf_t = tf_auc.T.copy()
tf_t.columns = [f"TF_{c}" for c in tf_t.columns]

# ── Block 4: Composite indices and BCI ────────────────────────────────────────
composite_sel = composite.set_index("sample_id")[["CDI", "NII", "IAI"]]
bci_sel       = bci.set_index("sample_id")[["bci_scaled"]]

# ── Merge all blocks ──────────────────────────────────────────────────────────
# Use left join from the metadata sample list so that samples with missing
# features get NaN (filled with 0 below) rather than being silently dropped.
feat = (
    meta[["sample_id"]].set_index("sample_id")
    .join(top_rre,        how="left")
    .join(ct_pivot,       how="left")
    .join(tf_t,           how="left")
    .join(composite_sel,  how="left")
    .join(bci_sel,        how="left")
)
# Fill NaN with 0 — missing features (e.g. TF BED not yet downloaded) are treated
# as uninformative rather than causing model fit failure.
feat = feat.fillna(0)

# Align feature matrix with metadata labels
common = [s for s in feat.index if s in meta["sample_id"].values]
feat   = feat.loc[common]
labels = meta.set_index("sample_id").loc[common, "group"]

print(f"[Module 10] Feature matrix: {feat.shape[0]} samples × {feat.shape[1]} features")


# ── Helper: cross-validated XGBoost with Optuna tuning ───────────────────────

def tune_and_cv(X, y, cv, seed):
    """
    Runs Optuna hyperparameter search and then cross-validated prediction.

    Steps:
      1. Encode string labels to integers with LabelEncoder.
      2. Define an Optuna objective function that trains XGBoost on the
         CV training folds and evaluates macro-averaged ROC-AUC on the
         validation fold. Return the mean validation AUC across folds.
      3. Run Optuna for n_trials evaluations (TPE sampler, maximise AUC).
      4. With the best hyperparameters, re-run CV collecting (true, proba)
         pairs for all samples — these form the out-of-fold predictions that
         are used for the confusion matrix and MEAS score.
      5. Refit the model on the entire training set for SHAP analysis.

    Returns: (y_true_array, y_proba_matrix, final_model, label_encoder)
    """
    le    = LabelEncoder()
    y_enc = le.fit_transform(y)

    # ── Optuna objective ──────────────────────────────────────────────────────
    def objective(trial):
        params = {
            "max_depth":        trial.suggest_int("max_depth", 2, 6),
            "learning_rate":    trial.suggest_float("learning_rate", 0.01, 0.3, log=True),
            "n_estimators":     trial.suggest_int("n_estimators", 50, 500),
            "subsample":        trial.suggest_float("subsample", 0.5, 1.0),
            "colsample_bytree": trial.suggest_float("colsample_bytree", 0.5, 1.0),
            "min_child_weight": trial.suggest_int("min_child_weight", 1, 10),
            "use_label_encoder": False,
            "eval_metric": "mlogloss",
            "random_state": seed,
        }
        clf    = xgb.XGBClassifier(**params)
        scores = []
        for tr_idx, val_idx in cv.split(X, y_enc):
            clf.fit(X.iloc[tr_idx], y_enc[tr_idx])
            prob = clf.predict_proba(X.iloc[val_idx])
            try:
                # Macro-averaged ROC-AUC handles class imbalance between groups
                auc = roc_auc_score(
                    y_enc[val_idx], prob, multi_class="ovr", average="macro"
                )
            except ValueError:
                auc = 0.0
            scores.append(auc)
        return np.mean(scores)

    # ── Run Optuna study ──────────────────────────────────────────────────────
    study = optuna.create_study(
        direction = "maximize",
        sampler   = optuna.samplers.TPESampler(seed=seed)  # reproducible tuning
    )
    study.optimize(objective, n_trials=n_trials, n_jobs=1, show_progress_bar=False)

    # ── Re-run CV with best parameters to collect out-of-fold predictions ─────
    best_params = study.best_params
    best_params.update({"use_label_encoder": False, "eval_metric": "mlogloss",
                        "random_state": seed})
    best_model = xgb.XGBClassifier(**best_params)

    all_true, all_proba = [], []
    for tr_idx, val_idx in cv.split(X, y_enc):
        best_model.fit(X.iloc[tr_idx], y_enc[tr_idx])
        prob = best_model.predict_proba(X.iloc[val_idx])
        all_true.extend(y_enc[val_idx])
        all_proba.extend(prob)

    # ── Final refit on full data for SHAP and feature importance ─────────────
    best_model.fit(X, y_enc)
    return np.array(all_true), np.array(all_proba), best_model, le


def get_cv(y_enc, strategy, folds):
    """
    Select the cross-validation scheme.
    Use LOOCV when the smallest class has < 20 samples (maximises training data).
    Otherwise use stratified k-fold (more efficient for larger cohorts).
    """
    counts = np.bincount(y_enc)
    if strategy == "loocv" or counts.min() < 20:
        return LeaveOneOut()
    return StratifiedKFold(n_splits=folds, shuffle=True, random_state=42)


# ── Classifier A: 4-class disease staging ────────────────────────────────────
print("[Module 10] Classifier A: 4-class staging")
mask_A  = labels.isin(GROUPS)
X_A, y_A = feat[mask_A], labels[mask_A]
le_A     = LabelEncoder()
y_A_enc  = le_A.fit_transform(y_A)
cv_A     = get_cv(y_A_enc, cv_strategy, cv_folds)

true_A, proba_A, model_A, le_A = tune_and_cv(X_A, y_A, cv_A, snakemake.params.random_seed)

# ── MEAS score derivation ──────────────────────────────────────────────────────
# MEAS = P(MS-Rituximab-Progressive) − P(Ctrl)
# This is a single continuous number per patient that represents where they
# sit on the disease severity spectrum, independent of their assigned group.
# Higher MEAS → more "Progressive-like" epigenome; lower MEAS → more Ctrl-like.
# MEAS is expected to correlate with EDSS, NfL, and CDI across groups.
prog_idx = list(le_A.classes_).index("MS-Rituximab-Progressive") if \
           "MS-Rituximab-Progressive" in le_A.classes_ else 0
ctrl_idx = list(le_A.classes_).index("Ctrl") if "Ctrl" in le_A.classes_ else 1

meas_df = pd.DataFrame({
    "sample_id":  X_A.index.tolist(),
    "group":      y_A.tolist(),
    "meas_score": proba_A[:, prog_idx] - proba_A[:, ctrl_idx],
})
# Merge BCI for downstream correlation analysis (MEAS vs BCI)
meas_df = meas_df.merge(
    bci[["sample_id", "bci_scaled"]], on="sample_id", how="left"
)

# ── Classifier B: binary Rituximab response (Stable vs Progressive) ────────
print("[Module 10] Classifier B: binary response")
rtx_groups = ["MS-Rituximab-Stable", "MS-Rituximab-Progressive"]
mask_B     = labels.isin(rtx_groups)
X_B, y_B   = feat[mask_B], labels[mask_B]

if len(y_B.unique()) == 2 and len(y_B) >= 4:
    # Only run if both classes are present and there are at least 4 samples total
    le_B    = LabelEncoder()
    y_B_enc = le_B.fit_transform(y_B)
    cv_B    = get_cv(y_B_enc, cv_strategy, cv_folds)
    true_B, proba_B, model_B, le_B = tune_and_cv(
        X_B, y_B, cv_B, snakemake.params.random_seed
    )
    # Binary ROC-AUC: probability of the "Progressive" class (index 1 after LE)
    roc_b = roc_auc_score(true_B, proba_B[:, 1])
    print(f"[Module 10] Classifier B ROC-AUC: {roc_b:.3f}")
else:
    model_B, proba_B, true_B, le_B = None, None, None, None
    print("[Module 10] Classifier B skipped (insufficient samples)")

# ── Step 3: Save tabular outputs ─────────────────────────────────────────────
Path(snakemake.output.meas).parent.mkdir(parents=True, exist_ok=True)
meas_df.to_csv(snakemake.output.meas, sep="\t", index=False)

# XGBoost built-in feature importance (gain-based)
imp_df = pd.Series(
    model_A.feature_importances_, index=X_A.columns
).sort_values(ascending=False).reset_index()
imp_df.columns = ["feature", "importance"]
imp_df.to_csv(snakemake.output.importance, sep="\t", index=False)

# ── Step 4: Confusion matrices ────────────────────────────────────────────────
# The confusion matrix shows how often each class is predicted correctly
# and which classes are confused. For LOOCV, this is the aggregated matrix
# across all held-out samples.
with pdf_backend.PdfPages(snakemake.output.confusion) as pp:
    # Classifier A: 4×4 matrix
    fig, ax = plt.subplots(figsize=(6, 5))
    pred_A        = le_A.inverse_transform(np.argmax(proba_A, axis=1))
    true_A_labels = le_A.inverse_transform(true_A)
    cm  = confusion_matrix(true_A_labels, pred_A, labels=le_A.classes_)
    disp = ConfusionMatrixDisplay(cm, display_labels=le_A.classes_)
    disp.plot(ax=ax, colorbar=False, cmap="Blues")
    ax.set_title("Classifier A — 4-class staging")
    plt.tight_layout()
    pp.savefig(fig)
    plt.close(fig)

    if model_B is not None:
        # Classifier B: 2×2 matrix with ROC-AUC in title
        fig, ax = plt.subplots(figsize=(4, 4))
        pred_B        = le_B.inverse_transform((proba_B[:, 1] > 0.5).astype(int))
        true_B_labels = le_B.inverse_transform(true_B)
        cm  = confusion_matrix(true_B_labels, pred_B, labels=le_B.classes_)
        disp = ConfusionMatrixDisplay(cm, display_labels=le_B.classes_)
        disp.plot(ax=ax, colorbar=False, cmap="Reds")
        ax.set_title(f"Classifier B — Rituximab response (AUC={roc_b:.2f})")
        plt.tight_layout()
        pp.savefig(fig)
        plt.close(fig)

# ── Step 5: SHAP summary plot ────────────────────────────────────────────────
# SHAP values decompose each prediction into the contribution of each feature.
# For multi-class XGBoost, shap_values is a list with one array per class.
# We take the mean absolute SHAP across classes to get an overall feature
# importance that accounts for all four disease stages simultaneously.
# The beeswarm/summary plot shows direction (positive = drives towards a class)
# and magnitude of each feature's effect.
with pdf_backend.PdfPages(snakemake.output.shap) as pp:
    explainer = shap.TreeExplainer(model_A)
    shap_vals = explainer.shap_values(X_A)

    if isinstance(shap_vals, list):
        # Multi-class: shap_vals is a list of n_classes arrays (samples × features)
        shap_mean = np.abs(np.stack(shap_vals)).mean(axis=0)
    else:
        shap_mean = np.abs(shap_vals)

    fig = plt.figure(figsize=(8, 10))
    shap.summary_plot(
        shap_vals if not isinstance(shap_vals, list) else shap_mean,
        X_A, show=False, max_display=30   # show top 30 most important features
    )
    plt.title("SHAP Feature Importance — Classifier A")
    plt.tight_layout()
    pp.savefig(fig)
    plt.close(fig)

print("[Module 10] Done.")
