"""
Module 7c — PyDESeq2 differential analysis.
Snakemake script: called via script: directive.

PURPOSE
Identify regulatory elements (RREs or 5 kb bins) with significantly
different H3K27ac cfDNA signal between two groups. This script runs a
single contrast on a single count matrix; Snakemake calls it once per
combination of (contrast, feature_set).

STATISTICAL APPROACH
We use the DESeq2 negative-binomial model because:
  1. Raw counts at regulatory elements are overdispersed integers —
     the NB model is the standard for this type of data.
  2. DESeq2 computes size factors internally from the count matrix
     itself, so we do NOT pre-scale the input.
  3. The design formula includes:
       W_1, W_2   — RUVg latent factors (from Module 7b) capturing
                    technical batch variation data-adaptively.
       sex        — biological confound.
       age_group  — binned age (young ≤35 / middle 36–50 / older >50).
       rituximab_treated — binary flag separating treated from naive/ctrl.
       bci_scaled — continuous B Cell Index / Ctrl mean (from Module 3),
                    removing residual B-cell reconstitution variance so
                    that the group contrast reflects disease biology, not
                    B-cell cfDNA contamination.
       group      — the factor of interest; reference level = Ctrl.

  If any group has < 10 samples, the design is reduced to:
    ~ W_1 + W_2 + rituximab_treated + bci_scaled + group
  to avoid over-parameterisation.

FOUR CONTRASTS (one Snakemake call per contrast):
  1. MS-Rituximab-Progressive vs MS-Rituximab-Stable  (non-response mechanism)
  2. NEW vs Ctrl                                       (pure disease signal)
  3. MS-Rituximab-Stable vs NEW                        (drug effect)
  4. MS-Rituximab-Progressive vs Ctrl                  (cumulative burden)

SIGNIFICANCE THRESHOLD: FDR (Benjamini-Hochberg) < 0.05, |log2FC| > 1.
"""

import warnings
warnings.filterwarnings("ignore")

import pandas as pd
import numpy as np
from pathlib import Path
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from pydeseq2.dds import DeseqDataSet
from pydeseq2.ds  import DeseqStats


# ── Step 1: Load count matrix, metadata, RUVg factors, and BCI ───────────────
# counts_df: raw integer count matrix (regions × samples) from Module 5 or 6.
# meta:      sample metadata with clinical covariates.
# ruv:       RUVg W matrix (samples × k latent factors) from Module 7b.
# bci:       per-sample scaled BCI values from Module 3.
counts_df = pd.read_csv(snakemake.input.counts, sep="\t", index_col=0)
meta      = pd.read_csv(snakemake.input.meta,   sep="\t", comment="#")
ruv       = pd.read_csv(snakemake.input.ruv,    sep="\t")
bci       = pd.read_csv(snakemake.input.bci,    sep="\t")

# ── Step 2: Build the clinical design frame ───────────────────────────────────
# Start from QC-passing samples only, then merge RUVg factors and BCI scaled.
meta = meta[meta["qc_include"].astype(str).str.upper() == "TRUE"].copy()

# Merge RUVg latent factors (columns named W_1, W_2, …)
meta = meta.merge(
    ruv[["sample_id"] + [c for c in ruv.columns if c.startswith("W_")]],
    on="sample_id", how="left",
)
# Merge the scaled BCI covariate computed in Module 3
meta = meta.merge(
    bci[["sample_id", "bci_scaled"]], on="sample_id", how="left",
)

# ── Step 3: Derive computed covariates ───────────────────────────────────────
# rituximab_treated: binary indicator (1 = either Rituximab group)
# Used to partially capture the shared effect of Rituximab across both treated
# groups before estimating the Stable vs Progressive difference.
meta["rituximab_treated"] = meta["group"].isin(
    ["MS-Rituximab-Stable", "MS-Rituximab-Progressive"]
).astype(int)

# age_group: continuous age binned to 3 levels so DESeq2 can treat it as a
# categorical covariate without assuming a linear effect.
def age_group(x):
    try:
        v = float(x)
        if v <= 35:   return "young"
        if v <= 50:   return "middle"
        return "older"
    except (ValueError, TypeError):
        return "NA"

meta["age_group"] = meta["age"].apply(age_group)

# ── Step 4: Align samples between count matrix and metadata ──────────────────
# Only samples present in BOTH the count matrix and the metadata are kept.
# This handles cases where a sample failed QC after counts were generated.
common    = [s for s in counts_df.columns if s in meta["sample_id"].values]
counts_df = counts_df[common]
meta      = meta[meta["sample_id"].isin(common)].set_index("sample_id")
meta      = meta.loc[common]   # enforce same ordering as count matrix columns

# ── Step 4b: Pre-filter low-count features ────────────────────────────────────
# Features with very low mean counts across all samples have essentially zero
# power to detect differential signal but consume multiple-testing budget.
# Keeping them inflates the BH correction, making the volcano appear flat.
# We filter: keep features where at least min_samples have count > 0 AND
# the row mean >= min_mean_count. This matches DESeq2's independent filtering
# spirit but is applied before model fitting for a cleaner multiple-test space.
min_mean_count = getattr(snakemake.params, "min_mean_count", 5)
n_samples = counts_df.shape[1]
row_means  = counts_df.mean(axis=1)
row_nonzero = (counts_df > 0).sum(axis=1)
keep = (row_means >= min_mean_count) & (row_nonzero >= max(2, n_samples * 0.1))
n_before = counts_df.shape[0]
counts_df = counts_df[keep]
print(f"[Module 7c] Pre-filter: kept {counts_df.shape[0]:,}/{n_before:,} features "
      f"(mean >= {min_mean_count}, nonzero in ≥10% of samples)")

# ── Step 5: Choose full or reduced design ────────────────────────────────────
# If any group has < 10 samples, the full design with sex and age_group would
# have too many parameters relative to sample size. The reduced design uses
# fewer covariates and is more appropriate for small cohorts.
group_counts = meta["group"].value_counts()
use_reduced  = (group_counts < 10).any()

w_cols = sorted([c for c in meta.columns if c.startswith("W_")])

# rituximab_treated is intentionally excluded: with all 4 groups present it is
# a perfect linear combination of the group dummies (Stable + Progressive),
# making the design matrix singular.  group already captures this information.
if use_reduced or not w_cols:
    design_factors = w_cols + ["bci_scaled", "group"]
    print("[Module 7c] Using reduced design (n < 10 in some group or no W factors)")
else:
    design_factors = w_cols + ["sex", "age_group", "bci_scaled", "group"]
    print("[Module 7c] Using full design")

# ── Step 6: Validate and clean design covariates ─────────────────────────────
# Drop any covariate that is missing from metadata or is constant across samples
# (a constant covariate provides no information and breaks DESeq2 model fitting).
valid_factors = []
for f in design_factors:
    if f not in meta.columns:
        print(f"[Module 7c] WARNING: covariate '{f}' not in metadata, skipping")
        continue
    col = meta[f].replace("NA", np.nan).dropna()
    if col.nunique() > 1:
        valid_factors.append(f)
    else:
        print(f"[Module 7c] Dropping constant covariate: {f}")

design_factors = valid_factors
# group must always be in the design (it is the factor being tested)
if "group" not in design_factors:
    design_factors.append("group")

# Fill NA values in numeric covariates with column mean so DESeq2 does not
# crash on missing clinical data. Imputing with the mean is conservative —
# it effectively treats missing samples as average for that covariate.
for f in design_factors:
    if f in meta.columns:
        meta[f] = meta[f].replace("NA", np.nan)
        if meta[f].dtype in (float, int) or pd.to_numeric(meta[f], errors="coerce").notna().any():
            meta[f] = pd.to_numeric(meta[f], errors="coerce")
            meta[f] = meta[f].fillna(meta[f].mean())

# ── Step 7: Set factor levels and types expected by pydeseq2 ─────────────────
# group must be a pd.Categorical with Ctrl as the reference (leftmost) level.
# sex and age_group must be string-typed for DESeq2's categorical encoding.
meta["group"] = pd.Categorical(
    meta["group"],
    categories=["Ctrl", "NEW", "MS-Rituximab-Stable", "MS-Rituximab-Progressive"],
)
meta["sex"]       = meta["sex"].astype(str)
meta["age_group"] = meta["age_group"].astype(str)

# pydeseq2 expects counts as a samples × features integer matrix.
# Transpose: counts_df is regions × samples → counts_t is samples × regions.
counts_t = counts_df.T.astype(int)
counts_t = counts_t.loc[meta.index]

# ── Step 8: Fit the DESeq2 model ──────────────────────────────────────────────
# DeseqDataSet runs:
#   1. Estimation of size factors (internal library-size normalisation).
#   2. Estimation of feature-wise dispersion (variance as a function of mean).
#   3. GLM fitting of the negative-binomial model with the design formula.
#   4. Cook's distance outlier replacement (refit_cooks=True) to handle single
#      extreme samples that would otherwise dominate a feature's LFC estimate.
print(f"[Module 7c] Running DESeq2: {counts_t.shape[0]} samples × "
      f"{counts_t.shape[1]} features, design: {design_factors}")

dds = DeseqDataSet(
    counts         = counts_t,
    metadata       = meta,
    design_factors = design_factors,
    ref_level      = snakemake.params.ref_level,   # ["group", "Ctrl"]
    n_cpus         = snakemake.threads,
    refit_cooks    = True,
)
dds.deseq2()

# ── Step 9: Extract the requested contrast ────────────────────────────────────
# DeseqStats computes Wald-test statistics and Benjamini-Hochberg adjusted
# p-values for the specified contrast: group1 vs group2.
# alpha sets the FDR threshold used by pydeseq2 for its internal reporting
# (we apply our own threshold again when writing the significant results file).
contrast = snakemake.params.contrast   # e.g. ["MS-Rituximab-Progressive", "Ctrl"]
print(f"[Module 7c] Contrast: {contrast[0]} vs {contrast[1]}")

stat_res = DeseqStats(
    dds,
    contrast = ["group", contrast[0], contrast[1]],
    alpha    = snakemake.params.fdr,
    n_cpus   = snakemake.threads,
)
stat_res.summary()
results  = stat_res.results_df.copy()
results.index.name = "region_id"

# ── Step 10: Save all results and the significant subset ──────────────────────
# All features: baseMean, log2FoldChange, lfcSE, stat, pvalue, padj for every
#   region tested (including those that failed to converge — padj = NaN).
# Significant: subset where padj < fdr AND |log2FC| > lfc (both criteria must
#   hold; LFC filters out features that are statistically significant but too
#   small to be biologically meaningful).
Path(snakemake.output.all_res).parent.mkdir(parents=True, exist_ok=True)
results.to_csv(snakemake.output.all_res, sep="\t")

sig = results[
    (results["padj"] < snakemake.params.fdr) &
    (results["log2FoldChange"].abs() > snakemake.params.lfc)
]
sig.to_csv(snakemake.output.sig_res, sep="\t")
print(f"[Module 7c] Significant: {len(sig)} features "
      f"(FDR<{snakemake.params.fdr}, |LFC|>{snakemake.params.lfc})")

# ── Step 11: Volcano plot ─────────────────────────────────────────────────────
# x-axis: log2FoldChange.  y-axis: -log10(padj) (higher = more significant).
# Significant features (red) are those passing both FDR and LFC thresholds.
# Dashed lines mark the significance and LFC cut-offs.
# rasterized=True keeps file size manageable when there are hundreds of thousands
# of points (full RRE universe).
fig, ax = plt.subplots(figsize=(6, 5))
plot_df = results.dropna(subset=["log2FoldChange", "padj"])
neg_log = -np.log10(plot_df["padj"].clip(lower=1e-300))

is_sig = (plot_df["padj"] < snakemake.params.fdr) & \
         (plot_df["log2FoldChange"].abs() > snakemake.params.lfc)

ax.scatter(plot_df.loc[~is_sig, "log2FoldChange"], neg_log[~is_sig],
           s=2, alpha=0.4, color="grey", rasterized=True)
ax.scatter(plot_df.loc[is_sig,  "log2FoldChange"], neg_log[is_sig],
           s=4, alpha=0.8, color="#D65F5F", rasterized=True)
ax.axhline(-np.log10(snakemake.params.fdr), ls="--", lw=0.8, color="navy",
           label=f"FDR={snakemake.params.fdr}")
ax.axvline(-snakemake.params.lfc, ls="--", lw=0.8, color="green",
           label=f"|LFC|={snakemake.params.lfc}")
ax.axvline( snakemake.params.lfc, ls="--", lw=0.8, color="green")
ax.set_xlabel("log2 Fold Change")
ax.set_ylabel("-log10(padj)")
ax.set_title(f"{contrast[0]} vs {contrast[1]}  (n={len(sig)} sig)")
plt.tight_layout()
fig.savefig(snakemake.output.volcano, dpi=150)
plt.close(fig)
