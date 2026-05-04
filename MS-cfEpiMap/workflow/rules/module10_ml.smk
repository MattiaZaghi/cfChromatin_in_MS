"""
Module 10 — Machine Learning Scoring
======================================
XGBoost classifiers with Optuna hyperparameter search and SHAP.
Produces the MEAS (MS Epigenomic Activity Score) per patient.
"""

rule run_ml:
    """
    Build two classifiers:
      A — 4-class disease staging (Ctrl / NEW / Stable / Progressive)
      B — binary Rituximab response (Stable vs Progressive)

    Features: top DESeq2 RRE hits, cell-type z-scores (Module 4),
    TF AUC scores (Module 8), composite indices (BCI, CDI, NII, IAI).
    Derive MEAS score from Classifier A probabilities.
    """
    input:
        # Top differential features
        deseq2_full = expand(
            "results/differential/full/{contrast}_significant.tsv",
            contrast=CONTRASTS,
        ),
        rre_counts   = "results/counts/rre_full_counts.tsv",
        # Cell-type scores
        composite    = "results/deconvolution/composite_indices.tsv",
        deconv       = "results/deconvolution/signature_scores.tsv",
        # TF activity
        tf_auc       = "results/tf_activity/tf_auc_scores.tsv",
        # BCI
        bci          = "results/bcell_qc/bci_scores.tsv",
        # Metadata
        meta         = config["sample_metadata"],
    output:
        meas         = "results/ml/meas_scores.tsv",
        confusion    = "results/ml/confusion_matrices.pdf",
        shap         = "results/ml/shap_summary.pdf",
        importance   = "results/ml/feature_importance.tsv",
    params:
        random_seed     = config["ml"]["random_seed"],
        n_optuna_trials = config["ml"]["n_optuna_trials"],
        cv_strategy     = config["ml"]["cv_strategy"],
        cv_folds        = config["ml"]["cv_folds"],
        max_features    = config["ml"]["max_features"],
    conda:  "../../envs/python_analysis.yaml"
    threads: config["threads"]
    script: "../../scripts/module10_ml.py"

rule generate_report:
    """Aggregate all results into a single HTML QC + analysis report."""
    input:
        sf          = "results/normalization/constitutive_scaling_factors.tsv",
        anchor_qc   = "results/normalization/anchor_qc.pdf",
        bci         = "results/bcell_qc/bci_scores.tsv",
        batch_table = "results/batch_qc/batch_group_table.tsv",
        rle_plot    = "results/batch_qc/rle_before_after.pdf",
        deconv      = "results/deconvolution/composite_indices.tsv",
        diff_sig    = expand(
            "results/differential/full/{contrast}_significant.tsv",
            contrast=CONTRASTS,
        ),
        meas        = "results/ml/meas_scores.tsv",
        meta        = config["sample_metadata"],
    output:
        html        = "results/reports/ms_cfchipmaps_report.html",
        qc_summary  = "results/qc_summary.tsv",
    params:
        contrasts   = CONTRASTS,
    conda:  "../../envs/python_analysis.yaml"
    script: "../../scripts/generate_report.py"
