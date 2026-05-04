"""
Module 7 — Differential Enrichment Analysis
=============================================
7a  batch balance check
7b  RUVg latent factor estimation (R/RUVSeq)
7c  PyDESeq2 on RRE subsets and bins
7d  glmmTMB ZINB sensitivity (bins only)
"""

# ── 7a — Batch balance ────────────────────────────────────────────────────────

rule batch_balance_check:
    """Cross-tabulate batch × group; warn on confounding."""
    input:
        meta = config["sample_metadata"],
    output:
        table = "results/batch_qc/batch_group_table.tsv",
    conda:  "../../envs/python_analysis.yaml"
    script: "../../scripts/module7a_batch_check.py"

# ── 7b — RUVg ────────────────────────────────────────────────────────────────

rule ruvg_estimation:
    """
    Estimate k latent technical factors from the constitutive anchor
    count matrix (biologically invariant negative controls).
    Generates W matrix and RLE plots.
    """
    input:
        anchor_counts = "results/normalization/anchor_counts.tsv",
        meta          = config["sample_metadata"],
    output:
        factors       = "results/batch_qc/ruv_factors.tsv",
        rle_plot      = "results/batch_qc/rle_before_after.pdf",
    params:
        k = config["batch"]["ruv"]["n_factors"],
        # Pass any RRE count matrix for RLE visualisation
        rre_counts = "results/counts/rre_full_counts.tsv",
    conda:  "../../envs/r_analysis.yaml"
    script: "../../scripts/module7b_ruvg.R"

# ── 7c — DESeq2 ───────────────────────────────────────────────────────────────

DESEQ2_COUNT_FILES = {
    "full":          "results/counts/rre_full_counts.tsv",
    "cns":           "results/counts/rre_cns_counts.tsv",
    "immune":        "results/counts/rre_immune_counts.tsv",
    "gwas_proximal": "results/counts/rre_gwas_proximal_counts.tsv",
    "bins":          "results/counts/bin_counts_5kb.tsv",
    # Per-cell-type CNS subsets — added dynamically so module7 runs DESeq2 on each
    **{f"cns_{ct}": f"results/counts/rre_cns_{ct}_counts.tsv"
       for ct in config["rre"].get("cns_cell_types", [])},
}

rule run_deseq2:
    """
    PyDESeq2 differential analysis for one feature set × one contrast.
    Design includes RUVg W factors, BCI covariate, and group.
    """
    input:
        counts  = lambda wc: DESEQ2_COUNT_FILES[wc.fset],
        meta    = config["sample_metadata"],
        ruv     = "results/batch_qc/ruv_factors.tsv",
        bci     = "results/bcell_qc/bci_scores.tsv",
    output:
        all_res = "results/differential/{fset}/{contrast}_all.tsv",
        sig_res = "results/differential/{fset}/{contrast}_significant.tsv",
        volcano = "results/differential/{fset}/{contrast}_volcano.pdf",
    params:
        contrast        = lambda wc: wc.contrast.split("_vs_"),
        fdr             = config["deseq2"]["fdr"],
        lfc             = config["deseq2"]["lfc"],
        ref_level       = config["deseq2"]["ref_level"],
        min_mean_count  = config["deseq2"].get("min_mean_count", 5),
    conda:  "../../envs/python_analysis.yaml"
    threads: config["threads"]
    script: "../../scripts/module7c_deseq2.py"

# ── 7d — glmmTMB (bins only) ─────────────────────────────────────────────────

rule run_glmmtmb:
    """
    Zero-inflated negative binomial sensitivity check on bin counts.
    Adds zinb_confirmed column to the bin DESeq2 results.
    """
    input:
        bin_counts = "results/counts/bin_counts_5kb.tsv",
        meta       = config["sample_metadata"],
        ruv        = "results/batch_qc/ruv_factors.tsv",
        bci        = "results/bcell_qc/bci_scores.tsv",
        deseq2_res = expand(
            "results/differential/bins/{contrast}_all.tsv",
            contrast=CONTRASTS,
        ),
    output:
        zinb_res   = expand(
            "results/differential/bins/{contrast}_zinb.tsv",
            contrast=CONTRASTS,
        ),
    params:
        contrasts          = CONTRASTS,
        family             = config["glmmtmb"]["family"],
        ziformula          = config["glmmtmb"]["ziformula"],
        zinb_confirmed_fdr = config["glmmtmb"]["zinb_confirmed_fdr"],
    conda:  "../../envs/r_analysis.yaml"
    threads: config["threads"]
    script: "../../scripts/module7d_glmmtmb.R"
