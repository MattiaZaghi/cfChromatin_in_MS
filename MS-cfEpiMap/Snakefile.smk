"""
MS-cfEpiMap — H3K27ac cfChIP-seq analysis pipeline
====================================================
Run from the MS-cfEpiMap/ directory:

    snakemake --cores 8 --use-conda

Module 1 (RRE universe) must be built once before the main pipeline:

    snakemake --cores 4 --use-conda build_rre_universe
"""

configfile: "config/config.yaml"

# All relative paths in every rule file resolve under outdir.
# Change outdir in config/config.yaml to redirect all output.
workdir: config["outdir"]

import pandas as pd

# ── Load sample sheet ─────────────────────────────────────────────────────────
meta = pd.read_csv(config["sample_metadata"], sep="\t", comment="#")
if "qc_include" in meta.columns:
    meta = meta[meta["qc_include"].astype(str).str.upper() == "TRUE"]

SAMPLES   = meta["sample_id"].tolist()
SAMPLES_P = meta.loc[meta["sample_type"] == "plasma", "sample_id"].tolist()

GROUPS = ["Ctrl", "NEW", "MS-Rituximab-Stable", "MS-Rituximab-Progressive"]

CONTRASTS = [
    f"{c[0]}_vs_{c[1]}"
    for c in config["deseq2"]["contrasts"]
]

CNS_CELL_TYPES = config["rre"]["cns_cell_types"]

# Only include RRE subsets whose BED file is non-empty
def _rre_bed(subset):
    """Return the absolute path for a given RRE subset key."""
    outdir = config["outdir"]
    mapping = {
        "full":          config["rre"]["universe_bed"],
        "cns":           config["rre"]["cns_bed"],
        "immune":        config["rre"]["immune_bed"],
        "gwas_proximal": config["rre"]["gwas_proximal_bed"],
    }
    for ct in CNS_CELL_TYPES:
        mapping[f"cns_{ct}"] = config["rre"]["cns_celltype_beds"][ct]
    rel = mapping.get(subset, "")
    return os.path.join(outdir, rel) if rel else ""

_ALL_RRE = ["full", "cns", "immune", "gwas_proximal"] + [f"cns_{ct}" for ct in CNS_CELL_TYPES]
RRE_SUBSETS  = [s for s in _ALL_RRE
                if _rre_bed(s) and os.path.exists(_rre_bed(s))
                and os.path.getsize(_rre_bed(s)) > 0]
FEATURE_SETS = RRE_SUBSETS + ["bins"]

if set(_ALL_RRE) - set(RRE_SUBSETS):
    print(f"Skipping empty RRE subsets: {sorted(set(_ALL_RRE) - set(RRE_SUBSETS))}")

# ── Include rule modules ──────────────────────────────────────────────────────
include: "workflow/rules/module1_rre_universe.smk"
include: "workflow/rules/module2_normalization.smk"
include: "workflow/rules/module3_bcell_index.smk"
include: "workflow/rules/module4_deconvolution.smk"
include: "workflow/rules/module5_rre_counts.smk"
include: "workflow/rules/module6_bin_counts.smk"
include: "workflow/rules/module7_differential.smk"
include: "workflow/rules/module8_tf_activity.smk"
include: "workflow/rules/module9_gsea.smk"
include: "workflow/rules/module10_ml.smk"


# ── Default target ────────────────────────────────────────────────────────────
rule all:
    input:
        # Module 2 — normalization
        "results/normalization/constitutive_scaling_factors.tsv",
        "results/normalization/anchor_counts.tsv",
        expand("results/normalization/bigwigs/{sample}_anchnorm.bw", sample=SAMPLES),
        "results/normalization/anchor_qc.pdf",
        # Module 3 — B cell index
        "results/bcell_qc/bci_scores.tsv",
        "results/bcell_qc/flagged_samples.tsv",
        "results/bcell_qc/bci_reconstitution.pdf",
        # Module 4 — deconvolution
        "results/deconvolution/signature_scores.tsv",
        "results/deconvolution/composite_indices.tsv",
        "results/deconvolution/deconvolution_heatmap.pdf",
        # Module 5 — RRE count matrices
        expand("results/counts/rre_{subset}_counts.tsv", subset=RRE_SUBSETS),
        # Module 6 — bin count matrix
        "results/counts/bin_counts_5kb.tsv",
        "results/counts/bin_universe.bed",
        # Module 7 — differential analysis
        "results/batch_qc/batch_group_table.tsv",
        "results/batch_qc/ruv_factors.tsv",
        "results/batch_qc/rle_before_after.pdf",
        expand("results/differential/{fset}/{contrast}_all.tsv",
               fset=FEATURE_SETS, contrast=CONTRASTS),
        expand("results/differential/{fset}/{contrast}_significant.tsv",
               fset=FEATURE_SETS, contrast=CONTRASTS),
        # Module 8 — TF activity
        "results/tf_activity/tf_auc_scores.tsv",
        "results/tf_activity/tf_activity_heatmap.pdf",
        "results/tf_activity/bcell_tf_reconstitution.pdf",
        # Module 9 — GSEA
        expand("results/gsea/{contrast}_dotplot.pdf", contrast=CONTRASTS),
        # Module 10 — ML
        "results/ml/meas_scores.tsv",
        "results/ml/confusion_matrices.pdf",
        "results/ml/shap_summary.pdf",
        "results/ml/feature_importance.tsv",
        # QC summary + HTML report
        "results/qc_summary.tsv",
        "results/reports/ms_cfchipmaps_report.html",
