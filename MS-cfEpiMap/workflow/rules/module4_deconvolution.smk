"""
Module 4 — Cell-Type Deconvolution
====================================
Estimates per-sample cell-type contributions via Poisson z-scores at
cell-type-specific H3K27ac signature regions.  Derives CDI, NII, IAI.
"""

rule compute_deconvolution:
    """
    For each cell type, count midpoints at cell-type-specific signature
    regions (present in that type, absent in all others).  Compute
    Poisson z-scores vs Ctrl distribution.  Derive composite indices.
    Run within-sample NB-GLM enrichment vs Sadeh-style local background
    built from a curated off-target BED.
    """
    input:
        frags_list     = expand(
            "{frags_dir}/{sample}.bed",
            frags_dir=config["data"]["frags_dir"],
            sample=SAMPLES,
        ),
        rre_universe   = config["rre"]["universe_bed"],
        anchors        = config["normalization"]["constitutive_anchors_bed"],
        off_target_bed = config["deconvolution"]["off_target_bed"],   # NEW
        sf             = "results/normalization/constitutive_scaling_factors.tsv",
        anchor_matrix  = "results/normalization/anchor_counts.tsv",
        meta           = config["sample_metadata"],
    output:
        scores                     = "results/deconvolution/signature_scores.tsv",
        composite                  = "results/deconvolution/composite_indices.tsv",
        heatmap                    = "results/deconvolution/deconvolution_heatmap.pdf",
        violin                     = "results/deconvolution/deconvolution_violin.pdf",
        within_sample_scores       = "results/deconvolution/scores_within_sample.tsv",
        within_sample_qc           = "results/deconvolution/qc_within_sample.tsv",
        within_sample_log2fe       = "results/deconvolution/fig_within_sample_log2FE.pdf",
        within_sample_contribution = "results/deconvolution/fig_within_sample_contribution.pdf",
        within_sample_qc_png       = "results/deconvolution/qc_within_sample.png",
        within_sample_windows_bed  = "results/deconvolution/background_signal_windows.bed.gz",
    params:
        samples      = SAMPLES,
        frags_dir    = config["data"]["frags_dir"],
        genome_fasta = config["genome"].get("fasta", ""),
    conda:   "../../envs/python_analysis.yaml"
    threads: config["threads"]
    script:  "../../scripts/module4_deconvolution.py"
