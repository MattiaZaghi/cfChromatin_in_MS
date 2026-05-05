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
    """
    input:
        frags_list   = expand(
            "{frags_dir}/{sample}.bed",
            frags_dir=config["data"]["frags_dir"],
            sample=SAMPLES,
        ),
        rre_universe = config["rre"]["universe_bed"],
        anchors      = config["normalization"]["constitutive_anchors_bed"],
        sf             = "results/normalization/constitutive_scaling_factors.tsv",
        anchor_matrix  = "results/normalization/anchor_counts.tsv",
        meta           = config["sample_metadata"],
    output:
        scores       = "results/deconvolution/signature_scores.tsv",
        composite    = "results/deconvolution/composite_indices.tsv",
        heatmap      = "results/deconvolution/deconvolution_heatmap.pdf",
        violin       = "results/deconvolution/deconvolution_violin.pdf",
    params:
        samples   = SAMPLES,
        frags_dir = config["data"]["frags_dir"],
    conda:  "../../envs/python_analysis.yaml"
    threads: config["threads"]
    script: "../../scripts/module4_deconvolution.py"
