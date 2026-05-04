"""
Module 8 — Transcription Factor Activity
==========================================
AUC-based TF activity scoring (Baca et al. 2023 approach).
TF ChIP-seq BEDs must be present in reference/tf_sites/ (hg19).
"""

rule compute_tf_activity:
    """
    For each TF, compute AUC of fragment midpoints at TF sites vs
    constitutive anchor background.  Run Kruskal-Wallis + Dunn tests.
    """
    input:
        frags_list  = expand(
            "{frags_dir}/{sample}.bed",
            frags_dir=config["data"]["frags_dir"],
            sample=SAMPLES,
        ),
        anchors     = config["normalization"]["constitutive_anchors_bed"],
        sf          = "results/normalization/constitutive_scaling_factors.tsv",
        tf_sites_dir = config["tf_activity"]["tf_sites_dir"],
        meta        = config["sample_metadata"],
    output:
        auc_matrix  = "results/tf_activity/tf_auc_scores.tsv",
        heatmap     = "results/tf_activity/tf_activity_heatmap.pdf",
        bcell_plot  = "results/tf_activity/bcell_tf_reconstitution.pdf",
    params:
        samples   = SAMPLES,
        frags_dir = config["data"]["frags_dir"],
        min_sites = config["tf_activity"]["min_sites"],
    conda:  "../../envs/python_analysis.yaml"
    threads: config["threads"]
    script: "../../scripts/module8_tf_activity.py"
