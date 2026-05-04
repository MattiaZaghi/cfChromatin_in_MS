"""
Module 3 — B Cell Index and Reconstitution Analysis
====================================================
Computes per-sample BCI, pairwise group tests, correlations with CD19+
and months_since_rtx, and reconstitution figures.
"""

rule compute_bci:
    """
    Count fragment midpoints at bcell_rre and constitutive anchors,
    compute BCI = (frags/kb at bcell_rre) / (frags/kb at constitutive),
    run statistical tests, generate figures.
    """
    input:
        frags_list  = expand(
            "{frags_dir}/{sample}.bed",
            frags_dir=config["data"]["frags_dir"],
            sample=SAMPLES,
        ),
        bcell_rre   = config["rre"]["bcell_bed"],
        anchors     = config["normalization"]["constitutive_anchors_bed"],
        meta        = config["sample_metadata"],
    output:
        scores      = "results/bcell_qc/bci_scores.tsv",
        flagged     = "results/bcell_qc/flagged_samples.tsv",
        plots       = "results/bcell_qc/bci_reconstitution.pdf",
    params:
        samples             = SAMPLES,
        frags_dir           = config["data"]["frags_dir"],
        flag_fold           = config["bcell"]["flag_above_ctrl_fold"],
    conda:  "../../envs/python_analysis.yaml"
    threads: config["threads"]
    script: "../../scripts/module3_bcell_index.py"
