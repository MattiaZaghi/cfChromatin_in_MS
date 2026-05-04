"""
Module 2 — Constitutive Anchor Normalization
=============================================
Counts fragment midpoints at constitutive anchor regions per sample,
computes scaling factors, generates anchor-normalised bigWigs and QC plot.
"""

rule count_anchor_midpoints:
    """Count fragment midpoints overlapping constitutive anchor regions."""
    input:
        frags   = lambda wc: f"{config['data']['frags_dir']}/{wc.sample}.bed",
        anchors = config["normalization"]["constitutive_anchors_bed"],
    output:
        temp("results/normalization/anchor_counts_per_sample/{sample}.counts"),
    params:
        script = workflow.basedir + "/scripts/count_midpoints.py",
    threads: 1
    conda:  "../../envs/python_analysis.yaml"
    shell:
        """
        python {params.script} \
            --frags   {input.frags} \
            --regions {input.anchors} \
            --output  {output}
        """

rule compute_scaling_factors:
    """Assemble anchor count matrix and compute per-sample scaling factors."""
    input:
        counts  = expand(
            "results/normalization/anchor_counts_per_sample/{sample}.counts",
            sample=SAMPLES,
        ),
        anchors = config["normalization"]["constitutive_anchors_bed"],
        meta    = config["sample_metadata"],
    output:
        sf      = "results/normalization/constitutive_scaling_factors.tsv",
        matrix  = "results/normalization/anchor_counts.tsv",
    params:
        samples = SAMPLES,
        flag_sd = config["normalization"]["flag_sf_sd_threshold"],
    conda:  "../../envs/python_analysis.yaml"
    script: "../../scripts/module2_normalization.py"

rule generate_bigwig:
    """Anchor-normalised bigWig for one sample via bamCoverage."""
    input:
        bam = lambda wc: (
            f"{config['data']['bam_dir']}/{wc.sample}"
            f"{config['data']['bam_suffix']}"
        ),
        bai = lambda wc: (
            f"{config['data']['bam_dir']}/{wc.sample}"
            f"{config['data']['bam_suffix']}.bai"
        ),
        sf  = "results/normalization/constitutive_scaling_factors.tsv",
    output:
        "results/normalization/bigwigs/{sample}_anchnorm.bw",
    threads: 4
    conda:  "../../envs/python_analysis.yaml"
    shell:
        """
        sf=$(awk -v s={wildcards.sample} 'NR>1 && $1==s {{print $3}}' {input.sf})
        bamCoverage \
            -b {input.bam} \
            -o {output} \
            --scaleFactor  $sf \
            --binSize      10 \
            --smoothLength 75 \
            --numberOfProcessors {threads} \
            --normalizeUsing None
        """

rule anchor_qc_plot:
    """Barplot of anchor read counts and scaling factors coloured by group."""
    input:
        sf   = "results/normalization/constitutive_scaling_factors.tsv",
        meta = config["sample_metadata"],
    output:
        "results/normalization/anchor_qc.pdf",
    conda:  "../../envs/r_analysis.yaml"
    script: "../../scripts/anchor_qc_plot.R"
