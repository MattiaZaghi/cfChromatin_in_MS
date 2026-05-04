"""
Module 6 — Genome-wide Bin Count Matrix
=========================================
Tiles the genome into 5 kb bins, optionally restricts to RRE-overlapping
bins, counts fragment midpoints per bin per sample, filters by mean ≥ 5.
"""

rule make_bin_universe:
    """Generate 5 kb genome-tiling BED, remove DAC regions."""
    input:
        chrom_sizes = config["genome"]["chrom_sizes"],
        rre         = config["rre"]["universe_bed"],
    output:
        bins = "results/counts/bin_universe.bed",
    params:
        bin_size      = config["bins"]["size_bp"],
        restrict_rre  = config["bins"]["restrict_to_rre"],
        dac           = config["genome"]["dac_regions"],
        script        = workflow.basedir + "/scripts/module6_make_bins.py",
    conda:  "../../envs/python_analysis.yaml"
    shell:
        """
        python {params.script} \
            --chrom_sizes  {input.chrom_sizes} \
            --dac          {params.dac} \
            --rre          {input.rre} \
            --bin_size     {params.bin_size} \
            --restrict_rre {params.restrict_rre} \
            --output       {output.bins}
        """

rule count_bin_midpoints_per_sample:
    """Count midpoints overlapping genome bins for one sample."""
    input:
        frags = lambda wc: f"{config['data']['frags_dir']}/{wc.sample}.bed",
        bins  = "results/counts/bin_universe.bed",
    output:
        temp("results/counts/bins_per_sample/{sample}.counts"),
    params:
        script = workflow.basedir + "/scripts/count_midpoints.py",
    threads: 1
    conda:  "../../envs/python_analysis.yaml"
    shell:
        """
        python {params.script} \
            --frags   {input.frags} \
            --regions {input.bins} \
            --output  {output}
        """

rule build_bin_count_matrix:
    """Assemble bin counts matrix and apply mean-count filter."""
    input:
        counts = expand(
            "results/counts/bins_per_sample/{sample}.counts",
            sample=SAMPLES,
        ),
        bins   = "results/counts/bin_universe.bed",
    output:
        matrix = "results/counts/bin_counts_5kb.tsv",
    params:
        samples        = SAMPLES,
        min_mean_count = config["bins"]["min_mean_count"],
    conda:  "../../envs/python_analysis.yaml"
    script: "../../scripts/module6_bin_counts.py"
