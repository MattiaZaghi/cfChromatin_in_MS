# deeptools bamCoverage — the only bigWig generation step in the pipeline.
# The bedgraph → bedGraphToBigWig path has been removed; deeptools handles
# normalisation, binning, and smoothing in a single step directly from the BAM.

rule coverage_deeptools:
    conda: "envs/common.yaml"
    input:
        bam=lambda wildcards: (
            config['outputFolder'] + f"/align/dac/{wildcards.sample}.dac_filtered.dedup.unique.sorted.bam"
            if config.get('exclude_dac_regions', False) else
            config['outputFolder'] + f"/align/dedup/{wildcards.sample}.dedup.unique.sorted.bam"
        ),
        bai=lambda wildcards: (
            config['outputFolder'] + f"/align/dac/{wildcards.sample}.dac_filtered.dedup.unique.sorted.bam.bai"
            if config.get('exclude_dac_regions', False) else
            config['outputFolder'] + f"/align/dedup/{wildcards.sample}.dedup.unique.sorted.bam.bai"
        )
    output:
        bw=config['outputFolder'] + "/bigwig/{sample}.bw"
    params:
        binsize=config.get('binsize', '10'),
        norm=config.get('norm_method', 'RPKM'),
        smooth=config.get('smooth_length', '300')
    threads: 20
    shell:
        """
        mkdir -p $(dirname {output.bw})
        bamCoverage -b {input.bam} \
            --outFileName {output.bw} \
            --normalizeUsing {params.norm} \
            --binSize {params.binsize} \
            --smoothLength {params.smooth} \
            --numberOfProcessors {threads} \
            --exactScaling || true
        """
