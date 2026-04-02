rule bedgraph_to_bigwig:
    conda: "envs/align.yaml"
    input:
        bedgraph=config['outputFolder'] + "/bedgraph/{sample}.bedgraph",
        chromsizes=config.get('chrom_sizes', '')
    output:
        bw=config['outputFolder'] + "/bigwig/{sample}.bw"
    params:
        dirname=config['outputFolder'] + "/bigwig/"
    shell:
        """
        mkdir -p {params.dirname};
        bedGraphToBigWig {input.bedgraph} {input.chromsizes} {output.bw} || true
        """
rule coverage_deeptools:
    conda: "envs/common.yaml"
    input:
        dedup=lambda wildcards: (
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
        bw_deeptools=config['outputFolder'] + "/bigwig/deeptools/{sample}.bw"
    params:
        dirname=config['outputFolder'] + "/bigwig/deeptools/",
        mapping_qual_bw = config['binsize'],
        norm= config['norm_method'],
        smooth= config['smooth_length']
    threads: 20
    shell:
        """
        mkdir -p {params.dirname};
        bamCoverage -b {input.dedup} --outFileName {output.bw_deeptools} --normalizeUsing {params.norm} --binSize {params.mapping_qual_bw} --smoothLength {params.smooth} --numberOfProcessors {threads} --exactScaling   || true
        """
