rule bam_to_bedgraph:
    conda: "envs/align.yaml"
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
        bedgraph=config['outputFolder'] + "/bedgraph/{sample}.bedgraph"
    params:
        dirname=config['outputFolder'] + "/bedgraph/",
        read_method=config.get('read_method', 'PE')
    threads: 4
    shell:
        """
        mkdir -p {params.dirname}
        if [ "{params.read_method}" = "PE" ]; then
            bedtools genomecov -ibam {input.bam} -bg -pc > {output.bedgraph} || true
        else
            bedtools genomecov -ibam {input.bam} -bg > {output.bedgraph} || true
        fi
        """
