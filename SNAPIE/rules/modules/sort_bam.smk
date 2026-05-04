rule sort_bam:
    conda: "envs/align.yaml"
    input:
           bam=config['outputFolder'] + "/align/raw/{sample}.bam"
    output:
           sorted=temp(config['outputFolder'] + "/align/sorted/{sample}.sorted.bam")
    params:
        dirname=config['outputFolder'] + "/align/"
    threads: 10
    shell:
        """
        mkdir -p {params.dirname};
        samtools sort -@ {threads} {input.bam} -o {output.sorted} || true
        """

rule sort_readname_bam:
    conda: "envs/align.yaml"
    input:
        bam=lambda wildcards: (
            config['outputFolder'] + f"/align/dac/{wildcards.sample}.dac_filtered.dedup.unique.sorted.bam"
            if config.get('exclude_dac_regions', False) else
            config['outputFolder'] + f"/align/dedup/{wildcards.sample}.dedup.unique.sorted.bam"
        )
    output:
        bam=temp(config['outputFolder'] + "/align/namesorted/{sample}.n_sorted.bam")
    params:
        read_method=config.get('read_method', 'PE')
    threads: 10
    shell:
        """
        mkdir -p $(dirname {output.bam})
        if [ "{params.read_method}" = "PE" ]; then
            samtools sort -@ {threads} -n {input.bam} -o {output.bam}
        else
            touch {output.bam}
        fi
        """
