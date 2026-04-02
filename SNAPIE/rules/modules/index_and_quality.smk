# Index the DAC-filtered BAM (used when exclude_dac_regions=True)
rule index_dac_bam:
    conda: "envs/align.yaml"
    input:
        bam=config['outputFolder'] + "/align/dac/{sample}.dac_filtered.dedup.unique.sorted.bam"
    output:
        bai=config['outputFolder'] + "/align/dac/{sample}.dac_filtered.dedup.unique.sorted.bam.bai"
    threads: 4
    shell:
        """
        samtools index -@ {threads} {input.bam}
        """

# Index the plain dedup BAM (used when exclude_dac_regions=False)
rule index_dedup_bam:
    conda: "envs/align.yaml"
    input:
        bam=config['outputFolder'] + "/align/dedup/{sample}.dedup.unique.sorted.bam"
    output:
        bai=config['outputFolder'] + "/align/dedup/{sample}.dedup.unique.sorted.bam.bai"
    threads: 4
    shell:
        """
        samtools index -@ {threads} {input.bam}
        """
