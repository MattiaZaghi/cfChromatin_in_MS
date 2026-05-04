rule unique_sam:
    conda: "envs/align.yaml"
    input:
        bam=config['outputFolder'] + "/align/pp/{sample}.pp.sorted.bam"
    output:
        bam=temp(config['outputFolder'] + "/align/unique/{sample}.unique.sorted.bam")
    threads: 4
    shell:
        """
        mkdir -p $(dirname {output.bam})
        samtools view --threads {threads} -b -q 1 {input.bam} -o {output.bam}
        """
