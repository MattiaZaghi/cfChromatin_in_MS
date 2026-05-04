rule align_bwa:
    conda: "envs/align.yaml"
    input:
        reads1=config['outputFolder'] + "/trim/{sample}_R1_trimmed.fq.gz",
        reads2=config['outputFolder'] + "/trim/{sample}_R2_trimmed.fq.gz"
    output:
        sam=temp(config['outputFolder'] + "/align/raw/{sample}.sam")
    params:
        threads=10,
        time="20:00:00",
        bwa_params=config.get('bwa_params', ''),
        dirname=config['outputFolder'] + "/align/raw/",
        genome=config.get('genome_fa', '')
    threads: 10
    resources:
        walltime=72000
    shell:
        """
        mkdir -p {params.dirname};
        bwa mem -M -t {params.threads} {params.genome} {input.reads1} {input.reads2}  {params.bwa_params} > {output.sam}; 
        """

rule sam_to_bam:
    conda: "envs/align.yaml"
    input:
        sam=config['outputFolder'] + "/align/raw/{sample}.sam"
    output:
        bam=temp(config['outputFolder'] + "/align/raw/{sample}.bam")
    params:
        threads=10,
        time="20:00:00",
        dirname=config['outputFolder'] + "/align/raw/"
    threads: 10
    resources:
        walltime=72000
    shell:
        """
        mkdir -p {params.dirname};
        samtools view {input.sam} -@ {params.threads} -Sb -o {output.bam};
        """