rule filter_properly_paired:
    conda: "envs/align.yaml"
    input:
        bam=config['outputFolder'] + "/align/sorted/{sample}.sorted.bam"
    output:
        bam=temp(config['outputFolder'] + "/align/pp/{sample}.pp.sorted.bam")
    params:
        read_method=config.get('read_method', 'PE')
    threads: 4
    shell:
        """
        mkdir -p $(dirname {output.bam})
        if [ "{params.read_method}" = "PE" ]; then
            samtools view -b -f 2 -@ {threads} {input.bam} -o {output.bam}
        else
            cp {input.bam} {output.bam}
        fi
        """
