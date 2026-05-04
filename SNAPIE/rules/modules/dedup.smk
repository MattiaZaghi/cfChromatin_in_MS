rule dedup_bam:
    conda: "envs/align.yaml"
    input:
        bam=config['outputFolder'] + "/align/filtered/{sample}.filtered.unique.sorted.bam"
    output:
        RG=temp(config['outputFolder'] + "/align/dedup/{sample}.dedup.RG.bam"),
        dedup=config['outputFolder'] + "/align/dedup/{sample}.dedup.unique.sorted.bam",
        metrics=config['outputFolder'] + "/align/dedup/{sample}-MarkDuplicates.metrics.txt"
    params:
        dirname=config['outputFolder'] + "/align/dedup/"
    resources:
        mem_mb=300000
    shell:
        """
        mkdir -p {params.dirname};

        export _JAVA_OPTIONS="-Xmx16g -Xms2g"

        picard AddOrReplaceReadGroups I={input.bam} O={output.RG} RGID=1 RGLB=lib1 RGPL=illumina RGPU=unit1 RGSM=sample1

        picard MarkDuplicates I={output.RG} O={output.dedup} REMOVE_DUPLICATES=true ASSUME_SORT_ORDER=coordinate VALIDATION_STRINGENCY=LENIENT METRICS_FILE={output.metrics}

        samtools index {output.dedup}

        """
