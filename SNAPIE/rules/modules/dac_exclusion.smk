import os as _os

rule dac_exclusion:
    conda: "envs/align.yaml"
    input:
        bam=config['outputFolder'] + "/align/dedup/{sample}.dedup.unique.sorted.bam",
        dac=_os.path.join(
            config.get('ref_dir', 'ref_files'),
            config.get('genome', 'hg19') + '.DAC.bed'
        )
    output:
        bam=config['outputFolder'] + "/align/dac/{sample}.dac_filtered.dedup.unique.sorted.bam"
    shell:
        """
        mkdir -p $(dirname {output.bam})
        if [ -s {input.dac} ]; then
            bedtools intersect -v -abam {input.bam} -b {input.dac} > {output.bam}
        else
            cp {input.bam} {output.bam}
        fi
        """
