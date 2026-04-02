import os as _os

rule filter_bam_fragle:
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
        bam=config['outputFolder'] + "/align/fragle/{sample}.filtered.fragle.bam",
        bai=config['outputFolder'] + "/align/fragle/{sample}.filtered.fragle.bam.bai"
    params:
        fragle_sites=config.get('fragle_sites_ref', ''),
        enrichment_mark=config.get('enrichment_mark', '')
    shell:
        """
        mkdir -p $(dirname {output.bam})
        SITES_BED="{params.fragle_sites}/{params.enrichment_mark}/sites.bed"
        if [ -f "$SITES_BED" ]; then
            samtools view -b -L "$SITES_BED" {input.bam} > {output.bam}
        else
            cp {input.bam} {output.bam}
            cp {input.bai} {output.bai} 2>/dev/null || samtools index {output.bam}
            exit 0
        fi
        samtools index {output.bam}
        """
