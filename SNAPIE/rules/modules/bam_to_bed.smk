rule bam_to_bed:
    conda: "envs/align.yaml"
    input:
        bam=config['outputFolder'] + "/align/namesorted/{sample}.n_sorted.bam"
    output:
        bed=config['outputFolder'] + "/frags/{sample}.bed"
    params:
        dirname=config['outputFolder'] + "/frags/",
        read_method=config.get('read_method', 'PE')
    shell:
        """
        mkdir -p {params.dirname}
        if [ "{params.read_method}" = "PE" ] && [ -s {input.bam} ]; then
            bedtools bamtobed -i {input.bam} -bedpe | \
                awk 'BEGIN{{OFS="\\t";FS="\\t"}} ($1==$4){{print $1, $2, $6}}' > {output.bed} || \
                bedtools bamtobed -i {input.bam} > {output.bed}
        else
            bedtools bamtobed -i {input.bam} | \
                awk 'BEGIN{{OFS="\\t"}} {{print $1, $2, $3}}' > {output.bed} || touch {output.bed}
        fi
        """
