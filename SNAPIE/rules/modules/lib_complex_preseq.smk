rule lib_complex_preseq:
    conda: "envs/common.yaml"
    input:
        sorted_bam=config['outputFolder'] + "/align/pp/{sample}.pp.sorted.bam"
    output:
        lc=config['outputFolder'] + "/align/{sample}/{sample}.lc_extrap.txt",
        versions=config['outputFolder'] + "/align/{sample}/preseq_mqc_versions.yml"
    params:
        outdir=config['outputFolder'] + "/align/{sample}/"
    threads: 1
    shell:
        """
        mkdir -p {params.outdir}
        (preseq lc_extrap -B {input.sorted_bam} > {output.lc} || \
         preseq lc_extrap -D -B {input.sorted_bam} > {output.lc}) || touch {output.lc}
        cat <<-END_VERSIONS > {output.versions}
        "lib_complex_preseq":
          preseq: $(preseq 2>&1 | grep Version | sed 's/.*Version: //' || echo 'not-detected')
        END_VERSIONS
        """

