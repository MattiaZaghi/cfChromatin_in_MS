rule calcFragsLength:
    conda: "envs/common.yaml"
    input:
        sorted_bam=lambda wildcards: (
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
        frag_txt=config['outputFolder'] + "/frags/{sample}/{sample}.fragment_sizes.txt",
        versions=config['outputFolder'] + "/frags/{sample}/bamPEFragmentSize_mqc_versions.yml"
    params:
        read_method=config.get('read_method', 'PE')
    threads: 1
    shell:
        """
        mkdir -p $(dirname {output.frag_txt})
        if [ "{params.read_method}" = "PE" ] && [ -f {input.sorted_bam} ]; then \
            bamPEFragmentSize -b {input.sorted_bam} --outRawFragmentLengths {output.frag_txt} || true ; \
        else \
            touch {output.frag_txt} ; \
        fi
        cat <<-END_VERSIONS > {output.versions}
        "calcFragsLength":
          bamPEFragmentSize: $(bamPEFragmentSize --version 2>&1 | head -1 || echo 'not-installed')
        END_VERSIONS
        """
