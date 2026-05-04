rule quality_filter:
    conda: "envs/common.yaml"
    input:
        bam=config['outputFolder'] + "/align/unique/{sample}.unique.sorted.bam"
    output:
        filtered=temp(config['outputFolder'] + "/align/filtered/{sample}.filtered.unique.sorted.bam"),
        versions=config['outputFolder'] + "/align/{sample}/samtools_QualityFilter_mqc_versions.yml"
    params:
        pe_flags=config.get('filter_samtools_pe_params', '-f 3 -F 3844 -q 30'),
        se_flags=config.get('filter_samtools_se_params', '-F 3844 -q 30'),
        read_method=config.get('read_method', 'PE')
    threads: 4
    shell:
        """
        mkdir -p $(dirname {output.filtered}) $(dirname {output.versions})
        if [ "{params.read_method}" = "PE" ]; then \
            FLAGS="{params.pe_flags}"; \
        else \
            FLAGS="{params.se_flags}"; \
        fi
        echo "Running samtools quality filter for sample {wildcards.sample} in {params.read_method} mode"
        samtools view -bh $FLAGS --threads {threads} {input.bam} > {output.filtered} || true

        cat <<-END_VERSIONS > {output.versions}
        "quality_filter":
          samtools: $(samtools --version | sed 's/^.*samtools //; s/Using.*$//')
        END_VERSIONS
        """
