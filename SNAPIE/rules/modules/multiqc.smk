rule multiqc_report:
    conda: "envs/qc.yaml"
    input:
        files_dir=config['outputFolder'] + "/fastqc/{sample}/"
    output:
        report=config['outputFolder'] + "/reports/multiqc/multiqc_report.html"
    params:
        mqc_config=config.get('multiqc_config_no_peak_annotation', '')
    shell:
        """
        mkdir -p $(dirname {output.report});
        if [ -n "{params.mqc_config}" ]; then \
            multiqc {input.files_dir} -c {params.mqc_config} -o $(dirname {output.report}) || true ; \
        else \
            multiqc {input.files_dir} -o $(dirname {output.report}) || true ; \
        fi
        """
