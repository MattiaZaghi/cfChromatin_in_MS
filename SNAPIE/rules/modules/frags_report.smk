rule frags_report:
    conda: "envs/reports.yaml"
    input:
        script=config.get('pathReportFrags', 'auxiliar_programs/report_frags.py'),
        unique_frags=expand(
            config['outputFolder'] + "/frags/{sample}/{sample}_unique_frags.csv",
            sample=config.get('_samples_', [])
        )
    params:
        frags_dir=config.get('outputFolder', '') + "/frags",
        out_dir=config.get('outputFolder', '') + "/reports/multiqc"
    output:
        report=config['outputFolder'] + "/reports/multiqc/frags_mqc.csv"
    threads: 1
    shell:
        """
        mkdir -p {params.out_dir}
        if [ -f "{input.script}" ] && [ -d "{params.frags_dir}" ]; then
            pushd {params.frags_dir} >/dev/null
            python {input.script} || true
            popd >/dev/null
            if [ -f "{params.frags_dir}/frags_mqc.csv" ]; then
                cp "{params.frags_dir}/frags_mqc.csv" {params.out_dir}/frags_mqc.csv || true
            fi
        fi
        touch {output.report}
        """
