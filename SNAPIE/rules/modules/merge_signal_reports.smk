rule merge_signal_reports:
    conda: "envs/reports.yaml"
    input:
        script=config.get('pathMergeReportSignal', 'auxiliar_programs/merge_signal_reports.py')
    params:
        out_dir=config.get('outputFolder', '') + "/reports/multiqc"
    output:
        merged=config['outputFolder'] + "/reports/multiqc/merged_signal.csv"
    threads: 1
    shell:
        """
        mkdir -p {params.out_dir}
        if [ -f "{input.script}" ]; then
            python {input.script} || true
        fi
        touch {output.merged}
        """
