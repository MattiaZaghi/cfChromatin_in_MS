rule peaks_report:
    conda: "envs/reports.yaml"
    input:
        script=config.get('pathReportPeaks', 'auxiliar_programs/report_peaks.py')
    params:
        peaks_dir=config.get('outputFolder', '') + "/peaks",
        out_dir=config.get('outputFolder', '') + "/reports/multiqc"
    output:
        report=config.get('outputFolder', '') + "/reports/multiqc/peaks_mqc.csv"
    threads: 1
    shell:
        """
        mkdir -p {params.out_dir}
        if [ -f {input.script} ] && [ -d {params.peaks_dir} ]; then
            pushd {params.peaks_dir} >/dev/null
            python {input.script} || true
            popd >/dev/null
            if [ -f {params.peaks_dir}/peaks_mqc.csv ]; then
                cp {params.peaks_dir}/peaks_mqc.csv {params.out_dir}/peaks_mqc.csv
            else
                echo "peaks_mqc.csv not produced" >&2
            fi
        else
            echo "Script {input.script} or peaks dir {params.peaks_dir} missing" >&2
        fi
        """
