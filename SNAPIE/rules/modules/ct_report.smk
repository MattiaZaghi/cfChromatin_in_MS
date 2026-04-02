rule ct_report:
    conda: "envs/reports.yaml"
    input:
        fragle=config['outputFolder'] + "/reports/fragle/Fragle.txt"
    output:
        report=config['outputFolder'] + "/reports/multiqc/ct_fragle_mqc.csv"
    params:
        outdir=config['outputFolder'] + "/reports/multiqc/",
        report_script=config.get('pathReportCT', 'auxiliar_programs/report_ct_fragle.py')
    shell:
        """
        mkdir -p {params.outdir}
        if [ -f "{params.report_script}" ]; then
            cp {input.fragle} {params.outdir}/Fragle.txt
            cd {params.outdir}
            python {params.report_script} || true
            if ls *_mqc.csv 1>/dev/null 2>&1; then
                cp $(ls *_mqc.csv | head -1) {output.report} 2>/dev/null || true
            fi
        fi
        touch {output.report}
        """
