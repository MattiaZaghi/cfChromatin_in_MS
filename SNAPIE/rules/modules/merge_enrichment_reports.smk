rule merge_enrichment_reports:
    conda: "envs/reports.yaml"
    input:
        script=config.get('pathMergeReportEnrichment', 'auxiliar_programs/merge_enrichment_reports.py')
    params:
        out_dir=config.get('outputFolder', '') + "/reports/multiqc"
    output:
        merged=config['outputFolder'] + "/reports/multiqc/merged_enrichment.csv"
    threads: 1
    shell:
        """
        mkdir -p {params.out_dir}
        if [ -f "{input.script}" ]; then
            python {input.script} || true
        fi
        touch {output.merged}
        """
