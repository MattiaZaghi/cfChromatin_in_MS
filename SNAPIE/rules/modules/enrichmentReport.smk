rule enrichmentReport_module:
    conda: "envs/common.yaml"
    params:
        report_dir=config['outputFolder'] + "/reports/multiqc"
    output:
        done=config['outputFolder'] + "/modules/enrichmentReport.done"
    run:
        import os, glob, subprocess

        outdir = os.path.dirname(output.done)
        os.makedirs(outdir, exist_ok=True)

        # Expect enrichment CSVs under peaks/*/*.csv
        csv_pattern = os.path.join(config['outputFolder'],'peaks','*','*_enrichment_states.csv')
        csvs = sorted(glob.glob(csv_pattern))

        # locate report script
        script = config.get('chReportEnrichment','')
        if not script:
            possible = glob.glob('**/*report*enrichment*', recursive=True)
            script = possible[0] if possible else ''

        # For each csv, run the report script (if available), otherwise create empty report
        for csv in csvs:
            sample = os.path.basename(csv).replace('_enrichment_states.csv','')
            report_file = os.path.join(outdir, sample + '_report.csv')
            open(report_file, 'w').close()
            if script:
                cmd = ['python', script, '--mark', config.get('enrichment_mark',''), '--samplename', sample]
                subprocess.check_call(cmd, cwd=outdir)

        with open(output.done, 'w') as f:
            f.write('done\n')
