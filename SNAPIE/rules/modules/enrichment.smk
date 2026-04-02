rule enrichment_module:
    conda: "envs/common.yaml"
    params:
        peaks_dir=config['outputFolder'] + "/peaks",
        enrichment_states_ref=config.get('enrichment_states_ref','')
    output:
        done=config['outputFolder'] + "/modules/enrichment.done"
    run:
        import os, glob, subprocess

        peaks_dir = params.peaks_dir
        outdir = os.path.dirname(output.done)
        os.makedirs(outdir, exist_ok=True)

        # Find BAMs under align that match dedup unique sorted pattern
        bam_pattern = os.path.join(config['outputFolder'],'align','*.dedup.unique.sorted.bam')
        bams = sorted(glob.glob(bam_pattern))
        if not bams:
            with open(output.done, 'w') as f:
                f.write('no-bams-found\n')
            print('No BAMs found for enrichment step')
            return

        for bam in bams:
            sample = os.path.basename(bam).replace('.dedup.unique.sorted.bam','')
            sample_dir = os.path.join(peaks_dir, sample)
            os.makedirs(sample_dir, exist_ok=True)

            # locate enrichment script(s)
            # Expect script path to be provided in config under 'chEnrichmentScript' or as a param
            script = config.get('chEnrichmentScript','')
            if not script:
                # fallback: look for any script named '*enrichment*' in repo modules
                possible = glob.glob('**/*enrichment*', recursive=True)
                script = possible[0] if possible else ''

            csv_out = os.path.join(sample_dir, sample + '_enrichment_states.csv')
            # If script missing, create empty CSV and continue
            if not script:
                open(csv_out, 'w').close()
                continue

            cmd = [script, bam, params.enrichment_states_ref, sample]
            subprocess.check_call(cmd, cwd=sample_dir)

        with open(output.done, 'w') as f:
            f.write('done\n')
