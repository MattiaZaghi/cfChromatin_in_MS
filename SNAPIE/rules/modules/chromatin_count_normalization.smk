rule chromatin_count_normalization_module:
    conda: "envs/common.yaml"
    params:
        frags_dir=config['outputFolder'] + "/frags",
        target_sites=config.get('target_sites',''),
        reference_sites=config.get('reference_sites','')
    input:
        target=config.get('target_sites',''),
        reference=config.get('reference_sites',''),
        frags=lambda wildcards: sorted(__import__('glob').glob(config['outputFolder'] + "/frags/*.bed"))
    output:
        done=config['outputFolder'] + "/modules/chromatin_count_normalization.done",
        versions=config['outputFolder'] + "/modules/chromatin_count_normalization_versions.yml"
    run:
        import os, glob, subprocess

        outdir = os.path.dirname(output.done)
        os.makedirs(outdir, exist_ok=True)

        target = input.target
        ref = input.reference
        frags = list(input.frags)

        # Guard: if target sites file is missing/empty, create marker and exit cleanly
        if not target or (not os.path.exists(target)) or os.path.getsize(target) == 0:
            open(os.path.join(outdir, 'EMPTY_TARGET_SITES'), 'w').close()
            with open(output.done, 'w') as f:
                f.write('EMPTY_TARGET_SITES\n')
            with open(output.versions, 'w') as vf:
                vf.write('chromatin_count_normalization: skipped_empty_target_sites\n')
            print('Skipped chromatin_count_normalization: target sites missing or empty')
            return

        # Create samplesheet
        sample_file = os.path.join(outdir, 'sample_name')
        with open(sample_file, 'w') as sf:
            sf.write('sample_name\n')
            for fpath in frags:
                bn = os.path.basename(fpath)
                name = os.path.splitext(bn)[0]
                sf.write(name + '\n')

        # Run the R script from a working directory under outdir
        cmd = ['Rscript', '/workspace/chromatin_count_norm_v2.R', '--samplesheet', sample_file, '--target-sites', target, '--frags-dir', '.', '--verbose']
        if ref:
            cmd = ['Rscript', '/workspace/chromatin_count_norm_v2.R', '--samplesheet', sample_file, '--target-sites', target, '--reference-sites', ref, '--frags-dir', '.', '--verbose']

        subprocess.check_call(cmd, cwd=outdir)

        # Record versions (Rscript version)
        try:
            rv = subprocess.check_output(['Rscript', '--version'], stderr=subprocess.STDOUT).decode().strip()
        except Exception:
            rv = 'Rscript: not-detected'
        with open(output.versions, 'w') as vf:
            vf.write('"chromatin_count_normalization":\n')
            vf.write('  Rscript: {}\n'.format(rv))

        with open(output.done, 'w') as f:
            f.write('done\n')
