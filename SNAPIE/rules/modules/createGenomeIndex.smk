rule createGenomeIndex_module:
    conda: "envs/common.yaml"
    params:
        ref_dir=config.get('referenceDir','reference')
    output:
        done=config['outputFolder'] + "/modules/createGenomeIndex.done",
        versions=config['outputFolder'] + "/modules/index_creation_mqc_versions.yml"
    run:
        import os, glob, subprocess

        refdir = params.ref_dir
        outdir = os.path.dirname(output.done)
        os.makedirs(outdir, exist_ok=True)

        # Find all .fa files in reference directory
        fasta_files = sorted(glob.glob(os.path.join(refdir, '*.fa')))
        if not fasta_files:
            with open(output.done, 'w') as f:
                f.write('no-fasta-found\n')
            with open(output.versions, 'w') as vf:
                vf.write('index_creation: no-fasta-found\n')
            print('No fasta files found in', refdir)
            return

        for fa in fasta_files:
            genome = os.path.splitext(os.path.basename(fa))[0]
            pac = os.path.join(refdir, genome + '.fa.pac')
            if not os.path.exists(pac):
                print('Creating index for', genome)
                subprocess.check_call(['bwa', 'index', os.path.join(refdir, genome + '.fa')])
            else:
                print('Index already exists for', genome)

        # Write versions (bwa)
        try:
            bwa_v = subprocess.check_output(['bwa', '2>&1'], shell=True, stderr=subprocess.STDOUT).decode()
            # extract Version line if present
            import re
            m = re.search(r'Version[: ]*([\d\.\-a-zA-Z]+)', bwa_v)
            bwa_ver = m.group(1) if m else bwa_v.splitlines()[0]
        except Exception:
            bwa_ver = 'bwa: not-detected'
        with open(output.versions, 'w') as vf:
            vf.write('"index_creation":\n')
            vf.write('  bwa: {}\n'.format(bwa_ver))

        with open(output.done, 'w') as f:
            f.write('done\n')
