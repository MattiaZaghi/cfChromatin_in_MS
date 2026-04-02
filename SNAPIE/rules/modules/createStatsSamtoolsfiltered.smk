rule createStatsSamtoolsfiltered_module:
    conda: "envs/common.yaml"
    params:
        align_dir=config['outputFolder'] + "/align"
    output:
        done=config['outputFolder'] + "/modules/createStatsSamtoolsfiltered.done"
    run:
        import os, glob, subprocess

        align_dir = params.align_dir
        outdir = os.path.dirname(output.done)
        os.makedirs(outdir, exist_ok=True)

        bam_pattern = os.path.join(align_dir, '*.dedup.unique.sorted.bam')
        bams = sorted(glob.glob(bam_pattern))
        if not bams:
            with open(output.done, 'w') as f:
                f.write('no-bams-found\n')
            print('No BAM files found matching', bam_pattern)
            return

        for bam in bams:
            sample = os.path.basename(bam).replace('.dedup.unique.sorted.bam','')
            sample_dir = os.path.join(align_dir, sample)
            os.makedirs(sample_dir, exist_ok=True)

            stats_file = os.path.join(sample_dir, sample + '.AfterFilter.stats')
            idxstats_file = os.path.join(sample_dir, sample + '.AfterFilter.idxstats')
            flagstat_file = os.path.join(sample_dir, sample + '.AfterFilter.flagstat')

            subprocess.check_call(['samtools', 'stats', bam], stdout=open(stats_file,'w'))
            subprocess.check_call(['samtools', 'idxstats', bam], stdout=open(idxstats_file,'w'))
            subprocess.check_call(['samtools', 'flagstat', bam], stdout=open(flagstat_file,'w'))

            # versions file per sample
            versions_path = os.path.join(sample_dir, 'samtools_stats_filtered_mqc_versions.yml')
            try:
                st = subprocess.check_output(['samtools', '--version'], stderr=subprocess.STDOUT).decode().strip()
                # Simplify grab of version token
                st = st.splitlines()[0]
            except Exception:
                st = 'samtools: not-detected'
            with open(versions_path, 'w') as vf:
                vf.write('"createStatsSamtoolsfiltered":\n')
                vf.write('  samtools: {}\n'.format(st))

        with open(output.done, 'w') as f:
            f.write('done\n')
