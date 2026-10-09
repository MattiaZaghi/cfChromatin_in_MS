
rule call_peaks:
    conda: "envs/peaks.yaml"
    input:
        bam=config['outputFolder'] + "/align/dedup/{sample}.dedup.unique.sorted.bam"
    output:
        narrowPeak=config['outputFolder'] + "/peaks/{sample}.narrowPeak",
        xls=config['outputFolder'] + "/peaks/{sample}_peaks.xls"
    params:
        genome_size="hs",
        dirname=config['outputFolder'] + "/peaks/"
    shell:
        """
        set -euo pipefail

        # MACS2 2.2.7.1 ships Cython extensions compiled against numpy 1.x. The
        # user site-packages (~/.local/lib/python3.9) precedes the conda env on
        # sys.path and holds numpy 2.x, which shadows the env's pinned 1.26.4 and
        # makes every `import MACS2` die with "numpy.dtype size changed".
        # Isolating user site keeps the conda env hermetic.
        export PYTHONNOUSERSITE=1

        echo "Running MACS2 for {wildcards.sample}"
        mkdir -p {params.dirname}

        # No `|| true`: a MACS2 crash must fail the job rather than leave an
        # empty output behind that Snakemake records as success.
        macs2 callpeak --SPMR -B -q 0.01 --keep-dup 1 \
            -g {params.genome_size} -f BAMPE \
            -t {input.bam} -n {wildcards.sample} \
            --bdg --outdir {params.dirname}

        # MACS2 names its output SAMPLE_peaks.narrowPeak; the declared output is
        # SAMPLE.narrowPeak, so it has to be renamed. The _peaks.xls it writes
        # already matches the declared name and needs no move.
        # (No braces in these comments: Snakemake str.format()s the whole block.)
        mv {params.dirname}{wildcards.sample}_peaks.narrowPeak {output.narrowPeak}

        echo "{wildcards.sample}: $(wc -l < {output.narrowPeak}) peaks"
        """
