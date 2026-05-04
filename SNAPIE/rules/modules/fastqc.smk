# ── FastQC on raw reads ────────────────────────────────────────────────────────
#
# FastQC names its output files after the *input fastq basename*, not the
# sample name, and produces one HTML + one ZIP per fastq file (so 2 pairs for
# PE data).  Because the fastq filenames are read dynamically from the
# samplesheet we cannot predict the exact output filenames at DAG-build time.
#
# Solution: direct FastQC output into a per-sample directory with -o, then
# touch a sentinel file (.done) when it finishes.  MultiQC reads the whole
# fastqc/ tree recursively, so it picks up all per-sample directories
# regardless of the individual filenames inside.

import csv, os

def _fastqc_reads(wildcards):
    """Return (r1, r2) paths for a given sample from the samplesheet."""
    ss = config.get('samplesheet', '')
    if ss and os.path.exists(ss):
        with open(ss) as fh:
            for row in csv.DictReader(fh):
                sid = row.get('sampleId') or row.get('sample') or row.get('sample_id')
                if sid == wildcards.sample:
                    r1 = row.get('read1') or row.get('read1_path') or row.get('Read1')
                    r2 = row.get('read2') or row.get('read2_path') or row.get('Read2')
                    if r1 and r2:
                        return r1, r2
    # Fallback: conventional filename pattern in reads_dir
    reads_dir = config.get('reads_dir', 'reads')
    r1 = os.path.join(reads_dir, f"{wildcards.sample}_R1_001.fastq.gz")
    r2 = os.path.join(reads_dir, f"{wildcards.sample}_R2_001.fastq.gz")
    return r1, r2


rule fastqc:
    """Run FastQC on both R1 and R2 for a single sample.

    Output is a sentinel file; the actual HTML/ZIP reports land in the same
    directory and are collected by MultiQC in the next step.
    """
    conda: "envs/qc.yaml"
    input:
        reads=lambda w: list(_fastqc_reads(w))
    output:
        done=touch(config['outputFolder'] + "/fastqc/{sample}/.done")
    params:
        outdir=config['outputFolder'] + "/fastqc/{sample}",
        threads=config.get('threads', 4)
    shell:
        """
        mkdir -p {params.outdir}
        fastqc --threads {params.threads} -o {params.outdir} {input.reads}
        """
