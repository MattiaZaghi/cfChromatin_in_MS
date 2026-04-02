rule fastqc:
    conda: "envs/qc.yaml"
    input:
        reads1=lambda w: __reads_for_fastqc(w)['reads1'],
        reads2=lambda w: __reads_for_fastqc(w)['reads2']
    output:
        html=config['outputFolder'] + "/fastqc/{sample}/{sample}_fastqc.html",
        zip=config['outputFolder'] + "/fastqc/{sample}/{sample}_fastqc.zip"
    params:
        threads=2
    shell:
        """
        mkdir -p $(dirname {output.html});
        fastqc --threads {params.threads} {input.reads1} {input.reads2} || true
        """

def __reads_for_fastqc(wildcards):
    # Reuse samplesheet lookup from trim module if available, otherwise fallback
    import csv, os
    ss = config.get('samplesheet', '')
    if ss and os.path.exists(ss):
        with open(ss) as fh:
            reader = csv.DictReader(fh)
            for row in reader:
                sid = row.get('sampleId') or row.get('sample') or row.get('sample_id')
                if sid == wildcards.sample:
                    r1 = row.get('read1') or row.get('read1_path') or row.get('Read1')
                    r2 = row.get('read2') or row.get('read2_path') or row.get('Read2')
                    if r1 and r2:
                        return dict(reads1=r1, reads2=r2)
    reads_dir = config.get('reads_dir', 'reads')
    r1 = os.path.join(reads_dir, f"{wildcards.sample}_R1_001.fastq.gz")
    r2 = os.path.join(reads_dir, f"{wildcards.sample}_R2_001.fastq.gz")
    # prefer existing files
    if not os.path.exists(r1):
        r1 = os.path.join(reads_dir, f"{wildcards.sample}_R1.fastq.gz")
    if not os.path.exists(r2):
        r2 = os.path.join(reads_dir, f"{wildcards.sample}_R2.fastq.gz")
    return dict(reads1=r1, reads2=r2)
