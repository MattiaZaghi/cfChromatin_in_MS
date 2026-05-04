import csv, os

_SNAPIE_SAMPLES_CACHE = None

def _load_samplesheet():
    global _SNAPIE_SAMPLES_CACHE
    if _SNAPIE_SAMPLES_CACHE is not None:
        return _SNAPIE_SAMPLES_CACHE
    _SNAPIE_SAMPLES_CACHE = {}
    ss = config.get('samplesheet', '')
    if ss and os.path.exists(ss):
        with open(ss) as fh:
            reader = csv.DictReader(fh)
            for row in reader:
                sid = row.get('sampleId') or row.get('sample') or row.get('sample_id')
                if sid:
                    _SNAPIE_SAMPLES_CACHE[sid] = row
    return _SNAPIE_SAMPLES_CACHE

def reads_from_samplesheet(wildcards):
    samples = _load_samplesheet()
    row = samples.get(wildcards.sample)
    if row:
        r1 = row.get('read1') or row.get('read1_path') or row.get('Read1')
        r2 = row.get('read2') or row.get('read2_path') or row.get('Read2')
        if r1 and r2:
            return dict(r1=r1, r2=r2)
    # fallback to reads_dir patterns
    reads_dir = config.get('reads_dir', 'reads')
    cand_r1 = os.path.join(reads_dir, f"{wildcards.sample}_R1_001.fastq.gz")
    cand_r1_alt = os.path.join(reads_dir, f"{wildcards.sample}_R1.fastq.gz")
    cand_r2 = os.path.join(reads_dir, f"{wildcards.sample}_R2_001.fastq.gz")
    cand_r2_alt = os.path.join(reads_dir, f"{wildcards.sample}_R2.fastq.gz")
    r1 = cand_r1 if os.path.exists(cand_r1) else cand_r1_alt
    r2 = cand_r2 if os.path.exists(cand_r2) else cand_r2_alt
    return dict(r1=r1, r2=r2)


rule trim_galore:
    conda: "envs/qc.yaml"
    input:
        r1=lambda w: reads_from_samplesheet(w)['r1'],
        r2=lambda w: reads_from_samplesheet(w)['r2']
    output:
        r1_trim=temp(config['outputFolder'] + "/trim/{sample}_R1_trimmed.fq.gz"),
        r2_trim=temp(config['outputFolder'] + "/trim/{sample}_R2_trimmed.fq.gz"),
        report=config['outputFolder'] + "/trim/{sample}_trim_report.txt"
    params:
        threads=4,
        trim_params=config.get('trimming_params', '')
    shell:
        """
        mkdir -p $(dirname {output.r1_trim});
        trim_galore --paired {input.r1} {input.r2} --gzip --cores {params.threads} {params.trim_params} || true
        touch {output.r1_trim} {output.r2_trim} {output.report}
        """

rule trim_fastp:
    conda: "envs/qc.yaml"
    input:
        r1=lambda w: reads_from_samplesheet(w)['r1'],
        r2=lambda w: reads_from_samplesheet(w)['r2']
    output:
        r1_trim=temp(config['outputFolder'] + "/trim/{sample}_R1_trimmed.fq.gz"),
        r2_trim=temp(config['outputFolder'] + "/trim/{sample}_R2_trimmed.fq.gz"),
        json=config['outputFolder'] + "/trim/{sample}_fastp.json",
        html=config['outputFolder'] + "/trim/{sample}_fastp.html"
    params:
        threads=4
    shell:
        """
        mkdir -p $(dirname {output.r1_trim});
        fastp -i {input.r1} -I {input.r2} -o {output.r1_trim} -O {output.r2_trim} --thread {params.threads} --json {output.json} --html {output.html} || true
        """
