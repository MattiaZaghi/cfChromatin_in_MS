# Generate comprehensive quality control reports
# Aggregates various QC metrics into HTML reports and CSV summaries

configfile: "config_Cut_Tag.yaml"

import json
FILES = json.load(open(config['SAMPLES_JSON']))

SAMPLES = sorted(FILES.keys())

# List all samples by sample_name, sample_type, and assay
MARK_SAMPLES = []
for sample in SAMPLES:
    for sample_type in FILES[sample].keys():
        for assay in FILES[sample][sample_type].keys():
            MARK_SAMPLES.append(sample + "_" + sample_type+ "_" + assay)

CUT_TAG = config["c_t"]
CHIP = config["chip"]
CUT_TAGS  = [sample for sample in MARK_SAMPLES if CUT_TAG in sample]
CHIPS = [sample for sample in MARK_SAMPLES if CHIP in sample]
ALL_SAMPLES =  CHIPS + CUT_TAGS

RUNID = config["RUN_ID"]

QC_SUMMARY = "{myrun}/qc_reports/quality_summary.csv"
MULTIQC_HTML = "{myrun}/qc_reports/multiqc_report.html"

rule all:
    input: QC_SUMMARY, MULTIQC_HTML

rule aggregate_qc_metrics:
    """
    Aggregate QC metrics from all samples into a summary CSV
    Includes: mapping stats, duplication rates, fragment length, enrichment, peak counts
    """
    input:
        flagstats = expand("{myrun}/filter/samtools/{sample}.flagstat", 
                          sample=ALL_SAMPLES, myrun=RUNID),
        metrics = expand("{myrun}/dedup/picard/{sample}.bam_metrics.txt", 
                        sample=ALL_SAMPLES, myrun=RUNID),
        insert_metrics = expand("{myrun}/filter/samtools/{sample}_insert.txt", 
                               sample=ALL_SAMPLES, myrun=RUNID),
        peaks = expand("{myrun}/peaks/macs3/{sample}_peaks.narrowPeak", 
                      sample=ALL_SAMPLES, myrun=RUNID)
    output:
        summary = "{myrun}/qc_reports/quality_summary.csv"
    params:
        dir = "{myrun}/qc_reports"
    resources:
        mem_mb = 32000
    threads: 1
    conda:
        "/home/mattia/miniconda3/envs/python.yml"
    script:
        """
        import pandas as pd
        import re
        
        params_dir = "{params.dir}"
        
        # Create output directory
        import os
        os.makedirs(params_dir, exist_ok=True)
        
        metrics_data = []
        
        # Process each sample
        sample_list = {str(sorted(set([f.split('/')[-1].split('.')[0] for f in {input.flagstats}])))!r}
        
        for sample in sample_list:
            sample_metrics = {{'Sample': sample}}
            
            # Read flagstat file
            flagstat_file = f"{RUNID}/filter/samtools/{{sample}}.flagstat"
            if os.path.exists(flagstat_file):
                with open(flagstat_file) as f:
                    lines = f.readlines()
                    # Parse flagstat format
                    for line in lines:
                        if 'mapped' in line and '%' in line:
                            match = re.search(r'(\d+) \+ \d+ mapped', line)
                            if match:
                                sample_metrics['mapped_reads'] = int(match.group(1))
                        elif 'properly paired' in line:
                            match = re.search(r'(\d+) \+ \d+ properly paired', line)
                            if match:
                                sample_metrics['properly_paired'] = int(match.group(1))
            
            # Count peaks
            peak_file = f"{RUNID}/peaks/macs3/{{sample}}_peaks.narrowPeak"
            if os.path.exists(peak_file):
                with open(peak_file) as f:
                    peak_count = sum(1 for line in f if line.strip())
                    sample_metrics['peak_count'] = peak_count
            
            metrics_data.append(sample_metrics)
        
        # Create DataFrame and save
        df = pd.DataFrame(metrics_data)
        df.to_csv("{output.summary}", index=False)
        print(f"QC summary written to {output.summary}")
        """

rule multiqc:
    """
    Run MultiQC to generate comprehensive HTML quality report
    Aggregates FastQC, alignment, and peak-calling statistics
    """
    input:
        fastqc_dirs = expand("{myrun}/fastqc/{sample}", 
                            sample=ALL_SAMPLES, myrun=RUNID),
        flagstats = expand("{myrun}/filter/samtools/{sample}.flagstat", 
                          sample=ALL_SAMPLES, myrun=RUNID),
        insert_pdfs = expand("{myrun}/filter/samtools/{sample}_insert.pdf", 
                            sample=ALL_SAMPLES, myrun=RUNID)
    output:
        html = "{myrun}/qc_reports/multiqc_report.html"
    params:
        dir = "{myrun}/qc_reports",
        input_dir = "{myrun}"
    resources:
        mem_mb = 32000
    threads: config.get('THREADS', 4)
    conda:
        "/home/mattia/miniconda3/envs/multiqc.yml"
    shell:
        """
        mkdir -p {params.dir}
        
        # Run MultiQC on the entire run directory
        # It will automatically detect and process quality control files
        multiqc {params.input_dir} -o {params.dir} -n multiqc_report.html -f --quiet
        """
