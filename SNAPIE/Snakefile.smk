configfile: "config.yaml"

import csv
import os
from pathlib import Path

# ─── Sample list ──────────────────────────────────────────────────────────────
SAMPLES = []
if config.get('samplesheet'):
    samplesheet_path = Path(config['samplesheet'])
    if samplesheet_path.exists():
        with open(samplesheet_path, newline='') as fh:
            reader = csv.DictReader(fh)
            for row in reader:
                sid = row.get('sampleId') or row.get('sample') or row.get('sample_id')
                if sid:
                    SAMPLES.append(sid)
    else:
        print(f"Warning: samplesheet {samplesheet_path} not found; SAMPLES empty")
else:
    print("No samplesheet provided in config.yaml; define 'samplesheet'")

# Expose sample list to rules that need it (e.g. fragle_ct_estimation)
config['_samples_'] = SAMPLES

# ─── Include rule modules ──────────────────────────────────────────────────────
# Preprocessing
include: "rules/modules/fastqc.smk"
include: "rules/modules/trim.smk"

# Alignment
include: "rules/modules/align.smk"

# BAM processing chain
include: "rules/modules/sort_bam.smk"             # sort_bam, sort_readname_bam
include: "rules/modules/filter_properly_paired.smk"
include: "rules/modules/unique_sam.smk"
include: "rules/modules/quality_filter.smk"
include: "rules/modules/dedup.smk"
include: "rules/modules/dac_exclusion.smk"
include: "rules/modules/index_and_quality.smk"     # index_sam
include: "rules/modules/lib_complex_preseq.smk"
include: "rules/modules/createStatsSamtoolsfiltered.smk"
include: "rules/modules/snp_smash_fingerprint.smk"

# Fragments processing
include: "rules/modules/end_motif_gc.smk"
include: "rules/modules/calcFragsLength.smk"
include: "rules/modules/bam_to_bed.smk"
include: "rules/modules/unique_frags.smk"
include: "rules/modules/frags_report.smk"

# Fragle CT estimation (optional)
include: "rules/modules/filter_bam_fragle.smk"
include: "rules/modules/fragle_ct_estimation.smk"
include: "rules/modules/ct_report.smk"

# Signal processing
include: "rules/modules/bam_to_bedgraph.smk"
include: "rules/modules/bedgraph_to_bigwig.smk"
include: "rules/modules/call_peaks.smk"
include: "rules/modules/peaks_annotations.smk"
include: "rules/modules/peaks_report.smk"
include: "rules/modules/enrichment.smk"
include: "rules/modules/enrichmentReport.smk"
include: "rules/modules/merge_enrichment_reports.smk"
include: "rules/modules/chromatin_count_normalization.smk"
include: "rules/modules/igv_reports.smk"

# Meta plots
include: "rules/modules/meta_plot_housekeeping.smk"

# Reports
include: "rules/modules/quality_report_lite.smk"
include: "rules/modules/signal_report_lite.smk"
include: "rules/modules/merge_signal_reports.smk"
include: "rules/modules/multiqc.smk"

# Reference/download rules (used when ref files need to be fetched)
include: "rules/modules/download.smk"
include: "rules/modules/createGenomeIndex.smk"
include: "rules/modules/fetch_chrom_sizes.smk"

# ─── Rule ordering ────────────────────────────────────────────────────────────
ruleorder: trim_fastp > trim_galore

# ─── Helper: final processed BAM path (dac-filtered or dedup) ─────────────────
def get_final_bam(wildcards):
    base = config['outputFolder']
    if config.get('exclude_dac_regions', False):
        return f"{base}/align/dac/{wildcards.sample}.dac_filtered.dedup.unique.sorted.bam"
    return f"{base}/align/dedup/{wildcards.sample}.dedup.unique.sorted.bam"

def get_final_bai(wildcards):
    base = config['outputFolder']
    if config.get('exclude_dac_regions', False):
        return f"{base}/align/dac/{wildcards.sample}.dac_filtered.dedup.unique.sorted.bam.bai"
    return f"{base}/align/dedup/{wildcards.sample}.dedup.unique.sorted.bam.bai"

# ─── Targets ──────────────────────────────────────────────────────────────────
OUT = config['outputFolder']

DEDUP = expand(OUT + "/align/dedup/{sample}.dedup.unique.sorted.bam", sample=SAMPLES)

if config.get('exclude_dac_regions', False):
    FINAL_BAM = expand(OUT + "/align/dac/{sample}.dac_filtered.dedup.unique.sorted.bam", sample=SAMPLES)
    FINAL_BAI = expand(OUT + "/align/dac/{sample}.dac_filtered.dedup.unique.sorted.bam.bai", sample=SAMPLES)
else:
    FINAL_BAM = DEDUP
    FINAL_BAI = expand(OUT + "/align/dedup/{sample}.dedup.unique.sorted.bam.bai", sample=SAMPLES)

BED          = expand(OUT + "/frags/{sample}.bed",               sample=SAMPLES)
BEDGRAPH     = expand(OUT + "/bedgraph/{sample}.bedgraph",        sample=SAMPLES)
BIGWIG       = expand(OUT + "/bigwig/{sample}.bw",               sample=SAMPLES)
DEEPTOOLS_BW = expand(OUT + "/bigwig/deeptools/{sample}.bw",     sample=SAMPLES)
PEAKS        = expand(OUT + "/peaks/{sample}.narrowPeak",        sample=SAMPLES)
FRAG_SIZES   = expand(OUT + "/frags/{sample}/{sample}.fragment_sizes.txt", sample=SAMPLES)
MOTIFS       = expand(OUT + "/motifs/{sample}/{sample}_4NMER_bp_motif.bed", sample=SAMPLES)
LIB_COMPLEX  = expand(OUT + "/align/{sample}/{sample}.lc_extrap.txt", sample=SAMPLES)
UNIQUE_FRAGS = expand(OUT + "/frags/{sample}/{sample}_unique_frags.csv", sample=SAMPLES)

META_PLOTS    = expand(OUT + "/reports/meta_plots/{sample}_housekeeping_meta.pdf", sample=SAMPLES)

TARGETS = []
TARGETS.extend(FINAL_BAM)
TARGETS.extend(FINAL_BAI)
TARGETS.extend(BED)
TARGETS.extend(BEDGRAPH)
TARGETS.extend(BIGWIG)
TARGETS.extend(DEEPTOOLS_BW)
TARGETS.extend(PEAKS)
TARGETS.extend(FRAG_SIZES)
TARGETS.extend(UNIQUE_FRAGS)
TARGETS.extend(LIB_COMPLEX)
TARGETS.extend(META_PLOTS)

if config.get('read_method', 'PE') == 'PE':
    TARGETS.extend(MOTIFS)

if config.get('fragle_ct_estimation', False) and SAMPLES:
    TARGETS.append(OUT + "/reports/fragle/Fragle.txt")
    TARGETS.append(OUT + "/reports/multiqc/ct_fragle_mqc.csv")

rule all:
    input: TARGETS
