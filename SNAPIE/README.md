# SNAPIE — Snakemake Edition
### Streamlined Analysis Pipeline for Immunoprecipitation-Based Epigenomic Profiling of Circulating Chromatin

![Pipeline Workflow](workflow-SNAP.png)
![DAG](dag.png)

> **SNAPIE** was originally developed at the BacaLab (Dana-Farber Cancer Institute) as a Nextflow pipeline. This repository contains a fully converted **Snakemake** implementation that preserves all functionality of the original pipeline while benefiting from Snakemake's file-based dependency graph, native conda environment management, and cluster integration.

---

## Table of Contents

1. [Overview](#overview)
2. [Directory Structure](#directory-structure)
3. [Installation & Dependencies](#installation--dependencies)
4. [Configuration (`config.yaml`)](#configuration-configyaml)
5. [Input Data](#input-data)
6. [Running the Pipeline](#running-the-pipeline)
7. [Pipeline Steps in Detail](#pipeline-steps-in-detail)
8. [Output Structure](#output-structure)
9. [Conda Environments](#conda-environments)
10. [Advanced Options](#advanced-options)
11. [Troubleshooting](#troubleshooting)

---

## Overview

SNAPIE processes cell-free chromatin data from raw FASTQ reads (or pre-aligned BAMs) through the following major stages:

```
FASTQ reads
    │
    ▼
[Quality Control]   FastQC
    │
    ▼
[Trimming]          fastp / Trim-Galore
    │
    ▼
[Alignment]         BWA-MEM → SAM → BAM
    │
    ▼
[BAM Processing]    sort → filter properly-paired → remove non-unique
                    → quality filter → deduplicate (Picard)
                    → DAC/blacklist exclusion → index
    │
    ├──► Library complexity (preseq)
    │
    ▼
[Fragment Processing]
    │
    ├──► Name-sorted BAM → BED fragments → unique fragment count
    ├──► Fragment length distribution (bamPEFragmentSize)
    ├──► End-motif / GC content analysis (4-mer)
    └──► Fragle CT estimation (optional)
    │
    ▼
[Signal Processing]
    │
    ├──► BAM → bedgraph → bigWig (bedtools + bedGraphToBigWig)
    ├──► deepTools bamCoverage (RPKM-normalised bigWig)
    ├──► Peak calling (MACS2)
    ├──► Peak annotations (genomic annotation R script)
    ├──► Enrichment scoring (on/off-target ratio)
    ├──► Chromatin count normalisation
    └──► Meta-coverage plots around housekeeping gene TSSs (R/ggplot2 PDF)
    │
    ▼
[Reports]
    ├──► Fragments report (MultiQC-compatible CSV)
    ├──► Peaks report (MultiQC-compatible CSV)
    ├──► CT report from Fragle (MultiQC-compatible CSV)
    ├──► IGV consolidated report (HTML)
    └──► MultiQC
```

---

## Directory Structure

```
SNAPIE/
├── Snakefile.smk                  # Main workflow entry point
├── config.yaml                    # All runtime parameters
├── dag.png                        # Auto-generated pipeline DAG
├── samplesheet.csv                # Example samplesheet
├── auxiliar_programs/             # Helper scripts (Python / R / Bash)
│   ├── enrichment.sh
│   ├── report_frags.py
│   ├── report_peaks.py
│   ├── report_ct_fragle.py
│   ├── genomic_annotation.R
│   ├── meta_plot_housekeeping.R   # Meta-coverage plots at housekeeping TSSs
│   └── ...
├── config/                        # Additional config files / samplesheets
├── envs/                          # (optional) top-level conda envs
├── ref_files/                     # Reference data (genomes, BED files, etc.)
│   ├── enrichment_states/         # On/off-target BED files per mark
│   ├── fragle_sites/              # Fragle site BED files per mark
│   ├── pileup_report/             # Housekeeping regions BED
│   └── multiqc/                   # MultiQC config / header files
├── rules/
│   └── modules/                   # One .smk file per processing step
│       ├── envs/                  # Conda environment YAML files
│       │   ├── align.yaml         # BWA, samtools, bedtools, picard, R
│       │   ├── common.yaml        # deeptools, preseq, python utilities
│       │   ├── peaks.yaml         # MACS2, R
│       │   ├── qc.yaml            # FastQC, fastp, Trim-Galore
│       │   ├── reports.yaml       # Python, pandas
│       │   └── r_plots.yaml      # R 4.2+, ggplot2, cowplot, rtracklayer, GenomicRanges
│       ├── align.smk
│       ├── bam_to_bed.smk
│       ├── bam_to_bedgraph.smk
│       ├── bedgraph_to_bigwig.smk
│       ├── calcFragsLength.smk
│       ├── call_peaks.smk
│       ├── chromatin_count_normalization.smk
│       ├── createGenomeIndex.smk
│       ├── createStatsSamtoolsfiltered.smk
│       ├── ct_report.smk
│       ├── dac_exclusion.smk
│       ├── dedup.smk
│       ├── download.smk
│       ├── end_motif_gc.smk
│       ├── enrichment.smk
│       ├── enrichmentReport.smk
│       ├── fastqc.smk
│       ├── fetch_chrom_sizes.smk
│       ├── filter_bam_fragle.smk
│       ├── filter_properly_paired.smk
│       ├── fragle_ct_estimation.smk
│       ├── frags_report.smk
│       ├── igv_reports.smk
│       ├── index_and_quality.smk
│       ├── lib_complex_preseq.smk
│       ├── merge_enrichment_reports.smk
│       ├── merge_signal_reports.smk
│       ├── meta_plot_housekeeping.smk
│       ├── multiqc.smk
│       ├── peaks_annotations.smk
│       ├── peaks_report.smk
│       ├── quality_filter.smk
│       ├── quality_report_lite.smk
│       ├── signal_report_lite.smk
│       ├── snp_smash_fingerprint.smk
│       ├── sort_bam.smk           # contains sort_bam + sort_readname_bam
│       ├── trim.smk
│       ├── unique_frags.smk
│       └── unique_sam.smk
└── test_datasets/                 # Small test data for validation
```

---

## Installation & Dependencies

### 1. Install Snakemake

```bash
conda create -n snakemake -c bioconda -c conda-forge snakemake graphviz
conda activate snakemake
```

### 2. Clone / set up the repository

```bash
git clone <repo_url>
cd SNAPIE
```

### 3. Conda environments

All per-rule conda environments are defined in `rules/modules/envs/`. Snakemake will automatically create them on first use when you pass `--use-conda`. The five environments are:

| Environment | Key tools |
|---|---|
| `qc` | FastQC 0.12, fastp 1.3, Trim-Galore 0.6 |
| `align` | BWA 0.7.19, samtools 1.23, bedtools 2.31, Picard 3.4, R 4.5, bedGraphToBigWig |
| `common` | deeptools 3.5, preseq 3.2, Python 3.12, bamPEFragmentSize |
| `peaks` | MACS2 2.2.7, R 4.5 |
| `reports` | Python 3, pandas |
| `r_plots` | R 4.2+, ggplot2, cowplot, reshape2, rtracklayer, GenomicRanges |

> The Fragle CT estimation step requires a separate container/environment (`fragle`) that wraps the Fragle tool. Set up the container image separately and ensure it is accessible to the compute environment.

---

## Configuration (`config.yaml`)

Edit `config.yaml` before running. Key parameters:

```yaml
# ── Reference ────────────────────────────────────────────────────────────
genome: "hg19"           # Reference genome name (hg19 or hg38)
ref_dir: "/path/to/genomes/hg19/"   # Directory where genome FASTA will be found/stored
genome_fa: "/path/to/genomes/hg19/fa/bwa/hg19"  # BWA index prefix
chrom_sizes: "/path/to/genomes/hg19/hg19.chrom.sizes"

# ── Input ─────────────────────────────────────────────────────────────────
samplesheet: "config/samplesheet.csv"   # Path to CSV samplesheet

# ── Output ────────────────────────────────────────────────────────────────
outputFolder: "/path/to/output"

# ── Processing flags ──────────────────────────────────────────────────────
read_method: "PE"                # PE (paired-end) or SE (single-end)
trim_method: "FASTP"             # FASTP or TRIM_GALORE
exclude_dac_regions: true        # Remove ENCODE DAC blacklist regions
deduped_bam: false               # Set true to skip deduplication (input already deduped)
fragle_ct_estimation: true       # Run Fragle CT burden estimation
report_peak_genomic_annotation: false  # Run R genomic annotation on peaks

# ── Tool parameters ───────────────────────────────────────────────────────
filter_samtools_pe_params: "-f 3 -F 3844 -q 30"
filter_samtools_se_params: "-F 3844 -q 30"
trimming_params_fastp: "--adapter_sequence=AGATCGGAAGAGCACACGTCTGAACTCCAGTCA --adapter_sequence_r2=AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT"
bwa_params: "-I 250,200,2000,10"

# ── deepTools bigWig ──────────────────────────────────────────────────────
binsize: "10"
norm_method: "RPKM"
smooth_length: "300"

# ── Chromatin count normalisation ─────────────────────────────────────────
chromatin_count_mode: "single"   # single or batch
chromatin_count_target_sites: null
chromatin_count_reference: null

# ── Script paths (relative to SNAPIE/) ───────────────────────────────────
pathEnrichmentScript: "auxiliar_programs/enrichment.sh"
pathReportFrags: "auxiliar_programs/report_frags.py"
pathReportPeaks: "auxiliar_programs/report_peaks.py"
pathReportCT: "auxiliar_programs/report_ct_fragle.py"
pathRGenomicAnnotation: "auxiliar_programs/genomic_annotation.R"

# ── Meta-coverage plots ───────────────────────────────────────────────────
housekeeping_bed: "housekeeping_genes_comprehensive.bed"  # BED with housekeeping gene coordinates
pathMetaPlotScript: "auxiliar_programs/meta_plot_housekeeping.R"
meta_plot_window: 10000    # Total window in bp centred on TSS (default ±5 kb)
meta_plot_binsize: 25      # Coverage bin size in bp
meta_plot_color: "steelblue"  # Heatmap fill colour
```

---

## Input Data

### Samplesheet format

The samplesheet is a CSV file with the following columns:

**FASTQ input:**
```
sampleId,enrichment_mark,read1,read2,control
sample1,H3K4me3,/data/sample1_R1.fq.gz,/data/sample1_R2.fq.gz,ctrl
ctrl,H3K4me3,/data/ctrl_R1.fq.gz,/data/ctrl_R2.fq.gz,
```

- `read2` can be left blank for single-end samples.
- `control` should be the `sampleId` of the control sample (typically IgG or input). Leave blank if no control is available — MACS2 will run without a control file.

**BAM input** (set `deduped_bam: true` if BAMs are already processed):
```
sampleId,enrichment_mark,bam,control,read_method
sample1,H3K4me3,/data/sample1.bam,ctrl,PE
ctrl,H3K4me3,/data/ctrl.bam,,PE
```

### Supported enrichment marks

Pre-configured reference files are included for:
- `H3K4me3`
- `H3K27ac`
- `MeDIP`
- `no_enrichment_mark` (skip enrichment scoring)

Custom marks: place `on.target.filt.bed` and `off.target.filt.bed` in `ref_files/enrichment_states/<mark_name>/` and set `enrichment_mark` in the samplesheet accordingly.

---

## Running the Pipeline

### Dry-run (check that the DAG builds correctly)

```bash
snakemake -s Snakefile.smk --dry-run --use-conda
```

### Local execution

```bash
snakemake -s Snakefile.smk \
    --use-conda \
    --cores 16 \
    --rerun-incomplete
```

### SLURM cluster (recommended for large cohorts)

```bash
snakemake -s Snakefile.smk \
    --use-conda \
    --executor slurm \
    --default-resources slurm_partition=short mem_mb=8000 \
    --jobs 200
```

### Regenerate the DAG image

```bash
snakemake -s Snakefile.smk --dag | dot -Tpng -o dag.png
```

---

## Pipeline Steps in Detail

### 1. Quality Control (`fastqc.smk`)
- **Tool:** FastQC 0.12
- **Input:** Raw FASTQ reads (from samplesheet)
- **Output:** `{outputFolder}/fastqc/{sample}/`
- Per-sample quality metrics used later in MultiQC

### 2. Read Trimming (`trim.smk`)
- **Tool:** fastp (default) or Trim-Galore (set `trim_method: TRIM_GALORE`)
- **Input:** Raw FASTQ
- **Output:** `{outputFolder}/trim/{sample}_R[12]_trimmed.fq.gz`
- fastp removes Illumina TruSeq adapters and low-quality bases
- Rule ordering ensures `trim_fastp` takes precedence over `trim_galore` for samples with identical output names

### 3. Alignment (`align.smk`)
- **Tool:** BWA-MEM + samtools view
- **Input:** Trimmed FASTQ, BWA index
- **Output:** `{outputFolder}/align/raw/{sample}.bam`
- Alignments are streamed directly to BAM (no intermediate SAM on disk)

### 4. BAM Processing (`sort_bam.smk`, `filter_properly_paired.smk`, `unique_sam.smk`, `quality_filter.smk`, `dedup.smk`, `dac_exclusion.smk`)

The full BAM processing chain is:

```
align/raw/{sample}.bam
    │ sort_bam (samtools sort)
    ▼
align/sorted/{sample}.sorted.bam
    │ filter_properly_paired
    │   PE: samtools view -b -f 2  (keep properly paired)
    │   SE: copy as-is
    ▼
align/pp/{sample}.pp.sorted.bam
    │ unique_sam (samtools view -q 1)
    ▼
align/unique/{sample}.unique.sorted.bam
    │ quality_filter (samtools view with -f/-F/-q flags from config)
    ▼
align/filtered/{sample}.filtered.unique.sorted.bam
    │ dedup_bam (Picard AddOrReplaceReadGroups + MarkDuplicates)
    ▼
align/dedup/{sample}.dedup.unique.sorted.bam
    │ [if exclude_dac_regions: true]
    │ dac_exclusion (bedtools intersect -v against {genome}.DAC.bed)
    ▼
align/dac/{sample}.dac_filtered.dedup.unique.sorted.bam  ← FINAL BAM
    │
    └──► index_dac_bam / index_dedup_bam (samtools index)
```

> All downstream steps (bedgraph, peaks, fragments, etc.) use the **final BAM**: `align/dac/` when `exclude_dac_regions: true`, otherwise `align/dedup/`.

**Library complexity** (`lib_complex_preseq.smk`): preseq lc_extrap is run on the `pp` BAM (before deduplication) to estimate library complexity. Output: `align/{sample}/{sample}.lc_extrap.txt`.

**Alignment statistics** (`createStatsSamtoolsfiltered.smk`): samtools stats/idxstats/flagstat run on the final BAM.

### 5. Fragment Processing

#### Name-sorted BAM (`sort_bam.smk :: sort_readname_bam`)
- **Tool:** samtools sort -n
- **Input:** Final BAM
- **Output:** `align/namesorted/{sample}.n_sorted.bam`
- Required for bedtools bamtobed -bedpe (paired-end BED)

#### BED fragments (`bam_to_bed.smk`)
- **Tool:** bedtools bamtobed
- **PE mode:** `-bedpe` → extract fragment start (read1 start) to end (read2 end), keeping only read-pairs on the same chromosome
- **SE mode:** single-end BED intervals
- **Output:** `frags/{sample}.bed`

#### Unique fragment count (`unique_frags.smk`)
- Counts lines in the BED file per sample
- **Output:** `frags/{sample}/{sample}_unique_frags.csv`

#### Fragment length distribution (`calcFragsLength.smk`)
- **Tool:** bamPEFragmentSize (deeptools) — PE only
- **Output:** `frags/{sample}/{sample}.fragment_sizes.txt`

#### End-motif & GC content (`end_motif_gc.smk`)
- **Tools:** bedtools bamtobed, bedtools nuc, bedtools getfasta
- **PE only** — skipped for SE samples (empty output file created)
- Generates BEDPE → filters for GC content → extracts 4-mer sequences from both read ends
- **Output:** `motifs/{sample}/{sample}_4NMER_bp_motif.bed`
- `nmer` (default 4) can be changed in `config.yaml`

#### Fragle CT estimation (`filter_bam_fragle.smk`, `fragle_ct_estimation.smk`, `ct_report.smk`)
- **Activated by:** `fragle_ct_estimation: true` in `config.yaml`
- `filter_bam_fragle`: subsets each BAM to Fragle-specific genomic sites (`ref_files/fragle_sites/{enrichment_mark}/sites.bed`)
- `fragle_ct_estimation`: runs the Fragle tool across all filtered BAMs jointly → `reports/fragle/Fragle.txt`
- `ct_report`: converts Fragle output to MultiQC-compatible CSV → `reports/multiqc/ct_fragle_mqc.csv`

### 6. Signal Processing

#### bedgraph & bigWig (`bam_to_bedgraph.smk`, `bedgraph_to_bigwig.smk`)
- **Tools:** bedtools genomecov, bedGraphToBigWig
- PE: `-pc` flag (pair-concordant coverage)
- **Outputs:** `bedgraph/{sample}.bedgraph`, `bigwig/{sample}.bw`

#### deeptools bigWig (`bedgraph_to_bigwig.smk :: coverage_deeptools`)
- **Tool:** bamCoverage (deeptools)
- RPKM-normalised, binsize 10 bp, smoothing 300 bp (all configurable)
- **Output:** `bigwig/deeptools/{sample}.bw`

#### Peak calling (`call_peaks.smk`)
- **Tool:** MACS2 callpeak
- PE format: `BAMPE`; SE format: `BAM`
- Uses control BAM when specified in samplesheet; otherwise runs without control
- **Output:** `peaks/{sample}.narrowPeak`, `peaks/{sample}_peaks.xls`

#### Peak annotations (`peaks_annotations.smk`)
- **Tool:** R genomic annotation script (`auxiliar_programs/genomic_annotation.R`)
- Activated by `report_peak_genomic_annotation: true`
- **Output:** annotation plots in `reports/multiqc/`

#### Enrichment scoring (`enrichment.smk`, `enrichmentReport.smk`, `merge_enrichment_reports.smk`)
- Calculates on-target vs off-target read enrichment using the reference BED files in `ref_files/enrichment_states/{enrichment_mark}/`
- Per-sample CSV → per-sample HTML report → merged report
- **Outputs:** `peaks/{sample}/{sample}_enrichment_states.csv`, `reports/multiqc/merged_enrichment.csv`

#### Chromatin count normalisation (`chromatin_count_normalization.smk`)
- Two modes: `single` (per-sample) or `batch` (all samples jointly)
- Requires `chromatin_count_target_sites` to be set; `chromatin_count_reference` is optional
- Based on the [chromatin-frags-normalization](https://github.com/chhetribsurya/chromatin-frags-normalization) library
- **Output:** `{outputFolder}/modules/chromatin_count_normalization.done`

#### IGV report (`igv_reports.smk`)
- Consolidates per-sample IGV HTML reports into a linked index page
- **Output:** `reports/multiqc/igv_housekeeping_genes_mqc.html`

#### Meta-coverage plots at housekeeping TSSs (`meta_plot_housekeeping.smk`)
- **Tool:** R (`rtracklayer`, `GenomicRanges`, `ggplot2`, `cowplot`)
- **Input:** deeptools-normalised bigWig (`bigwig/deeptools/{sample}.bw`) + housekeeping gene BED (`housekeeping_genes_comprehensive.bed`)
- **What it does:**
  1. Imports the bigWig as an `RleList` coverage track via `rtracklayer::import()`
  2. Converts each housekeeping gene's start coordinate to a single-bp TSS anchor
  3. Calls `CollectMeta()` (from `MetaPlot.R` logic) to build a coverage matrix of binned signal in a configurable window (default ±5 kb, 25 bp bins) centred on each TSS, with strand-aware reversal for minus-strand genes
  4. Produces two complementary panels per sample:
     - **Average meta-profile** (line plot): mean coverage across all housekeeping TSSs
     - **Heatmap**: per-gene coverage ranked by total signal
  5. Combines both panels with `cowplot::plot_grid()` into a single two-panel PDF
- **Output:** `reports/meta_plots/{sample}_housekeeping_meta.pdf`
- **Config options:**

| Key | Default | Description |
|---|---|---|
| `housekeeping_bed` | `housekeeping_genes_comprehensive.bed` | BED file with housekeeping gene coordinates |
| `pathMetaPlotScript` | `auxiliar_programs/meta_plot_housekeeping.R` | Path to the meta-plot R script |
| `meta_plot_window` | `10000` | Total window size in bp centred on TSS |
| `meta_plot_binsize` | `25` | Coverage aggregation bin size (bp) |
| `meta_plot_color` | `steelblue` | Heatmap fill colour |

The included `housekeeping_genes_comprehensive.bed` contains 15 constitutively-expressed genes (ACTB, GAPDH, TUBB, RPL13, HPRT1, TBP, B2M, RPLP0, GUSB, PPIA, YWHAZ, TFRC, UBC, NONO, PGK1) whose nucleosome-depleted TSS regions serve as positive controls for open-chromatin enrichment.

---

## Output Structure

```
{outputFolder}/
├── trim/                          # Trimmed FASTQ files
├── fastqc/                        # FastQC HTML + ZIP per sample
├── align/
│   ├── raw/                       # Initial BAM from BWA
│   ├── sorted/                    # Coordinate-sorted BAM
│   ├── pp/                        # Properly-paired filtered BAM
│   ├── unique/                    # Multi-mapper removed BAM
│   ├── filtered/                  # Quality-filtered BAM
│   ├── dedup/                     # Deduplicated BAM + metrics
│   ├── dac/                       # DAC blacklist-excluded BAM (final BAM)
│   ├── namesorted/                # Name-sorted BAM (for BED conversion)
│   ├── fragle/                    # Fragle-filtered BAMs (if enabled)
│   └── {sample}/                  # Per-sample stats and preseq output
├── frags/
│   ├── {sample}.bed               # Fragment BED file
│   └── {sample}/                  # Fragment length txt + unique frag count
├── motifs/
│   └── {sample}/                  # 4-mer end-motif BED files
├── bedgraph/                      # Per-sample bedgraph coverage files
├── bigwig/
│   ├── {sample}.bw                # bedGraphToBigWig output
│   └── deeptools/{sample}.bw      # RPKM-normalised deeptools output
├── peaks/
│   └── {sample}/                  # narrowPeak, XLS, MACS2 bedgraph
├── chromatin_count_normalization/ # CCN output matrices
├── reports/
│   ├── fragle/Fragle.txt          # CT burden estimates
│   ├── meta_plots/                # Per-sample meta-coverage PDFs
│   │   └── {sample}_housekeeping_meta.pdf
│   └── multiqc/
│       ├── frags_mqc.csv          # Fragment count summary
│       ├── peaks_mqc.csv          # Peak count summary
│       ├── ct_fragle_mqc.csv      # CT burden (MultiQC table)
│       ├── merged_enrichment.csv  # Merged enrichment report
│       └── igv_housekeeping_genes_mqc.html
└── modules/                       # Sentinel .done files for module-level rules
```

---

## Conda Environments

All environments are pre-pinned to exact package versions for reproducibility. They live in `rules/modules/envs/` and are activated automatically by Snakemake per rule. To pre-build all environments:

```bash
snakemake -s Snakefile.smk --use-conda --conda-create-envs-only
```

---

## Advanced Options

### Skip deduplication (pre-processed BAMs)

Set in `config.yaml`:
```yaml
deduped_bam: true
```
The pipeline will use the BAM supplied in the samplesheet directly, bypassing sorting, filtering, and deduplication.

### Single-end reads

```yaml
read_method: "SE"
```
Affects: trimming, alignment flags, BAM filtering, bedtools genomecov, bam_to_bed, calcFragsLength (disabled for SE), end_motif_gc (disabled for SE), sort_readname_bam (disabled for SE).

### Custom genome

1. Add the genome to `ref_files/genome/genome_paths.csv` following the existing format.
2. Place (or let the pipeline download) the FASTA at `{ref_dir}/{genome}.fa`.
3. Set `genome: "custom_name"` in `config.yaml`.

### Blacklist exclusion

The pipeline downloads the ENCODE DAC blacklist BED (`{genome}.DAC.bed`) into `ref_dir` via `rules/modules/download.smk`. To use a local file, pre-place it at `{ref_dir}/{genome}.DAC.bed` before running.

### Changing the 4-mer size for end-motif analysis

```yaml
nmer: 6   # default is 4
```

### deeptools bigWig parameters

```yaml
binsize: "10"       # bin size in bp
norm_method: "RPKM" # RPKM, CPM, BPM, RPGC, or None
smooth_length: "300"
```

### Chromatin count normalisation modes

```yaml
chromatin_count_mode: "batch"                             # joint across all samples
chromatin_count_target_sites: "ref_files/tss_peaks.bed"  # required
chromatin_count_reference: "ref_files/reference.bed"     # optional
```

---

## Troubleshooting

| Issue | Solution |
|---|---|
| `MissingInputException` for genome FASTA | Run `download_genome` rule first, or pre-place `{ref_dir}/{genome}.fa` |
| `MissingInputException` for DAC BED | Run `download_dac` rule first, or pre-place `{ref_dir}/{genome}.DAC.bed` |
| Fragle step fails silently | Ensure the Fragle container is available; the rule falls back gracefully to an empty output file |
| Empty BED fragments for PE sample | Check that name-sorted BAM is non-empty; verify that `read_method: PE` is set |
| `dot: graph is too large` warning | Expected for large cohorts; DAG image is auto-scaled. Use `--dag` + `dot -Tsvg` for vector output |
| `400 jobs have missing provenance/metadata` | Normal warning for pre-existing output files that pre-date Snakemake tracking; use `--rerun-incomplete` |

To force re-run of specific rules:
```bash
snakemake -s Snakefile.smk --use-conda --cores 16 \
    --forcerun dedup_bam dac_exclusion
```

---

## Citation

If you use SNAPIE in your research, please cite the original pipeline:

> BacaLab — Dana-Farber Cancer Institute. SNAPIE: Streamlined Analysis Pipeline for Immunoprecipitation-Based Epigenomic Profiling of Circulating Chromatin. https://github.com/prc992/SNAPIE
