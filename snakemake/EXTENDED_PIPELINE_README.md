# Master Snakemake Pipeline for cfChromatin Analysis with SNAPIE-inspired Features
# This file documents how to use the extended pipeline with additional QC and analysis features

## Quick Start

The cfChromatin_in_MS pipeline now includes additional optional analysis modules inspired by SNAPIE.
The core pipeline (trimming, alignment, filtering, peak calling) remains unchanged.

### Running the main pipeline only:
```
snakemake -s snakefile.smk --configfile config_Cut_Tag.yaml -j 8
```

### Running with additional quality control:

#### 1. FastQC (raw read quality assessment)
```
snakemake -s fastqc.smk --configfile config_Cut_Tag.yaml -j 8
```

#### 2. Filter properly paired reads
```
snakemake -s filter_properly_paired.smk --configfile config_Cut_Tag.yaml -j 8
```

#### 3. Library complexity analysis using preseq
```
snakemake -s lib_complex_preseq.smk --configfile config_Cut_Tag.yaml -j 8
```

#### 4. Fragment length distribution analysis
```
snakemake -s frag_length_distribution.smk --configfile config_Cut_Tag.yaml -j 8
```

#### 5. Peak annotation with genomic features
```
snakemake -s peaks_annotation.smk --configfile config_Cut_Tag.yaml -j 8
```
*Requires: GTF annotation file configured in config_Cut_Tag.yaml as 'gtf_annotation'*

#### 6. Motif and GC content analysis at fragment ends
```
snakemake -s motif_gc_content.smk --configfile config_Cut_Tag.yaml -j 8
```
*Requires: Genome FASTA file configured in config_Cut_Tag.yaml as 'genome_fasta'*

#### 7. Enrichment analysis for histone marks
```
snakemake -s enrichment_analysis.smk --configfile config_Cut_Tag.yaml -j 8
```
*Requires: On-target and off-target BED files for enrichment marks*

#### 8. SNP-based sample fingerprinting
```
snakemake -s snp_fingerprint.smk --configfile config_Cut_Tag.yaml -j 8
```
*Optional: SNP reference file for sample identity verification*

#### 9. Comprehensive Quality Reports
```
snakemake -s quality_reports.smk --configfile config_Cut_Tag.yaml -j 8
```
*Generates CSV summary and MultiQC HTML report*

## Configuration Requirements

Add the following to your config_Cut_Tag.yaml if using advanced features:

```yaml
# For peak annotation
gtf_annotation: "path/to/reference.gtf"

# For motif and GC content analysis
genome_fasta: "path/to/genome.fa"
nmer_motif: 3  # k-mer length for motif analysis (default 3 bp)

# For enrichment analysis
enrichment_marks: ['H3K4me3', 'H3K27ac', 'MeDIP']  # Histone marks to analyze
enrichment_on_target_dir: "ref_files/enrichment_states"

# For SNP fingerprinting
snp_reference: "path/to/snps.vcf.gz"
smash_script: "auxiliar_programs/SMaSH.py"
```

## Pipeline Features Summary

### Original cfChromatin Features (UNCHANGED):
- Trimmomatic adapter/quality trimming
- Bowtie2/BWA alignment
- SAMtools filtering and sorting
- Picard deduplication
- Chromosome filtering
- DeepTools coverage/BigWig generation
- MACS3 peak calling
- Insert size distribution

### New Features (Optional):

| Module | Purpose | Input | Output |
|--------|---------|-------|--------|
| fastqc.smk | Quality control on raw reads | FASTQ files | HTML QC reports |
| filter_properly_paired.smk | Keep only properly paired reads | BAM | Filtered BAM |
| lib_complex_preseq.smk | Estimate library complexity | BAM | Learning curve |
| frag_length_distribution.smk | Analyze fragment size distribution | BAM | Fragment length file & plot |
| peaks_annotation.smk | Annotate peaks with genomic features | narrowPeak + GTF | Annotated peaks |
| motif_gc_content.smk | Analyze nucleotide composition at fragment ends | BAM + FASTA | Motif BED file |
| enrichment_analysis.smk | Calculate histone mark enrichment | BED + on/off-target BED | Enrichment ratios |
| snp_fingerprint.smk | Sample identity verification using SNPs | BAM + SNP refs | Fingerprint report |
| quality_reports.smk | Aggregate all QC metrics | All QC files | CSV summary + HTML report |

## Running Combined Analyses

To run multiple analyses in sequence:

```bash
# Run main pipeline + basic QC
snakemake -s snakefile.smk --configfile config_Cut_Tag.yaml -j 8 && \
snakemake -s fastqc.smk --configfile config_Cut_Tag.yaml -j 8 && \
snakemake -s lib_complex_preseq.smk --configfile config_Cut_Tag.yaml -j 8 && \
snakemake -s frag_length_distribution.smk --configfile config_Cut_Tag.yaml -j 8 && \
snakemake -s quality_reports.smk --configfile config_Cut_Tag.yaml -j 8
```

## Key Differences from SNAPIE

### What's the Same:
- Overall analysis philosophy: raw reads → aligned BAM → peaks → annotation
- Support for paired-end Cut&Tag and ChIP-seq data
- Quality control and enrichment analysis
- Fragment-level analysis

### What's Different:
- **Alignment**: Uses Bowtie2/BWA (proven tools) instead of custom SNAPIE aligner
- **Format**: Snakemake instead of Nextflow (Snakemake runs locally/HPC, Nextflow designed for clouds)
- **Deduplication**: Uses Picard (well-tested) instead of custom SAM/BAM handling
- **Configuration**: Uses JSON sample sheets instead of CSV/directory-based inputs
- **No external containers**: Uses conda environments already on your system
- **Modular design**: Run only what you need; optional modules don't affect main pipeline

## Customization Tips

### Adding conda environments:
Most rules reference conda environments at standard paths. Update these paths if your conda envs are located elsewhere:
```yaml
# In snakemake files, check paths like:
conda: "/home/mattia/miniconda3/envs/fastqc.yml"
```

### Adjusting resource requirements:
Each rule has `resources` and `threads` parameters. Modify for your HPC system:
```python
resources:
    mem_mb = 64000  # Adjust RAM in MB
threads: config['THREADS']  # Adjust based on your system
```

### Disabling specific rules:
Comment out targets in the `rule all:` input to skip modules you don't need.

## Troubleshooting

### Error: "Environment not found"
Install missing conda environment or update path to existing one

### Error: "File not found" for reference files
Ensure GTF, FASTA, or BED files are provided in config and paths are correct

### Memory errors
Increase `mem_mb` in resource definitions or reduce `THREADS`

### Missing output files
Check the snakefile targets in `rule all:` input list - add paths to target files you need

## Output Directory Structure

```
{RUNID}/
├── trimmed/trimmomatic/          # Trimmed FASTQ files
├── mapped/bwa/ or /bowtie2/       # Aligned SAM files
├── sorted/samtools/               # Sorted BAM files
├── dedup/picard/                  # Deduplicated BAM files
├── filter/samtools/               # Final filtered BAM files
├── coverage/deeptools/            # BigWig files for visualization
├── peaks/macs3/                   # Called peaks (narrowPeak, summits)
├── peaks/annotated/               # Annotated peaks (if annotation enabled)
├── bed/                           # BED format conversions
├── fastqc/                        # FastQC reports (if enabled)
├── preseq/                        # Library complexity data (if enabled)
├── fragments/                     # Fragment length analysis (if enabled)
├── motifs/                        # Motif/GC analysis (if enabled)
├── enrichment/                    # Enrichment calculations (if enabled)
├── snp_fingerprint/               # SNP fingerprinting results (if enabled)
└── qc_reports/                    # Aggregated QC summary and MultiQC
```

## Citation & References

**cfChromatin_in_MS**: Original pipeline for circulating cfChromatin analysis
**SNAPIE**: Streamlined Nextflow Analysis Pipeline for Immunoprecipitation-Based Epigenomic Profiling
**Tools Used**:
- Bowtie2: Langmead & Salzberg (2012)
- BWA: Li & Durbin (2009)
- SAMtools: Danecek et al. (2021)
- Picard: Broad Institute
- DeepTools: Ramírez et al. (2016)
- MACS3: Feng et al. (2011) / Zang et al. (2021)
- FastQC: Babraham Institute
- Preseq: Daley & Smith (2013)
- MultiQC: Ewels et al. (2016)

---
**Last Updated**: January 2026
**Contact**: See original cfChromatin_in_MS repository
