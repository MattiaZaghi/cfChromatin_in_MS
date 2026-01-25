# Quick Start Guide: Extended cfChromatin Pipeline

## Overview

The cfChromatin pipeline has been extended with SNAPIE-inspired features while maintaining your original alignment and peak-calling workflow. All new features are **optional** and **modular** - your existing pipeline runs unchanged.

## What's New

9 new analysis modules have been added:

1. **FastQC** - Raw read quality assessment
2. **Filter Properly Paired** - QC filtering for paired-end data
3. **Library Complexity (Preseq)** - Estimate sequencing depth requirements
4. **Fragment Length Distribution** - Analyze fragment size patterns
5. **Peak Annotation** - Annotate peaks with genomic features
6. **Motif & GC Content** - Analyze nucleotide composition at fragment ends
7. **Enrichment Analysis** - Calculate histone mark enrichment ratios
8. **SNP Fingerprinting** - Verify sample identities and detect cross-contamination
9. **Quality Reports** - Aggregate all QC metrics into summaries

## Running the Pipeline

### Basic: Original pipeline only (no changes)
```bash
cd cfChromatin_in_MS/snakemake
snakemake -s snakefile.smk --configfile config_Cut_Tag.yaml -j 8
```

### With specific optional modules
```bash
# FastQC only
snakemake -s fastqc.smk --configfile config_Cut_Tag.yaml -j 8

# Library complexity analysis
snakemake -s lib_complex_preseq.smk --configfile config_Cut_Tag.yaml -j 8

# Fragment length analysis
snakemake -s frag_length_distribution.smk --configfile config_Cut_Tag.yaml -j 8

# Quality report generation
snakemake -s quality_reports.smk --configfile config_Cut_Tag.yaml -j 8
```

### Run all analyses (main + all optional)
```bash
# Main pipeline
snakemake -s snakefile.smk --configfile config_Cut_Tag.yaml -j 8

# Then run optional modules
for module in fastqc lib_complex_preseq frag_length_distribution quality_reports; do
    snakemake -s ${module}.smk --configfile config_Cut_Tag.yaml -j 8
done
```

## File Descriptions

### Core Pipeline (UNCHANGED)
- **snakefile.smk** - Your original pipeline with trimming, alignment, filtering, peaks
- **config_Cut_Tag.yaml** - Original configuration

### Extended Pipeline
- **snakefile_master.smk** - Optional: Combined pipeline with conditional module inclusion
- **config_extended.yaml** - Template showing new configuration parameters

### New Standalone Snakemake Files
Each can be run independently:
- **fastqc.smk** - QC on raw reads
- **filter_properly_paired.smk** - Filter paired-end reads
- **lib_complex_preseq.smk** - Library complexity estimation
- **frag_length_distribution.smk** - Fragment size distribution analysis
- **peaks_annotation.smk** - Annotate peaks with genomic features
- **motif_gc_content.smk** - Fragment end motif analysis
- **enrichment_analysis.smk** - Histone mark enrichment calculation
- **snp_fingerprint.smk** - Sample identity verification
- **quality_reports.smk** - Aggregate QC metrics

### Documentation
- **EXTENDED_PIPELINE_README.md** - Full documentation of all features
- **config_extended.yaml** - Configuration template with all new parameters
- **QUICKSTART.md** - This file

## Configuration Changes (if using extended features)

Add to your `config_Cut_Tag.yaml`:

```yaml
# For peak annotation (optional)
gtf_annotation: "ref_files/hg38/hg38_genes.gtf"

# For motif analysis (optional)  
genome_fasta: "ref_files/hg38/hg38.fa"
nmer_motif: 3

# For enrichment analysis (optional)
enrichment_marks:
  - H3K4me3
  - H3K27ac

# For SNP fingerprinting (optional)
snp_reference: "ref_files/snps/common_snps.vcf.gz"
```

## Common Use Cases

### Case 1: Run main pipeline + QC reports
```bash
snakemake -s snakefile.smk --configfile config_Cut_Tag.yaml -j 8 && \
snakemake -s fastqc.smk --configfile config_Cut_Tag.yaml -j 8 && \
snakemake -s quality_reports.smk --configfile config_Cut_Tag.yaml -j 8
```

### Case 2: Check library quality and complexity
```bash
snakemake -s fastqc.smk --configfile config_Cut_Tag.yaml -j 8 && \
snakemake -s lib_complex_preseq.smk --configfile config_Cut_Tag.yaml -j 8 && \
snakemake -s frag_length_distribution.smk --configfile config_Cut_Tag.yaml -j 8
```

### Case 3: Full analysis with annotations and enrichment
```bash
snakemake -s snakefile.smk --configfile config_Cut_Tag.yaml -j 8 && \
snakemake -s peaks_annotation.smk --configfile config_Cut_Tag.yaml -j 8 && \
snakemake -s enrichment_analysis.smk --configfile config_Cut_Tag.yaml -j 8
```

### Case 4: Verify sample identities
```bash
snakemake -s snakefile.smk --configfile config_Cut_Tag.yaml -j 8 && \
snakemake -s snp_fingerprint.smk --configfile config_Cut_Tag.yaml -j 8
```

## Output Locations

All new outputs go to subdirectories under your existing RUNID folder:

```
{RUNID}/
├── fastqc/                    # FastQC reports (if enabled)
├── preseq/                    # Library complexity data
├── fragments/                 # Fragment length analysis
├── peaks/annotated/           # Annotated peak files
├── motifs/                    # Fragment end motifs
├── enrichment/                # Enrichment results
├── snp_fingerprint/           # Fingerprinting results
└── qc_reports/                # Aggregated reports
```

## Troubleshooting

### "Environment not found" error
- Ensure conda environments exist: check paths in snakemake files
- Install missing tools: `conda env create -f environment.yml`

### Missing reference files
- Check gtf_annotation, genome_fasta, snp_reference paths in config
- Leave blank (`""`) to skip those analyses

### Memory errors
- Reduce THREADS in config
- Increase RAM in rule `resources: mem_mb = ...`

### Slow performance
- Lower THREADS if hitting system limits
- Run modules in parallel on separate HPC nodes
- Use `--dry-run` first to check file dependencies

## Key Differences from Original Pipeline

| Aspect | Original | Extended |
|--------|----------|----------|
| Trimming | Trimmomatic | Same |
| Alignment | Bowtie2/BWA | Same |
| Filtering | SAMtools | Same |
| Dedup | Picard | Same |
| Peak calling | MACS3 | Same |
| QC | Manual | FastQC + MultiQC |
| Library complexity | Not assessed | Preseq |
| Fragment analysis | Insert size only | Full distribution + motif |
| Peak annotation | Not available | Available (optional) |
| Enrichment | Not assessed | Available (optional) |

## Support & Questions

- **Full documentation**: See `EXTENDED_PIPELINE_README.md`
- **Original pipeline**: See original cfChromatin_in_MS README.md
- **SNAPIE info**: See `/cfs/klemming/home/m/matzag/SNAPIE/README.md`

## Citation

Please cite the original tools used:
- Trimmomatic
- Bowtie2/BWA
- SAMtools
- Picard
- MACS3
- DeepTools
- Preseq
- FastQC
- MultiQC

---

Last Updated: January 2026
