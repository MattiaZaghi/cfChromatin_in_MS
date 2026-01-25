# Quick Reference Card - Extended cfChromatin Pipeline

## 📝 Files at a Glance

**Location**: `/cfs/klemming/home/m/matzag/cfChromatin_in_MS/snakemake/`

### 🎯 Start with These
- `INDEX.md` - Navigation guide
- `QUICKSTART.md` - Quick start (5 min)
- `DELIVERY_SUMMARY.md` - Project overview

### 📊 Analysis Modules (9)
1. `fastqc.smk` - Raw read QC
2. `lib_complex_preseq.smk` - Library complexity
3. `frag_length_distribution.smk` - Fragment sizes
4. `filter_properly_paired.smk` - Paired-end filtering
5. `peaks_annotation.smk` - Peak annotation
6. `motif_gc_content.smk` - Motif/GC analysis
7. `enrichment_analysis.smk` - Enrichment ratios
8. `snp_fingerprint.smk` - Sample verification
9. `quality_reports.smk` - QC aggregation

### 📚 Documentation (5)
- `INDEX.md` - File navigation
- `QUICKSTART.md` - Quick examples
- `EXTENDED_PIPELINE_README.md` - Full docs
- `FILE_OVERVIEW.md` - File organization
- `IMPLEMENTATION_SUMMARY.md` - Technical details

### ⚙️ Configuration
- `config_extended.yaml` - Configuration template
- `snakefile_master.smk` - Combined pipeline

---

## 🚀 Usage

### Original Pipeline (No Changes)
```bash
snakemake -s snakefile.smk --configfile config_Cut_Tag.yaml -j 8
```

### Run One Module
```bash
snakemake -s fastqc.smk --configfile config_Cut_Tag.yaml -j 8
```

### Run Multiple Modules
```bash
snakemake -s snakefile.smk --configfile config_Cut_Tag.yaml -j 8 && \
snakemake -s fastqc.smk --configfile config_Cut_Tag.yaml -j 8 && \
snakemake -s lib_complex_preseq.smk --configfile config_Cut_Tag.yaml -j 8 && \
snakemake -s quality_reports.smk --configfile config_Cut_Tag.yaml -j 8
```

### Test First (Dry Run)
```bash
snakemake -s fastqc.smk --configfile config_Cut_Tag.yaml --dry-run
```

---

## ⚙️ Configuration Template

Add to `config_Cut_Tag.yaml`:

```yaml
fastqc_enabled: true
preseq_enabled: true
frag_length_enabled: true
peaks_annotation_enabled: true
gtf_annotation: "path/to/hg38_genes.gtf"
motif_gc_enabled: true
genome_fasta: "path/to/hg38.fa"
nmer_motif: 3
enrichment_enabled: true
enrichment_marks: [H3K4me3, H3K27ac]
snp_fingerprint_enabled: true
snp_reference: "path/to/snps.vcf.gz"
generate_qc_reports: true
```

---

## 📊 What Each Module Does

| Module | Input | Output | Use Case |
|--------|-------|--------|----------|
| fastqc | FASTQ | HTML report | Check raw read quality |
| lib_complex_preseq | BAM | Learning curve | Estimate sequencing depth needs |
| frag_length_distribution | BAM | Fragment data | Analyze fragment sizes |
| filter_properly_paired | BAM | Filtered BAM | Quality filtering |
| peaks_annotation | narrowPeak | Annotated peaks | Add genomic feature info |
| motif_gc_content | BAM | Motif BED | Analyze cutting patterns |
| enrichment_analysis | BED | Enrichment ratios | Quantify enrichment |
| snp_fingerprint | BAM | Fingerprints | Verify sample identity |
| quality_reports | All QC | CSV + HTML | Summary reports |

---

## 📂 Output Locations

All go under your `{RUNID}/` folder:

```
{RUNID}/
├── fastqc/              (FastQC reports)
├── preseq/              (Library complexity)
├── fragments/           (Fragment analysis)
├── peaks/annotated/     (Annotated peaks)
├── motifs/              (Motif/GC data)
├── enrichment/          (Enrichment results)
├── snp_fingerprint/     (Fingerprints)
└── qc_reports/          (QC summaries)
```

---

## 💡 Common Scenarios

### "I just want to run my normal pipeline"
```bash
snakemake -s snakefile.smk --configfile config_Cut_Tag.yaml -j 8
```
✅ No changes needed, works as before

### "I want FastQC + basic QC"
```bash
snakemake -s snakefile.smk --configfile config_Cut_Tag.yaml -j 8
snakemake -s fastqc.smk --configfile config_Cut_Tag.yaml -j 8
snakemake -s quality_reports.smk --configfile config_Cut_Tag.yaml -j 8
```

### "I want library quality assessment"
```bash
snakemake -s lib_complex_preseq.smk --configfile config_Cut_Tag.yaml -j 8
snakemake -s frag_length_distribution.smk --configfile config_Cut_Tag.yaml -j 8
```

### "I want peak annotation + enrichment"
```bash
# Ensure config has:
# gtf_annotation: "path/to/genes.gtf"
# enrichment_enabled: true

snakemake -s peaks_annotation.smk --configfile config_Cut_Tag.yaml -j 8
snakemake -s enrichment_analysis.smk --configfile config_Cut_Tag.yaml -j 8
```

### "I want everything"
```bash
for mod in snakefile fastqc lib_complex_preseq frag_length_distribution \
           peaks_annotation enrichment_analysis quality_reports; do
    snakemake -s ${mod}.smk --configfile config_Cut_Tag.yaml -j 8
done
```

---

## ❓ Troubleshooting

| Error | Solution |
|-------|----------|
| "Environment not found" | Check conda paths in .smk file |
| "File not found" | Ensure GTF/FASTA paths in config |
| "Missing reference file" | Leave config blank to skip that analysis |
| "Out of memory" | Reduce THREADS or increase RAM |
| "Module not found" | Check file is in snakemake/ directory |

---

## 📖 Documentation Index

| Document | Time | Purpose |
|----------|------|---------|
| QUICKSTART.md | 5 min | Get running fast |
| INDEX.md | 5 min | Navigate files |
| EXTENDED_PIPELINE_README.md | 20 min | Deep dive |
| config_extended.yaml | 5 min | See all options |
| FILE_OVERVIEW.md | 10 min | Understand structure |
| IMPLEMENTATION_SUMMARY.md | 15 min | Learn design |

**Start with**: QUICKSTART.md or INDEX.md

---

## ✅ Checklist Before Running

- [ ] Read QUICKSTART.md or INDEX.md
- [ ] Check conda environments exist
- [ ] Update config_Cut_Tag.yaml with reference files (if using advanced features)
- [ ] Run --dry-run first: `snakemake -s {module}.smk --dry-run`
- [ ] Check that outputs look reasonable
- [ ] Run for real with `-j` threads
- [ ] Monitor logs if errors occur

---

## 🔑 Key Facts

✅ **Backward Compatible** - Original pipeline works unchanged  
✅ **Optional** - Use only what you need  
✅ **Independent** - Each module runs alone  
✅ **Configurable** - Single config file controls all  
✅ **Documented** - 5 comprehensive guides  
✅ **Ready to Use** - Production quality  

---

## 🎯 Next Steps

1. **Read**: QUICKSTART.md (5 minutes)
2. **Configure**: Add references to config (5 minutes)
3. **Test**: Run `--dry-run` (1 minute)
4. **Execute**: Run snakemake -j 8 (varies)

---

## 📞 Need Help?

| Question | See This |
|----------|----------|
| Quick reference | This card |
| How to run | QUICKSTART.md |
| File locations | INDEX.md or FILE_OVERVIEW.md |
| All options | EXTENDED_PIPELINE_README.md |
| Configuration | config_extended.yaml |
| Design details | IMPLEMENTATION_SUMMARY.md |

---

## Summary

**16 new files** = **9 modules** + **5 docs** + **1 config** + **1 master**

All fully documented, production-ready, 100% backward compatible.

**Status**: ✅ Ready to use  
**Compatibility**: ✅ 100% backward compatible  
**Documentation**: ✅ Complete  

---

**Version**: 1.0 | **Date**: January 2026 | **Status**: Production Ready
