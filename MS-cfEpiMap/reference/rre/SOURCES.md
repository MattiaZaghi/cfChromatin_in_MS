# RRE Universe — Reference Data Sources

Fill in accession numbers and download URLs in `scripts/module1_build_rre.py`
before running `snakemake build_rre_universe`.

## Required datasets (all hg19)

| Cell type | Category | Source | Accession | Status |
|---|---|---|---|---|
| B_naive | bcell/immune | BLUEPRINT | TODO | ⬜ |
| B_memory | bcell/immune | BLUEPRINT | TODO | ⬜ |
| CD4_Th1 | immune | BLUEPRINT | TODO | ⬜ |
| CD4_Th17 | immune | BLUEPRINT | TODO | ⬜ |
| CD4_Treg | immune | BLUEPRINT | TODO | ⬜ |
| CD8_T | immune | BLUEPRINT | TODO | ⬜ |
| Monocyte_CD14 | immune | BLUEPRINT | TODO | ⬜ |
| Monocyte_CD16 | immune | BLUEPRINT | TODO | ⬜ |
| NK | immune | BLUEPRINT | TODO | ⬜ |
| Neutrophil | immune | BLUEPRINT | TODO | ⬜ |
| Megakaryocyte | other | BLUEPRINT | TODO | ⬜ |
| Oligodendrocyte | cns | ENCODE GRCh37 | TODO | ⬜ |
| OPC | cns | ENCODE GRCh37 | TODO | ⬜ |
| Neuron | cns | ENCODE GRCh37 | TODO | ⬜ |
| Astrocyte | cns | ENCODE GRCh37 | TODO | ⬜ |
| Microglia | cns | Roadmap/neuro-epigenomics | TODO | ⬜ |

## MS GWAS SNPs (for gwas_proximal_rre.bed)
- Source: IMSGC (International MS Genetics Consortium)
- File to place at: `reference/rre/ms_gwas_snps_hg19.bed`
- Coordinates must be hg19

## Housekeeping seed (already present)
- `reference/rre/housekeeping_seed.bed` — 14 loci from SNAPIE pipeline
- Used to initialise constitutive anchor set in Module 1
