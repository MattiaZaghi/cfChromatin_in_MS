#!/usr/bin/env Rscript
# qc_flag_samples.R
#
# QC flagging for cfChIP-seq samples, coloured by disease group.
#
# Cut-offs come from config['qc_thresholds'] (passed in as a TSV) and are
# applied per enrichment mark: a sample is only ever compared against others
# sharing its mark, because on/off-target region sets — and therefore the score
# scale — differ between marks. One PDF page per mark.
#
# All shared machinery lives in qc_flag_common.R; this file only supplies the
# grouping used for colouring.
#
# Args:
#   1  enrichment_tsv   — all_enrichment_scores.tsv  (sample_id, mark, …, enrichment_score)
#   2  frag_counts_tsv  — fragment_counts.tsv        (sample_id, n_fragments)
#   3  output_tsv       — qc_summary.tsv
#   4  output_pdf       — qc_distributions.pdf
#   5  thresholds_tsv   — qc_thresholds.tsv          (mark, enrich_thresh, frag_thresh)

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 4)
  stop("Usage: qc_flag_samples.R <enrichment_tsv> <frag_counts_tsv> ",
       "<output_tsv> <output_pdf> [thresholds_tsv]")

enrich_tsv <- args[1]
frag_tsv   <- args[2]
out_tsv    <- args[3]
out_pdf    <- args[4]
thresh_tsv <- if (length(args) >= 5) args[5] else NA

# Locate qc_flag_common.R next to this script, whatever the working directory
.this <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])
source(file.path(dirname(.this), "qc_flag_common.R"))

# ── Disease group, from the sample ID ────────────────────────────────────────
parse_group <- function(id) {
  if      (grepl("MS-Rituximab-Progressive", id)) "MS-Rit-Progressive"
  else if (grepl("MS-Rituximab-Stable",      id)) "MS-Rit-Stable"
  else if (grepl("MS-New-PPMS",              id)) "MS-New-PPMS"
  else if (grepl("MS-New",                   id)) "MS-New-RR"
  else if (grepl("Ctrl",                     id)) "Ctrl"
  else                                            "Other"
}
group_pal <- c("Ctrl"               = "#4DAF4A",
               "MS-New-RR"          = "#377EB8",
               "MS-New-PPMS"        = "#984EA3",
               "MS-Rit-Stable"      = "#FF7F00",
               "MS-Rit-Progressive" = "#E41A1C",
               "Other"              = "#999999")

dat <- qc_load(enrich_tsv, frag_tsv, thresh_tsv)
qc_report_marks(dat)
qc_write(dat, out_tsv)
qc_render(dat, out_pdf, parse_group, group_pal, group_name = "Group")
