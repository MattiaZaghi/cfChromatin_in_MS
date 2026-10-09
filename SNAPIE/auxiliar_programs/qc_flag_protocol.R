#!/usr/bin/env Rscript
# qc_flag_protocol.R
#
# As qc_flag_samples.R, but coloured by library protocol version (V1–V4)
# instead of disease group. Same per-mark cut-offs, same one-page-per-mark
# layout; all shared machinery is in qc_flag_common.R.
#
# Args:
#   1  enrichment_tsv   — all_enrichment_scores.tsv
#   2  frag_counts_tsv  — fragment_counts.tsv
#   3  output_tsv       — qc_summary_protocol.tsv
#   4  output_pdf       — qc_distributions_protocol.pdf
#   5  thresholds_tsv   — qc_thresholds.tsv  (mark, enrich_thresh, frag_thresh)

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 4)
  stop("Usage: qc_flag_protocol.R <enrichment_tsv> <frag_counts_tsv> ",
       "<output_tsv> <output_pdf> [thresholds_tsv]")

enrich_tsv <- args[1]
frag_tsv   <- args[2]
out_tsv    <- args[3]
out_pdf    <- args[4]
thresh_tsv <- if (length(args) >= 5) args[5] else NA

.this <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])
source(file.path(dirname(.this), "qc_flag_common.R"))

# ── Protocol version, from the sample ID ────────────────────────────────────
parse_group <- function(id) {
  if      (grepl("V1", id)) "V1"
  else if (grepl("V2", id)) "V2"
  else if (grepl("V3", id)) "V3"
  else if (grepl("V4", id)) "V4"
  else                      "Other"
}
group_pal <- c("V1"    = "#4DAF4A",
               "V2"    = "#377EB8",
               "V3"    = "#984EA3",
               "V4"    = "#FF7F00",
               "Other" = "#999999")

dat <- qc_load(enrich_tsv, frag_tsv, thresh_tsv)
qc_report_marks(dat)
qc_write(dat, out_tsv)
qc_render(dat, out_pdf, parse_group, group_pal, group_name = "Protocol")
