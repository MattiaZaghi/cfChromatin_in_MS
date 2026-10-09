#!/usr/bin/env Rscript
# enrichment_score_comparison_plots.R
#
# Three comparison boxplots of enrichment scores across the cohort, on one page
# per enrichment mark:
#   1. By protocol version  (V1 / V2 / V3 / V4 / 1D — parsed from sample IDs)
#   2. By sample type       (Plasma / CSF — parsed from sample IDs)
#   3. By disease group     (Ctrl / MS-New / MS-Rit-Stable / … — parsed from IDs)
#
# Protocol version and sample type are inferred from the sample ID naming
# convention: {PatientID}-{P|C}-{Disease…}-{V1|V2|V3|V4}-{Amount}
#
# A Kruskal-Wallis p-value is annotated on each plot, computed within a mark:
# the on/off-target region sets differ per mark, so a test pooling marks would
# mostly measure the region sets rather than the biology.
#
# Args:
#   1  aggregate_tsv  — TSV with columns: sample_id, mark, enrichment_score (+ others)
#   2  output_pdf     — output PDF path

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(cowplot)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) {
  stop("Usage: enrichment_score_comparison_plots.R <aggregate_tsv> <output_pdf>")
}
agg_tsv <- args[1]
out_pdf  <- args[2]

# ── Load ──────────────────────────────────────────────────────────────────────
dat <- read.table(agg_tsv, header = TRUE, sep = "\t",
                  stringsAsFactors = FALSE, comment.char = "") %>%
  mutate(enrichment_score = suppressWarnings(as.numeric(enrichment_score))) %>%
  filter(!is.na(enrichment_score))

# Aggregates written before the score became mark-aware have no 'mark' column.
if (!"mark" %in% names(dat)) dat$mark <- "H3K27ac"

if (nrow(dat) == 0) {
  message("No valid enrichment scores. Writing blank PDF.")
  dir.create(dirname(out_pdf), recursive = TRUE, showWarnings = FALSE)
  pdf(out_pdf, width = 14, height = 5); plot.new(); dev.off()
  quit(status = 0)
}

# ── Parse sample metadata from IDs ───────────────────────────────────────────
# Sample ID convention:
#   {PatientID}-{P|C?}-{Disease…}-{V1|V2|V3|V4|1D?}-{Amount?}
# e.g. "04061-C-MS-Rituximab-Stable-V1-100"
#      "H1-P-Ctrl-V2-1000"
#      "15032-MS-Rituximab-Stable-V4-1000"   (no P/C separator)
#      "H6-P-Ctrl-1D-1000"                   (1D protocol)

parse_meta <- function(id) {
  # Protocol version: first V[0-9] or digit followed by uppercase letter (e.g. 1D)
  proto_match <- regmatches(id, regexpr("V[0-9]|(?<=[^A-Za-z])[0-9][A-Z]", id, perl = TRUE))
  proto <- if (length(proto_match) > 0 && nchar(proto_match) > 0) proto_match else "Unknown"

  # Sample type: second field after splitting on '-', if exactly one letter P or C
  parts <- strsplit(id, "-")[[1]]
  stype <- if (length(parts) >= 2 && nchar(parts[2]) == 1 && parts[2] %in% c("P", "C")) {
    c("P" = "Plasma", "C" = "CSF")[parts[2]]
  } else {
    "Unknown"
  }

  # Disease group: matched from longest to shortest to avoid substring collisions
  disease <- if      (grepl("MS-Rituximab-Progressive", id)) "MS-Rit-Progressive"
             else if (grepl("MS-Rituximab-Stable",      id)) "MS-Rit-Stable"
             else if (grepl("MS-New-PPMS",              id)) "MS-New-PPMS"
             else if (grepl("MS-New",                   id)) "MS-New-RR"
             else if (grepl("Ctrl",                     id)) "Ctrl"
             else if (grepl("Mix",                      id)) "Mix"
             else                                            "Unknown"

  data.frame(sample_id   = id,
             protocol    = proto,
             sample_type = stype,
             disease     = disease,
             stringsAsFactors = FALSE)
}

meta <- do.call(rbind, lapply(dat$sample_id, parse_meta))
dat  <- merge(dat, meta, by = "sample_id", all.x = TRUE)

# ── Factor levels ─────────────────────────────────────────────────────────────
proto_levs <- c(sort(grep("^V", unique(dat$protocol), value = TRUE)),
                grep("^[0-9]", unique(dat$protocol), value = TRUE),
                "Unknown")
proto_levs <- proto_levs[proto_levs %in% unique(dat$protocol)]
dat$protocol    <- factor(dat$protocol,    levels = unique(proto_levs))

dat$sample_type <- factor(dat$sample_type,
                           levels = c("Plasma", "CSF", "Unknown"))

dis_levs <- c("Ctrl", "MS-New-RR", "MS-New-PPMS",
              "MS-Rit-Stable", "MS-Rit-Progressive", "Mix", "Unknown")
dat$disease <- factor(dat$disease,
                       levels = dis_levs[dis_levs %in% unique(dat$disease)])

# ── Shared theme ──────────────────────────────────────────────────────────────
theme_box <- theme_bw(base_size = 22) +
  theme(
    plot.title    = element_text(face = "bold", size = 22),
    plot.subtitle = element_text(size = 18, color = "grey45"),
    axis.title.y  = element_text(size = 20),
    axis.text.x   = element_text(angle = 35, hjust = 1, size = 18),
    legend.position  = "none",
    panel.grid.minor = element_blank()
  )

# Kruskal-Wallis subtitle helper
kw_subtitle <- function(data, group_col) {
  g <- data[[group_col]]
  if (length(unique(g[!is.na(g)])) < 2) return("")
  tryCatch({
    p <- kruskal.test(data$enrichment_score, g)$p.value
    sprintf("Kruskal-Wallis  p = %.3g", p)
  }, error = function(e) "")
}

dis_pal <- c(
  "Ctrl"               = "#4DAF4A",
  "MS-New-RR"          = "#377EB8",
  "MS-New-PPMS"        = "#984EA3",
  "MS-Rit-Stable"      = "#FF7F00",
  "MS-Rit-Progressive" = "#E41A1C",
  "Mix"                = "#999999",
  "Unknown"            = "#CCCCCC"
)

# ── One page of three boxplots for a single mark ──────────────────────────────
# Each page gets its own y-axis scale and its own Kruskal-Wallis tests, because
# scores are only comparable within a mark.
mark_page <- function(d, mark) {
  d$protocol    <- droplevels(d$protocol)
  d$sample_type <- droplevels(d$sample_type)
  d$disease     <- droplevels(d$disease)
  y_lab <- "Enrichment score\n(on-target / off-target density)"

  p1 <- ggplot(d, aes(x = protocol, y = enrichment_score, fill = protocol)) +
    geom_boxplot(outlier.shape = NA, alpha = 0.65, linewidth = 0.45, width = 0.55) +
    geom_jitter(width = 0.18, size = 1.7, alpha = 0.75, stroke = 0.3, color = "grey25") +
    scale_fill_brewer(palette = "Set1") +
    labs(title    = sprintf("%s - by protocol version", mark),
         subtitle = kw_subtitle(d, "protocol"),
         x = "Protocol", y = y_lab) +
    theme_box

  p2 <- ggplot(d, aes(x = sample_type, y = enrichment_score, fill = sample_type)) +
    geom_boxplot(outlier.shape = NA, alpha = 0.65, linewidth = 0.45, width = 0.45) +
    geom_jitter(width = 0.15, size = 1.7, alpha = 0.75, stroke = 0.3, color = "grey25") +
    scale_fill_brewer(palette = "Set2") +
    labs(title    = sprintf("%s - by sample type", mark),
         subtitle = kw_subtitle(d, "sample_type"),
         x = "Sample type", y = y_lab) +
    theme_box

  p3 <- ggplot(d, aes(x = disease, y = enrichment_score, fill = disease)) +
    geom_boxplot(outlier.shape = NA, alpha = 0.65, linewidth = 0.45, width = 0.55) +
    geom_jitter(width = 0.18, size = 1.7, alpha = 0.75, stroke = 0.3, color = "grey25") +
    scale_fill_manual(values = dis_pal[levels(d$disease)]) +
    labs(title    = sprintf("%s - by disease group  (n = %d)", mark, nrow(d)),
         subtitle = kw_subtitle(d, "disease"),
         x = "Disease group", y = y_lab) +
    theme_box

  plot_grid(p1, p2, p3, nrow = 1, align = "h", axis = "tb")
}

# ── Render one page per mark ───────────────────────────────────────────────────
dir.create(dirname(out_pdf), recursive = TRUE, showWarnings = FALSE)
pdf(out_pdf, width = 24, height = 10)
for (m in sort(unique(dat$mark))) {
  d <- dat[dat$mark == m, , drop = FALSE]
  message(sprintf("[%s] n=%d  median=%.2f  range=%.2f–%.2f",
                  m, nrow(d), median(d$enrichment_score),
                  min(d$enrichment_score), max(d$enrichment_score)))
  print(mark_page(d, m))
}
dev.off()
message("Saved: ", out_pdf)
