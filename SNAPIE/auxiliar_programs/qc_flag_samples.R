#!/usr/bin/env Rscript
# qc_flag_samples.R
#
# Data-driven QC flagging for cfChIP-seq samples.
# Thresholds are computed from the user's own distribution — no external
# reference values are hard-coded, so they adapt to whatever on/off-target
# regions were used for the enrichment score.
#
# Threshold method (Kneedle / maximum-distance-from-diagonal):
#   1. Sort values ascending; normalise to [0,1] on both axes.
#   2. Find the point of maximum perpendicular distance from the line
#      connecting (0,0) to (1,1). That is the natural "elbow/knee".
#   Both enrichment score and fragment count use this method independently.
#
# Args:
#   1  enrichment_tsv   — all_enrichment_scores.tsv
#   2  frag_counts_tsv  — fragment_counts.tsv  (sample_id, n_fragments)
#   3  output_tsv       — qc_summary.tsv
#   4  output_pdf       — qc_distributions.pdf

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(cowplot)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 4)
  stop("Usage: qc_flag_samples.R <enrichment_tsv> <frag_counts_tsv> <output_tsv> <output_pdf>")

enrich_tsv  <- args[1]
frag_tsv    <- args[2]
out_tsv     <- args[3]
out_pdf     <- args[4]

# ── Kneedle: threshold = value at max perpendicular distance from diagonal ────
kneedle_threshold <- function(x) {
  x   <- sort(x[!is.na(x)])
  n   <- length(x)
  if (n < 4) return(min(x))
  xs  <- (seq_along(x) - 1) / (n - 1)       # normalised rank [0,1]
  ys  <- (x - min(x)) / (max(x) - min(x))   # normalised value [0,1]
  # perpendicular distance from y = x line
  dist <- abs(ys - xs) / sqrt(2)
  x[which.max(dist)]
}

# ── Load data ─────────────────────────────────────────────────────────────────
enrich <- read.table(enrich_tsv, header = TRUE, sep = "\t",
                     stringsAsFactors = FALSE, comment.char = "") %>%
  mutate(enrichment_score = suppressWarnings(as.numeric(enrichment_score))) %>%
  filter(!is.na(enrichment_score))

frags <- read.table(frag_tsv, header = TRUE, sep = "\t",
                    stringsAsFactors = FALSE, comment.char = "")

dat <- merge(enrich, frags, by = "sample_id", all = TRUE)

# ── Compute thresholds ────────────────────────────────────────────────────────
enrich_thresh <- 7  # fixed: enriched above background; knee-point was biased by control-heavy distribution
frag_thresh   <- 1e6   # fixed at 1M fragments (hard floor for cfChIP-seq)

message(sprintf("Enrichment score threshold (knee-point): %.4f", enrich_thresh))
message(sprintf("Fragment count threshold   (fixed):      %s",
                format(frag_thresh, big.mark = ",")))

# ── Flag samples ──────────────────────────────────────────────────────────────
dat <- dat %>%
  mutate(
    enrich_pass = !is.na(enrichment_score) & enrichment_score >= enrich_thresh,
    frag_pass   = !is.na(n_fragments)      & n_fragments      >= frag_thresh,
    qc_pass     = enrich_pass & frag_pass
  )

n_pass <- sum(dat$qc_pass, na.rm = TRUE)
n_tot  <- nrow(dat)
message(sprintf("QC pass: %d / %d samples", n_pass, n_tot))

# ── Write QC summary ─────────────────────────────────────────────────────────
dir.create(dirname(out_tsv), recursive = TRUE, showWarnings = FALSE)
write.table(dat %>%
              select(sample_id, n_fragments, enrichment_score,
                     frag_pass, enrich_pass, qc_pass),
            out_tsv, sep = "\t", quote = FALSE, row.names = FALSE)
message("Written: ", out_tsv)

# ── Parse group for colouring ─────────────────────────────────────────────────
parse_group <- function(id) {
  if      (grepl("MS-Rituximab-Progressive", id)) "MS-Rit-Progressive"
  else if (grepl("MS-Rituximab-Stable",      id)) "MS-Rit-Stable"
  else if (grepl("MS-New-PPMS",              id)) "MS-New-PPMS"
  else if (grepl("MS-New",                   id)) "MS-New-RR"
  else if (grepl("Ctrl",                     id)) "Ctrl"
  else                                            "Other"
}
dat$group <- sapply(dat$sample_id, parse_group)
group_pal <- c("Ctrl"               = "#4DAF4A",
               "MS-New-RR"          = "#377EB8",
               "MS-New-PPMS"        = "#984EA3",
               "MS-Rit-Stable"      = "#FF7F00",
               "MS-Rit-Progressive" = "#E41A1C",
               "Other"              = "#999999")

# ── Plot 1: enrichment score distribution with threshold ──────────────────────
p_enrich <- ggplot(dat, aes(x = enrichment_score, fill = group, colour = group)) +
  geom_histogram(bins = 30, alpha = 0.7, linewidth = 0.2, position = "stack") +
  geom_vline(xintercept = enrich_thresh, colour = "black", linewidth = 0.9,
             linetype = "dashed") +
  annotate("text", x = enrich_thresh, y = Inf,
           label = sprintf(" threshold\n %.3f", enrich_thresh),
           hjust = 0, vjust = 1.4, size = 3, fontface = "bold") +
  scale_fill_manual(values   = group_pal[names(group_pal) %in% dat$group]) +
  scale_colour_manual(values = group_pal[names(group_pal) %in% dat$group]) +
  labs(
    title    = sprintf("Enrichment score distribution  (threshold = %.1f  [fixed])",
                       enrich_thresh),
    subtitle = sprintf("%d / %d samples pass  |  %d fail",
                       sum(dat$enrich_pass, na.rm=TRUE), n_tot,
                       sum(!dat$enrich_pass, na.rm=TRUE)),
    x = "Enrichment score  (on-target / off-target density)",
    y = "Number of samples",
    fill = "Group", colour = "Group"
  ) +
  theme_bw(base_size = 11) +
  theme(panel.grid.minor = element_blank(),
        plot.title    = element_text(face = "bold", size = 10),
        plot.subtitle = element_text(size = 8.5, colour = "grey40"))

# ── Plot 2: fragment count distribution with threshold ────────────────────────
frag_M <- dat$n_fragments / 1e6   # display in millions

p_frags <- ggplot(dat, aes(x = frag_M, fill = group, colour = group)) +
  geom_histogram(bins = 30, alpha = 0.7, linewidth = 0.2, position = "stack") +
  geom_vline(xintercept = frag_thresh / 1e6, colour = "black", linewidth = 0.9,
             linetype = "dashed") +
  annotate("text", x = frag_thresh / 1e6, y = Inf,
           label = sprintf(" threshold\n %.1f M", frag_thresh / 1e6),
           hjust = 0, vjust = 1.4, size = 3, fontface = "bold") +
  scale_fill_manual(values   = group_pal[names(group_pal) %in% dat$group]) +
  scale_colour_manual(values = group_pal[names(group_pal) %in% dat$group]) +
  labs(
    title    = sprintf("Fragment count distribution  (threshold = %.0f M  [fixed])",
                       frag_thresh / 1e6),
    subtitle = sprintf("%d / %d samples pass  |  %d fail",
                       sum(dat$frag_pass, na.rm=TRUE), n_tot,
                       sum(!dat$frag_pass, na.rm=TRUE)),
    x = "Fragment count (millions)",
    y = "Number of samples",
    fill = "Group", colour = "Group"
  ) +
  theme_bw(base_size = 11) +
  theme(panel.grid.minor = element_blank(),
        plot.title    = element_text(face = "bold", size = 10),
        plot.subtitle = element_text(size = 8.5, colour = "grey40"))

# ── Plot 3: scatter enrichment vs fragments, coloured by QC pass/fail ─────────
dat$qc_label <- ifelse(dat$qc_pass, "PASS", "FAIL")

p_scatter <- ggplot(dat, aes(x = n_fragments / 1e6, y = enrichment_score,
                              colour = group, shape = qc_label)) +
  geom_point(size = 2.5, alpha = 0.85, stroke = 0.5) +
  geom_vline(xintercept = frag_thresh   / 1e6, linetype = "dashed",
             colour = "grey30", linewidth = 0.6) +
  geom_hline(yintercept = enrich_thresh,        linetype = "dashed",
             colour = "grey30", linewidth = 0.6) +
  scale_colour_manual(values = group_pal[names(group_pal) %in% dat$group]) +
  scale_shape_manual(values = c("PASS" = 16, "FAIL" = 4)) +
  labs(
    title    = "QC: enrichment score vs. fragment count",
    subtitle = sprintf("%d PASS  |  %d FAIL  (dashed lines = knee-point thresholds)",
                       sum(dat$qc_pass, na.rm=TRUE),
                       sum(!dat$qc_pass, na.rm=TRUE)),
    x      = "Fragment count (millions)",
    y      = "Enrichment score",
    colour = "Group",
    shape  = "QC"
  ) +
  theme_bw(base_size = 11) +
  theme(panel.grid.minor = element_blank(),
        plot.title    = element_text(face = "bold", size = 10),
        plot.subtitle = element_text(size = 8.5, colour = "grey40"))

# ── Plot 4: QC pass/fail summary per disease group ────────────────────────────
summary_df <- dat %>%
  mutate(qc_label = ifelse(qc_pass, "PASS", "FAIL")) %>%
  count(group, qc_label) %>%
  group_by(group) %>%
  mutate(total     = sum(n),
         pct       = round(100 * n / total),
         group_lab = sprintf("%s\n(n=%d)", group, total[1])) %>%
  ungroup()

group_order_sum <- summary_df %>%
  filter(qc_label == "PASS") %>%
  arrange(desc(n)) %>%
  pull(group)
summary_df$group_lab <- factor(
  summary_df$group_lab,
  levels = unique(summary_df$group_lab[match(group_order_sum,
                                              summary_df$group[summary_df$qc_label == "PASS"])])
)

p_summary <- ggplot(summary_df,
                    aes(x = group_lab, y = n, fill = qc_label)) +
  geom_col(width = 0.65, colour = "white", linewidth = 0.3) +
  geom_text(aes(label = sprintf("%d\n(%d%%)", n, pct)),
            position = position_stack(vjust = 0.5),
            size = 3, fontface = "bold", colour = "white") +
  scale_fill_manual(values = c("PASS" = "#388E3C", "FAIL" = "#C62828"),
                    name = "QC") +
  labs(
    title    = sprintf("QC outcome by disease group  (%d PASS / %d total)",
                       n_pass, n_tot),
    subtitle = sprintf("Enrichment score >= %.2f  &  fragments >= %.0f M",
                       enrich_thresh, frag_thresh / 1e6),
    x = NULL,
    y = "Number of samples"
  ) +
  theme_bw(base_size = 11) +
  theme(panel.grid.minor  = element_blank(),
        panel.grid.major.x = element_blank(),
        plot.title    = element_text(face = "bold", size = 10),
        plot.subtitle = element_text(size = 8.5, colour = "grey40"),
        axis.text.x   = element_text(size = 9))

# ── Save ──────────────────────────────────────────────────────────────────────
dir.create(dirname(out_pdf), recursive = TRUE, showWarnings = FALSE)
pdf(out_pdf, width = 10, height = 14)
print(plot_grid(p_enrich, p_frags, p_scatter, p_summary,
                ncol = 1, align = "v"))
dev.off()
message("Written: ", out_pdf)
