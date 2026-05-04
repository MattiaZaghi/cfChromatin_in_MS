#!/usr/bin/env Rscript
# enrichment_score_plot.R
#
# Per-sample enrichment score QC plot.
# Shows the sorted distribution of enrichment scores across all samples
# (empirical CDF), with the current sample highlighted, quantile reference
# lines, and the knee-point of the distribution marked.
#
# X-axis: enrichment score value (used as cutoff threshold)
# Y-axis: cumulative fraction of samples at or below that score
#
# Args:
#   1  aggregate_tsv  — TSV with columns sample_id, enrichment_score (all samples)
#   2  sample_id      — ID of the sample to highlight
#   3  output_pdf     — Path for the output PDF

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 3) {
  stop("Usage: enrichment_score_plot.R <aggregate_tsv> <sample_id> <output_pdf>")
}

agg_tsv   <- args[1]
sample_id <- args[2]
out_pdf   <- args[3]

# ── Load data ─────────────────────────────────────────────────────────────────
dat <- read.table(agg_tsv, header = TRUE, sep = "\t",
                  stringsAsFactors = FALSE, comment.char = "")

dat <- dat %>%
  filter(!is.na(enrichment_score), enrichment_score != "NA") %>%
  mutate(enrichment_score = as.numeric(enrichment_score)) %>%
  filter(!is.na(enrichment_score)) %>%
  arrange(enrichment_score) %>%
  mutate(
    rank     = seq_len(n()),
    fraction = rank / n()
  )

if (nrow(dat) == 0) {
  message("No valid enrichment scores found. Writing blank PDF.")
  pdf(out_pdf, width = 7, height = 5)
  plot.new()
  title("No enrichment score data available")
  dev.off()
  quit(status = 0)
}

# ── Quantile reference lines ─────────────────────────────────────────────────
q_probs  <- c(0.10, 0.25, 0.50, 0.75, 0.90)
q_vals   <- quantile(dat$enrichment_score, q_probs)
q_df     <- data.frame(
  label = paste0("Q", q_probs * 100),
  val   = q_vals,
  stringsAsFactors = FALSE
)

# ── Knee-point (Kneedle algorithm) ───────────────────────────────────────────
# Normalize both axes to [0, 1], find point with max perpendicular distance
# from the line connecting first and last point.
x <- dat$enrichment_score
y <- dat$fraction

x_range <- max(x) - min(x)
y_range <- max(y) - min(y)

if (x_range == 0) {
  knee_score <- x[1]
} else {
  x_n <- (x - min(x)) / x_range
  y_n <- (y - min(y)) / y_range

  # Direction vector of the chord
  dx <- x_n[length(x_n)] - x_n[1]
  dy <- y_n[length(y_n)] - y_n[1]

  # Signed perpendicular distance
  denom <- sqrt(dx^2 + dy^2)
  dist  <- abs(dy * x_n - dx * y_n +
                 x_n[length(x_n)] * y_n[1] -
                 y_n[length(y_n)] * x_n[1]) / denom

  knee_idx   <- which.max(dist)
  knee_score <- x[knee_idx]
}

knee_frac <- dat$fraction[which.max(abs(dat$enrichment_score - knee_score) ==
                                      min(abs(dat$enrichment_score - knee_score)))[1]]

# ── Current sample ────────────────────────────────────────────────────────────
this_row <- dat %>% filter(sample_id == !!sample_id)
has_score <- nrow(this_row) == 1

# ── Plot ──────────────────────────────────────────────────────────────────────
lt_vals  <- setNames(rep("dotted", length(q_probs)), q_df$label)
lt_vals["Q50"] <- "dashed"

p <- ggplot(dat, aes(x = enrichment_score, y = fraction)) +
  geom_line(color = "grey50", linewidth = 0.9) +

  # Quantile vertical lines
  geom_vline(
    data       = q_df,
    aes(xintercept = val, linetype = label),
    color      = "#1565C0",
    linewidth  = 0.55
  ) +
  scale_linetype_manual(values = lt_vals, name = "Quantile") +

  # Knee-point vertical line
  geom_vline(
    xintercept = knee_score,
    color      = "#2E7D32",
    linetype   = "solid",
    linewidth  = 0.8
  ) +
  annotate(
    "label",
    x      = knee_score,
    y      = 0.10,
    label  = sprintf("Knee\n%.3f", knee_score),
    color  = "#2E7D32",
    fill   = "white",
    label.size = 0.3,
    hjust  = -0.05,
    size   = 3.0,
    fontface = "bold"
  ) +

  labs(
    title   = sprintf("H3K27ac Enrichment Score — %s", sample_id),
    x       = "Enrichment Score  (on-target density / off-target density)",
    y       = "Cumulative fraction of samples"
  ) +
  theme_bw(base_size = 12) +
  theme(
    plot.title    = element_text(face = "bold", size = 12),
    legend.position = c(0.82, 0.25),
    legend.background = element_rect(fill = "white", color = "grey80")
  )

# Overlay current sample point + label
if (has_score) {
  p <- p +
    geom_point(
      data  = this_row,
      aes(x = enrichment_score, y = fraction),
      color = "#E53935", fill = "#E53935",
      shape = 21, size = 4.5, stroke = 1.2
    ) +
    annotate(
      "label",
      x      = this_row$enrichment_score,
      y      = this_row$fraction,
      label  = sprintf("%s\n%.3f", sample_id, this_row$enrichment_score),
      color  = "#E53935",
      fill   = "white",
      label.size = 0.3,
      hjust  = -0.1,
      size   = 3.0
    )
}

dir.create(dirname(out_pdf), recursive = TRUE, showWarnings = FALSE)
ggsave(out_pdf, p, width = 7, height = 5, device = "pdf")
message(sprintf("Saved: %s", out_pdf))
