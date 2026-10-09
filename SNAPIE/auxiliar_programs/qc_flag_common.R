#!/usr/bin/env Rscript
# qc_flag_common.R
#
# Shared machinery for qc_flag_samples.R and qc_flag_protocol.R. The two
# differ only in how they group samples for colouring (disease group vs
# protocol version); everything else — loading, per-mark threshold joining,
# flagging and the 4-panel page layout — lives here.
#
# Cut-offs are NOT derived from the data. They are read from
# config['qc_thresholds'] via the thresholds TSV written by
# rule qc_thresholds_table, and applied independently within each enrichment
# mark. Scores from different marks are computed against different chromHMM
# region sets and are only comparable within a mark, which is why every
# distribution, threshold line and pass/fail tally below is per-mark.

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(cowplot)
})

DEFAULT_MARK          <- "H3K27ac"
FALLBACK_ENRICH_THRESH <- 7
FALLBACK_FRAG_THRESH   <- 1e6

# ── Load enrichment scores, fragment counts and per-mark cut-offs ────────────
qc_load <- function(enrich_tsv, frag_tsv, thresh_tsv = NA) {

  enrich <- read.table(enrich_tsv, header = TRUE, sep = "\t",
                       stringsAsFactors = FALSE, comment.char = "") %>%
    mutate(enrichment_score = suppressWarnings(as.numeric(enrichment_score))) %>%
    filter(!is.na(enrichment_score))

  # Aggregates written before the enrichment score became mark-aware have no
  # 'mark' column; treat those as H3K27ac so old outputs still render.
  if (!"mark" %in% names(enrich)) enrich$mark <- DEFAULT_MARK

  frags <- read.table(frag_tsv, header = TRUE, sep = "\t",
                      stringsAsFactors = FALSE, comment.char = "")

  dat <- merge(enrich, frags, by = "sample_id", all = TRUE)
  dat$mark[is.na(dat$mark)] <- DEFAULT_MARK

  # Per-mark cut-offs
  thr <- NULL
  if (!is.na(thresh_tsv) && nzchar(thresh_tsv) && file.exists(thresh_tsv)) {
    thr <- read.table(thresh_tsv, header = TRUE, sep = "\t",
                      stringsAsFactors = FALSE, comment.char = "")
  }
  if (is.null(thr) || nrow(thr) == 0) {
    message("No thresholds TSV; falling back to ", FALLBACK_ENRICH_THRESH,
            " / ", format(FALLBACK_FRAG_THRESH, big.mark = ","))
    thr <- data.frame(mark = unique(dat$mark),
                      enrich_thresh = FALLBACK_ENRICH_THRESH,
                      frag_thresh   = FALLBACK_FRAG_THRESH,
                      stringsAsFactors = FALSE)
  }

  dat <- merge(dat, thr, by = "mark", all.x = TRUE)
  dat$enrich_thresh[is.na(dat$enrich_thresh)] <- FALLBACK_ENRICH_THRESH
  dat$frag_thresh[is.na(dat$frag_thresh)]     <- FALLBACK_FRAG_THRESH

  dat %>%
    mutate(
      enrich_pass = !is.na(enrichment_score) & enrichment_score >= enrich_thresh,
      frag_pass   = !is.na(n_fragments)      & n_fragments      >= frag_thresh,
      qc_pass     = enrich_pass & frag_pass
    )
}

# ── Log per-mark distribution, so cut-offs can be re-tuned from the run log ──
qc_report_marks <- function(dat) {
  for (m in sort(unique(dat$mark))) {
    d <- dat[dat$mark == m, ]
    s <- d$enrichment_score[!is.na(d$enrichment_score)]
    message(sprintf(
      "[%s] n=%d  cut-offs: enrichment >= %.2f, fragments >= %s  |  pass %d/%d",
      m, nrow(d), d$enrich_thresh[1],
      format(d$frag_thresh[1], big.mark = ",", scientific = FALSE),
      sum(d$qc_pass, na.rm = TRUE), nrow(d)))
    if (length(s) > 0) {
      q <- quantile(s, c(0, .10, .25, .50, .75, .90, 1))
      message(sprintf(
        "[%s] score  min=%.2f  Q10=%.2f  Q25=%.2f  median=%.2f  Q75=%.2f  Q90=%.2f  max=%.2f",
        m, q[1], q[2], q[3], q[4], q[5], q[6], q[7]))
    }
  }
}

qc_write <- function(dat, out_tsv) {
  dir.create(dirname(out_tsv), recursive = TRUE, showWarnings = FALSE)
  write.table(dat %>%
                arrange(mark, sample_id) %>%
                # plain integer, not 1e+06, so the column is greppable
                mutate(frag_thresh = format(frag_thresh, scientific = FALSE,
                                            trim = TRUE)) %>%
                select(sample_id, mark, n_fragments, enrichment_score,
                       enrich_thresh, frag_thresh,
                       frag_pass, enrich_pass, qc_pass),
              out_tsv, sep = "\t", quote = FALSE, row.names = FALSE)
  message("Written: ", out_tsv)
}

# ── One 4-panel page for a single mark ───────────────────────────────────────
# d           subset of the QC table for this mark (already has $group)
# mark        mark name, used in every panel title
# group_pal   named colour vector for $group
# group_name  legend title ("Group" / "Protocol")
qc_page <- function(d, mark, group_pal, group_name = "Group") {

  n_tot         <- nrow(d)
  n_pass        <- sum(d$qc_pass, na.rm = TRUE)
  enrich_thresh <- d$enrich_thresh[1]
  frag_thresh   <- d$frag_thresh[1]
  pal           <- group_pal[names(group_pal) %in% d$group]

  base <- theme_bw(base_size = 22) +
    theme(panel.grid.minor = element_blank(),
          plot.title    = element_text(face = "bold", size = 20),
          plot.subtitle = element_text(size = 16, colour = "grey40"),
          legend.title  = element_text(size = 18),
          legend.text   = element_text(size = 16))

  # Panel 1 — enrichment score distribution
  p_enrich <- ggplot(d, aes(x = enrichment_score, fill = group, colour = group)) +
    geom_histogram(bins = 30, alpha = 0.7, linewidth = 0.2, position = "stack") +
    geom_vline(xintercept = enrich_thresh, colour = "black", linewidth = 0.9,
               linetype = "dashed") +
    annotate("text", x = enrich_thresh, y = Inf,
             label = sprintf(" threshold\n %.2f", enrich_thresh),
             hjust = 0, vjust = 1.4, size = 7, fontface = "bold") +
    scale_fill_manual(values = pal) + scale_colour_manual(values = pal) +
    labs(title    = sprintf("%s - enrichment score  (threshold = %.2f)",
                            mark, enrich_thresh),
         subtitle = sprintf("%d / %d pass  |  %d fail",
                            sum(d$enrich_pass, na.rm = TRUE), n_tot,
                            sum(!d$enrich_pass, na.rm = TRUE)),
         x = "Enrichment score  (on-target / off-target density)",
         y = "Number of samples", fill = group_name, colour = group_name) +
    base

  # Panel 2 — fragment count distribution
  p_frags <- ggplot(d, aes(x = n_fragments / 1e6, fill = group, colour = group)) +
    geom_histogram(bins = 30, alpha = 0.7, linewidth = 0.2, position = "stack") +
    geom_vline(xintercept = frag_thresh / 1e6, colour = "black", linewidth = 0.9,
               linetype = "dashed") +
    annotate("text", x = frag_thresh / 1e6, y = Inf,
             label = sprintf(" threshold\n %.1f M", frag_thresh / 1e6),
             hjust = 0, vjust = 1.4, size = 7, fontface = "bold") +
    scale_fill_manual(values = pal) + scale_colour_manual(values = pal) +
    labs(title    = sprintf("%s - fragment count  (threshold = %.0f M)",
                            mark, frag_thresh / 1e6),
         subtitle = sprintf("%d / %d pass  |  %d fail",
                            sum(d$frag_pass, na.rm = TRUE), n_tot,
                            sum(!d$frag_pass, na.rm = TRUE)),
         x = "Fragment count (millions)", y = "Number of samples",
         fill = group_name, colour = group_name) +
    base

  # Panel 3 — scatter, both cut-offs
  d$qc_label <- ifelse(d$qc_pass, "PASS", "FAIL")
  p_scatter <- ggplot(d, aes(x = n_fragments / 1e6, y = enrichment_score,
                             colour = group, shape = qc_label)) +
    geom_point(size = 2.5, alpha = 0.85, stroke = 0.5) +
    geom_vline(xintercept = frag_thresh / 1e6, linetype = "dashed",
               colour = "grey30", linewidth = 0.6) +
    geom_hline(yintercept = enrich_thresh, linetype = "dashed",
               colour = "grey30", linewidth = 0.6) +
    scale_colour_manual(values = pal) +
    scale_shape_manual(values = c("PASS" = 16, "FAIL" = 4)) +
    labs(title    = sprintf("%s - enrichment score vs. fragment count", mark),
         subtitle = sprintf("%d PASS  |  %d FAIL  (dashed = cut-offs)",
                            n_pass, n_tot - n_pass),
         x = "Fragment count (millions)", y = "Enrichment score",
         colour = group_name, shape = "QC") +
    base

  # Panel 4 — pass/fail per group
  totals <- d %>% count(group, name = "total")
  summary_df <- d %>%
    mutate(qc_label = ifelse(qc_pass, "PASS", "FAIL")) %>%
    count(group, qc_label) %>%
    left_join(totals, by = "group") %>%
    mutate(pct       = round(100 * n / total),
           group_lab = sprintf("%s\n(n=%d)", group, total))

  ord <- summary_df %>% filter(qc_label == "PASS") %>%
    arrange(desc(n)) %>% pull(group_lab)
  ord <- c(ord, setdiff(unique(summary_df$group_lab), ord))
  summary_df$group_lab <- factor(summary_df$group_lab, levels = ord)

  p_summary <- ggplot(summary_df, aes(x = group_lab, y = n, fill = qc_label)) +
    geom_col(width = 0.65, colour = "white", linewidth = 0.3) +
    geom_text(aes(label = sprintf("%d\n(%d%%)", n, pct)),
              position = position_stack(vjust = 0.5),
              size = 7, fontface = "bold", colour = "white") +
    scale_fill_manual(values = c("PASS" = "#388E3C", "FAIL" = "#C62828"),
                      name = "QC") +
    labs(title    = sprintf("%s - QC outcome by %s  (%d PASS / %d total)",
                            mark, tolower(group_name), n_pass, n_tot),
         subtitle = sprintf("Enrichment score >= %.2f  &  fragments >= %.0f M",
                            enrich_thresh, frag_thresh / 1e6),
         x = NULL, y = "Number of samples") +
    base +
    theme(panel.grid.major.x = element_blank(),
          axis.text.x = element_text(size = 18))

  plot_grid(p_enrich, p_frags, p_scatter, p_summary, ncol = 1, align = "v")
}

# ── Write one page per mark into a single PDF ────────────────────────────────
qc_render <- function(dat, out_pdf, parse_group, group_pal, group_name) {
  dat$group <- sapply(dat$sample_id, parse_group)

  dir.create(dirname(out_pdf), recursive = TRUE, showWarnings = FALSE)
  pdf(out_pdf, width = 16, height = 24)
  on.exit(dev.off(), add = TRUE)

  for (m in sort(unique(dat$mark))) {
    d <- dat[dat$mark == m, ]
    if (nrow(d) == 0 || all(is.na(d$enrichment_score))) {
      plot.new(); title(sprintf("%s - no enrichment score data", m))
      next
    }
    print(qc_page(d, m, group_pal, group_name))
  }
  message("Written: ", out_pdf)
}
