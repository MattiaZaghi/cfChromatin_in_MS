#!/usr/bin/env Rscript
# V4-protocol QC pass summary, H3K27ac only, split MS vs Ctrl.
#
# Standalone reporting script (not wired into the Snakefile) — regenerates the
# collaborator-facing figure for the V4 pass rates from the QC table the
# pipeline already writes.
#
# Usage: Rscript v4_qc_pass_summary.R <qc_summary.tsv> <out_prefix>

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(cowplot)
})

args <- commandArgs(trailingOnly = TRUE)
qc_file <- if (length(args) >= 1) args[1] else "/date/gcb/gcb_MZ/SNAPIE/qc/qc_summary.tsv"
out_pre <- if (length(args) >= 2) args[2] else "/date/gcb/gcb_MZ/SNAPIE/reports/qc/v4_qc_pass_MS_vs_Ctrl"

dat <- read.delim(qc_file, stringsAsFactors = FALSE)

# Protocol from the sample ID: V1/V2/V3/V4, or the older "1D"-style tag.
protocol <- function(x) {
  m <- regmatches(x, regexpr("V[0-9]|(?<=[^A-Za-z])[0-9][A-Z]", x, perl = TRUE))
  ifelse(lengths(m) == 0, NA_character_, unlist(lapply(m, function(z) if (length(z)) z else NA)))
}

# Disease from the sample ID. Deliberately broader than parse_group() in
# qc_flag_samples.R: the V4 Natalizumab and bare "MS-Unknown" IDs match no
# disease pattern there and fall through to "Unknown", but they are MS patients
# and must be tallied as such here.
disease <- function(x) {
  ifelse(grepl("Ctrl", x), "Ctrl",
  ifelse(grepl("MS-",  x), "MS", "Other"))
}

dat <- dat %>%
  mutate(
    protocol = vapply(sample_id, function(s) {
      m <- regmatches(s, regexpr("V[0-9]|(?<=[^A-Za-z])[0-9][A-Z]", s, perl = TRUE))
      if (length(m)) m else NA_character_
    }, character(1)),
    disease  = disease(sample_id),
    qc_pass     = toupper(trimws(as.character(qc_pass)))     == "TRUE",
    enrich_pass = toupper(trimws(as.character(enrich_pass)))  == "TRUE",
    frag_pass   = toupper(trimws(as.character(frag_pass)))    == "TRUE"
  )

v4 <- dat %>% filter(protocol == "V4", mark == "H3K27ac", disease %in% c("Ctrl", "MS"))
stopifnot(nrow(v4) > 0)

grp_lab <- v4 %>%
  group_by(disease) %>%
  summarise(n = n(), pass = sum(qc_pass), .groups = "drop") %>%
  mutate(lab = sprintf("%s\nn = %d", disease, n),
         pct = sprintf("%d / %d\n(%.0f%%)", pass, n, 100 * pass / n))

lab_map <- setNames(grp_lab$lab, grp_lab$disease)

# --- panel a: pass / fail per group ------------------------------------------
pf <- v4 %>%
  mutate(status = factor(ifelse(qc_pass, "Pass", "Fail"), levels = c("Fail", "Pass"))) %>%
  count(disease, status) %>%
  mutate(group = lab_map[disease])

p_pass <- ggplot(pf, aes(x = group, y = n, fill = status)) +
  geom_col(width = 0.62, colour = "grey25", linewidth = 0.3) +
  geom_text(aes(label = n), position = position_stack(vjust = 0.5),
            size = 4.4, fontface = "bold", colour = "white") +
  geom_text(data = grp_lab, aes(x = lab, y = n + 2, label = pct),
            inherit.aes = FALSE, size = 4, lineheight = 0.95, vjust = 0) +
  scale_fill_manual(values = c(Pass = "#2A7F62", Fail = "#B0562E"), name = NULL) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.22))) +
  labs(title = "a  Samples passing QC", x = NULL, y = "Samples") +
  theme_cowplot(13) +
  # legends below the panel: with a legend on top, panels a/b would push their
  # titles up a row relative to c, which has no legend.
  theme(legend.position = "bottom")

# --- panel b: why the failures failed ---------------------------------------
reason <- v4 %>%
  filter(!qc_pass) %>%
  mutate(reason = case_when(
    !enrich_pass &  frag_pass ~ "Enrichment only",
     enrich_pass & !frag_pass ~ "Fragments only",
    TRUE                      ~ "Both"
  )) %>%
  count(disease, reason) %>%
  mutate(reason = factor(reason, levels = c("Enrichment only", "Fragments only", "Both")),
         group  = lab_map[disease])

p_reason <- ggplot(reason, aes(x = group, y = n, fill = reason)) +
  geom_col(width = 0.62, colour = "grey25", linewidth = 0.3) +
  geom_text(aes(label = n), position = position_stack(vjust = 0.5),
            size = 4.2, fontface = "bold", colour = "white") +
  scale_fill_manual(values = c("Enrichment only" = "#B0562E",
                               "Fragments only"  = "#C9A227",
                               "Both"            = "#6B3A2E"), name = NULL) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.12))) +
  labs(title = "b  Which cut-off the failures missed",
       x = NULL, y = "Samples failing QC") +
  theme_cowplot(13) +
  theme(legend.position = "bottom")

# --- panel c: the underlying score distribution -----------------------------
ymax <- max(v4$enrichment_score)
med <- v4 %>% group_by(disease) %>%
  summarise(m = median(enrichment_score), .groups = "drop") %>%
  mutate(group = lab_map[disease])

p_score <- ggplot(v4, aes(x = lab_map[disease], y = enrichment_score, fill = disease)) +
  geom_hline(yintercept = 7, linetype = "dashed", colour = "#B0562E", linewidth = 0.6) +
  geom_boxplot(width = 0.45, outlier.shape = NA, alpha = 0.55, linewidth = 0.4) +
  geom_jitter(width = 0.14, height = 0, size = 1.5, alpha = 0.7, shape = 21,
              colour = "grey20", stroke = 0.2) +
  # medians above each box rather than beside it, so they clear the points
  geom_text(data = med, aes(x = group, y = ymax * 1.9, label = sprintf("median %.1f", m)),
            inherit.aes = FALSE, size = 3.8) +
  scale_fill_manual(values = c(Ctrl = "#3B6EA5", MS = "#A53B5B"), guide = "none") +
  scale_y_continuous(trans = "log10", expand = expansion(mult = c(0.05, 0.22))) +
  annotate("text", x = 1.5, y = 7, label = "cut-off 7", vjust = -0.6,
           size = 3.5, colour = "#B0562E") +
  labs(title = "c  Enrichment score distribution", x = NULL, y = "Enrichment score") +
  theme_cowplot(13)

cap <- paste0(
  "V4-protocol H3K27ac libraries (n = 108). QC pass requires enrichment score >= 7 AND >= 1,000,000 unique fragments.\n",
  "Sequencing depth is equivalent between groups (median 2.84M vs 2.85M fragments), so the gap is enrichment, not depth. ",
  "MS includes the Natalizumab and MS-Unknown V4 IDs."
)

fig <- plot_grid(p_pass, p_reason, p_score, nrow = 1,
                 rel_widths = c(1, 1, 1.1), align = "h", axis = "tb")
fig <- plot_grid(
  fig,
  ggdraw() + draw_label(cap, x = 0.008, hjust = 0, vjust = 0.5, size = 9.5,
                            colour = "grey30", lineheight = 1.15),
  ncol = 1, rel_heights = c(1, 0.11)
)

dir.create(dirname(out_pre), recursive = TRUE, showWarnings = FALSE)
ggsave(paste0(out_pre, ".pdf"), fig, width = 15, height = 5.4, device = "pdf")
ggsave(paste0(out_pre, ".png"), fig, width = 15, height = 5.4, dpi = 200)

message("wrote ", out_pre, ".pdf / .png")
print(grp_lab)
print(reason)
