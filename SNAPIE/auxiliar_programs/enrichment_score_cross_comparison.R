#!/usr/bin/env Rscript
# enrichment_score_cross_comparison.R
#
# Cross-cohort H3K27ac cfChIP-seq enrichment score comparison.
#
# Panel 1 — Protocol comparison (primary):
#   Internal samples broken down by protocol version (V1/V2/V3/V4/1D) vs
#   all external published samples collapsed into one "Literature (Baca et al.)"
#   group. Useful for assessing whether our protocols match published quality.
#
# Panel 2 — Detailed cohort breakdown (secondary):
#   Internal groups (disease) vs external cohort subtypes
#   (cancer type, healthy plasma, post-surgery). Groups ordered by median.
#
# Args:
#   1  internal_tsv  — all_enrichment_scores.tsv  (internal MS cohort)
#   2  external_tsv  — external_enrichment_scores.tsv (GSM samples)
#   3  output_pdf    — output PDF path

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(cowplot)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 3)
  stop("Usage: enrichment_score_cross_comparison.R <internal_tsv> <external_tsv> <output_pdf>")

int_tsv <- args[1]
ext_tsv <- args[2]
out_pdf  <- args[3]

# ── Load both datasets ────────────────────────────────────────────────────────
read_scores <- function(path, dataset_label) {
  read.table(path, header = TRUE, sep = "\t",
             stringsAsFactors = FALSE, comment.char = "") %>%
    mutate(enrichment_score = suppressWarnings(as.numeric(enrichment_score))) %>%
    filter(!is.na(enrichment_score)) %>%
    mutate(dataset = dataset_label)
}

int_dat <- read_scores(int_tsv, "Internal")
ext_dat <- read_scores(ext_tsv, "External")

# ── Parse internal metadata ───────────────────────────────────────────────────
parse_internal_protocol <- function(id) {
  m <- regmatches(id, regexpr("V[0-9]|(?<=[^A-Za-z])[0-9][A-Z]", id, perl = TRUE))
  if (length(m) > 0 && nchar(m) > 0) m else "Unknown"
}

parse_internal_group <- function(id) {
  if      (grepl("MS-Rituximab-Progressive", id)) "MS-Rit-Progressive"
  else if (grepl("MS-Rituximab-Stable",      id)) "MS-Rit-Stable"
  else if (grepl("MS-New-PPMS",              id)) "MS-New-PPMS"
  else if (grepl("MS-New",                   id)) "MS-New-RR"
  else if (grepl("Ctrl",                     id)) "Ctrl"
  else if (grepl("Mix",                      id)) "Mix"
  else                                            "Unknown"
}

int_dat$protocol <- sapply(int_dat$sample_id, parse_internal_protocol)
int_dat$group    <- sapply(int_dat$sample_id, parse_internal_group)

# ── Parse external cohort labels ──────────────────────────────────────────────
parse_external_cohort <- function(id) {
  body <- sub("^GSM[0-9]+_", "", id)
  if      (grepl("^HP",      body)) "Healthy Plasma"
  else if (grepl("^gLC",     body)) "Cancer - Gastric (gLC)"
  else if (grepl("^gPC",     body)) "Cancer - Gastric (gPC)"
  else if (grepl("^mLC",     body)) "Cancer - Lung (mLC)"
  else if (grepl("^mPC",     body)) "Cancer - Pancreatic (mPC)"
  else if (grepl("AMP",      body)) "Cancer - Ampullary"
  else if (grepl("^Merk",    body)) "Cancer - Merkel"
  else if (grepl("^reS",     body)) "Post-surgery"
  else if (grepl("BCDC",     body)) "Cancer - Biliary"
  else if (grepl("TRHC",     body)) "Cancer - HCC"
  else if (grepl("_TRC",     body)) "Cancer - CRC"
  else if (grepl("^K27[A-Z]",body)) "Cancer - Other"
  else                              "External (unclassified)"
}

ext_dat$cohort <- sapply(ext_dat$sample_id, parse_external_cohort)

# ── Shared theme ──────────────────────────────────────────────────────────────
theme_box <- theme_bw(base_size = 11) +
  theme(
    plot.title       = element_text(face = "bold", size = 11),
    plot.subtitle    = element_text(size = 8.5, color = "grey45"),
    axis.title.y     = element_text(size = 9),
    axis.text.x      = element_text(angle = 35, hjust = 1, size = 9),
    legend.position  = "none",
    panel.grid.minor = element_blank()
  )

kw_subtitle <- function(scores, groups) {
  if (length(unique(groups)) < 2) return("")
  tryCatch({
    p <- kruskal.test(scores, as.factor(groups))$p.value
    sprintf("Kruskal-Wallis  p = %.3g", p)
  }, error = function(e) "")
}

# helper: add (n=X) to x-axis labels and return labelled factor + named palette
label_groups <- function(df, group_col, pal_named) {
  levs <- levels(df[[group_col]])
  n_tab <- df %>% count(.data[[group_col]])
  gl_map <- setNames(
    paste0(levs, "\n(n=", n_tab$n[match(levs, n_tab[[group_col]])], ")"),
    levs
  )
  df$group_label <- factor(gl_map[as.character(df[[group_col]])],
                            levels = gl_map)
  pal_out <- setNames(pal_named[levs], gl_map[levs])
  list(df = df, pal = pal_out)
}

# ══════════════════════════════════════════════════════════════════════════════
# PANEL 1 — Protocol comparison: internal protocols vs Literature as one group
# ══════════════════════════════════════════════════════════════════════════════

proto_df <- int_dat %>%
  select(sample_id, enrichment_score, protocol) %>%
  rename(group = protocol) %>%
  mutate(source = "Internal")

lit_df <- ext_dat %>%
  select(sample_id, enrichment_score) %>%
  mutate(group  = "Literature\n(Baca et al.)",
         source = "Literature")

proto_combined <- bind_rows(proto_df, lit_df)

int_proto_levs <- c(
  sort(grep("^V",  unique(int_dat$protocol), value = TRUE)),
  grep("^[0-9]",   unique(int_dat$protocol), value = TRUE),
  if ("Unknown" %in% unique(int_dat$protocol)) "Unknown"
)
int_proto_levs <- int_proto_levs[!is.na(int_proto_levs)]
all_proto_levs <- c(int_proto_levs, "Literature\n(Baca et al.)")
proto_combined$group <- factor(
  proto_combined$group,
  levels = all_proto_levs[all_proto_levs %in% unique(proto_combined$group)]
)

# Colour: Blues for internal protocols, dark green for Literature
blues_ramp <- c("#BBDEFB","#90CAF9","#64B5F6","#42A5F5","#2196F3","#1E88E5","#1565C0","#0D47A1")
n_int      <- sum(int_proto_levs %in% unique(proto_combined$group))
int_pal    <- blues_ramp[seq_len(n_int) + max(0, length(blues_ramp) - n_int)]
proto_pal  <- c(
  setNames(int_pal, int_proto_levs[int_proto_levs %in% unique(proto_combined$group)]),
  "Literature\n(Baca et al.)" = "#2E7D32"
)

res1 <- label_groups(proto_combined, "group", proto_pal)
proto_combined <- res1$df

p1 <- ggplot(proto_combined,
             aes(x = group_label, y = enrichment_score, fill = group_label)) +
  geom_violin(trim = TRUE, scale = "width", alpha = 0.65,
              linewidth = 0.35, width = 0.75) +
  geom_boxplot(width = 0.18, outlier.shape = NA, alpha = 0.9,
               linewidth = 0.4, fill = "white") +
  geom_jitter(width = 0.08, size = 1.0, alpha = 0.55,
              stroke = 0, color = "grey20") +
  stat_summary(fun = mean, geom = "point", shape = 23,
               size = 3, fill = "gold", color = "black", stroke = 0.5) +
  stat_summary(fun = mean, geom = "text",
               aes(label = sprintf("%.2f", after_stat(y))),
               vjust = -1, size = 2.8, fontface = "bold", color = "black") +
  scale_fill_manual(values = res1$pal) +
  labs(
    title    = "Enrichment score: internal protocols vs. Literature (Baca et al.)",
    subtitle = kw_subtitle(proto_combined$enrichment_score,
                           as.character(proto_combined$group)),
    x        = NULL,
    y        = "Enrichment score  (on-target / off-target density)"
  ) +
  theme_box

# ══════════════════════════════════════════════════════════════════════════════
# PANEL 2 — Detailed cohort breakdown (disease groups vs cancer subtypes)
# ══════════════════════════════════════════════════════════════════════════════

int_det <- int_dat %>%
  select(sample_id, enrichment_score, dataset) %>%
  mutate(group = sapply(sample_id, parse_internal_group))

ext_det <- ext_dat %>%
  select(sample_id, enrichment_score, dataset, group = cohort)

detail_combined <- bind_rows(int_det, ext_det)

group_order <- detail_combined %>%
  group_by(group) %>%
  summarise(med = median(enrichment_score, na.rm = TRUE), .groups = "drop") %>%
  arrange(med) %>%
  pull(group)
detail_combined$group <- factor(detail_combined$group, levels = group_order)
all_groups <- levels(detail_combined$group)

int_groups  <- c("Ctrl","MS-New-RR","MS-New-PPMS",
                 "MS-Rit-Stable","MS-Rit-Progressive","Mix","Unknown")
ext_healthy <- "Healthy Plasma"
ext_cancer  <- setdiff(unique(ext_det$group), ext_healthy)

pal <- setNames(rep("grey60", length(all_groups)), all_groups)
int_blues <- c("#BBDEFB","#90CAF9","#64B5F6","#42A5F5","#2196F3","#1E88E5","#1565C0")
int_present <- intersect(int_groups, all_groups)
for (i in seq_along(int_present)) pal[int_present[i]] <- int_blues[min(i, length(int_blues))]
for (g in intersect(ext_healthy, all_groups)) pal[g] <- "#388E3C"
ext_reds <- c("#FFCCBC","#FFAB91","#FF8A65","#FF7043","#F4511E",
              "#E64A19","#D84315","#BF360C","#8D6E63","#795548","#6D4C41")
ext_c_present <- intersect(ext_cancer, all_groups)
for (i in seq_along(ext_c_present)) pal[ext_c_present[i]] <- ext_reds[min(i, length(ext_reds))]

res2 <- label_groups(detail_combined, "group", pal)
detail_combined <- res2$df

p2 <- ggplot(detail_combined,
             aes(x = group_label, y = enrichment_score, fill = group_label)) +
  geom_violin(trim = TRUE, scale = "width", alpha = 0.7,
              linewidth = 0.35, width = 0.8) +
  geom_boxplot(width = 0.18, outlier.shape = NA, alpha = 0.9,
               linewidth = 0.4, fill = "white") +
  geom_jitter(width = 0.08, size = 0.9, alpha = 0.55,
              stroke = 0, color = "grey20") +
  stat_summary(fun = mean, geom = "point", shape = 23,
               size = 2.5, fill = "gold", color = "black", stroke = 0.5) +
  stat_summary(fun = mean, geom = "text",
               aes(label = sprintf("%.2f", after_stat(y))),
               vjust = -1, size = 2.5, fontface = "bold", color = "black") +
  scale_fill_manual(values = res2$pal) +
  labs(
    title    = "Detailed breakdown: MS cohort vs. Baca et al. cohort subtypes",
    subtitle = sprintf("Internal: %d samples  |  External: %d samples",
                       nrow(int_dat), nrow(ext_dat)),
    x        = NULL,
    y        = "Enrichment score  (on-target / off-target density)"
  ) +
  theme_box +
  theme(axis.text.x = element_text(angle = 40, hjust = 1, size = 8.5))

# ── Legend strip (for panel 2) ────────────────────────────────────────────────
legend_df <- data.frame(
  label = c("Internal\n(MS cohort)", "External\nHealthy", "External\nCancer"),
  fill  = c("#1E88E5", "#388E3C", "#F4511E"),
  x = 1:3, y = 1
)
p_legend <- ggplot(legend_df, aes(x = x, y = y, fill = label)) +
  geom_tile(width = 0.8, height = 0.6, color = "white", linewidth = 1) +
  geom_text(aes(label = label), size = 3, fontface = "bold", color = "white") +
  scale_fill_manual(values = setNames(legend_df$fill, legend_df$label)) +
  scale_x_continuous(expand = c(0.3, 0.3)) +
  theme_void() + theme(legend.position = "none")

# ── Save ──────────────────────────────────────────────────────────────────────
dir.create(dirname(out_pdf), recursive = TRUE, showWarnings = FALSE)
pdf(out_pdf,
    width  = max(14, length(all_groups) * 1.1),
    height = 13)
print(plot_grid(
  p1,
  plot_grid(p2, p_legend, ncol = 1, rel_heights = c(10, 1)),
  ncol = 1,
  rel_heights = c(1, 1.2)
))
dev.off()
message(sprintf("Saved: %s  (%d internal, %d external samples)",
                out_pdf, nrow(int_dat), nrow(ext_dat)))
