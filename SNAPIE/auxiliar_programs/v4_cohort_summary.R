#!/usr/bin/env Rscript
# V4 cohort composition: what was profiled in the latest round vs earlier V4 runs.
#
# Descriptive companion slide that precedes the QC-pass figure
# (see v4_qc_pass_summary.R). No pass/fail here by design.
#
# The "latest round" is defined from the samplesheet git history, not guessed
# from sample IDs: commit fc67794 ("samplesheet with new samples added",
# 2026-09-21) introduced 97 sample IDs, 95 of them V4. Regenerate the ID list
# with:
#   git show fc67794~1:SNAPIE/config/samplesheet.csv | tail -n +2 | cut -d, -f1 | sort -u > prev.txt
#   git show fc67794:SNAPIE/config/samplesheet.csv   | tail -n +2 | cut -d, -f1 | sort -u > curr.txt
#   comm -13 prev.txt curr.txt > new_ids.txt
#
# Usage: Rscript v4_cohort_summary.R <qc_summary.tsv> <new_ids.txt> <out_prefix>

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(cowplot)
})

args    <- commandArgs(trailingOnly = TRUE)
qc_file <- if (length(args) >= 1) args[1] else "/date/gcb/gcb_MZ/SNAPIE/qc/qc_summary.tsv"
new_file<- if (length(args) >= 2) args[2] else "/tmp/new_ids.txt"
out_pre <- if (length(args) >= 3) args[3] else "/date/gcb/gcb_MZ/SNAPIE/reports/qc/v4_cohort_summary"

dat      <- read.delim(qc_file, stringsAsFactors = FALSE)
new_ids  <- readLines(new_file)
new_ids  <- trimws(new_ids[nzchar(trimws(new_ids))])

# Clinical group from the sample ID, longest pattern first so that e.g.
# MS-Rituximab-Progressive is not captured by a shorter MS-Rituximab rule.
# "MS untyped" is the MS-Unknown block, which is MS with the subtype not yet
# annotated in the sheet -- not an unknown disease status.
group_of <- function(x) {
  out <- rep("Mix/other", length(x))
  pats <- c("Ctrl" = "Ctrl",
            "MS-Rituximab-Progressive" = "MS Rituximab-Progressive",
            "MS-Rituximab-Stable"      = "MS Rituximab-Stable",
            "MS-Natalizumab"           = "MS Natalizumab",
            "MS-New-PPMS"              = "MS New PPMS",
            "MS-New-RR"                = "MS New RR",
            "MS-Unknown"               = "MS untyped")
  for (p in names(pats)) out[out == "Mix/other" & grepl(p, x, fixed = TRUE)] <- pats[[p]]
  out
}

v4 <- dat %>%
  mutate(
    protocol = vapply(sample_id, function(s) {
      m <- regmatches(s, regexpr("V[0-9]|(?<=[^A-Za-z])[0-9][A-Z]", s, perl = TRUE))
      if (length(m)) m else NA_character_
    }, character(1)),
    group = group_of(sample_id),
    round = factor(ifelse(sample_id %in% new_ids, "Latest round", "Earlier V4"),
                   levels = c("Earlier V4", "Latest round")),
    mark  = factor(mark, levels = c("H3K27ac", "H3K4me3"))
  ) %>%
  filter(protocol == "V4")

stopifnot(nrow(v4) > 0)

n_tot   <- nrow(v4)
n_new   <- sum(v4$round == "Latest round")
n_old   <- sum(v4$round == "Earlier V4")
n_k27   <- sum(v4$mark == "H3K27ac")
n_k4    <- sum(v4$mark == "H3K4me3")
n_ms    <- sum(grepl("^MS", v4$group) & v4$mark == "H3K27ac")
n_ctrl  <- sum(v4$group == "Ctrl" & v4$mark == "H3K27ac")
# donors carrying both marks: the K4me3 libraries are suffixed on an otherwise
# identical ID, so stripping the suffix collapses the pair onto one donor
paired  <- v4 %>%
  mutate(donor = sub("-H3K4me3$", "", sample_id)) %>%
  count(donor) %>% filter(n > 1) %>% nrow()

# --- panel a: composition by clinical group ---------------------------------
comp <- v4 %>% count(group, round) %>%
  group_by(group) %>% mutate(tot = sum(n)) %>% ungroup() %>%
  mutate(group = reorder(group, tot))

p_comp <- ggplot(comp, aes(x = n, y = group, fill = round)) +
  geom_col(width = 0.68, colour = "grey25", linewidth = 0.3) +
  geom_text(aes(label = n), position = position_stack(vjust = 0.5),
            size = 3.9, fontface = "bold", colour = "white") +
  geom_text(data = distinct(comp, group, tot), aes(x = tot, y = group, label = tot),
            inherit.aes = FALSE, hjust = -0.45, size = 3.9, colour = "grey20") +
  scale_fill_manual(values = c("Earlier V4" = "#8C9BAB", "Latest round" = "#2C5F8A"),
                    name = NULL) +
  scale_x_continuous(expand = expansion(mult = c(0, 0.13))) +
  labs(title = "a  V4 libraries profiled, by clinical group",
       x = "Libraries", y = NULL) +
  theme_cowplot(13) +
  theme(legend.position = "bottom",
        # the longest group label ("MS Rituximab-Stable") is clipped otherwise
        plot.margin = margin(5.5, 5.5, 5.5, 22))

# --- panel b: mark breakdown per round --------------------------------------
mk <- v4 %>% count(round, mark)

p_mark <- ggplot(mk, aes(x = round, y = n, fill = mark)) +
  geom_col(width = 0.6, colour = "grey25", linewidth = 0.3) +
  geom_text(aes(label = n), position = position_stack(vjust = 0.5),
            size = 4.1, fontface = "bold", colour = "white") +
  scale_fill_manual(values = c("H3K27ac" = "#4C7A34", "H3K4me3" = "#C9A227"), name = NULL) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.1))) +
  labs(title = "b  Round and mark", x = NULL, y = "Libraries") +
  theme_cowplot(13) +
  theme(legend.position = "bottom")

# --- panel c: headline counts as text ---------------------------------------
lines <- c(
  sprintf("%d V4 libraries total", n_tot),
  sprintf("%d in the latest round, %d earlier", n_new, n_old),
  "",
  sprintf("%d H3K27ac  ·  %d H3K4me3", n_k27, n_k4),
  sprintf("%d MS  ·  %d control  (H3K27ac)", n_ms, n_ctrl),
  "",
  sprintf("%d donors profiled for both marks", paired),
  "All V4 libraries are plasma"
)
p_txt <- ggdraw() +
  draw_label("c  At a glance", x = 0.03, y = 0.995, hjust = 0, vjust = 1,
             size = 14, fontface = "bold") +
  draw_label(paste(lines, collapse = "\n"), x = 0.06, y = 0.46, hjust = 0, vjust = 0.5,
             size = 11.5, lineheight = 1.5, colour = "grey15")

cap <- paste0(
  "V4-protocol libraries only. \"Latest round\" = the ", length(new_ids),
  " sample IDs added to config/samplesheet.csv in commit fc67794 (2026-09-21), of which ", n_new, " are V4.\n",
  "\"MS untyped\" are MS patients whose subtype is not yet annotated in the samplesheet. ",
  "Counts are libraries, not donors: the ", paired, " dual-mark donors contribute two rows each."
)

fig <- plot_grid(p_comp, p_mark, p_txt, nrow = 1, rel_widths = c(1.6, 0.8, 1.05))
fig <- plot_grid(
  fig,
  ggdraw() + draw_label(cap, x = 0.008, hjust = 0, vjust = 0.5, size = 9.5,
                        colour = "grey30", lineheight = 1.15),
  ncol = 1, rel_heights = c(1, 0.12)
)

dir.create(dirname(out_pre), recursive = TRUE, showWarnings = FALSE)
ggsave(paste0(out_pre, ".pdf"), fig, width = 15.5, height = 5.2, device = "pdf")
ggsave(paste0(out_pre, ".png"), fig, width = 15.5, height = 5.2, dpi = 200)

# machine-readable version of the same table
tbl <- v4 %>% count(group, round, mark) %>% arrange(group, round, mark)
write.table(tbl, paste0(out_pre, ".tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

message("wrote ", out_pre, ".pdf / .png / .tsv")
print(as.data.frame(tbl))
