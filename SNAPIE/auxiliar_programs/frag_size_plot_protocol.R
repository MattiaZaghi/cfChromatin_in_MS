#!/usr/bin/env Rscript
# Fragment size distribution — group-level summary plot (mean ± SD per condition)
# Usage: Rscript frag_size_plot.R <frags_dir> <output.pdf> [samplesheet.csv]

suppressMessages({
  library(ggplot2)
  library(dplyr)
})

args        <- commandArgs(trailingOnly = TRUE)
frags_dir   <- args[1]
out_pdf     <- args[2]
samplesheet <- if (length(args) >= 3 && file.exists(args[3])) args[3] else NULL

# ── Helpers ───────────────────────────────────────────────────────────────────

infer_group <- function(s) {
  dplyr::case_when(
    grepl("V1",s) ~ "V1",
    grepl("V2",s) ~ "V2",
    grepl("V3",s) ~ "V3",
    grepl("V4",s) ~ "V4"
  )
}

read_frag_file <- function(f) {
  sample_name <- basename(dirname(f))
  dat <- tryCatch(
    read.table(f, skip = 2, sep = "\t", header = FALSE,
               col.names = c("Size", "Occurrences", "BAM"),
               stringsAsFactors = FALSE),
    error = function(e) NULL
  )
  if (is.null(dat) || nrow(dat) == 0) return(NULL)
  dat$sample <- sample_name
  dat[, c("sample", "Size", "Occurrences")]
}

# ── Load data ─────────────────────────────────────────────────────────────────

frag_files <- list.files(frags_dir, pattern = "fragment_sizes\\.txt$",
                         recursive = TRUE, full.names = TRUE)

if (length(frag_files) == 0) {
  pdf(out_pdf, width = 8, height = 5)
  plot.new(); title("No fragment size data found")
  dev.off()
  quit(save = "no")
}

all_data <- dplyr::bind_rows(lapply(frag_files, read_frag_file))

if (nrow(all_data) == 0) {
  pdf(out_pdf, width = 8, height = 5)
  plot.new(); title("Fragment size files are empty")
  dev.off()
  quit(save = "no")
}

# ── Annotate groups ───────────────────────────────────────────────────────────

all_data$group <- infer_group(all_data$sample)

if (!is.null(samplesheet)) {
  ss <- tryCatch(read.csv(samplesheet, stringsAsFactors = FALSE), error = function(e) NULL)
  if (!is.null(ss) && all(c("sampleId", "group") %in% colnames(ss))) {
    all_data <- dplyr::left_join(all_data, ss[, c("sampleId", "group")],
                                 by = c("sample" = "sampleId"), suffix = c("_inferred", ""))
    all_data$group <- dplyr::coalesce(all_data$group, all_data$group_inferred)
    all_data$group_inferred <- NULL
  }
}

# ── Normalise and filter ──────────────────────────────────────────────────────

all_data <- all_data %>%
  filter(Size >= 20, Size <= 800) %>%
  group_by(sample) %>%
  mutate(Frequency = Occurrences / sum(Occurrences)) %>%
  ungroup()

# ── Group summary (mean ± SD) ─────────────────────────────────────────────────

group_summary <- all_data %>%
  group_by(group, Size) %>%
  summarise(mean_freq = mean(Frequency),
            sd_freq   = sd(Frequency),
            n         = n(),
            .groups   = "drop")

# ── Colours ───────────────────────────────────────────────────────────────────

group_colors <- c(
  "V1"                     = "#2196F3",
  "V2"                      = "#F44336",
  "V3"                       = "#4CAF50",
  "V4" =                      "#FF9800"
)
present <- unique(all_data$group)
group_colors <- group_colors[names(group_colors) %in% present]

# ── Nucleosomal markers ───────────────────────────────────────────────────────

nucl_lines <- data.frame(
  xint  = c(147, 294, 441),
  label = c("Mono\n147 bp", "Di\n294 bp", "Tri\n441 bp")
)
y_top <- max(group_summary$mean_freq, na.rm = TRUE)

# ── Plot: group mean ± SD ribbons ─────────────────────────────────────────────

p <- ggplot(group_summary, aes(x = Size, color = group, fill = group)) +
  geom_ribbon(aes(ymin = pmax(0, mean_freq - sd_freq),
                  ymax = mean_freq + sd_freq),
              alpha = 0.15, color = NA) +
  geom_line(aes(y = mean_freq), linewidth = 0.9) +
  geom_vline(xintercept = nucl_lines$xint,
             linetype = "dashed", color = "grey50", linewidth = 0.4) +
  annotate("text", x = nucl_lines$xint + 4, y = y_top * 0.97,
           label = nucl_lines$label,
           hjust = 0, vjust = 1, size = 2.8, color = "grey40") +
  scale_color_manual(values = group_colors, name = "Group") +
  scale_fill_manual(values  = group_colors, name = "Group") +
  scale_x_continuous(breaks = c(50, 147, 200, 294, 400, 441, 600, 800),
                     minor_breaks = NULL) +
  scale_y_continuous(limits = c(0, 0.03)) +
  labs(title    = "Fragment length distribution by condition",
       subtitle = sprintf("Mean \u00b1 1 SD — %d samples across %d groups",
                          length(unique(all_data$sample)),
                          length(unique(all_data$group))),
       x = "Fragment size (bp)", y = "Normalized frequency") +
  theme_bw(base_size = 12) +
  theme(panel.grid.minor = element_blank(),
        axis.text.x      = element_text(angle = 45, hjust = 1),
        plot.title       = element_text(face = "bold"),
        legend.position  = "right")

pdf(out_pdf, width = 10, height = 6)
print(p)
dev.off()

message("Group fragment size plot written to: ", out_pdf)
