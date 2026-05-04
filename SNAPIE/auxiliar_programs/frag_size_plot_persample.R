#!/usr/bin/env Rscript
# Fragment size distribution — single sample plot
# Usage: Rscript frag_size_plot_persample.R <fragment_sizes.txt> <output.pdf> <sample_name>

suppressMessages({
  library(ggplot2)
})

args        <- commandArgs(trailingOnly = TRUE)
frag_file   <- args[1]
out_pdf     <- args[2]
sample_name <- args[3]

# ── Check input ───────────────────────────────────────────────────────────────

info <- file.info(frag_file)
if (is.na(info$size) || info$size == 0) {
  pdf(out_pdf, width = 7, height = 5)
  plot.new()
  title(main = sample_name, sub = "No fragment size data available")
  dev.off()
  quit(save = "no")
}

# ── Load data ─────────────────────────────────────────────────────────────────

dat <- tryCatch(
  read.table(frag_file, skip = 2, sep = "\t", header = FALSE,
             col.names = c("Size", "Occurrences", "BAM"),
             stringsAsFactors = FALSE),
  error = function(e) NULL
)

if (is.null(dat) || nrow(dat) == 0) {
  pdf(out_pdf, width = 7, height = 5)
  plot.new()
  title(main = sample_name, sub = "Could not parse fragment size file")
  dev.off()
  quit(save = "no")
}

dat <- dat[dat$Size >= 20 & dat$Size <= 800, ]
dat$Frequency <- dat$Occurrences / sum(dat$Occurrences)

# Weighted median
cum_freq    <- cumsum(dat$Frequency)
median_size <- dat$Size[which(cum_freq >= 0.5)[1]]

y_top <- max(dat$Frequency, na.rm = TRUE)

# ── Nucleosomal markers ───────────────────────────────────────────────────────

nucl_lines <- data.frame(
  xint  = c(147, 294, 441),
  label = c("Mono\n147 bp", "Di\n294 bp", "Tri\n441 bp")
)

# ── Plot ──────────────────────────────────────────────────────────────────────

p <- ggplot(dat, aes(x = Size, y = Frequency)) +
  geom_area(fill = "#2196F3", alpha = 0.15) +
  geom_line(color = "#2196F3", linewidth = 0.8) +
  geom_vline(xintercept = nucl_lines$xint,
             linetype = "dashed", color = "grey50", linewidth = 0.4) +
  annotate("text", x = nucl_lines$xint + 4, y = y_top * 0.97,
           label = nucl_lines$label,
           hjust = 0, vjust = 1, size = 3, color = "grey40") +
  geom_vline(xintercept = median_size,
             color = "#E53935", linewidth = 0.7, linetype = "dotted") +
  annotate("text", x = median_size + 4, y = y_top * 0.82,
           label = sprintf("median = %d bp", median_size),
           hjust = 0, size = 3.2, color = "#E53935") +
  scale_x_continuous(breaks = c(50, 147, 200, 294, 400, 441, 600, 800),
                     minor_breaks = NULL) +
  labs(title    = sample_name,
       subtitle = "Fragment length distribution",
       x = "Fragment size (bp)", y = "Normalized frequency") +
  theme_bw(base_size = 11) +
  theme(panel.grid.minor = element_blank(),
        axis.text.x      = element_text(angle = 45, hjust = 1),
        plot.title       = element_text(face = "bold"))

pdf(out_pdf, width = 7, height = 5)
print(p)
dev.off()

message("Fragment size plot written to: ", out_pdf)
