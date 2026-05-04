#!/usr/bin/env Rscript
# meta_plot_hk_promoters.R
#
# Per-sample meta-coverage plot at housekeeping gene promoters.
# Window: -5 kb to +20 kb around the TSS (TSS = centre of input BED regions).
# Top panel : mean coverage — single red line, no ribbon.
# Bottom panel: heatmap sorted by total signal (black → red).
# Matches the Sadeh cfChIP-seq pipeline aesthetic.
#
# Args:
#   --bw       bigWig file (normalised, e.g. RPKM deeptools output)
#   --bed      BED file with ~2 kb promoter windows (centre used as TSS)
#   --output   output PDF path
#   --sample   sample name for plot title
#   --binsize  bin size in bp (default 25)

suppressPackageStartupMessages({
  library(optparse)
  library(rtracklayer)
  library(GenomicRanges)
  library(IRanges)
  library(ggplot2)
  library(cowplot)
})

option_list <- list(
  make_option("--bw",      type = "character", help = "Input bigWig file"),
  make_option("--bed",     type = "character", help = "Housekeeping promoter BED"),
  make_option("--output",  type = "character", help = "Output PDF path"),
  make_option("--sample",  type = "character", default = "sample"),
  make_option("--binsize", type = "integer",   default = 25L)
)
opt <- parse_args(OptionParser(option_list = option_list))

if (is.null(opt$bw) || is.null(opt$bed) || is.null(opt$output))
  stop("--bw, --bed, and --output are required.")

# ── Import bigWig ─────────────────────────────────────────────────────────────
message("Importing bigWig: ", opt$bw)
bw_gr <- import(opt$bw)
Cov   <- coverage(bw_gr, weight = bw_gr$score)

# ── Import BED and derive TSS anchor points ───────────────────────────────────
message("Loading promoter BED: ", opt$bed)
bed <- import(opt$bed)
# TSS = centre of each BED window (BED regions are ~2 kb pre-centred on TSS)
tss <- resize(bed, width = 1L, fix = "center")

# ── Filter to chromosomes present in the bigWig coverage ─────────────────────
cov_chroms     <- names(Cov)
tss            <- tss[as.character(seqnames(tss)) %in% cov_chroms]
seqlevels(tss) <- as.character(unique(seqnames(tss)))
message(sprintf("TSS sites on covered chromosomes: %d", length(tss)))

if (length(tss) == 0) {
  dir.create(dirname(opt$output), recursive = TRUE, showWarnings = FALSE)
  pdf(opt$output, width = 8, height = 10); plot.new()
  title(paste0(opt$sample, " - no TSS on covered chromosomes")); dev.off()
  quit(status = 0)
}

# ── Create -5 kb / +20 kb windows ─────────────────────────────────────────────
upstream   <- 5000L
downstream <- 20000L
win_width  <- upstream + downstream    # 25,000 bp

windows <- GRanges(
  seqnames = seqnames(tss),
  ranges   = IRanges(start = start(tss) - upstream,
                     end   = start(tss) + downstream - 1L),
  strand   = "*"
)
seqlevels(windows) <- seqlevels(tss)

# Drop windows that extend past chromosome boundaries
chrom_lens <- lengths(Cov)
win_chroms <- as.character(seqnames(windows))
in_bounds  <- (start(windows) >= 1L) &
              (end(windows)   <= chrom_lens[win_chroms])
n_dropped  <- sum(!in_bounds)
if (n_dropped > 0)
  message(sprintf("Dropping %d windows that extend past chromosome end", n_dropped))
windows <- windows[in_bounds]

if (length(windows) == 0) {
  dir.create(dirname(opt$output), recursive = TRUE, showWarnings = FALSE)
  pdf(opt$output, width = 8, height = 10); plot.new()
  title(paste0(opt$sample, " - no valid windows")); dev.off()
  quit(status = 0)
}
message(sprintf("Computing coverage for %d windows x %d bp ...",
                length(windows), win_width))

# ── Build coverage matrix (windows x bins) ────────────────────────────────────
n_bins     <- floor(win_width / opt$binsize)
bin_ranges <- IRanges(start = seq(1L, by = opt$binsize, length.out = n_bins),
                      width = opt$binsize)

sub_covs <- Cov[windows]   # RleList: one element per window

mat <- do.call(rbind, lapply(sub_covs, function(x) {
  if (length(x) < win_width) return(rep(NA_real_, n_bins))
  as.numeric(aggregate(x[seq_len(win_width)], bin_ranges, sum)) / opt$binsize
}))

# Position axis: 0 = TSS; negative = upstream, positive = downstream
pos <- (seq_len(n_bins) - 0.5) * opt$binsize - upstream

# ── X-axis labels ─────────────────────────────────────────────────────────────
x_breaks <- c(-5000, 0, 5000, 10000, 15000, 20000)
x_labels <- c("-5 kb", "TSS", "+5 kb", "+10 kb", "+15 kb", "+20 kb")
x_lim    <- c(min(pos) - opt$binsize / 2, max(pos) + opt$binsize / 2)

# ── Profile panel (top) — mean line only ──────────────────────────────────────
mn      <- colMeans(mat, na.rm = TRUE)
prof_df <- data.frame(pos = pos, mean = mn)

p_prof <- ggplot(prof_df, aes(x = pos, y = mean)) +
  geom_line(color = "#cc2222", linewidth = 0.9) +
  geom_vline(xintercept = 0, color = "grey60",
             linetype = "dashed", linewidth = 0.5) +
  scale_x_continuous(breaks = x_breaks, labels = x_labels,
                     expand = c(0, 0), limits = x_lim) +
  labs(
    title = sprintf("%s  |  Housekeeping promoters  (n = %d)", opt$sample, nrow(mat)),
    y     = "Mean coverage (RPKM/bp)",
    x     = NULL
  ) +
  theme_bw(base_size = 10) +
  theme(
    plot.title       = element_text(face = "bold", size = 9),
    panel.grid.minor = element_blank(),
    axis.text.x      = element_blank(),
    axis.ticks.x     = element_blank(),
    axis.title.y     = element_text(size = 9)
  )

# ── Heatmap panel (bottom) ────────────────────────────────────────────────────
row_order  <- order(rowSums(mat, na.rm = TRUE), decreasing = TRUE)
mat_sorted <- mat[row_order, , drop = FALSE]

clip_hi            <- quantile(mat, 0.99, na.rm = TRUE)
mat_clipped        <- mat_sorted
mat_clipped[is.na(mat_clipped)]    <- 0
mat_clipped[mat_clipped > clip_hi] <- clip_hi
mat_clipped[mat_clipped < 0]       <- 0

heat_df <- data.frame(
  region = rep(seq_len(nrow(mat_clipped)), times = n_bins),
  pos    = rep(pos, each = nrow(mat_clipped)),
  value  = as.vector(mat_clipped)
)

p_heat <- ggplot(heat_df, aes(x = pos, y = region, fill = value)) +
  geom_raster(interpolate = FALSE) +
  geom_vline(xintercept = 0, color = "grey75", linewidth = 0.35) +
  scale_fill_gradient(low = "black", high = "#cc2222",
                      limits = c(0, clip_hi), na.value = "black",
                      name = "Coverage") +
  scale_x_continuous(breaks = x_breaks, labels = x_labels,
                     expand = c(0, 0), limits = x_lim) +
  scale_y_reverse(expand = c(0, 0)) +
  labs(x = "Position relative to TSS", y = "Promoters (sorted by signal)") +
  theme_minimal(base_size = 10) +
  theme(
    panel.background = element_rect(fill = "black", color = NA),
    plot.background  = element_rect(fill = "white", color = NA),
    panel.grid       = element_blank(),
    axis.text.y      = element_blank(),
    axis.ticks.y     = element_blank(),
    axis.title       = element_text(size = 9),
    legend.title     = element_text(size = 8),
    legend.text      = element_text(size = 8)
  )

# ── Combine and save ──────────────────────────────────────────────────────────
dir.create(dirname(opt$output), recursive = TRUE, showWarnings = FALSE)
pdf(opt$output, width = 8, height = 10)
print(plot_grid(p_prof, p_heat,
                ncol = 1, align = "v", axis = "lr",
                rel_heights = c(1, 3)))
dev.off()
message("Done: ", opt$output)
