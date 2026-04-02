#!/usr/bin/env Rscript
# meta_plot_housekeeping.R
# Generates meta-coverage plots around housekeeping gene TSSs from a bigWig file.
# Uses MetaPlot.R functions (CollectMeta, ggPlotCovergeGroups, ggPlotHeatMap).
#
# Usage:
#   Rscript meta_plot_housekeeping.R \
#     --bw        /path/to/sample.bw \
#     --regions   /path/to/housekeeping_genes.bed \
#     --output    /path/to/sample_housekeeping_meta.pdf \
#     --sample    SAMPLE_NAME \
#     --metaplot  /path/to/MetaPlot.R \
#     --window    10000 \
#     --binsize   25

suppressPackageStartupMessages({
  library(optparse)
  library(rtracklayer)
  library(GenomicRanges)
  library(ggplot2)
  library(cowplot)
  library(reshape2)
})

# ─── Argument parsing ─────────────────────────────────────────────────────────
option_list <- list(
  make_option("--bw",       type="character", help="Input bigWig file"),
  make_option("--regions",  type="character", help="Housekeeping genes BED file (4-col: chr start end name)"),
  make_option("--output",   type="character", help="Output PDF file"),
  make_option("--sample",   type="character", default="sample", help="Sample name for plot title"),
  make_option("--metaplot", type="character", default=NULL,
              help="Path to MetaPlot.R (if not provided, functions are defined internally)"),
  make_option("--window",   type="integer",   default=10000,
              help="Total window size in bp around TSS [default: 10000 = ±5kb]"),
  make_option("--binsize",  type="integer",   default=25,
              help="Bin size in bp for coverage aggregation [default: 25]"),
  make_option("--color",    type="character", default="steelblue",
              help="Heatmap color [default: steelblue]")
)
opt <- parse_args(OptionParser(option_list=option_list))

if (is.null(opt$bw) || is.null(opt$regions) || is.null(opt$output)) {
  stop("--bw, --regions, and --output are required.")
}

# ─── Load MetaPlot functions ──────────────────────────────────────────────────
if (!is.null(opt$metaplot) && file.exists(opt$metaplot)) {
  source(opt$metaplot)
} else {
  # Embed CollectMeta and plotting functions (from MetaPlot.R)
  CollectMeta <- function(Cov, Regions, Width, Offset, WindowSize=10, Norm=1) {
    Regions <- resize(Regions, width=Offset + 0.5 * width(Regions), fix="end")
    Regions <- resize(Regions, width=Width, fix="start")
    ws <- IRanges(start=seq(1, Width - 1, by=WindowSize), width=WindowSize)
    n  <- length(Regions)
    m  <- ceiling(Width / WindowSize)
    X  <- matrix(0, nr=n, nc=m, dimnames=list(Regions$name, start(ws)))
    I  <- as.logical(strand(Regions) == "-")
    X[,] <- t(sapply(Cov[Regions], function(x) aggregate(x, ws, sum)))
    if (any(I))
      X[I, 1:m] <- X[I, m:1]
    X <- X * Norm / WindowSize
    X <- X[order(rowSums(X)), ]
    X
  }

  ggPlotCovergeGroups <- function(A, pos=5000, label="TSS", Kbytes=100,
                                  xtickSpace=5, extraSpace=4, main="", ylim=NULL) {
    if (is.null(ylim)) {
      M <- max(0, max(sapply(A, max), na.rm=TRUE) * 1.1)
      m <- min(0, min(sapply(A, min), na.rm=TRUE) * 1.1)
      ylim <- c(m, M)
    }
    w <- length(A[[1]])
    p <- ggplot()
    p <- p + theme(legend.position=c(0.95, 0.95),
                   legend.justification=c("right", "top"))
    p <- p + scale_color_brewer(palette="Set1",
                                guide=guide_legend(title=NULL, reverse=FALSE))
    p <- p + geom_vline(xintercept=0, colour="gray", linewidth=1)
    for (i in seq_along(A)) {
      p <- p + geom_line(
        data=data.frame(x=(1:w) - Kbytes * pos / 1000,
                        y=A[[i]],
                        color=rep(names(A)[[i]], w)),
        aes(x=x, y=y, colour=color), linewidth=0.5, show.legend=TRUE)
    }
    ipos   <- as.integer(pos / (xtickSpace * 1000)) * Kbytes * xtickSpace
    iwidth <- as.integer(w / Kbytes) * Kbytes
    labelat  <- seq(-ipos, iwidth, xtickSpace * Kbytes)
    labelval <- paste0(labelat / Kbytes, "Kb")
    labelval[labelval == "0Kb"] <- label
    xlim <- c(-Kbytes * pos / 1000, iwidth - Kbytes * (pos / 1000 + extraSpace))
    p <- p + scale_x_continuous(name="", breaks=labelat, labels=labelval,
                                 limits=xlim, expand=c(0, 0))
    p <- p + scale_y_continuous(name="Coverage", expand=c(0, 0), limits=ylim)
    p <- p + labs(title=main)
    p <- p + theme(text=element_text(colour="black", size=8),
                   axis.text=element_text(colour="black", size=8))
    p
  }

  ggPlotHeatMap <- function(A, xlab="", ylab="Gene", zlab="Coverage", title="",
                            zlim=NULL, offset=5000, label="TSS", Kbytes=100,
                            color="red", bg.color="white", xtickSpace=5) {
    m   <- do.call(rbind, A)
    n   <- ncol(m)
    rownames(m) <- seq_len(nrow(m))
    colnames(m) <- seq_len(n)
    ioff   <- Kbytes * offset / 1000
    breaks <- seq(0, n, xtickSpace * Kbytes) - ioff
    labels <- paste0(breaks / Kbytes, "Kb")
    labels[labels == "0Kb"] <- label
    ybreaks <- cumsum(sapply(A, nrow))
    df      <- melt(m, varnames=c("y", "x"))
    df$x    <- df$x - ioff
    t1 <- quantile(m, probs=0.99, na.rm=TRUE)
    t2 <- quantile(m, probs=0.05, na.rm=TRUE)
    if (t2 > 0) t2 <- 0
    if (is.null(zlim)) zlim <- c(t2, t1)
    df$value[df$value < zlim[[1]]] <- zlim[[1]]
    df$value[df$value > zlim[[2]]] <- zlim[[2]]
    p <- ggplot(df, aes(x=x, y=y, fill=value))
    p <- p + geom_raster()
    p <- p + geom_vline(xintercept=0, colour="gray", linewidth=1)
    for (h in ybreaks[seq_len(length(ybreaks) - 1)])
      p <- p + geom_hline(yintercept=h, colour="gray", linewidth=1)
    p <- p + scale_fill_gradient(low=bg.color, high=color,
                                  na.value="black", limits=zlim)
    p <- p + scale_x_continuous(breaks=breaks, labels=labels,
                                 limits=c(-ioff, n - ioff), expand=c(0, 0))
    p <- p + labs(x=xlab, y=ylab, fill=zlab, title=title)
    p <- p + theme_minimal()
    p <- p + theme(text=element_text(colour="black", size=8),
                   axis.text=element_text(colour="black", size=8))
    p
  }
}

# ─── Import bigWig ─────────────────────────────────────────────────────────────
message("Importing bigWig: ", opt$bw)
bw_gr  <- import(opt$bw, as="GRanges")
# Build RleList coverage from the bigWig signal
Cov    <- coverage(bw_gr, weight=bw_gr$score)

# ─── Load housekeeping gene regions → TSS GRanges ─────────────────────────────
message("Loading housekeeping regions: ", opt$regions)
bed    <- read.table(opt$regions, header=FALSE, sep="\t", stringsAsFactors=FALSE,
                     col.names=c("chr", "start", "end", "name"))

# Use gene start as TSS (1-based), single-bp anchor
tss_gr <- GRanges(
  seqnames = bed$chr,
  ranges   = IRanges(start=bed$start + 1L, width=1L),   # BED is 0-based
  strand   = "*",
  name     = bed$name
)
# Keep only chromosomes present in the bigWig
tss_gr <- tss_gr[as.character(seqnames(tss_gr)) %in% names(Cov)]

if (length(tss_gr) == 0) {
  warning("No housekeeping regions overlap with bigWig chromosomes. Writing empty output.")
  pdf(opt$output, width=8, height=6)
  plot.new()
  title(paste0(opt$sample, " — no overlapping regions found"))
  dev.off()
  quit(status=0)
}

# ─── Compute coverage matrix ───────────────────────────────────────────────────
Width    <- opt$window      # total bp
Offset   <- Width / 2       # TSS at window centre
BinSize  <- opt$binsize

message(sprintf("Computing coverage matrix (%d regions, ±%d bp, %d bp bins)",
                length(tss_gr), Offset, BinSize))

# Group all housekeeping genes as one group named by sample
region_list <- list()
region_list[[opt$sample]] <- tss_gr

Ms <- lapply(region_list, function(g) {
  CollectMeta(Cov, g, Width=Width, Offset=Offset, WindowSize=BinSize, Norm=1)
})

# Average profiles
Ls <- lapply(Ms, colMeans)

# Colour limits
Ls_max <- max(unlist(Ls))
if (Ls_max > 3) {
  clim <- c(0, 5 * ceiling(Ls_max / 5))
} else {
  clim <- c(0, max(ceiling(Ls_max), 1))
}

Kbytes   <- 1000 / BinSize   # windows per kb
TickStep <- 2000             # x-axis tick every 2 kb

# ─── Build plots ──────────────────────────────────────────────────────────────
p_meta <- ggPlotCovergeGroups(
  Ls,
  pos       = Offset,
  label     = "TSS",
  Kbytes    = Kbytes,
  xtickSpace= TickStep / 1000,
  extraSpace= 0,
  main      = paste0(opt$sample, " — Housekeeping gene coverage"),
  ylim      = clim
)

p_heat <- ggPlotHeatMap(
  Ms,
  offset    = Offset,
  label     = "TSS",
  Kbytes    = Kbytes,
  ylab      = "Genes (ranked by coverage)",
  zlab      = "Coverage",
  title     = "",
  color     = opt$color,
  bg.color  = "white",
  zlim      = clim,
  xtickSpace= TickStep / 1000
)
p_heat <- p_heat +
  guides(fill="none") +
  theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.title.x=element_blank())

# ─── Combine and save ─────────────────────────────────────────────────────────
message("Saving plot: ", opt$output)
dir.create(dirname(opt$output), showWarnings=FALSE, recursive=TRUE)

p_combined <- plot_grid(p_heat, p_meta, ncol=1, align="v",
                        rel_heights=c(2, 1), axis="lrtb")

pdf(opt$output, width=8, height=10)
print(p_combined)
dev.off()

message("Done.")
