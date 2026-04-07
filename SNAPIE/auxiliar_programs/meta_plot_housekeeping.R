#!/usr/bin/env Rscript
# meta_plot_housekeeping.R
# Generates meta-coverage plots matching the cfChIP-seq Sadeh pipeline style:
#   - Left panel:  TSS meta-plot (Meta-genes.bed, groups: CpG.High/Low, NonCpG.High/Low, NotExpressed)
#   - Right panel: Enhancer meta-plot (Meta-enhancers.bed, groups: Many/Middle/none)
#   - Black heatmap background with red color, line plot above
#
# Usage:
#   Rscript meta_plot_housekeeping.R \
#     --bw            /path/to/sample.bw \
#     --meta_genes    /path/to/Meta-genes.bed \
#     --meta_enhancers /path/to/Meta-enhancers.bed \
#     --output        /path/to/sample_meta.pdf \
#     --sample        SAMPLE_NAME \
#     --metaplot_r    /path/to/MetaPlot.R   # optional, functions embedded if absent

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
  make_option("--bw",             type="character", help="Input bigWig file"),
  make_option("--meta_genes",     type="character", default=NULL,
              help="Meta-genes BED file (6-col with name groups: a.CpG.High, b.CpG.Low, ...)"),
  make_option("--meta_enhancers", type="character", default=NULL,
              help="Meta-enhancers BED file (6-col with name groups: a.Many, b.Middle, z.none)"),
  make_option("--output",         type="character", help="Output PDF file"),
  make_option("--sample",         type="character", default="sample", help="Sample name for plot title"),
  make_option("--metaplot_r",     type="character", default=NULL,
              help="Path to MetaPlot.R (functions embedded if not provided)"),
  make_option("--binsize",        type="integer",   default=25,
              help="Bin size in bp [default: 25]")
)
opt <- parse_args(OptionParser(option_list=option_list))

if (is.null(opt$bw) || is.null(opt$output)) {
  stop("--bw and --output are required.")
}
if (is.null(opt$meta_genes) && is.null(opt$meta_enhancers)) {
  stop("At least one of --meta_genes or --meta_enhancers is required.")
}

# ─── Load MetaPlot functions ──────────────────────────────────────────────────
if (!is.null(opt$metaplot_r) && file.exists(opt$metaplot_r)) {
  message("Sourcing MetaPlot functions from: ", opt$metaplot_r)
  source(opt$metaplot_r)
} else {
  message("Using embedded MetaPlot functions")

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
                                  xtickSpace=5, extraSpace=0, main="", ylim=NULL) {
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
    p <- p + scale_y_continuous(name="", expand=c(0, 0), limits=ylim)
    p <- p + labs(title=main)
    p <- p + theme(text=element_text(colour="black", size=8),
                   axis.text=element_text(colour="black", size=8))
    p
  }

  ggPlotHeatMap <- function(A, xlab="", ylab="", zlab="Coverage", title="",
                            zlim=NULL, offset=5000, label="TSS", Kbytes=100,
                            color="red", bg.color="black", xtickSpace=5) {
    m   <- do.call(rbind, A)
    n   <- ncol(m)
    rownames(m) <- seq_len(nrow(m))
    colnames(m) <- seq_len(n)
    ioff   <- Kbytes * offset / 1000
    breaks <- seq(0, n, xtickSpace * Kbytes) - ioff
    labels <- paste0(breaks / Kbytes, "Kb")
    labels[labels == "0Kb"] <- label
    ybreaks <- cumsum(sapply(A, nrow))
    ylabels <- names(A)
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

  PlotMeta <- function(Cov, PlotList, WindowSize=25, Norm=1) {
    # fixCoverage: pad chromosomes to seqlengths from the coverage object itself
    sl <- lengths(Cov)
    fixCoverage <- function(Cv) {
      for (chr in names(sl)) {
        l <- length(Cv[[chr]])
        m <- sl[[chr]]
        if (!is.na(m) && l < m)
          Cv[[chr]] <- append(Cv[[chr]], Rle(0, m - l))
      }
      Cv
    }
    Cov <- fixCoverage(Cov)
    message("Norm: ", Norm)

    p.list <- lapply(PlotList, function(pl) {
      RegionGroups <- split(pl$BED, pl$BED$name)
      RegionGroups <- RegionGroups[order(names(RegionGroups), decreasing=TRUE)]
      message("Groups: ", paste(names(RegionGroups), collapse=", "))
      names(RegionGroups) <- sub("^[[:lower:]]\\.", "", names(RegionGroups))

      Ms <- lapply(RegionGroups, function(g) {
        CollectMeta(Cov, g, pl$Width, pl$Offset, WindowSize, Norm)
      })
      Ls <- lapply(Ms, colMeans)

      if (!is.null(pl$Max) && as.numeric(pl$Max) > 0) {
        clim <- c(0, as.numeric(pl$Max))
      } else {
        Ls_max <- max(max(unlist(Ls)))
        if (Ls_max > 3) {
          clim <- c(0, 5 * ceiling(Ls_max / 5))
        } else {
          clim <- c(0, ceiling(Ls_max))
        }
      }
      message("clim: ", paste(clim, collapse=", "))

      Kb   <- 1000 / WindowSize
      p_meta <- ggPlotCovergeGroups(Ls, pos=pl$Offset, label=pl$Label,
                                    xtickSpace=pl$Tick / 1000,
                                    Kbytes=Kb, ylim=clim, extraSpace=0)
      p_heat <- ggPlotHeatMap(Ms, offset=pl$Offset, label=pl$Label,
                              Kbytes=Kb, ylab="", color=pl$Color,
                              bg.color=pl$BGColor, zlim=clim,
                              xtickSpace=pl$Tick / 1000)
      p_heat <- p_heat + guides(fill="none") +
        theme(axis.text.x=element_blank(),
              axis.ticks.x=element_blank(),
              axis.title.x=element_blank())
      list(p_meta, p_heat)
    })

    p.list.flat <- c(lapply(p.list, function(l) l[[1]]),
                     lapply(p.list, function(l) l[[2]]))
    plot_grid(plotlist=p.list.flat, ncol=length(PlotList), align="vh",
              rel_heights=rep(c(1, 2), length(PlotList)), axis="lrtb")
  }
}

# ─── Import bigWig as RleList coverage ────────────────────────────────────────
message("Importing bigWig: ", opt$bw)
bw_gr <- import(opt$bw, as="GRanges")
Cov   <- coverage(bw_gr, weight=bw_gr$score)

# ─── Build PlotList ───────────────────────────────────────────────────────────
mp_list <- list()

if (!is.null(opt$meta_genes) && file.exists(opt$meta_genes)) {
  message("Loading Meta-genes BED: ", opt$meta_genes)
  meta_genes <- import(opt$meta_genes)
  # Keep only chromosomes present in the bigWig
  meta_genes <- meta_genes[as.character(seqnames(meta_genes)) %in% names(Cov)]
  if (length(meta_genes) > 0) {
    mp_list[["Gene"]] <- list(
      BED     = meta_genes,
      Offset  = 5000,
      Width   = 25000,
      Tick    = 5000,
      Label   = "TSS",
      Max     = -1,
      Color   = "red",
      BGColor = "black"
    )
  } else {
    warning("No Meta-genes regions overlap with bigWig chromosomes — skipping TSS panel.")
  }
}

if (!is.null(opt$meta_enhancers) && file.exists(opt$meta_enhancers)) {
  message("Loading Meta-enhancers BED: ", opt$meta_enhancers)
  meta_enh <- import(opt$meta_enhancers)
  meta_enh <- meta_enh[as.character(seqnames(meta_enh)) %in% names(Cov)]
  if (length(meta_enh) > 0) {
    mp_list[["Enhancer"]] <- list(
      BED     = meta_enh,
      Offset  = 25000,
      Width   = 50000,
      Tick    = 10000,
      Label   = "Enhancer",
      Max     = -1,
      Color   = "red",
      BGColor = "black"
    )
  } else {
    warning("No Meta-enhancer regions overlap with bigWig chromosomes — skipping Enhancer panel.")
  }
}

if (length(mp_list) == 0) {
  warning("No valid regions found. Writing empty output.")
  dir.create(dirname(opt$output), showWarnings=FALSE, recursive=TRUE)
  pdf(opt$output, width=8, height=6)
  plot.new()
  title(paste0(opt$sample, " — no overlapping regions found"))
  dev.off()
  quit(status=0)
}

# ─── Plot and save ────────────────────────────────────────────────────────────
message(sprintf("Generating meta-plot for %s (%d panel(s))", opt$sample, length(mp_list)))
dir.create(dirname(opt$output), showWarnings=FALSE, recursive=TRUE)

p <- PlotMeta(Cov=Cov, PlotList=mp_list, WindowSize=opt$binsize, Norm=1)
p <- p + ggtitle(opt$sample)

pdf(opt$output, width=8 * length(mp_list), height=11)
print(p)
dev.off()

message("Done: ", opt$output)
