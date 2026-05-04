# Module 7b — RUVg batch factor estimation
# Snakemake R script: snakemake object is injected automatically.
#
# Input:  anchor_counts.tsv   (regions × samples, raw integers)
#         sample_metadata.tsv
# Output: ruv_factors.tsv     (samples × W_1 … W_k)
#         rle_before_after.pdf

suppressPackageStartupMessages({
    library(RUVSeq)
    library(EDASeq)
    library(ggplot2)
    library(reshape2)
})

# ── Load inputs ───────────────────────────────────────────────────────────────
anchor_mat <- read.table(
    snakemake@input$anchor_counts,
    header = TRUE, row.names = 1, sep = "\t", check.names = FALSE
)
anchor_mat <- as.matrix(anchor_mat)
mode(anchor_mat) <- "integer"

meta <- read.table(
    snakemake@input$meta,
    header = TRUE, sep = "\t", comment.char = "#", stringsAsFactors = FALSE
)
meta <- meta[meta$qc_include == "TRUE", ]

# Align columns to metadata order
common_samples <- intersect(colnames(anchor_mat), meta$sample_id)
anchor_mat     <- anchor_mat[, common_samples, drop = FALSE]
meta           <- meta[meta$sample_id %in% common_samples, ]
meta           <- meta[match(common_samples, meta$sample_id), ]

k <- as.integer(snakemake@params$k)
cat(sprintf("[Module 7b] RUVg: %d samples, %d anchor regions, k=%d\n",
            ncol(anchor_mat), nrow(anchor_mat), k))

# ── RLE before correction ─────────────────────────────────────────────────────
rle_matrix <- function(mat, title) {
    log_mat   <- log1p(mat)
    row_meds  <- apply(log_mat, 1, median)
    rle_vals  <- sweep(log_mat, 1, row_meds, "-")
    df <- reshape2::melt(rle_vals, varnames = c("region", "sample"),
                         value.name = "RLE")
    df$group <- meta$group[match(df$sample, meta$sample_id)]
    ggplot(df, aes(x = sample, y = RLE, fill = group)) +
        geom_boxplot(outlier.size = 0.3, lwd = 0.3) +
        geom_hline(yintercept = 0, linetype = "dashed", colour = "grey40") +
        theme_bw(base_size = 8) +
        theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 4)) +
        labs(title = title, x = NULL, y = "RLE") +
        scale_fill_brewer(palette = "Set2")
}

# ── Run RUVg ─────────────────────────────────────────────────────────────────
# All anchor regions are negative controls (biologically invariant)
set   <- EDASeq::newSeqExpressionSet(
    counts = anchor_mat,
    phenoData = data.frame(
        row.names  = common_samples,
        group      = meta$group
    )
)
ruv_fit <- RUVSeq::RUVg(set, cIdx = rownames(anchor_mat), k = k)

W <- as.data.frame(pData(ruv_fit))
W <- W[, grepl("^W_", colnames(W)), drop = FALSE]
W$sample_id <- rownames(W)
W <- W[, c("sample_id", grep("^W_", colnames(W), value = TRUE))]

# ── RLE after correction ──────────────────────────────────────────────────────
# Load any RRE count matrix for visualisation if available
rre_path <- snakemake@params$rre_counts
if (file.exists(rre_path)) {
    rre_mat <- as.matrix(read.table(rre_path, header = TRUE, row.names = 1,
                                    sep = "\t", check.names = FALSE))
    rre_mat <- rre_mat[, common_samples, drop = FALSE]
    mode(rre_mat) <- "integer"
    p_before <- rle_matrix(rre_mat, "RLE before RUVg (RRE full)")
    # Corrected matrix: divide by RUVg expected (approx via anchor fit)
    ruv_rre <- EDASeq::newSeqExpressionSet(
        counts = rre_mat,
        phenoData = data.frame(row.names = common_samples, group = meta$group)
    )
    ruv_rre_fit <- RUVSeq::RUVg(ruv_rre, cIdx = rownames(rre_mat)[seq_len(min(200, nrow(rre_mat)))], k = k)
    p_after  <- rle_matrix(
        normCounts(ruv_rre_fit), "RLE after RUVg correction (RRE full)"
    )
} else {
    p_before <- rle_matrix(anchor_mat, "RLE before RUVg (anchor regions)")
    p_after  <- ggplot() + labs(title = "RLE after — RRE matrix not yet available")
}

# ── Save outputs ──────────────────────────────────────────────────────────────
dir.create(dirname(snakemake@output$factors), recursive = TRUE, showWarnings = FALSE)
write.table(W, snakemake@output$factors, sep = "\t", quote = FALSE, row.names = FALSE)

pdf(snakemake@output$rle_plot, width = 14, height = 5)
print(p_before)
print(p_after)
dev.off()

cat(sprintf("[Module 7b] W matrix saved: %d samples × %d factors\n",
            nrow(W), k))
