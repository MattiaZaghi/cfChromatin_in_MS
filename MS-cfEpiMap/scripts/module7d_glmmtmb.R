# Module 7d — glmmTMB ZINB sensitivity analysis (bins only).
# Snakemake R script: snakemake object is injected automatically.
#
# Runs a zero-inflated negative binomial model on the bin count matrix
# for each contrast and flags features significant in both DESeq2 and ZINB.

suppressPackageStartupMessages({
    library(glmmTMB)
    library(data.table)
    library(dplyr)
})

# ── Load inputs ───────────────────────────────────────────────────────────────
bin_counts <- as.matrix(
    read.table(snakemake@input$bin_counts, header = TRUE,
               row.names = 1, sep = "\t", check.names = FALSE)
)

meta <- read.table(
    snakemake@input$meta, header = TRUE, sep = "\t",
    comment.char = "#", stringsAsFactors = FALSE
)
meta <- meta[meta$qc_include == "TRUE", ]

ruv <- read.table(snakemake@input$ruv, header = TRUE, sep = "\t",
                  stringsAsFactors = FALSE)
bci <- read.table(snakemake@input$bci, header = TRUE, sep = "\t",
                  stringsAsFactors = FALSE)

meta <- merge(meta, ruv, by = "sample_id", all.x = TRUE)
meta <- merge(meta, bci[, c("sample_id", "bci_scaled")],
              by = "sample_id", all.x = TRUE)

meta$rituximab_treated <- as.integer(
    meta$group %in% c("MS-Rituximab-Stable", "MS-Rituximab-Progressive")
)
meta$group <- factor(
    meta$group,
    levels = c("Ctrl", "NEW", "MS-Rituximab-Stable", "MS-Rituximab-Progressive")
)
meta$total_frags <- colSums(bin_counts)[meta$sample_id]

# Align
common <- intersect(colnames(bin_counts), meta$sample_id)
bin_counts <- bin_counts[, common]
meta       <- meta[meta$sample_id %in% common, ]
meta       <- meta[match(common, meta$sample_id), ]

w_cols     <- grep("^W_", colnames(meta), value = TRUE)
n_factors  <- length(w_cols)

fdr_thr    <- snakemake@params$zinb_confirmed_fdr
contrasts  <- snakemake@params$contrasts
family     <- snakemake@params$family

cat(sprintf("[Module 7d] %d bins × %d samples\n",
            nrow(bin_counts), ncol(bin_counts)))

# ── Helper: fit one ZINB model for one bin ───────────────────────────────────
fit_zinb <- function(y, meta_df) {
    df <- data.frame(count = as.integer(y), meta_df,
                     stringsAsFactors = FALSE)
    fixed_parts <- c(w_cols, "rituximab_treated", "bci_scaled", "group")
    fixed_parts <- fixed_parts[fixed_parts %in% colnames(df)]
    formula_str <- paste("count ~", paste(fixed_parts, collapse = " + "),
                         "+ offset(log(total_frags + 1))")
    tryCatch({
        fit <- glmmTMB(
            as.formula(formula_str),
            data      = df,
            family    = nbinom2,
            ziformula = ~ 1,
        )
        coefs <- summary(fit)$coefficients$cond
        grp_rows <- rownames(coefs)[grepl("^group", rownames(coefs))]
        list(coefs = coefs[grp_rows, , drop = FALSE], converged = TRUE)
    }, error = function(e) {
        list(coefs = NULL, converged = FALSE)
    })
}

# ── Run per contrast ──────────────────────────────────────────────────────────
for (i in seq_along(contrasts)) {
    contrast <- contrasts[[i]]   # e.g. "MS-Rituximab-Progressive_vs_MS-Rituximab-Stable"
    parts    <- strsplit(contrast, "_vs_")[[1]]
    g1 <- parts[1]; g2 <- parts[2]

    # Restrict to samples in the two groups
    sub_meta <- meta[meta$group %in% c(g1, g2), ]
    sub_meta$group <- droplevels(sub_meta$group)
    sub_mat  <- bin_counts[, sub_meta$sample_id]

    cat(sprintf("[Module 7d] Contrast %s: %d samples\n",
                contrast, nrow(sub_meta)))

    results <- vector("list", nrow(sub_mat))
    for (j in seq_len(nrow(sub_mat))) {
        if (j %% 500 == 0)
            cat(sprintf("  bin %d / %d\n", j, nrow(sub_mat)))
        res <- fit_zinb(sub_mat[j, ], sub_meta)
        if (!is.null(res$coefs) && nrow(res$coefs) > 0) {
            # Take the coefficient for g1 vs g2
            row_name <- paste0("group", g1)
            if (row_name %in% rownames(res$coefs)) {
                r <- res$coefs[row_name, ]
                results[[j]] <- data.frame(
                    region_id = rownames(sub_mat)[j],
                    log2FC_zinb = r["Estimate"] / log(2),
                    pval_zinb   = r["Pr(>|z|)"],
                    stringsAsFactors = FALSE
                )
            }
        }
    }

    zinb_df <- rbindlist(results[!sapply(results, is.null)])

    # BH correction
    zinb_df$padj_zinb <- p.adjust(zinb_df$pval_zinb, method = "BH")

    # Load DESeq2 results for this contrast and feature set
    deseq_file <- snakemake@input$deseq2_res[i]
    if (file.exists(deseq_file)) {
        deseq_df <- read.table(deseq_file, header = TRUE, sep = "\t",
                               stringsAsFactors = FALSE)
        merged <- merge(deseq_df, zinb_df, by = "region_id", all.x = TRUE)
        merged$zinb_confirmed <- !is.na(merged$padj_zinb) &
            merged$padj_zinb < fdr_thr &
            !is.na(merged$padj) &
            merged$padj < fdr_thr
    } else {
        merged <- zinb_df
        merged$zinb_confirmed <- merged$padj_zinb < fdr_thr
    }

    out_path <- snakemake@output$zinb_res[i]
    dir.create(dirname(out_path), recursive = TRUE, showWarnings = FALSE)
    write.table(merged, out_path, sep = "\t", quote = FALSE, row.names = FALSE)
    cat(sprintf("[Module 7d] %s: %d ZINB-confirmed features\n",
                contrast, sum(merged$zinb_confirmed, na.rm = TRUE)))
}
