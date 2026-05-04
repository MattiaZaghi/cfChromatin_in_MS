# Module 2 — Anchor QC barplot.
# Snakemake R script.

suppressPackageStartupMessages({
    library(ggplot2)
    library(dplyr)
})

sf   <- read.table(snakemake@input$sf,   header = TRUE, sep = "\t", stringsAsFactors = FALSE)
meta <- read.table(snakemake@input$meta, header = TRUE, sep = "\t",
                   comment.char = "#",   stringsAsFactors = FALSE)

df <- sf
if (!("group" %in% colnames(df))) {
    extra <- meta[, c("sample_id", "group", "sample_type")]
    df <- merge(sf, extra, by = "sample_id", all.x = TRUE)
}

group_order  <- c("Ctrl", "NEW", "MS-Rituximab-Stable", "MS-Rituximab-Progressive")
group_colors <- c(
    Ctrl                         = "#4878CF",
    NEW                          = "#6ACC65",
    `MS-Rituximab-Stable`        = "#D65F5F",
    `MS-Rituximab-Progressive`   = "#B47CC7"
)

df$sf_flagged <- df$sf_flagged %in% c(TRUE, "True", "TRUE")
df$group <- factor(df$group, levels = group_order)
df <- df[order(df$group, df$constitutive_sf), ]
df$sample_id <- factor(df$sample_id, levels = df$sample_id)

p1 <- ggplot(df, aes(x = sample_id, y = anchor_reads / 1e6, fill = group)) +
    geom_col() +
    geom_point(data = df[df$sf_flagged, ], aes(y = anchor_reads / 1e6),
               shape = 8, colour = "black", size = 2) +
    scale_fill_manual(values = group_colors, na.value = "grey70") +
    theme_bw(base_size = 9) +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 5)) +
    labs(title = "Constitutive anchor read counts",
         x = NULL, y = "Anchor reads (M)", fill = "Group")

p2 <- ggplot(df, aes(x = sample_id, y = constitutive_sf, fill = group)) +
    geom_col() +
    geom_hline(yintercept = 1, linetype = "dashed", colour = "grey40") +
    geom_point(data = df[df$sf_flagged, ], aes(y = constitutive_sf),
               shape = 8, colour = "black", size = 2) +
    scale_fill_manual(values = group_colors, na.value = "grey70") +
    theme_bw(base_size = 9) +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 5)) +
    labs(title = "Constitutive scaling factors  (★ = flagged >3 SD)",
         x = NULL, y = "Scaling factor", fill = "Group")

dir.create(dirname(snakemake@output[[1]]), recursive = TRUE, showWarnings = FALSE)
pdf(snakemake@output[[1]], width = 14, height = 8)
print(p1)
print(p2)
dev.off()

cat(sprintf("[Module 2] %d samples plotted, %d flagged\n",
            nrow(df), sum(df$sf_flagged, na.rm = TRUE)))
