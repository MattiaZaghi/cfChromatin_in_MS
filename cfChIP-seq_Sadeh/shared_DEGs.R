library(dplyr)
library(readr)
library(tidyr)
library(ggplot2)
library(VennDiagram)
library(UpSetR)
library(clusterProfiler)
library(org.Hs.eg.db)

# Set the directory containing your CSV files
setwd("/date/gcb/gcb_MZ/Analysis/Output/H3K27ac_hg38/DiffGenes/")

# List all CSV files in the directory
files <- list.files(pattern = "*.csv")

# Specify the files you want to exclude (if any)
#exclude_files <- c("H5-P_H3K27ac_ChIP.csv", "sample2.csv")  # Replace with your actual file names

# Exclude the specified files
#files_to_use <- setdiff(files, exclude_files)

# Read the remaining CSV files into a list of data frames
data_list <- lapply(files, read.csv)

# Extract the gene names from each data frame
gene_lists <- lapply(data_list, function(df) df[, 1])

# Combine all gene lists into one vector
all_genes <- unlist(gene_lists)

# Count the occurrences of each gene
gene_counts <- table(all_genes)

# Find genes that appear in at least 2 samples
common_genes <- names(gene_counts[gene_counts >= 2])

# Print the common genes
print(common_genes)


# Create a named list of gene sets for the Venn diagram
gene_sets <- setNames(gene_lists, files)

# Convert gene_sets to a data frame suitable for UpSetR
gene_sets_df <- fromList(gene_sets)

# Save the UpSet plot as a PDF
pdf("upset_plot.pdf")
upset(gene_sets_df, sets = names(gene_sets), order.by = "freq")
dev.off()

# Save the UpSet plot as a PNG
png("upset_plot.png")
upset(gene_sets_df, sets = names(gene_sets), order.by = "freq")
dev.off()

# Specify the three files you want to compare
files_to_use <- c("18070-P-SP_H3K27ac_ChIP.csv", "12-179-P-RR_H3K27ac_ChIP.csv", "14-229-P-RR_H3K27ac_ChIP.csv")  # Replace with your actual file names

# Read the specified CSV files into a list of data frames
data_list <- lapply(files_to_use, read.csv)

# Extract the gene names from each data frame
gene_lists <- lapply(data_list, function(df) df[, 1])

# Find the intersection of genes present in all three datasets
common_genes <- Reduce(intersect, gene_lists)

# Print the common genes
print(common_genes)

# Save the common genes to a CSV file
write.csv(data.frame(Gene = common_genes), "common_genes.csv", row.names = FALSE)

# Read the gene list from a CSV file
gene_list <- read.csv("common_genes.csv", stringsAsFactors = FALSE)$Gene

# Convert gene symbols to Entrez IDs
gene_entrez <- bitr(gene_list, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)

# Perform KEGG pathway enrichment analysis
kegg_enrich <- enrichKEGG(gene = gene_entrez$ENTREZID, organism = 'hsa', pvalueCutoff = 0.05)

# Convert the results to a data frame
kegg_results <- as.data.frame(kegg_enrich)

# Select the top 10 enriched pathways
top_pathways <- kegg_results %>% top_n(10, wt = -p.adjust)

# Create a bar plot
ggplot(top_pathways, aes(x = reorder(Description, -p.adjust), y = -log10(p.adjust))) +
  geom_bar(stat = "identity") +
  coord_flip() +
  labs(title = "Top 10 Enriched KEGG Pathways", x = "Pathway", y = "-log10(p.adjust)") +
  theme_minimal()

# Save the plot
ggsave("top_enriched_kegg_pathways.png")