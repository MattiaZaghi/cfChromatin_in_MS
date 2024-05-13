suppressMessages(library(rtracklayer))
suppressMessages(library(tools))
suppressMessages(library(ggplot2))

# Set the file path
file_paths <- "/date/gcb/gcb_MZ/Analysis/Samples/H3K4me3"

# List all files in the directory that match the pattern
files <- list.files(path = file_paths, pattern = ".*H3K4me3.*_ChIP.rdata$", full.names = TRUE)

# Use lapply to load all files
data_list <- lapply(files, readRDS)

# Extract the base names from the file paths
base_names <- basename(files)

# Remove the pattern "_H3K4me3_ChIP.rdata" from the base names
names <- sub("_H3K4me3_ChIP.rdata$", "", base_names)

# Assign names to the list elements
names(data_list) <- names




# Assuming H1, H2, H3, and H4 are your objects with gene counts
gene_counts_list <- lapply(data_list, function(x) x$GeneCounts)

# Combine them into a data frame
gene_counts_df <- do.call(cbind, gene_counts_list)

# Calculate row means (average gene counts)
average_gene_counts <- rowMeans(gene_counts_df, na.rm = TRUE)

#load common Genes (HouseKeeping based on H3K4me3 calculation Sadeh et al.)

Common_genes<-readRDS("/date/gcb/gcb_MZ/Analysis/cfChIP-seq/SetupFiles/H3K4me3/CommonGenes.rds")

#filter the gene counts matrix keeping only the housekeeping genes
filtered_data <- average_gene_counts[names(average_gene_counts) %in% Common_genes]

# Print the result
print(filtered_gene_counts_df)

saveRDS(filtered_data,"/date/gcb/gcb_MZ/Analysis/cfChIP-seq/SetupFiles/H3K4me3/HealthyRef.rds")

