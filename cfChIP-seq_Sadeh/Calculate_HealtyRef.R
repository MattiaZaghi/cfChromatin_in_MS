suppressMessages(library(rtracklayer))
suppressMessages(library(tools))
suppressMessages(library(ggplot2))

file_paths <- "/date/gcb/gcb_MZ/Analysis/Samples/H3K27ac_ref"
files <- list.files(path = file_paths, pattern = "^Human.*rdata$", full.names = TRUE)


# Use lapply to load all files
data_list <- lapply(files, readRDS)

# Extract the base names from the file paths
base_names <- basename(files)

# Remove the pattern "__H3K27ac_ChIP-SE.rdata" from the base names
names <- sub("_H3K27ac_ChIP-SE.rdata$", "", base_names)

# Assign names to the list elements
names(data_list) <- object_names




# Assuming H1, H2, H3, and H4 are your objects with gene counts
gene_counts_list <- lapply(data_list, function(x) x$GeneCounts)

# Combine them into a data frame
gene_counts_df <- do.call(cbind, gene_counts_list)

# Calculate row means (average gene counts)
average_gene_counts <- rowMeans(gene_counts_df, na.rm = TRUE)

# Print the result
print(average_gene_counts)

saveRDS(average_gene_counts,"/date/gcb/gcb_MZ/Analysis/cfChIP-seq/SetupFiles/H3K27ac_ref/HealthyRef.rds")
