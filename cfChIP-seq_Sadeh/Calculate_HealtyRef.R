suppressMessages(library(rtracklayer))
suppressMessages(library(tools))
suppressMessages(library(ggplot2))

# Create a vector of file paths
file_paths <- c("/date/gcb/gcb_MZ/Analysis/Samples/H3K27ac/H1-P_H3K27ac_ChIP.rdata",
                "/date/gcb/gcb_MZ/Analysis/Samples/H3K27ac/H2-P_H3K27ac_ChIP.rdata",
                "/date/gcb/gcb_MZ/Analysis/Samples/H3K27ac/H3-P_H3K27ac_ChIP.rdata",
                "/date/gcb/gcb_MZ/Analysis/Samples/H3K27ac/H4-P_H3K27ac_ChIP.rdata",
                "/date/gcb/gcb_MZ/Analysis/Samples/H3K27ac/19019-P_H3K27ac_ChIP.rdata")

# Use lapply to load all files
data_list <- lapply(file_paths, readRDS)

# Create a vector of object names
object_names <- c("H1", "H2", "H3", "H4","19019")

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

saveRDS(average_gene_counts,"/date/gcb/gcb_MZ/Analysis/cfChIP-seq/SetupFiles/H3K27ac/HealthyRef.rds")
