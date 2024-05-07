# Define the directory
directory <- "/date/gcb/gcb_MZ/Analysis/Output/H3K27ac"

# Get a list of all CSV files in the directory that end with "_ChIP.csv"
file_paths <- list.files(directory, pattern = "_ChIP\\.csv$", full.names = TRUE)

# Initialize an empty data frame
merged_data <- data.frame()

# Loop through the file paths
for (file_path in file_paths) {
  # Read the CSV file
  data <- read.csv(file_path)
  
  # Merge the data
  merged_data <- rbind(merged_data, data)
}


# Rename the first column as "sample"
colnames(merged_data)[1] <- "sample"
# Rename the first column as "sample"
colnames(merged_data)[7] <- "Total Signal"
# Rename the first column as "sample"
colnames(merged_data)[8] <- "Background Signal"

# Save the merged data to a new CSV file
write.csv(merged_data, "/date/gcb/gcb_MZ/Analysis/Output/H3K27ac/QC_all_samples.csv", row.names = FALSE)

# Print a success message
print("Merged data has been saved to '/path/to/merged_file.csv'")


#plot Yield



ggplot(merged_data) +
  aes(x = sample, y = Frip, fill = sample, color = sample) +
  geom_violin(trim = TRUE, color = "#000000") +
  geom_boxplot(width = 0.1, color = "#000000", fill = "#ffffff", outlier.shape = NA) +
  scale_fill_manual(values = c("#FF0000", "#FF0000","#3399FF", "#3399FF", "#3399FF","#00CC33", "#00CC33","#00CC33")) +
  scale_color_manual(values = c("#FF0000", "#FF0000","#3399FF","#3399FF","#3399FF",  "#00CC33", "#00CC33","#00CC33")) +
  ggthemes::theme_base() +
  xlab("") +
  ylab("Fraction of reads into peaks") +
  scale_x_discrete(labels = c("nanoCT_ATAC","nanoCTAR_ATAC", "Droplet_Pair-Tag_H3K27ac","nanoCT_H3K27ac","nanoCTAR_H3K27ac","Droplet_Pair-Tag_H3K27me3","nanoCT_H3K27me3", "nanoCTAR_H3K27me3"))+
  #scale_y_continuous(breaks = seq(0, 25000, by = 5000)) +
  theme(panel.border = element_rect()) +
  theme_classic() + theme(legend.position = "none") +
  theme(axis.text.x = element_text(size = 18, family = "Arial", angle = 45, hjust = 1),
        axis.text.y = element_text(size = 18, family = "Arial"),
        axis.title.y = element_text(size = 18, family = "Arial"),
        axis.line = element_line(size = 1))
