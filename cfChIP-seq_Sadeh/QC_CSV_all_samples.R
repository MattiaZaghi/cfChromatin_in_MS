suppressMessages(library(rtracklayer))
suppressMessages(library(tools))
suppressMessages(library(ggplot2))
suppressMessages(library(tidyverse))
suppressMessages(library(dplyr))
suppressMessages(library(extrafont))

# Define the directory
directory <- "/proj/user/mattia/Analysis/Output/H3K27ac_hg38/QC/"

# Get a list of all CSV files in the directory
all_files <- list.files(directory, pattern = "\\.csv$", full.names = TRUE)

# Remove the files that contain "_H3K4me3_ChIP"
file_paths <- all_files


# Initialize an empty data frame
merged_data <- data.frame()

# Loop through the file paths
for (file_path in file_paths) {
  # Read the CSV file
  data <- read.csv(file_path)
  
  # Merge the data
  merged_data <- rbind(merged_data, data)
}
merged_data$Total<-NULL
merged_data$TSS<-NULL
# Rename the first column as "sample"
colnames(merged_data)[1] <- "sample"
# Rename the first column as "sample"
colnames(merged_data)[5] <- "Total"
# Rename the first column as "sample"
colnames(merged_data)[6] <- "Background"
# Rename the first column as "sample"
colnames(merged_data)[7] <- "TSS"
# Rename the first column as "sample"
colnames(merged_data)[8] <- "Background in TSS"
# Rename the first column as "sample"
colnames(merged_data)[9] <- "Enhancers"
# Rename the first column as "sample"
colnames(merged_data)[10] <- "Background in Enhancers"
# Rename the first column as "sample"
colnames(merged_data)[15] <- "Global SNR"

# Save the merged data to a new CSV file
write.csv(merged_data, "/date/gcb/gcb_MZ/Analysis/Output/H3K4me3/QC_all_Sadeh_samples.csv", row.names = FALSE)



#plot Yield


# Exclude specified samples
exclude_samples <- c("12179-P-MS-Rituximab-Progressive_H3K27ac_ChIP",
                     "14131-P-MS-Rituximab-Stable_H3K27ac_ChIP",
                     "14-229-P-MS-Rituximab-Stable_H3K27ac_ChIP",
                     "18070-P-MS-Rituximab-Progressive_H3K27ac_ChIP",
                     "18075-P-MS-Rituximab-Progressive_H3K27ac_ChIP")
data <- merged_data[!merged_data$sample %in% exclude_samples, ]

# Group samples
data <- data %>%
  dplyr::mutate(Group = case_when(
    grepl("S3-nanoCT", sample) ~ "S3-nanoCT",
    grepl("nanoCT", sample) ~ "nanoCT",
    grepl("^H", sample) ~ "ctrl-pA",
    grepl("P20", sample) ~ "ctrl-pA-old",
    grepl("^GSM", sample) ~ "baca et al.",
    grepl("New-RR", sample) ~ "New-pA-MS",
    TRUE ~ "Rituximab-pA-MS"
  ))
exclude_samples <- c("New-pA-MS",
                     "ctrl-pA-old",
                     "Rituximab-pA-MS")
data <- data[!data$Group %in% exclude_samples, ]

# Reshape data to long format
to_plot <- data %>%
  pivot_longer(cols = c("Global SNR"), names_to = "Signal_Type", values_to = "value")


# Reorder Signal_Type levels
to_plot$Signal_Type <- factor(to_plot$Signal_Type, levels = c("Global SNR"))


# Plotting
ggplot(to_plot, aes(x = Group, y = value, fill = Group)) +
  geom_boxplot(width = 0.1, color = "#000000", outlier.shape = NA) +
  geom_jitter(width = 0.1, size = 1, color = "#000000") +
  ggthemes::theme_base() +
  xlab("") +
  ylab("Signal/Noise") +
  ylim(0,80)+
  facet_wrap(~ Signal_Type, scales = "free_y") +
  theme(panel.border = element_rect()) +
  theme_classic() + theme(legend.position = "none") +
  theme(axis.text.x = element_text(size = 9,  angle = 45, hjust = 1),
        axis.text.y = element_text(size = 12, ),
        axis.title.y = element_text(size = 12, ),
        axis.line = element_line(size = 1))
ggsave("/proj/user/mattia/Analysis/Output/H3K27ac_hg38/Signal_SNR_nanoCT.png", plot = last_plot(), device = NULL, path = NULL, width = 130, height = 115, units = "mm", dpi = 300, limitsize = TRUE)  

