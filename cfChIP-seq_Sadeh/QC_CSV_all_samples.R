suppressMessages(library(rtracklayer))
suppressMessages(library(tools))
suppressMessages(library(ggplot2))
suppressMessages(library(tidyverse))
suppressMessages(library(dplyr))
suppressMessages(library(extrafont))

# Define the directory
directory <- "/date/gcb/gcb_MZ/Analysis/Output/H3K27ac/QC/"

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

# Save the merged data to a new CSV file
write.csv(merged_data, "/date/gcb/gcb_MZ/Analysis/Output/H3K27ac/QC/QC_all_samples.csv", row.names = FALSE)

merged_data$Total<-NULL
merged_data$Enhancer<-NULL
# Rename the first column as "sample"
colnames(merged_data)[1] <- "sample"
# Rename the first column as "sample"
colnames(merged_data)[2] <- "Total reads"
# Rename the first column as "sample"
colnames(merged_data)[7] <- "% of Signal"
# Rename the first column as "sample"
colnames(merged_data)[8] <- "% of Background"
# Rename the first column as "sample"
colnames(merged_data)[9] <- "% of Signal in TSS"
# Rename the first column as "sample"
colnames(merged_data)[10] <- "% of Background in TSS"
# Rename the first column as "sample"
colnames(merged_data)[11] <- "% of Signal in Enhancers"
# Rename the first column as "sample"
colnames(merged_data)[12] <- "% of Background in Enhancers"
# Rename the first column as "sample"
colnames(merged_data)[15] <- "Global SNR"


merged_data$ON_target_reads<-round(merged_data$`Total reads`*(merged_data$`% of Signal`/100))
merged_data$OFF_target_reads<-merged_data$`Total reads`-merged_data$ON_target_reads
merged_data$ON_OFF_ratio<-merged_data$ON_target_reads/merged_data$OFF_target_reads


# Save the merged data to a new CSV file
write.csv(merged_data, "/date/gcb/gcb_MZ/Analysis/Output/H3K27ac/QC/QC_all_samples.csv", row.names = FALSE)



#plot Yield


# Exclude specified samples
#exclude_samples <- c("12179-P-MS-Rituximab-Progressive_H3K27ac_ChIP",
                     #"14131-P-MS-Rituximab-Stable_H3K27ac_ChIP",
                     #"14-229-P-MS-Rituximab-Stable_H3K27ac_ChIP",
                     #"18070-P-MS-Rituximab-Progressive_H3K27ac_ChIP",
                     #"18075-P-MS-Rituximab-Progressive_H3K27ac_ChIP")
data <- data[!data$sample %in% exclude_samples, ]

# Group samples
data <- merged_data %>%
  dplyr::mutate(Type = case_when(
    grepl("-C-.*V1", sample) ~ "C-V1",
    grepl("-C-.*V2", sample) ~ "C-V2",
    grepl("-C-.*V3", sample) ~ "C-V3",
    grepl("-P-.*V1", sample) ~ "P-V1",
    grepl("-P-.*V2", sample) ~ "P-V2",
    grepl("-P-.*V3", sample) ~ "P-V3",
    grepl("GSM.*HP", sample) ~ "Baca healthy",
    grepl("GSM", sample) ~ "Baca cancer",
    TRUE ~ NA_character_
  ))

print(data)


print(data)

    
exclude_samples <- c("NA")
data <- data[!data$Group %in% exclude_samples, ]

# Reshape data to long format
to_plot <- data %>%
  pivot_longer(cols = c("% of Signal","% of Background"), names_to = "Signal_Type", values_to = "value")


# Reorder Signal_Type levels
to_plot$Signal_Type <- factor(to_plot$Signal_Type, levels = c("% of Signal","% of Background"))


# Reorder Group levels
to_plot$Group <- factor(to_plot$Type, levels = c("Baca cancerl", "Baca healthy", "P-V1","P-V2","P-V3",
                                                  "C-V1","C-V2","C-V3"))# Specify your desired order here


# Plotting
ggplot(to_plot, aes(x =Type, y = value, fill = Type)) +
  geom_boxplot(width = 0.1, color = "#000000", outlier.shape = NA) +
  geom_jitter(width = 0.1, size = 1, color = "#000000") +
  ggthemes::theme_base() +
  xlab("") +
  ylab("% of Reads") +
  #ylim(0,80)+
  facet_wrap(~ Signal_Type, scales = "free_y") +
  theme(panel.border = element_rect()) +
  theme_classic() + theme(legend.position = "none") +
  theme(axis.text.x = element_text(size = 9,  angle = 45, hjust = 1),
        axis.text.y = element_text(size = 12, ),
        axis.title.y = element_text(size = 12, ),
        axis.line = element_line(size = 1))
ggsave("/date/gcb/gcb_MZ/Analysis/QC/Signal_percentage_frag.png", plot = last_plot(), device = NULL, path = NULL, width = 150, height = 115, units = "mm", dpi = 300, limitsize = TRUE)  


#select sample with signal enrichment over 65%


data_65 <- data[!grepl("GSM", data$sample), ]

data_65<-data_65 %>% dplyr::filter(`% of Signal`>=50)







