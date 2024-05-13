suppressMessages(library(rtracklayer))
suppressMessages(library(tools))
suppressMessages(library(ggplot2))
suppressMessages(library(tidyverse))
suppressMessages(library(dplyr))

# Define the directory
directory <- "/date/gcb/gcb_MZ/Analysis/Output/H3K4me3"

# Get a list of all CSV files in the directory
all_files <- list.files(directory, pattern = "\\.csv$", full.names = TRUE)

# Remove the files that contain "_H3K4me3_ChIP"
file_paths <- all_files[!grepl("_H3K4me3_ChIP", all_files)]


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
colnames(merged_data)[6] <- "Total Signal"
# Rename the first column as "sample"
colnames(merged_data)[7] <- "Background Signal"

# Save the merged data to a new CSV file
write.csv(merged_data, "/date/gcb/gcb_MZ/Analysis/Output/H3K4me3/QC_all_Sadeh_samples.csv", row.names = FALSE)



#plot Yield

# Assuming the two columns you want to plot are `Total Signal` and `Another Column`
to_plot <-merged_data %>% 
  dplyr::select(`Total Signal`,`Background Signal`) %>% 
  gather(key=Group, value=value, "Total Signal","Background Signal")

# Define the order
levels_order <- c("Total Signal","Background Signal")  # replace with your actual labels

# Convert the column to a factor and specify the order of levels
to_plot$Category <- factor(to_plot$Group, levels = levels_order)

ggplot(to_plot, aes(x=Group, y=value, fill=Group)) +
  geom_boxplot(width = 0.1, color = "#000000", outlier.shape = NA) +
  geom_jitter(width = 0.1, size = 1, color = "#000000") +
  ggthemes::theme_base() +
  xlab("") +
  ylab("% of Signal") +
  scale_x_discrete(labels = c("H3K4me3 CfChromatin Background","H3K4me3 CfChromatin in TSS")) +
  theme(panel.border = element_rect()) +
  theme_classic() + theme(legend.position = "none") +
  theme(axis.text.x = element_text(size = 7, family = "Arial", angle = 45, hjust = 1),
        axis.text.y = element_text(size = 10, family = "Arial"),
        axis.title.y = element_text(size = 10, family = "Arial"),
        axis.line = element_line(size = 1))


