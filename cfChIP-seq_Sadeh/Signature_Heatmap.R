suppressMessages(library(rtracklayer))
suppressMessages(library(tools))
suppressMessages(library(ggplot2))
suppressMessages(library(tidyverse))
suppressMessages(library(dplyr))
suppressMessages(library(reshape2))
# Directory path
directory_path <- "/proj/user/mattia/Analysis/Output/H3K27ac_hg38/"

# List to store dataframes
dataframes <- list()

# Get list of files in the directory
files <- list.files(directory_path, pattern = "Zscores.csv$", full.names = TRUE)

# Iterate over files and read them into dataframes
for (file in files) {
  df <- read.csv(file)
  
  # Extract the part between PLOTSIGNATURES_ and -Zscores
  column_name <- sub(".*PLOTSIGNATURES_(.*)-Zscores.*", "\\1", basename(file))
  
  # Add the extracted part as a new column
  df$Signature <- column_name
  
  dataframes[[file]] <- df
}

# Merge all dataframes into a unique dataframe
merged_df <- do.call(rbind, dataframes)
# List to store dataframes
dataframes <- list()
# Get list of files in the directory
files <- list.files(directory_path, pattern = "pValues.csv$", full.names = TRUE)

# Iterate over files and read them into dataframes
for (file in files) {
  df <- read.csv(file)
  
  # Extract the part between PLOTSIGNATURES_ and -Zscores
  column_name <- sub(".*PLOTSIGNATURES_(.*)-pValues.*", "\\1", basename(file))
  
  # Add the extracted part as a new column
  df$Signature <- column_name
  
  dataframes[[file]] <- df
}
# Merge all dataframes into a unique dataframe
merged_df2 <- do.call(rbind, dataframes)
# Print the merged dataframe
print(merged_df)

merged_df<-merged_df %>% dplyr::rename(sample=1)
merged_df$Signature<-NULL
merged_df2<-merged_df2 %>% dplyr::rename(sample=1)
merged_df2$Signature<-NULL
# Load necessary libraries
library(ggplot2)

# Exclude specified samples
exclude_samples <- c("P20027-P-Ctrl_H3K27ac_ChIP",
                     "P20040-P-Ctrl_H3K27ac_ChIP",
                     "P20030-P-Ctrl_H3K27ac_ChIP",
                     "P20015-P-Ctrl_H3K27ac_ChIP",
                     
                     )
merged_df <- merged_df[!merged_df$sample %in% exclude_samples, ]
merged_df2 <- merged_df2[!merged_df2$sample %in% exclude_samples, ]
# Melt the dataframe to long format for ggplot2
data_long <- reshape2::melt(merged_df, id.vars = "sample")
# Melt the dataframe to long format for ggplot2
data_long2 <- reshape2::melt(merged_df2, id.vars = "sample") 
data_long2 <- data_long2 %>% dplyr::rename(pval=3)
data_long_merge<- full_join(data_long,data_long2)
# Add a group column based on the variable pattern
data <- data_long_merge %>%
  mutate(group = case_when(
    grepl("S3-nanoCT", sample) ~ "S3-nanoCT",
    grepl("nanoCT", sample) ~ "nanoCT",
    grepl("Ctrl", sample) ~ "Control-pA",
    grepl("Stable", sample) ~ "Rituximab-Stable",
    grepl("Progressive|Prog", sample) ~ "Rituximab-Progressive",
    grepl("GSM", sample) ~ "Baca et al.",
    TRUE ~ "New-MS"
  ))
# Exclude specified samples
exclude_samples <- c("Rituximab-Progressive",
                     "Rituximab-Stable",
                     "New-MS",
                     "Control-pA")
data <- data[!data$group %in% exclude_samples, ]
ggplot(data, aes(x = variable, y = sample, size = round(pval), color = value)) +
  geom_point() +
  scale_color_viridis_c(option = "plasma") +
  scale_size_continuous(range = c(1, 10)) +
  theme_minimal() +
  theme(
    # Increase size of all text elements
    text = element_text(size = 14),
    # Adjust specific elements
    axis.text.x = element_text(angle = 45, hjust = 1, size = 12),
    axis.text.y = element_text(size = 12),
    axis.title = element_text(size = 16),
    plot.title = element_text(size = 18),
    legend.text = element_text(size = 12),
    legend.title = element_text(size = 14),
    strip.text = element_text(size = 14, face = "bold"),  # Facet labels
    panel.grid.major = element_line(color = "gray90"),
    panel.grid.minor = element_blank(),
    legend.position = "right"
  ) +
  labs(
    title = "H3K27ac ChIP-seq Values by Sample and Tissue",
    x = "Tissue",
    y = "Sample",
    size = "Pvalue",
    color = "Zscore"
  ) +
  facet_grid(group ~ ., scales = "free_y", space = "free")


# Save the plot as a PDF (vector format for best quality)
ggsave(
  filename = "/proj/user/mattia/Analysis/Output/H3K27ac_hg38/Signature.pdf",
  plot = last_plot(),
  width = 15,
  height = 12,
  units = "in",
  dpi = 300
)

