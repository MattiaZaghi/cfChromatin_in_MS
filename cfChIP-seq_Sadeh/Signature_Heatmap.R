suppressMessages(library(rtracklayer))
suppressMessages(library(tools))
suppressMessages(library(ggplot2))
suppressMessages(library(tidyverse))
suppressMessages(library(dplyr))
suppressMessages(library(reshape2))
# Directory path
directory_path <- "/date/gcb/gcb_MZ/Analysis/Output/H3K27ac/"

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
dataframes2 <- list()
# Get list of files in the directory
files2 <- list.files(directory_path, pattern = "pValues.csv$", full.names = TRUE)

# Iterate over files and read them into dataframes
for (file in files2) {
  df <- read.csv(file)
  
  # Extract the part between PLOTSIGNATURES_ and -Zscores
  column_name <- sub(".*PLOTSIGNATURES_(.*)-pValues.*", "\\1", basename(file))
  
  # Add the extracted part as a new column
  df$Signature <- column_name
  
  dataframes2[[file]] <- df
}
# Merge all dataframes into a unique dataframe
merged_df2 <- do.call(rbind, dataframes2)
# Print the merged dataframe
print(merged_df)

merged_df<-merged_df %>% dplyr::rename(sample=1)
merged_df$Signature<-NULL
merged_df2<-merged_df2 %>% dplyr::rename(sample=1)
merged_df2$Signature<-NULL

merged_df<-merged_df %>%
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

merged_df2<-merged_df2 %>%
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



merged_df <- merged_df[merged_df$Type %in% c("P-V2"), ]
merged_df2 <- merged_df2[merged_df2$Type %in% c("P-V2"), ]


# Melt the dataframe to long format for ggplot2
data_long <- reshape2::melt(merged_df, id.vars = "sample")
# Melt the dataframe to long format for ggplot2
data_long2 <- reshape2::melt(merged_df2, id.vars = "sample") 
data_long2 <- data_long2 %>% dplyr::rename(pval=3)
data_long_merge<- full_join(data_long,data_long2)




data <- data_long_merge %>%
  mutate(group = case_when(
    grepl("-P-.*MS", sample) ~ "MS",
    grepl("-P-.*H", sample) ~ "Healthy donors",
    TRUE ~ NA_character_
  ))
# Add a group column based on the variable pattern
data <- data_long_merge %>%
  mutate(group = case_when(
    sample == "04061-C-100_H3K27ac_ChIP" ~ "CSF V1",
    sample == "04061-C-200_H3K27ac_ChIP" ~ "CSF V1",
    sample == "04061-P-100_H3K27ac_Cut-Tag" ~ "NA",
    sample == "04061-P-200_H3K27ac_Cut-Tag" ~ "NA",
    sample == "07068-P-MS-Rituximab-Stable_H3K27ac_ChIP" ~ "Plasma V3",
    sample == "12097-C_H3K27ac_ChIP" ~ "CSF V1",
    sample == "12097-P_H3K27ac_ChIP" ~ "Plasma V1",
    sample == "12179-P-MS-Rituximab-Progressive_H3K27ac_ChIP" ~ "NA",
    sample == "14-031-C-RR_H3K27ac_ChIP" ~ "CSF V1",
    sample == "14131-P-MS-Rituximab-Stable_H3K27ac_ChIP" ~ "Plasma V2",
    sample == "14223-C_H3K27ac_ChIP" ~ "CSF V1",
    sample == "14-223-C-RR_H3K27ac_ChIP" ~ "CSF V2",
    sample == "14223-P_H3K27ac_ChIP" ~ "Plasma V1",
    sample == "14-229-P-MS-Rituximab-Stable_H3K27ac_ChIP" ~ "NA",
    sample == "16057-P-100_H3K27ac_ChIP" ~ "Plasma V1",
    sample == "16057-P-200_H3K27ac_ChIP" ~ "Plasma V1",
    sample == "16170-C_H3K27ac_ChIP" ~ "CSF V1",
    sample == "16170-P_H3K27ac_ChIP" ~ "Plasma V1",
    sample == "18070-P-MS-Rituximab-Progressive_H3K27ac_ChIP" ~ "Plasma V2",
    sample == "18075-P-MS-Rituximab-Progressive_H3K27ac_ChIP" ~ "Plasma V2",
    sample == "19019-C_H3K27ac_ChIP" ~ "CSF V1",
    sample == "19019-P_H3K27ac_ChIP" ~ "Plasma V1",
    sample == "19-020-C-RR_H3K27ac_ChIP" ~ "CSF V2",
    sample == "20-014-C-RR_H3K27ac_ChIP" ~ "CSF V2",
    sample == "20-015-C-RR_H3K27ac_ChIP" ~ "CSF V2",
    sample == "20-027-C-RR_H3K27ac_ChIP" ~ "CSF V2",
    sample == "20-030-C-RR_H3K27ac_ChIP" ~ "CSF V2",
    sample == "20-031-C-RR_H3K27ac_ChIP" ~ "CSF V2",
    sample == "20-033-C-RR_H3K27ac_ChIP" ~ "CSF V2",
    sample == "20-034-C-RR_H3K27ac_ChIP" ~ "CSF V2",
    sample == "20-040-C-RR_H3K27ac_ChIP" ~ "CSF V2",
    sample == "20-042-C-RR_H3K27ac_ChIP" ~ "CSF V2",
    sample == "20-053-C-RR_H3K27ac_ChIP" ~ "CSF V2",
    sample == "20-056-C-RR_H3K27ac_ChIP" ~ "CSF V2",
    sample == "GSM7787973_HP030132_H3K27Ac" ~ "Baca Ctrl",
    sample == "GSM7787975_HP030642_H3K27Ac" ~ "Baca Ctrl",
    sample == "GSM7787978_HP031645_H3K27Ac" ~ "Baca Ctrl",
    sample == "GSM7787980_HP034881_H3K27Ac" ~ "Baca Ctrl",
    sample == "GSM7787982_HP035094_H3K27Ac" ~ "Baca Ctrl",
    sample == "GSM7787985_HP038748_H3K27Ac" ~ "Baca Ctrl",
    sample == "GSM7787993_HP056703_H3K27Ac" ~ "Baca Ctrl",
    sample == "GSM7788084_K27_BELI68_2" ~ "Baca Canc",
    sample == "GSM7788101_K27_BMS7_2" ~ "Baca Canc",
    sample == "GSM7788108_K27_HPC81_2" ~ "Baca Canc",
    sample == "GSM7788120_k27_JR6BR" ~ "Baca Canc",
    sample == "GSM7788134_K27_Merk165" ~ "Baca Canc",
    sample == "GSM7788168_K27_PS245_2" ~ "Baca Canc",
    sample == "GSM7788183_K27_SMAC05_2" ~ "Baca Canc",
    sample == "H10-P-Ctrl_H3K27ac_ChIP" ~ "Plasma V3 Fresh",
    sample == "H11_H3K27ac_S3-nanoCT" ~ "NA",
    sample == "H11-P-Ctrl_H3K27ac_ChIP" ~ "Plasma V3 Fresh",
    sample == "H12-P-Ctrl_H3K27ac_ChIP" ~ "Plasma V3 Fresh",
    sample == "H13-P-Ctrl_H3K27ac_ChIP" ~ "Plasma V3 Fresh",
    sample == "H14-P-Ctrl_H3K27ac_ChIP" ~ "Plasma V3 Fresh",
    sample == "H15_H3K27ac_nanoCT" ~ "NA",
    sample == "H15_H3K27ac_S3-nanoCT" ~ "NA",
    sample == "H15-P-Ctrl_H3K27ac_ChIP" ~ "Plasma V3 Fresh",
    sample == "H16-P-Ctrl_H3K27ac_ChIP" ~ "Plasma V3 Fresh",
    sample == "H17-P-Ctrl_H3K27ac_ChIP" ~ "Plasma V3 Fresh",
    sample == "H18-P-Ctrl_H3K27ac_ChIP" ~ "Plasma V3 Fresh",
    sample == "H19-P-Ctrl_H3K27ac_ChIP" ~ "Plasma V3 Fresh",
    sample == "H20-P-Ctrl_H3K27ac_ChIP" ~ "Plasma V3 Fresh",
    sample == "H21_H3K27ac_S3-nanoCT" ~ "Plasma V3 Fresh",
    sample == "H21-P-Ctrl_H3K27ac_ChIP" ~ "Plasma V3 Fresh",
    sample == "H22-P-Ctrl_H3K27ac_ChIP" ~ "Plasma V3 Fresh",
    sample == "H23-P-Ctrl_H3K27ac_ChIP" ~ "Plasma V3 Fresh",
    sample == "H24-P-Ctrl_H3K27ac_ChIP" ~ "Plasma V3 Fresh",
    sample == "H5-P-Ctrl_H3K27ac_ChIP" ~ "Plasma V3",
    sample == "H5-P_H3K27ac_ChIP" ~ "Plasma V2",
    sample == "H6-P-Ctrl_H3K27ac_ChIP" ~ "Plasma V3",
    sample == "H6-P_H3K27ac_ChIP" ~ "Plasma V2",
    sample == "H7_H3K27ac_nanoCT" ~ "NA",
    sample == "H7-P-Ctrl_H3K27ac_ChIP" ~ "Plasma V3",
    sample == "H7-P_H3K27ac_ChIP" ~ "Plasma V2",
    sample == "H8-P-1D_H3K36me3_ChIP" ~ "NA",
    sample == "H8-P_H3K27ac_ChIP" ~ "Plasma V2",
    sample == "H9-P-Ctrl_H3K27ac_ChIP" ~ "Plasma V3 Fresh",
    sample == "P12179-P-MS-Rituximab-Prog_H3K27ac_ChIP" ~ "Plasma V2",
    sample == "P12179-P-MS-Rituximab-Prog-pA_H3K27ac_ChIP" ~ "Plasma V3",
    sample == "P14020-P-MS-Rituximab-Stable_H3K27ac_ChIP" ~ "Plasma V3",
    sample == "P14131-P-MS-Rituximab-Stable_H3K27ac_ChIP" ~ "Plasma V2",
    sample == "P14131-P-MS-Rituximab-Stable-pA_H3K27ac_ChIP" ~ "Plasma V3",
    sample == "P14229-P-MS-Rituximab-Stable_H3K27ac_ChIP" ~ "Plasma V2",
    sample == "P14229-P-MS-Rituximab-Stable-pA_H3K27ac_ChIP" ~ "Plasma V3",
    sample == "P14245-P-MS-Rituximab-Stable_H3K27ac_ChIP" ~ "Plasma V3",
    sample == "P15024-P-MS-Rituximab-Stable_H3K27ac_ChIP" ~ "Plasma V3",
    sample == "P15041-P-MS-Rituximab-Stable_H3K27ac_ChIP" ~ "Plasma V3",
    sample == "P16186-P-MS-Rituximab-Stable_H3K27ac_ChIP" ~ "Plasma V3",
    sample == "P18070-P-MS-Rituximab-Prog_H3K27ac_ChIP" ~ "Plasma V2",
    sample == "P18070-P-MS-Rituximab-Prog-pA_H3K27ac_ChIP" ~ "Plasma V3",
    sample == "P18075-P-MS-Rituximab-Prog_H3K27ac_ChIP" ~ "Plasma V2",
    sample == "P18075-P-MS-Rituximab-Prog-pA_H3K27ac_ChIP" ~ "Plasma V3",
    sample == "P20015-P-Ctrl_H3K27ac_ChIP" ~ "Plasma V3 old",
    sample == "P20027-P-Ctrl_H3K27ac_ChIP" ~ "Plasma V3 old",
    sample == "P20030-P-Ctrl_H3K27ac_ChIP" ~ "Plasma V3 old",
    sample == "P20040-P-Ctrl_H3K27ac_ChIP" ~ "Plasma V3 old",
    sample == "P24116-P-MS-New-RR_H3K27ac_ChIP" ~ "Plasma V3",
    sample == "P24117-P-MS-New-RR_H3K27ac_ChIP" ~ "Plasma V3",
    sample == "P24118-P-MS-New-RR_H3K27ac_ChIP" ~ "Plasma V3",
    sample == "P24126-P-MS-New-RR_H3K27ac_ChIP" ~ "Plasma V3",
    sample == "P24132-P-MS-New-RR_H3K27ac_ChIP" ~ "Plasma V3",
    sample == "P24134-P-MS-New-RR_H3K27ac_ChIP" ~ "Plasma V3",
    sample == "P24136-P-MS-New-RR_H3K27ac_ChIP" ~ "Plasma V3",
    sample == "H1-P_H3K27ac_ChIP" ~ "Plasma V2",
    sample == "H2-P_H3K27ac_ChIP" ~ "Plasma V2",
    sample == "H3-P_H3K27ac_ChIP" ~ "Plasma V2",
    sample == "H4-P_H3K27ac_ChIP" ~ "Plasma V2",
    TRUE ~ sample  
    ))
# Exclude specified samples
exclude_samples <- c("S3-nanoCT",
                     "nanoCT",
                     "Rituximab-Progressive",
                     "Rituximab-Stable",
                     "New-MS",
                     "Control-pA",
                     "Baca et al."
  )
exclude_samples <- c("NA")
data <- data[!data$value %in% exclude_samples, ]


samples_to_keep <- c("H11", "H7", "H15", "H21")
pattern <- paste(samples_to_keep, collapse = "|")
data <- data[grepl(pattern, data$sample), ]

exclude_signatures <- c("Leukocytes","Lymphocytes","Monocytes","Placenta")
data <- data[!data$variable %in% exclude_signatures, ]

data<-data_long_merge

data$pval <- as.numeric(data$pval)
data$value <- as.numeric(data$value)

data <- na.omit(data)


data <- data[data$Type %in% c("Baca healthy", "Baca cancer","P-V2","P-V3"), ]

exclude_samples <- c("H5-P-Ctrl_H3K27ac_ChIP-V2",
                     "H6-P-Ctrl_H3K27ac_ChIP-V2",
                     "H8-P-Ctrl_H3K27ac_ChIP-V2",
                     "H1-P-Ctrl_H3K27ac_ChIP-V2",
                     "H4-P-Ctrl_H3K27ac_ChIP-V2",
                     "H7-P-Ctrl_H3K27ac_ChIP-V2"
)
data <- data[!data$sample %in% exclude_samples, ]

ggplot(data, aes(x = variable, y = sample, size = round(pval), color = value)) +
  geom_point() +
  scale_color_viridis_c(option = "plasma") +
  scale_size_continuous(range = c(1, 20)) +
  theme_minimal() + 
  theme(
    # Increase size of all text elements
    text = element_text(size = 20),
    # Adjust specific elements
    axis.text.x = element_text(angle = 45, hjust = 1, size = 20),
    axis.text.y = element_text(size = 20),
    axis.title = element_text(size = 24),
    plot.title = element_text(size = 26),
    legend.text = element_text(size = 20),
    legend.title = element_text(size = 22),
    strip.text = element_text(size = 22, face = "bold"),  # Facet labels
    panel.grid.major = element_line(color = "gray90"),
    panel.grid.minor = element_blank(),
    legend.position = "right",
    panel.spacing.y = unit(0.1, "lines")  # Reduce vertical spacing between panels
  ) +
  labs(
    title = "H3K27ac ChIP-seq Values by Sample and Tissue",
    x = "Tissue",
    y = "Sample",
    size = "Pvalue",
    color = "Zscore"
  ) +
  facet_grid(group ~ ., scales = "free_y", space = "free", switch = "y") +  # Move facet labels to the left
  coord_cartesian(clip="off")  # Adjust the coordinate system

# Save the plot as a PDF (vector format for best quality)
ggsave(
  filename = "/date/gcb/gcb_MZ/Analysis/Output/H3K27ac/plots/Signature_V2.png",
  plot = last_plot(),
  width =20,
  height = 25,
  units = "in",
  dpi = 300
)

exclude_signatures <- c("Leukocytes","Lymphocytes","Monocytes","Placenta")
data_no_immune <- data[!data$variable %in% exclude_signatures, ]

# Exclude specified samples
exclude_samples <- c("H5-P-Ctrl_H3K27ac_ChIP-V2",
                     "H6-P-Ctrl_H3K27ac_ChIP-V2",
                     "H8-P-Ctrl_H3K27ac_ChIP-V2",
                     "H1-P-Ctrl_H3K27ac_ChIP-V2",
                     "H4-P-Ctrl_H3K27ac_ChIP-V2",
                     "H7-P-Ctrl_H3K27ac_ChIP-V2"
)
data_no_immune<- data_no_immune[!data_no_immune$sample %in% exclude_samples, ]


ggplot(data_no_immune, aes(x = variable, y = sample, size = round(pval), color = value)) +
  geom_point() +
  scale_color_viridis_c(option = "plasma") +
  scale_size_continuous(range = c(1, 20)) +
  theme_minimal() + 
  theme(
    # Increase size of all text elements
    text = element_text(size = 20),
    # Adjust specific elements
    axis.text.x = element_text(angle = 45, hjust = 1, size = 20),
    axis.text.y = element_text(size = 20),
    axis.title = element_text(size = 24),
    plot.title = element_text(size = 26),
    legend.text = element_text(size = 20),
    legend.title = element_text(size = 22),
    strip.text = element_text(size = 22, face = "bold"),  # Facet labels
    panel.grid.major = element_line(color = "gray90"),
    panel.grid.minor = element_blank(),
    legend.position = "right",
    panel.spacing.y = unit(0.1, "lines")  # Reduce vertical spacing between panels
  ) +
  labs(
    title = "H3K27ac ChIP-seq Values by Sample and Tissue no Immune",
    x = "Tissue",
    y = "Sample",
    size = "Pvalue",
    color = "Zscore"
  ) +
  facet_grid(group ~ ., scales = "free_y", space = "free", switch = "y") +  # Move facet labels to the left
  coord_cartesian(clip="off")  # Adjust the coordinate system
ggsave(
  filename = "/date/gcb/gcb_MZ/Analysis/Output/H3K27ac/plots/Signature_V2_no_immune.png",
  plot = last_plot(),
  width =20,
  height = 25,
  units = "in",
  dpi = 300
)






