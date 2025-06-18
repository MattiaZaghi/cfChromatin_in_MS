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

merged_data <- data.frame(matrix(ncol = 17, nrow = 0))


# Loop through the file paths
for (file_path in file_paths) {
  # Read the CSV file
  data <- read.csv(file_path)
  
  # Merge the data
  merged_data <- rbind(merged_data, data)
} 

# Save the merged data to a new CSV file
write.csv(merged_data, "/proj/user/mattia/Analysis/Output/H3K27ac_hg38/QC/QC_all_samples.csv", row.names = FALSE)

merged_data$Total<-NULL
merged_data$Enhancer<-NULL
# Rename the first column as "sample"
colnames(data)[1] <- "sample"
# Rename the first column as "sample"
colnames(data)[2] <- "Total reads"
# Rename the first column as "sample"
colnames(data)[5] <- "Total"
# Rename the first column as "sample"
colnames(data)[6] <- "Background"
# Rename the first column as "sample"
colnames(data)[7] <- "% of Signal in TSS"
# Rename the first column as "sample"
colnames(data)[8] <- "% of Background in TSS"
# Rename the first column as "sample"
colnames(data)[9] <- "% of Signal in Enhancers"
# Rename the first column as "sample"
colnames(data)[10] <- "% of Background in Enhancers"
# Rename the first column as "sample"
colnames(data)[15] <- "Global SNR"

# Save the merged data to a new CSV file
write.csv(merged_data, "/proj/user/mattia/Analysis/Output/H3K27ac_hg38/QC/QC_all_samples.csv", row.names = FALSE)



#plot Yield


# Exclude specified samples
#exclude_samples <- c("12179-P-MS-Rituximab-Progressive_H3K27ac_ChIP",
                     #"14131-P-MS-Rituximab-Stable_H3K27ac_ChIP",
                     #"14-229-P-MS-Rituximab-Stable_H3K27ac_ChIP",
                     #"18070-P-MS-Rituximab-Progressive_H3K27ac_ChIP",
                     #"18075-P-MS-Rituximab-Progressive_H3K27ac_ChIP")
data <- data[!data$sample %in% exclude_samples, ]

# Group samples
data <- data %>%
  dplyr::mutate(Group = case_when(
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
    TRUE ~ sample # Keep the original name if no condition is met
    
      ))
    
exclude_samples <- c("NA")
data <- data[!data$Group %in% exclude_samples, ]

# Reshape data to long format
to_plot <- data %>%
  pivot_longer(cols = c("Total","Background"), names_to = "Signal_Type", values_to = "value")


# Reorder Signal_Type levels
to_plot$Signal_Type <- factor(to_plot$Signal_Type, levels = c("Total","Background"))


# Reorder Group levels
to_plot$Group <- factor(to_plot$Group, levels = c("Baca Ctrl", "Baca Canc", "Plasma V3 Fresh",
                                                  "Plasma V3","Plasma V3 old",
                                                  "Plasma V2","Plasma V1","CSF V2","CSF V1"))# Specify your desired order here


# Plotting
ggplot(to_plot, aes(x = Group, y = value, fill = Group)) +
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
ggsave("/proj/user/mattia/Analysis/Output/H3K27ac_hg38/SNR_ChIP_signal.png", plot = last_plot(), device = NULL, path = NULL, width = 130, height = 115, units = "mm", dpi = 300, limitsize = TRUE)  

