library(dplyr)
library(readr)
library(ggplot2)

# Define a list of directories
directories <- c("/date/gcb/gcb_MZ/Analysis/Output/H3K4me1/","/date/gcb/gcb_MZ/Analysis/Output/H3K4me3/", 
                 "/date/gcb/gcb_MZ/Analysis/Output/H3K27ac_hg38/","/date/gcb/gcb_MZ/Analysis/Output/H3K27ac/")

# Initialize an empty list to store file paths
file_list <- list()

# Loop through each directory and find files matching the pattern
for (dir in directories) {
  files <- list.files(path = dir, pattern = "^QC_.*\\.csv$", full.names = TRUE)
  file_list <- c(file_list, files)
}

merged_data <- file_list %>%
  lapply(read_csv) %>%
  bind_rows() %>%
  distinct()


merged_data<-merged_data %>%dplyr::rename(sample=1)


write_csv(merged_data, "all_samples_QC.csv")


# Load necessary libraries


# Select specific columns
QC_selected <-merged_data %>%
  select(sample, Total, Total.uniq, Total.uniq.est, X.Signal.Total, X.Background.Total,Global.signal.yield,Local.signal.yield,Global.SNR,Local.SNR,Seq.factor)

QC_selected <- QC_selected %>%
  dplyr::mutate(
    dataset = case_when(
      grepl("-C.*H3K27ac", sample) ~ "CSF H3K27ac Mattia",
      grepl("-P.*H3K27ac", sample) ~ "Plasma H3K27ac Mattia",
      grepl("-P.*H3K4me3", sample) ~ "Plasma H3K4me3 Mattia",
      grepl("^GSM", sample) ~ "H3K27ac baca",
      grepl("^H0", sample) ~ "H3K4me3 sadeh" # Keep the original value if none of the conditions match
    ))%>% 
  # Filter out rows with NA in the dataset column
  filter(!is.na(dataset))

# Rename the first four rows in the sample column
QC_selected$dataset[1:4]<-"H3K4me1 sadeh"

colors <- c(
  "CSF H3K27ac Mattia" = "#FF0000",
  "Plasma H3K27ac Mattia" = "#3399FF",
  "Plasma H3K4me3 Mattia" = "#33FF99",
  "H3K27ac baca" = "#FF9933",
  "H3K4me3 sadeh" = "#9933FF",
  "H3K4me1 sadeh" = "#FF33CC"
)

# Calculate mean values for each category
mean_values <- QC_selected %>%
  group_by(dataset) %>%
  summarize(mean_total = mean(Total))

p1 <- ggplot(QC_selected) +
  aes(x = dataset, y = Total, fill = dataset, color = dataset) +
  geom_violin(trim = TRUE, color = "#000000") +
  geom_boxplot(width = 0.1, color = "#000000", fill = "#ffffff", outlier.shape = NA) +
  scale_fill_manual(values = colors) +
  scale_color_manual(values = colors) +
  ggthemes::theme_base() +
  xlab("") +
  ylab("Number of reads") +
  theme(panel.border = element_rect()) +
  theme_classic() + theme(legend.position = "none") +
  theme(axis.text.x = element_text(size = 12, family = "sans", angle = 45, hjust = 1),
        axis.text.y = element_text(size = 12, family = "sans"),
        axis.title.y = element_text(size = 14, family = "sans"),
        axis.line = element_line(size = 1)) +
  geom_text(data = mean_values, aes(x = dataset, y = mean_total, label = round(mean_total, 2)), 
            color = "black", size = 5, vjust = -12)  # Adjusted vjust value


# Save the plot as a PDF
ggsave("/date/gcb/gcb_MZ/Analysis/QC/total_reads.pdf", plot = p1, device = "pdf", width = 10, height = 6)




# Calculate mean values for each category
mean_values <- QC_selected %>%
  group_by(dataset) %>%
  summarize(mean_total = mean(Total.uniq))


# Plot with default font family
p1 <- ggplot(QC_selected) +
  aes(x = dataset, y = Total.uniq, fill = dataset, color = dataset) +
  geom_violin(trim = TRUE, color = "#000000") +
  geom_boxplot(width = 0.1, color = "#000000", fill = "#ffffff", outlier.shape = NA) +
  scale_fill_manual(values = colors) +
  scale_color_manual(values = colors) +
  ggthemes::theme_base() +
  xlab("") +
  ylab("Number of unique reads") +
  theme(panel.border = element_rect()) +
  theme_classic() + theme(legend.position = "none") +
  theme(axis.text.x = element_text(size = 12, family = "sans", angle = 45, hjust = 1),
        axis.text.y = element_text(size = 12, family = "sans"),
        axis.title.y = element_text(size = 14, family = "sans"),
        axis.line = element_line(size = 1)) +
  geom_text(data = mean_values, aes(x = dataset, y = mean_total, label = round(mean_total, 2)), 
            color = "black", size = 5, vjust = -12)

# Save the plot as a PDF
ggsave("/date/gcb/gcb_MZ/Analysis/QC/total.uniq_reads.pdf", plot = p1, device = "pdf", width = 10, height = 6)


# Calculate mean values for each category
mean_values <- QC_selected %>%
  group_by(dataset) %>%
  summarize(mean_total = mean(Seq.factor))


# Plot with default font family
p1 <- ggplot(QC_selected) +
  aes(x = dataset, y = Seq.factor, fill = dataset, color = dataset) +
  geom_violin(trim = TRUE, color = "#000000") +
  geom_boxplot(width = 0.1, color = "#000000", fill = "#ffffff", outlier.shape = NA) +
  scale_fill_manual(values = colors) +
  scale_color_manual(values = colors) +
  ggthemes::theme_base() +
  xlab("") +
  ylab("Sequencing Factor") +
  theme(panel.border = element_rect()) +
  theme_classic() + theme(legend.position = "none") +
  theme(axis.text.x = element_text(size = 12, family = "sans", angle = 45, hjust = 1),
        axis.text.y = element_text(size = 12, family = "sans"),
        axis.title.y = element_text(size = 14, family = "sans"),
        axis.line = element_line(size = 1)) +
  geom_text(data = mean_values, aes(x = dataset, y = mean_total, label = round(mean_total, 2)), 
            color = "black", size = 5, vjust = -12)

# Save the plot as a PDF
ggsave("/date/gcb/gcb_MZ/Analysis/QC/Seq.factor.pdf", plot = p1, device = "pdf", width = 10, height = 6)



# Calculate mean values for each category
mean_values <- QC_selected %>%
  group_by(dataset) %>%
  summarize(mean_total = mean(X.Signal.Total))

# Plot with default font family
p1 <- ggplot(QC_selected) +
  aes(x = dataset, y = X.Signal.Total, fill = dataset, color = dataset) +
  geom_violin(trim = TRUE, color = "#000000") +
  geom_boxplot(width = 0.1, color = "#000000", fill = "#ffffff", outlier.shape = NA) +
  scale_fill_manual(values = colors) +
  scale_color_manual(values = colors) +
  ggthemes::theme_base() +
  xlab("") +
  ylab("% of signal in signal windows") +
  theme(panel.border = element_rect()) +
  theme_classic() + theme(legend.position = "none") +
  theme(axis.text.x = element_text(size = 12, family = "sans", angle = 45, hjust = 1),
        axis.text.y = element_text(size = 12, family = "sans"),
        axis.title.y = element_text(size = 14, family = "sans"),
        axis.line = element_line(size = 1)) +
  geom_text(data = mean_values, aes(x = dataset, y = mean_total, label = round(mean_total, 2)), 
            color = "black", size = 5, vjust =-3)
p1

# Save the plot as a PDF
ggsave("/date/gcb/gcb_MZ/Analysis/QC/Signal_percentage.pdf", plot = p1, device = "pdf", width = 10, height = 6)


# Calculate mean values for each category
mean_values <- QC_selected %>%
  group_by(dataset) %>%
  summarize(mean_total = mean(X.Background.Total))

# Plot with default font family
p1 <- ggplot(QC_selected) +
  aes(x = dataset, y =X.Background.Total, fill = dataset, color = dataset) +
  geom_violin(trim = TRUE, color = "#000000") +
  geom_boxplot(width = 0.1, color = "#000000", fill = "#ffffff", outlier.shape = NA) +
  scale_fill_manual(values = colors) +
  scale_color_manual(values = colors) +
  ggthemes::theme_base() +
  xlab("") +
  ylab("% of signal in background") +
  theme(panel.border = element_rect()) +
  theme_classic() + theme(legend.position = "none") +
  theme(axis.text.x = element_text(size = 12, family = "sans", angle = 45, hjust = 1),
        axis.text.y = element_text(size = 12, family = "sans"),
        axis.title.y = element_text(size = 14, family = "sans"),
        axis.line = element_line(size = 1)) +
  geom_text(data = mean_values, aes(x = dataset, y = mean_total, label = round(mean_total, 2)), 
            color = "black", size = 5, vjust = -3)
p1

# Save the plot as a PDF
ggsave("/date/gcb/gcb_MZ/Analysis/QC/background_percentage.pdf", plot = p1, device = "pdf", width = 10, height = 6)


# Calculate mean values for each category
mean_values <- QC_selected %>%
  group_by(dataset) %>%
  summarize(mean_total = mean(X.Background.Total))

# Plot with default font family
p1 <- ggplot(QC_selected) +
  aes(x = dataset, y =X.Background.Total, fill = dataset, color = dataset) +
  geom_violin(trim = TRUE, color = "#000000") +
  geom_boxplot(width = 0.1, color = "#000000", fill = "#ffffff", outlier.shape = NA) +
  scale_fill_manual(values = colors) +
  scale_color_manual(values = colors) +
  ggthemes::theme_base() +
  xlab("") +
  ylab("% of signal in background") +
  theme(panel.border = element_rect()) +
  theme_classic() + theme(legend.position = "none") +
  theme(axis.text.x = element_text(size = 12, family = "sans", angle = 45, hjust = 1),
        axis.text.y = element_text(size = 12, family = "sans"),
        axis.title.y = element_text(size = 14, family = "sans"),
        axis.line = element_line(size = 1)) +
  geom_text(data = mean_values, aes(x = dataset, y = mean_total, label = round(mean_total, 2)), 
            color = "black", size = 5, vjust = -0.5)
p1

# Save the plot as a PDF
ggsave("/date/gcb/gcb_MZ/Analysis/QC/background_percentage.pdf", plot = p1, device = "pdf", width = 10, height = 6)


# Calculate mean values for each category
mean_values <- QC_selected %>%
  group_by(dataset) %>%
  summarize(mean_total = mean(Global.signal.yield))

# Plot with default font family
p1 <- ggplot(QC_selected) +
  aes(x = dataset, y =Global.signal.yield, fill = dataset, color = dataset) +
  geom_violin(trim = TRUE, color = "#000000") +
  geom_boxplot(width = 0.1, color = "#000000", fill = "#ffffff", outlier.shape = NA) +
  scale_fill_manual(values = colors) +
  scale_color_manual(values = colors) +
  ggthemes::theme_base() +
  xlab("") +
  ylab("Global signal yield") +
  theme(panel.border = element_rect()) +
  theme_classic() + theme(legend.position = "none") +
  theme(axis.text.x = element_text(size = 12, family = "sans", angle = 45, hjust = 1),
        axis.text.y = element_text(size = 12, family = "sans"),
        axis.title.y = element_text(size = 14, family = "sans"),
        axis.line = element_line(size = 1)) #+
  #geom_text(data = mean_values, aes(x = dataset, y = mean_total, label = round(mean_total, 2)), 
            #color = "black", size = 5, vjust = -0.5)
p1

# Save the plot as a PDF
ggsave("/date/gcb/gcb_MZ/Analysis/QC/Global.signal.yield.pdf", plot = p1, device = "pdf", width = 10, height = 6)



# Calculate mean values for each category
mean_values <- QC_selected %>%
  group_by(dataset) %>%
  summarize(mean_total = mean(Local.signal.yield))

# Plot with default font family
p1 <- ggplot(QC_selected) +
  aes(x = dataset, y =Local.signal.yield, fill = dataset, color = dataset) +
  geom_violin(trim = TRUE, color = "#000000") +
  geom_boxplot(width = 0.1, color = "#000000", fill = "#ffffff", outlier.shape = NA) +
  scale_fill_manual(values = colors) +
  scale_color_manual(values = colors) +
  ggthemes::theme_base() +
  xlab("") +
  ylab("Local signal yield") +
  theme(panel.border = element_rect()) +
  theme_classic() + theme(legend.position = "none") +
  theme(axis.text.x = element_text(size = 12, family = "sans", angle = 45, hjust = 1),
        axis.text.y = element_text(size = 12, family = "sans"),
        axis.title.y = element_text(size = 14, family = "sans"),
        axis.line = element_line(size = 1)) #+
#geom_text(data = mean_values, aes(x = dataset, y = mean_total, label = round(mean_total, 2)), 
#color = "black", size = 5, vjust = -0.5)
p1

# Save the plot as a PDF
ggsave("/date/gcb/gcb_MZ/Analysis/QC/Local.signal.yield.pdf", plot = p1, device = "pdf", width = 10, height = 6)


# Calculate mean values for each category
mean_values <- QC_selected %>%
  group_by(dataset) %>%
  summarize(mean_total = mean(Global.SNR))

# Plot with default font family
p1 <- ggplot(QC_selected) +
  aes(x = dataset, y =Global.SNR, fill = dataset, color = dataset) +
  geom_violin(trim = TRUE, color = "#000000") +
  geom_boxplot(width = 0.1, color = "#000000", fill = "#ffffff", outlier.shape = NA) +
  scale_fill_manual(values = colors) +
  scale_color_manual(values = colors) +
  ggthemes::theme_base() +
  xlab("") +
  ylab("Global Signal to noise Ratio") +
  theme(panel.border = element_rect()) +
  theme_classic() + theme(legend.position = "none") +
  theme(axis.text.x = element_text(size = 12, family = "sans", angle = 45, hjust = 1),
        axis.text.y = element_text(size = 12, family = "sans"),
        axis.title.y = element_text(size = 14, family = "sans"),
        axis.line = element_line(size = 1)) +
geom_text(data = mean_values, aes(x = dataset, y = mean_total, label = round(mean_total, 2)), 
color = "black", size = 5, vjust = -2)
p1

# Save the plot as a PDF
ggsave("/date/gcb/gcb_MZ/Analysis/QC/Global.SNR.pdf", plot = p1, device = "pdf", width = 10, height = 6)


# Calculate mean values for each category
mean_values <- QC_selected %>%
  group_by(dataset) %>%
  summarize(mean_total = mean(Local.SNR))

# Plot with default font family
p1 <- ggplot(QC_selected) +
  aes(x = dataset, y =Local.SNR, fill = dataset, color = dataset) +
  geom_violin(trim = TRUE, color = "#000000") +
  geom_boxplot(width = 0.1, color = "#000000", fill = "#ffffff", outlier.shape = NA) +
  scale_fill_manual(values = colors) +
  scale_color_manual(values = colors) +
  ggthemes::theme_base() +
  xlab("") +
  ylab("Local Signal to noise Ratio") +
  theme(panel.border = element_rect()) +
  theme_classic() + theme(legend.position = "none") +
  theme(axis.text.x = element_text(size = 12, family = "sans", angle = 45, hjust = 1),
        axis.text.y = element_text(size = 12, family = "sans"),
        axis.title.y = element_text(size = 14, family = "sans"),
        axis.line = element_line(size = 1)) +
  geom_text(data = mean_values, aes(x = dataset, y = mean_total, label = round(mean_total, 2)), 
            color = "black", size = 5, vjust = -2)
p1

# Save the plot as a PDF
ggsave("/date/gcb/gcb_MZ/Analysis/QC/Local.SNR.pdf", plot = p1, device = "pdf", width = 10, height = 6)
