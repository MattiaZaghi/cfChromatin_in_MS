library(preprocessCore)
library(MASS)
library(parallel)
library(DESeq2)
library(dplyr)
library(biomaRt)
library(RColorBrewer)
library(pheatmap)
library(ggplot2)
library(ggrepel)
library(dplyr)
library(DT)
library(scales)
library(stringr)
library(ggpubr)
library(grid)
library(ggplotify)
library(plotrix)
library(GEOquery)
library(Biobase)
library(RColorBrewer)
library(colorspace)
library(EnhancedVolcano)
library(LSD)
library(graphics)
library(itsadug)
library(viridis)
library(M3C)
library(xlsx)
library(Rtsne)
library(affy)
library(gplots)
library(readxl)
library(VennDiagram)
library(pals)
library(gprofiler2)
library(openxlsx)
library(ggvenn)
library(tidyr)
library(data.table)
library(purrr)
library(Polychrome)
print("Setting up directories and loading utility functions...")
SourceDIR<-"/proj/user/mattia/Analysis/cfChIP-seq/"
SetupDIR<-"/proj/user/mattia/Analysis/cfChIP-seq/SetupFiles/H3K27ac_hg38/"
TSS.windows<-readRDS("/proj/user/mattia/Analysis/cfChIP-seq/SetupFiles/H3K27ac_hg38/Windows.rds")
TargetMod<-"H3K27ac_hg38"
source(paste0(SourceDIR, "Calculate_HealtyRef_fun.R")) # Utility functions for non-negative matrix factorization.
data_dir <- "/proj/user/mattia/Analysis/Samples/H3K27ac_hg38/"

print("Loading .rdata files...")
rdata_files <- list.files(path = data_dir, pattern = "\\.rdata$", full.names = TRUE)
data_list <- list()
for (rdata_file in rdata_files) {
  var_name <- tools::file_path_sans_ext(basename(rdata_file))
  loaded_data <- readRDS(rdata_file)
  data_list[[var_name]] <- loaded_data
}
# Your vector
Healthy <- c("H7-P-Ctrl_H3K27ac_ChIP","H7_H3K27ac_nanoCT",
             "H11-P-Ctrl_H3K27ac_ChIP","H11_H3K27ac_S3-nanoCT",
              "H15-P-Ctrl_H3K27ac_ChIP", "H15_H3K27ac_nanoCT",
              "H15_H3K27ac_S3-nanoCT",
              "H21-P-Ctrl_H3K27ac_ChIP", "H21_H3K27ac_S3-nanoCT"
)


# Get the names of the tissues that are in your datasets vector
data_list <- data_list[names(data_list) %in% Healthy]

# Print single_cells
print(names(data_list))




# Iterate over each element in data_list
for (name in names(data_list)) {
  # Calculate GeneCounts.Background
  gene_counts_background <- data_list[[name]][["GeneCounts"]] - data_list[[name]][["GeneBackground"]]
  
  # Replace negative values with 0
  gene_counts_background[gene_counts_background < 0] <- 0
  
  # Round the values to the nearest integer
  gene_counts_background <- round(gene_counts_background)
  
  # Assign the result back to the data_list
  data_list[[name]][["GeneCounts.Background"]] <- gene_counts_background
}

# Step 1: Extract Gene Counts and Filter Out Unwanted Samples
unwanted_samples <- c("P20015-P-Ctrl_H3K27ac_ChIP", "P20027-P-Ctrl_H3K27ac_ChIP", "P20030-P-Ctrl_H3K27ac_ChIP", "P20040-P-Ctrl_H3K27ac_ChIP",
                      "14131-P-MS-Rituximab-Stable_H3K27ac_ChIP","18070-P-MS-Rituximab-Progressive_H3K27ac_ChIP",
                      "12179-P-MS-Rituximab-Progressive_H3K27ac_ChIP","18075-P-MS-Rituximab-Progressive_H3K27ac_ChIP",
                      "14-229-P-MS-Rituximab-Stable_H3K27ac_ChIP")
filtered_data_list <- data_list[names(data_list) %in% Healthy]

gene_counts_list <- lapply(data_list ,function(x) x[["GeneCounts.Background"]])


# Step 2: Create Sample Names and Groups
sample_names <- names(data_list)
short_sample_names <- sapply(sample_names, function(x) {
  parts <- strsplit(x, "")[[1]]
  if (grepl("H", x)) {
    return(paste0(parts[1], ""))
  } else {
    return(parts[1])
  }
})
groups <- sapply(sample_names, function(x) {
  if (grepl("S3-nanoCT", x)) {
    return("S3")
  } else if (grepl("nanoCT", x)) {
    return("nanoCT")
  }  else if (grepl("ChIP", x)) {
    return("ChIP")
  }
})
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
# Step 3: Combine Gene Counts into a Matrix
cfChromatin_GeneCounts <- do.call(cbind, gene_counts_list)
colnames(cfChromatin_GeneCounts) <- sample_names

# Step 4: Create Metadata
metaData <- data.frame(Group = groups, replicates = sample_names)

# Ensure row names of countData are set correctly
rownames(cfChromatin_GeneCounts) <- rownames(cfChromatin_GeneCounts)

# Ensure column names of countData match row names of colData
colnames(cfChromatin_GeneCounts) <- rownames(metaData)

# Ensure row names of colData are set correctly
rownames(metaData) <- colnames(cfChromatin_GeneCounts)

cfChromatin_GeneCounts <- cfChromatin_GeneCounts[rowSums(cfChromatin_GeneCounts)>0,] #take out genes not expressed
cfChromatin_GeneCounts[cfChromatin_GeneCounts == 0] <- 1

# Step 5: Create DESeqDataSet
cfChromatin_GeneCountsDDS <- DESeqDataSetFromMatrix(countData = cfChromatin_GeneCounts, colData = metaData, design = ~ Group)
# removing rows of the DESeqDataSet that have no counts, or only a single count across all samples. 
cfChromatin_GeneCountsDDS <-cfChromatin_GeneCountsDDS[rowSums(counts(cfChromatin_GeneCountsDDS)) > 1, ]

#plot PCA
cfChromatin_GeneCountsDDS <-estimateSizeFactors(cfChromatin_GeneCountsDDS)
readcounts_norm <- as.data.frame(counts(cfChromatin_GeneCountsDDS, normalized = TRUE))
readcounts_norm$Gene_name <- rownames(readcounts_norm)
sizeFactors(cfChromatin_GeneCountsDDS)

se <- SummarizedExperiment(log2(counts(cfChromatin_GeneCountsDDS, normalized=TRUE)+1),
                           colData=colData(cfChromatin_GeneCountsDDS))

pca <- DESeqTransform( se )

p <- plotPCA(pca, intgroup = "Group") + 
  geom_point(aes(colour = metaData$Group), size = 0.1) +
  geom_jitter(width = 0.25) +
  #scale_fill_manual(values=okabe(8)) +
  #scale_color_manual(values=okabe(8)) +
  theme_classic() +
  geom_point(size=0.1) +
  theme(axis.text.x = element_text(size = 12, family = "Arial"),
        axis.text.y = element_text(size = 12, family = "Arial"),
        axis.title.y = element_text(size = 12, family = "Arial"),
        axis.title.x = element_text(size = 12, family = "Arial"),
        legend.text = element_text(size = 12, family = "Arial"),
        axis.line = element_line(size = 0.5))

# Save the plot with better parameters
ggsave(filename = "/proj/user/mattia/Analysis/Output/H3K27ac_hg38/PCA-nanoCT.png", 
       plot = p, 
       device = "png", 
       width = 6, 
       height = 5.75, 
       units = "in", 
       dpi = 300, 
       limitsize = TRUE)

vsd <- vst(cfChromatin_GeneCountsDDS, blind=F)
rld <- rlog(cfChromatin_GeneCountsDDS, blind=T)


### Sample to sample distance.
sampleDists <- dist(t(assay(rld)))
sampleDistMatrix <- as.matrix(sampleDists)
rownames(sampleDistMatrix) <- rownames(sampleDistMatrix)
colnames(sampleDistMatrix) <- rownames(sampleDistMatrix)



suppressMessages(library(stringr))

metadata = data.frame(samples = metaData$replicates,
                      condition = metaData$Group,
                      row.names = names(raw),
                      stringsAsFactors = T)

suppressMessages(library(RColorBrewer))
suppressMessages(library("viridis"))
library(pals)
library(pheatmap)
annotation_column <- metadata[,1:(dim(metadata)[2])]
mycolors_s <- c(polychrome(36), "#000000"); names(mycolors_s) = levels(annotation_column$samples)
mycolors_c <- kelly(22)[3:6]; names(mycolors_c) = levels(annotation_column$condition)
ann_colors = list(samples = mycolors_s, condition=mycolors_c)
colors <- colorRampPalette(rev(brewer.pal(9, "Blues")))(255)


sample_distace_plot <- pheatmap(sampleDistMatrix,
                                clustering_distance_rows=sampleDists,
                                clustering_distance_cols=sampleDists, border_color = "white",
                                col=colors, cellwidth = 20, cellheight = 20, annotation_col = annotation_column,
                                annotation_colors = ann_colors, filename = "/proj/user/mattia/Analysis/Output/H3K27ac_hg38/S_distance_Gene_nano.png") 


dds <- DESeq(cfChromatin_GeneCountsDDS)


#healthy vs New Diagnosis

## define contrasts (careful here!!) Control FIB
contrast=c("Ctrl","New-RR") # ratio numerator/denumerator is expressed always as c("numerator", "denumerator")
phenoData<-colData(dds)
phenoData<-phenoData[phenoData$sampletype %in% contrast,]
contast_samples_ordered=rownames(metaData[order(metaData$Group),])
fcT=1
fdr=0.99
res <- results(dds, contrast=c("Group",contrast), alpha = fdr, pAdjustMethod="BH")

sink("res_Contrast.txt")
summary(res)
sink()
res <- res %>% as.data.frame()
res$Gene_name <- rownames(res)

res <- left_join(res, readcounts_norm, "Gene_name")
rownames(res) <- res$Gene_name

# First create your MA plot
ma_plot<-ggmaplot(res, main = "DEGS Ctrl vs New Diagnosis",
         fdr = 0.05, fc = 0.5, size = 1,
         palette = c("#B31B21", "#1465AC", "gray"),
         genenames = as.vector(rownames(res)),
         legend = "top",
         font.label = c("bold", 10),
         font.legend = "bold",
         font.main = c("bold", 14), label.rectangle = TRUE,
         ggtheme = ggplot2::theme_minimal(),top = 0)


# Save the plot with optimized settings
ggsave(
  filename = "/proj/user/mattia/Analysis/Output/H3K27ac_hg38/MA_plot_Ctrl_vs_New.png",
  plot = ma_plot,  # Explicitly specify the plot object
  width = 8,       # Slightly wider for better readability
  height = 7,      # Keep height at 7
  dpi = 330,       # Good resolution for publication
  bg = "white",    # Ensure white background
  device = "png"   # Explicitly specify device
)
res <- res %>% filter(padj < 0.05) 
res$Gene_name<-rownames(res)
UP <- filter(res, res$log2FoldChange > 0)[,7]
DOWN <- filter(res, res$log2FoldChange < 0)[,7]

library(gprofiler2)
library(ggnewscale)

GO_list <- list(`Upregulated in New MS`=UP,
                `Upregulated in Ctrl`=DOWN)
GO_out <- list()
for (i in 1:2) {
  gostres <- gost(query = GO_list[[i]],
                  organism = "hsapiens",
                  evcodes = TRUE,
                  significant = TRUE,
                  correction_method = "fdr",
                  user_threshold = 0.05 , sources = c("GO:BP"))
  
  result <- as.data.frame(gostres$result)
  GO_out[[i]] <- result
  names(GO_out)[i] <- names(GO_list)[i]
}

for (i in 1:2) {
  GO_out[[i]]$Condition <- names(GO_out)[i]
}

GO_all <- bind_rows(GO_out)
GO_all <- GO_all[c(3, 4,6, 9:11,16:17)]
GO_all$`Percentage of enrichment` <- GO_all$intersection_size / GO_all$term_size *100
GO_all$`-log10 Pvalue` <- -log10(GO_all$p_value) 
write.xlsx(GO_all,"GO_Astro_RNA_seq.xlsx")
#healthy vs stable

## define contrasts (careful here!!) Control FIB
contrast=c("Ctrl","Stable") # ratio numerator/denumerator is expressed always as c("numerator", "denumerator")
phenoData<-colData(dds)
phenoData<-phenoData[phenoData$sampletype %in% contrast,]
contast_samples_ordered=rownames(metaData[order(metaData$Group),])
fcT=1
fdr=0.99
res <- results(dds, contrast=c("Group",contrast), alpha = fdr, pAdjustMethod="BH")

sink("res_Contrast.txt")
summary(res)
sink()
res <- res %>% as.data.frame()
res$Gene_name <- rownames(res)
res <- res %>% filter(padj < 0.05) 
res <- left_join(res, readcounts_norm, "Gene_name")
rownames(res) <- res$Gene_name

# First create your MA plot
ma_plot<-ggmaplot(res, main = "DEGS Ctrl vs Stable Rituximab",
                  fdr = 0.05, fc = 0.5, size = 1,
                  palette = c("#B31B21", "#1465AC", "gray"),
                  genenames = as.vector(rownames(res)),
                  legend = "top",
                  font.label = c("bold", 10),
                  font.legend = "bold",
                  font.main = c("bold", 14), label.rectangle = TRUE,
                  ggtheme = ggplot2::theme_minimal(),top = 0)


# Save the plot with optimized settings
ggsave(
  filename = "/proj/user/mattia/Analysis/Output/H3K27ac_hg38/MA_plot_Ctrl_vs_Stable.png",
  plot = ma_plot,  # Explicitly specify the plot object
  width = 8,       # Slightly wider for better readability
  height = 7,      # Keep height at 7
  dpi = 330,       # Good resolution for publication
  bg = "white",    # Ensure white background
  device = "png"   # Explicitly specify device
)

UP <- filter(res, res$log2FoldChange > 0.5)[,7]
DOWN <- filter(res, res$log2FoldChange < -0.5)[,7]

library(gprofiler2)
library(ggnewscale)

GO_list <- list(`Upregulated in Stable Rituximab`=UP,
                `Upregulated in Ctrl`=DOWN)
GO_out <- list()
for (i in 1:2) {
  gostres <- gost(query = GO_list[[i]],
                  organism = "hsapiens",
                  evcodes = TRUE,
                  significant = TRUE,
                  correction_method = "fdr",
                  user_threshold = 0.05 , sources = c("GO:BP"))
  
  result <- as.data.frame(gostres$result)
  GO_out[[i]] <- result
  names(GO_out)[i] <- names(GO_list)[i]
}

for (i in 1:2) {
  GO_out[[i]]$Condition <- names(GO_out)[i]
}

GO_all <- bind_rows(GO_out)

#healthy vs Progressive

## define contrasts (careful here!!) Control FIB
contrast=c("Ctrl","Progressive") # ratio numerator/denumerator is expressed always as c("numerator", "denumerator")
phenoData<-colData(dds)
phenoData<-phenoData[phenoData$sampletype %in% contrast,]
contast_samples_ordered=rownames(metaData[order(metaData$Group),])
fcT=1
fdr=0.99
res <- results(dds, contrast=c("Group",contrast), alpha = fdr, pAdjustMethod="BH")

sink("res_Contrast.txt")
summary(res)
sink()
res <- res %>% as.data.frame()
res$Gene_name <- rownames(res)
res <- res %>% filter(padj < 0.05) 
res <- left_join(res, readcounts_norm, "Gene_name")
rownames(res) <- res$Gene_name

# First create your MA plot
ma_plot<-ggmaplot(res, main = "DEGS Ctrl vs Progressive Rituximab",
                  fdr = 0.05, fc = 0.5, size = 1,
                  palette = c("#B31B21", "#1465AC", "gray"),
                  genenames = as.vector(rownames(res)),
                  legend = "top",
                  font.label = c("bold", 10),
                  font.legend = "bold",
                  font.main = c("bold", 14), label.rectangle = TRUE,
                  ggtheme = ggplot2::theme_minimal(),top = 0)


# Save the plot with optimized settings
ggsave(
  filename = "/proj/user/mattia/Analysis/Output/H3K27ac_hg38/MA_plot_Ctrl_vs_Progressive.png",
  plot = ma_plot,  # Explicitly specify the plot object
  width = 8,       # Slightly wider for better readability
  height = 7,      # Keep height at 7
  dpi = 330,       # Good resolution for publication
  bg = "white",    # Ensure white background
  device = "png"   # Explicitly specify device
)



#healthy vs Progressive

## define contrasts (careful here!!) Control FIB
contrast=c("Ctrl","Progressive") # ratio numerator/denumerator is expressed always as c("numerator", "denumerator")
phenoData<-colData(dds)
phenoData<-phenoData[phenoData$sampletype %in% contrast,]
contast_samples_ordered=rownames(metaData[order(metaData$Group),])
fcT=1
fdr=0.99
res <- results(dds, contrast=c("Group",contrast), alpha = fdr, pAdjustMethod="BH")

sink("res_Contrast.txt")
summary(res)
sink()
res <- res %>% as.data.frame()
res$Gene_name <- rownames(res)
res <- res %>% filter(padj < 0.05) 
res <- left_join(res, readcounts_norm, "Gene_name")
rownames(res) <- res$Gene_name

# First create your MA plot
ma_plot<-ggmaplot(res, main = "DEGS Ctrl vs New",
                  fdr = 0.05, fc = 0.5, size = 1,
                  palette = c("#B31B21", "#1465AC", "gray"),
                  genenames = as.vector(rownames(res)),
                  legend = "top",
                  font.label = c("bold", 10),
                  font.legend = "bold",
                  font.main = c("bold", 14), label.rectangle = TRUE,
                  ggtheme = ggplot2::theme_minimal(),top = 0)


# Save the plot with optimized settings
ggsave(
  filename = "/proj/user/mattia/Analysis/Output/H3K27ac_hg38/MA_plot_Ctrl_vs_Progressive.png",
  plot = ma_plot,  # Explicitly specify the plot object
  width = 8,       # Slightly wider for better readability
  height = 7,      # Keep height at 7
  dpi = 330,       # Good resolution for publication
  bg = "white",    # Ensure white background
  device = "png"   # Explicitly specify device
)

#healthy vs Progressive

## define contrasts (careful here!!) Control FIB
contrast=c("Stable","Progressive") # ratio numerator/denumerator is expressed always as c("numerator", "denumerator")
phenoData<-colData(dds)
phenoData<-phenoData[phenoData$sampletype %in% contrast,]
contast_samples_ordered=rownames(metaData[order(metaData$Group),])
fcT=1
fdr=0.99
res <- results(dds, contrast=c("Group",contrast), alpha = fdr, pAdjustMethod="BH")

sink("res_Contrast.txt")
summary(res)
sink()
res <- res %>% as.data.frame()
res$Gene_name <- rownames(res)
res <- res %>% filter(padj < 0.05) 
res <- left_join(res, readcounts_norm, "Gene_name")
rownames(res) <- res$Gene_name

# First create your MA plot
ma_plot<-ggmaplot(res, main = "DEGS Stable Rituximab vs Progressive Rituximab",
                  fdr = 0.05, fc = 0.5, size = 1,
                  palette = c("#B31B21", "#1465AC", "gray"),
                  genenames = as.vector(rownames(res)),
                  legend = "top",
                  font.label = c("bold", 10),
                  font.legend = "bold",
                  font.main = c("bold", 14), label.rectangle = TRUE,
                  ggtheme = ggplot2::theme_minimal(),top = 0)


# Save the plot with optimized settings
ggsave(
  filename = "/proj/user/mattia/Analysis/Output/H3K27ac_hg38/MA_plot_Stable_vs_Progressive.png",
  plot = ma_plot,  # Explicitly specify the plot object
  width = 8,       # Slightly wider for better readability
  height = 7,      # Keep height at 7
  dpi = 330,       # Good resolution for publication
  bg = "white",    # Ensure white background
  device = "png"   # Explicitly specify device
)

## define contrasts (careful here!!) Control FIB
contrast=c("Stable","New-RR") # ratio numerator/denumerator is expressed always as c("numerator", "denumerator")
phenoData<-colData(dds)
phenoData<-phenoData[phenoData$sampletype %in% contrast,]
contast_samples_ordered=rownames(metaData[order(metaData$Group),])
fcT=1
fdr=0.99
res <- results(dds, contrast=c("Group",contrast), alpha = fdr, pAdjustMethod="BH")

sink("res_Contrast.txt")
summary(res)
sink()
res <- res %>% as.data.frame()
res$Gene_name <- rownames(res)
res <- res %>% filter(padj < 0.05) 
res <- left_join(res, readcounts_norm, "Gene_name")
rownames(res) <- res$Gene_name

# First create your MA plot
ma_plot<-ggmaplot(res, main = "DEGS Stable Rituximab vs New Diagnosis",
                  fdr = 0.05, fc = 0.5, size = 1,
                  palette = c("#B31B21", "#1465AC", "gray"),
                  genenames = as.vector(rownames(res)),
                  legend = "top",
                  font.label = c("bold", 10),
                  font.legend = "bold",
                  font.main = c("bold", 14), label.rectangle = TRUE,
                  ggtheme = ggplot2::theme_minimal(),top = 0)


# Save the plot with optimized settings
ggsave(
  filename = "/proj/user/mattia/Analysis/Output/H3K27ac_hg38/MA_plot_Stable_vs_New.png",
  plot = ma_plot,  # Explicitly specify the plot object
  width = 8,       # Slightly wider for better readability
  height = 7,      # Keep height at 7
  dpi = 330,       # Good resolution for publication
  bg = "white",    # Ensure white background
  device = "png"   # Explicitly specify device
)



## define contrasts (careful here!!) Control FIB
contrast=c("Progressive","New-RR") # ratio numerator/denumerator is expressed always as c("numerator", "denumerator")
phenoData<-colData(dds)
phenoData<-phenoData[phenoData$sampletype %in% contrast,]
contast_samples_ordered=rownames(metaData[order(metaData$Group),])
fcT=1
fdr=0.99
res <- results(dds, contrast=c("Group",contrast), alpha = fdr, pAdjustMethod="BH")

sink("res_Contrast.txt")
summary(res)
sink()
res <- res %>% as.data.frame()
res$Gene_name <- rownames(res)
res <- res %>% filter(padj < 0.05) 
res <- left_join(res, readcounts_norm, "Gene_name")
rownames(res) <- res$Gene_name

# First create your MA plot
ma_plot<-ggmaplot(res, main = "DEGS Progressive Rituximab vs New Diagnosis",
                  fdr = 0.05, fc = 0.5, size = 1,
                  palette = c("#B31B21", "#1465AC", "gray"),
                  genenames = as.vector(rownames(res)),
                  legend = "top",
                  font.label = c("bold", 10),
                  font.legend = "bold",
                  font.main = c("bold", 14), label.rectangle = TRUE,
                  ggtheme = ggplot2::theme_minimal(),top = 0)


# Save the plot with optimized settings
ggsave(
  filename = "/proj/user/mattia/Analysis/Output/H3K27ac_hg38/MA_plot_Progressive_vs_New.png",
  plot = ma_plot,  # Explicitly specify the plot object
  width = 8,       # Slightly wider for better readability
  height = 7,      # Keep height at 7
  dpi = 330,       # Good resolution for publication
  bg = "white",    # Ensure white background
  device = "png"   # Explicitly specify device
)







