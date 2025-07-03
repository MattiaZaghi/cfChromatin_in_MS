library(preprocessCore)
library(MASS)
library(parallel)

print("Setting up directories and loading utility functions...")
SourceDIR<-"/home/mattia/cfChromatin_in_MS/cfChIP-seq_Sadeh/"
SetupDIR<-"/date/gcb/gcb_MZ/Analysis/cfChIP-seq/SetupFiles/H3K27ac/"
TSS.windows<-readRDS("/date/gcb/gcb_MZ/Analysis/cfChIP-seq/SetupFiles/H3K27ac/Windows.rds")
TargetMod<-"H3K27ac"
source(paste0(SourceDIR, "Calculate_HealtyRef_fun.R")) # Utility functions for non-negative matrix factorization.
data_dir <- "/date/gcb/gcb_MZ/Analysis/Samples/H3K27ac/"

print("Loading .rdata files...")
rdata_files <- list.files(path = data_dir, pattern = "\\.rdata$", full.names = TRUE)
data_list <- list()
for (rdata_file in rdata_files) {
  var_name <- tools::file_path_sans_ext(basename(rdata_file))
  loaded_data <- readRDS(rdata_file)
  data_list[[var_name]] <- loaded_data
}
# Your vector
Healthy <- c("H10-P-Ctrl_H3K27ac_ChIP-V3", "H19-P-Ctrl_H3K27ac_ChIP-V3", "H5-P-Ctrl_H3K27ac_ChIP-V2",
"H11-P-Ctrl_H3K27ac_ChIP-V3", "H1-P-Ctrl_H3K27ac_ChIP-V2", "H5-P-Ctrl_H3K27ac_ChIP-V3",
"H12-P-Ctrl_H3K27ac_ChIP-V3", "H20-P-Ctrl_H3K27ac_ChIP-V3", "H6-P-Ctrl_H3K27ac_ChIP-V2-1D",
"H13-P-Ctrl_H3K27ac_ChIP-V3", "H21-P-Ctrl_H3K27ac_ChIP-V3", "H6-P-Ctrl_H3K27ac_ChIP-V2",
"H14-P-Ctrl_H3K27ac_ChIP-V3", "H22-P-Ctrl_H3K27ac_ChIP-V3", "H6-P-Ctrl_H3K27ac_ChIP-V3",
"H15-P-Ctrl_H3K27ac_ChIP-V3", "H23-P-Ctrl_H3K27ac_ChIP-V3", "H7-P-Ctrl_H3K27ac_ChIP-V2",
"H16-P-Ctrl_H3K27ac_ChIP-V3", "H24-P-Ctrl_H3K27ac_ChIP-V3", "H7-P-Ctrl_H3K27ac_ChIP-V3",
"H17-P-Ctrl_H3K27ac_ChIP-V3", "H2-P-Ctrl_H3K27ac_ChIP-V2", "H8-P-Ctrl_H3K27ac_ChIP-V2-1D",
"H17-P-Ctrl_H3K4me3_ChIP-V3", "H3-P-Ctrl_H3K27ac_ChIP-V2-1D", "H8-P-Ctrl_H3K27ac_ChIP-V2",
"H18-P-Ctrl_H3K27ac_ChIP-V3", "H4-P-Ctrl_H3K27ac_ChIP-V2", "H9-P-Ctrl_H3K27ac_ChIP-V3"
 )


# Get the names of the tissues that are in your datasets vector
data_list <- data_list[names(data_list) %in% Healthy]

# Print single_cells
print(names(data_list))


print("Calculating signal-to-background ratio for each gene in each sample...")
ratios_list<-list()
for (i in seq_along(data_list)) {
  sample <- data_list[[i]]
  gene_names <- names(sample$GeneCounts)
  sample$GeneDiff = pmax(sample$GeneCounts - sample$GeneBackground, 0)
  data_list[[i]] <- sample
  sample_df <- data.frame(sample = names(data_list)[i], gene = gene_names, Ratio = sample$GeneDiff)
  ratios_list[[i]] <- sample_df
}
all_ratios <- do.call(rbind, ratios_list)
all_ratios <- all_ratios[!grepl("^ENST|^MIR|^UNKNOWN", all_ratios$gene), ]
mean_ratios <- aggregate(Ratio ~ gene, all_ratios, mean)
sorted_genes <- mean_ratios[order(-mean_ratios$Ratio), ]
threshold <- quantile(sorted_genes$Ratio, 0.75)
candidate_housekeeping_genes <- sorted_genes$gene[sorted_genes$Ratio > threshold]
sd_ratios <- aggregate(Ratio ~ gene, all_ratios[all_ratios$gene %in% candidate_housekeeping_genes, ], sd)
sd_threshold <- median(sd_ratios$Ratio)
final_housekeeping_genes <- sd_ratios$gene[sd_ratios$Ratio < sd_threshold]

print("Saving final housekeeping genes...")
saveRDS(final_housekeeping_genes,"/date/gcb/gcb_MZ/Analysis/cfChIP-seq/SetupFiles/H3K27ac/CommonGenes.rds")

print("Normalizing gene levels...")
GeneCounts.CommonGenes <- do.call(cbind, lapply(data_list, function(sample) sample$GeneDiff[names(sample$GeneDiff) %in% final_housekeeping_genes]))
GeneCounts.CommonGenes.qq <- normalize.quantiles(GeneCounts.CommonGenes)
for (i in seq_along(data_list)) {
  sample <- data_list[[i]]
  X <- GeneCounts.CommonGenes[, i]
  Y <- GeneCounts.CommonGenes.qq[, i]
  sample$normalization_factor <- tryCatch(coef(rlm(Y ~ X - 1)), error = function(e) coef(lm(Y ~ X - 1)))
  sample$GeneCounts.QQnorm <- sample$GeneDiff * sample$normalization_factor
  sample$QQNorm <- sample$normalization_factor * (1e6 / sum(sample$GeneCounts.QQnorm * sample$normalization_factor,na.rm = TRUE))
  sample$GeneCounts.QQnorm <- sample$GeneDiff * sample$QQNorm
  sample$WinBackground = getMultiBackgroundEstimate(sample$Background, 1:length(TSS.windows))
  sample$Counts.QQnorm = pmax(sample$Counts - sample$WinBackground, 0) * sample$QQNorm
  data_list[[i]] <- sample
}

print("Aggregating window counts and background data from all samples...")
WinCounts = do.call("cbind", lapply(data_list, function(sample) sample$Counts))
WinBackground = do.call("cbind", lapply(data_list, function(sample) sample$WinBackground))
QQnorm = sapply(data_list, function(sample) sample$QQNorm)
names(QQnorm) = colnames(WinCounts)

print("Estimating mean and variance for window data...")
win.est = cfChIP.EstimateMeanVarianceBasis(WinCounts, WinBackground, QQnorm)

print("Aggregating gene counts and background data from all samples...")
GeneCounts = do.call("cbind", lapply(data_list, function(sample) sample$GeneCounts))
GeneBackground = do.call("cbind", lapply(data_list, function(sample) sample$GeneBackground))

print("Estimating mean and variance for gene data...")
gene.est = cfChIP.EstimateMeanVarianceBasis(GeneCounts, GeneBackground, QQnorm)

print("Compiling consensus data and saving it as an RDS file...")
consensus = list(Win.avg = win.est$avg, Win.var = win.est$var, Gene.avg = gene.est$avg, Gene.var = gene.est$var)
saveRDS(consensus,"/date/gcb/gcb_MZ/Analysis/cfChIP-seq/SetupFiles/H3K27ac/HealthyRef.rds")

print("All tasks completed.")

