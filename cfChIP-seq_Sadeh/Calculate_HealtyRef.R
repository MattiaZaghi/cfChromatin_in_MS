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


# Define your healthy reference samples
Healthy <- c(
  "GSM7787973_HP030132_H3K27Ac", "GSM7787975_HP030642_H3K27Ac", 
  "GSM7787978_HP031645_H3K27Ac", "GSM7787980_HP034881_H3K27Ac", 
  "GSM7787982_HP035094_H3K27Ac", "GSM7787985_HP038748_H3K27Ac", 
  "GSM7787987_HP041556_H3K27Ac", "GSM7787993_HP056703_H3K27Ac", 
  "GSM7788000_HP098228_H3K27Ac"
)


print("Loading only required .rdata files...")
# Get all .rdata files in the directory
all_rdata_files <- list.files(path = data_dir, pattern = "\\.rdata$", full.names = TRUE)

# Filter files to only those matching your healthy samples
target_files <- character(0)
for (healthy_sample in Healthy) {
  # Look for files that contain the healthy sample name
  matching_files <- all_rdata_files[grepl(paste0(healthy_sample, "\\.rdata$"), all_rdata_files)]
  target_files <- c(target_files, matching_files)
}

# Load only the filtered files
data_list <- list()
for (rdata_file in target_files) {
  var_name <- tools::file_path_sans_ext(basename(rdata_file))
  # Only load if the variable name is in your Healthy vector
  if (var_name %in% Healthy) {
    loaded_data <- readRDS(rdata_file)
    data_list[[var_name]] <- loaded_data
    print(paste("Loaded:", var_name))
  }
}

# Print loaded samples
print("Loaded samples:")
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

