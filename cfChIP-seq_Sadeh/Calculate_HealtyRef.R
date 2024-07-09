library(preprocessCore)
library(MASS)
SetupDIR<-"/date/gcb/gcb_MZ/Analysis/cfChIP-seq/SetupFiles/H3K27ac_hg38/"
TSS.windows<-readRDS("/date/gcb/gcb_MZ/Analysis/cfChIP-seq/SetupFiles/H3K27ac_hg38/Windows_TSS_enh_bastien.rds")
TargetMod<-"H3K27ac_hg38"
source("/date/gcb/gcb_MZ/Analysis/cfChIP-seq/Background.R")
# Specify the directory containing the .rdata files
data_dir <- "/date/gcb/gcb_MZ/Analysis/Samples/H3K27ac_hg38/"

# Get a list of all .rdata files in the directory
rdata_files <- list.files(path = data_dir, pattern = "\\.rdata$", full.names = TRUE)

# Initialize an empty list to store the loaded data
data_list <- list()

# Load each .rdata file and add it to the list
for (rdata_file in rdata_files) {
  # Get the base name of the file (without the extension)
  var_name <- tools::file_path_sans_ext(basename(rdata_file))
  
  # Load the .rdata file
  loaded_data <- readRDS(rdata_file)
  
  # Add the loaded data to the list with the same name
  data_list[[var_name]] <- loaded_data
}

# Calculate the signal-to-background ratio for each gene in each sample
ratios_list<-list()
# Calculate the signal-to-background ratio for each gene in each sample
for (i in seq_along(data_list)) {
  sample <- data_list[[i]]
  gene_names <- names(sample$GeneCounts)
  sample$GeneDiff = pmax(sample$GeneCounts - sample$GeneBackground, 0)
  data_list[[i]] <- sample
  # Create a data frame with the sample identifier, gene names, and ratios
  sample_df <- data.frame(sample = names(data_list)[i], gene = gene_names, Ratio = sample$GeneDiff)
  
  ratios_list[[i]] <- sample_df
}

# Combine all the ratios into a single data frame
all_ratios <- do.call(rbind, ratios_list)

# Exclude genes that start with 'ENST' or 'MIR'
all_ratios <- all_ratios[!grepl("^ENST|^MIR|^UNKNOWN", all_ratios$gene), ]

# Calculate the mean ratio for each gene across all samples
mean_ratios <- aggregate(Ratio ~ gene, all_ratios, mean)

# Sort the genes by their median ratio in descending order
sorted_genes <- mean_ratios[order(-mean_ratios$Ratio), ]

# Define the threshold as the value at the top 25% (75th percentile)
threshold <- quantile(sorted_genes$Ratio, 0.75)

# The housekeeping genes are the ones with a ratio above the threshold
candidate_housekeeping_genes <- sorted_genes$gene[sorted_genes$Ratio > threshold]

# Calculate the standard deviation of the ratio for each candidate housekeeping gene
sd_ratios <- aggregate(Ratio ~ gene, all_ratios[all_ratios$gene %in% candidate_housekeeping_genes, ], sd)

# Define the threshold as the median standard deviation
sd_threshold <- median(sd_ratios$Ratio)

# The final housekeeping genes are the ones with a standard deviation below the threshold
final_housekeeping_genes <- sd_ratios$gene[sd_ratios$Ratio < sd_threshold]



# Now, data_list is a list where each element is the data loaded from an .rdata file,
# and the name of each element is the base name of the corresponding .rdata file.

saveRDS(final_housekeeping_genes,"/date/gcb/gcb_MZ/Analysis/cfChIP-seq/SetupFiles/H3K27ac_hg38/CommonGenes.rds")
# Assuming you have a data frame 'df' with columns 'window', 'sample', 'read_counts', 'background_rate', 'gene'



# Assuming data_list is a list of matrices where each matrix represents a sample
# and final_housekeeping_genes is a vector of common genes

# First, create a matrix of gene counts for the common genes across all samples
GeneCounts.CommonGenes <- do.call(cbind, lapply(data_list, function(sample) sample$GeneDiff[names(sample$GeneDiff) %in% final_housekeeping_genes]))

# Apply quantile normalization across all samples at once
GeneCounts.CommonGenes.qq <- normalize.quantiles(GeneCounts.CommonGenes)

# Process each sample in the data_list
for (i in seq_along(data_list)) {
  # Get the sample data
  sample <- data_list[[i]]
  
  
  # Get the original and quantile-normalized expression levels for the common genes
  X <- GeneCounts.CommonGenes[, i]
  Y <- GeneCounts.CommonGenes.qq[, i]
  
  # Fit a robust linear model, falling back to a standard linear model if it fails
  sample$normalization_factor <- tryCatch(coef(rlm(Y ~ X - 1)), error = function(e) coef(lm(Y ~ X - 1)))
  
  # Compute the normalized gene levels for each sample
  sample$GeneCounts.QQnorm <- sample$GeneDiff * sample$normalization_factor
  
  # Rescale the normalization factor so that the total normalized signal at the set of reference healthy samples will be a million on average
  sample$QQNorm <- sample$normalization_factor * (1e6 / sum(sample$GeneCounts.QQnorm * sample$normalization_factor,na.rm = TRUE))
  
  # Compute the normalized gene levels for each sample
  sample$GeneCounts.QQnorm <- sample$GeneDiff * sample$QQNorm
  
  sample$WinBackground = getMultiBackgroundEstimate(sample$Background, 1:length(TSS.windows))
  sample$Counts.QQnorm = pmax(sample$Counts - sample$WinBackground, 0) * sample$QQNorm
  
  # Update the sample data in the data_list
  data_list[[i]] <- sample
}


# Aggregate window counts and background data from all samples.
WinCounts = do.call("cbind", lapply(data_list, function(sample) sample$Counts))
WinBackground = do.call("cbind", lapply(data_list, function(sample) sample$WinBackground))

# Aggregate normalization factors from all samples.
QQnorm = sapply(data_list, function(sample) sample$QQNorm)
names(QQnorm) = colnames(WinCounts)

# Estimate mean and variance for window data.
win.est = cfChIP.EstimateMeanVarianceBasis(WinCounts, WinBackground, QQnorm)

# Aggregate gene counts and background data from all samples.
GeneCounts = do.call("cbind", lapply(data_list, function(sample) sample$GeneCounts))
GeneBackground = do.call("cbind", lapply(data_list, function(sample) sample$GeneBackground))

# Estimate mean and variance for gene data.
gene.est = cfChIP.EstimateMeanVarianceBasis(GeneCounts, GeneBackground, QQnorm)

# Compile consensus data into a list and save it as an RDS file.
consensus = list(Win.avg = win.est$avg, Win.var = win.est$var, Gene.avg = gene.est$avg, Gene.var = gene.est$var)
saveRDS(consensus,"/date/gcb/gcb_MZ/Analysis/cfChIP-seq/SetupFiles/H3K27ac_hg38/HealthyRef.rds")
