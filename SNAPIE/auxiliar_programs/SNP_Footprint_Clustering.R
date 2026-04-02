library(vcfR)
library(adegenet)
library(factoextra)

directory <- "."
vcf_files <- list.files(path = directory, pattern = "\\.(vcf|vcf\\.gz)$", full.names = TRUE, recursive = TRUE)

# Initialize an empty list to store genotype data and SNP positions for each sample
genotype_list <- list()
snp_positions <- list()

# Loop through each VCF file and extract genotype information along with SNP positions
for (vcf_file in vcf_files) {
  vcf <- read.vcfR(vcf_file)
  sample_name <- gsub(".vcf.gz", "", basename(vcf_file))  # Get sample name
  
  # Extract genotype matrix (GT field) and convert to numeric
  genotype_matrix <- extract.gt(vcf, element = "GT", as.numeric = TRUE)
  
  # Extract SNP positions (chromosome and position)
  snp_pos <- paste(vcf@fix[, "CHROM"], vcf@fix[, "POS"], sep = "_")
  
  # Add genotype data and SNP positions to the lists
  genotype_list[[sample_name]] <- genotype_matrix
  snp_positions[[sample_name]] <- snp_pos
}

# Find common SNP positions across all samples
common_snps <- Reduce(intersect, snp_positions)

# Align the genotype matrices based on common SNP positions
aligned_genotype_list <- lapply(names(genotype_list), function(sample_name) {
  sample_genotypes <- genotype_list[[sample_name]]
  sample_snps <- snp_positions[[sample_name]]
  
  # Subset the genotype matrix to only include common SNPs
  match_idx <- match(common_snps, sample_snps)
  aligned_genotypes <- sample_genotypes[match_idx, ]
  
  return(aligned_genotypes)
})

# Combine the aligned genotype matrices into a single matrix (samples as columns, SNPs as rows)
combined_genotype_matrix <- do.call(cbind, aligned_genotype_list)

colnames(combined_genotype_matrix) <- names(genotype_list)  # Sample names as columns
colnames(combined_genotype_matrix)<-gsub("_unique.sorted.dedup", "", colnames(combined_genotype_matrix))

# Transpose the matrix to get samples as rows and SNPs as columns
genotype_matrix_samples <- t(combined_genotype_matrix)

# Replace any missing values with 0 or an appropriate value
genotype_matrix_samples[is.na(genotype_matrix_samples)] <- 0

pca_result <- prcomp(genotype_matrix_samples)

qty_samples <- length(pca_result$x[, 1])#

if (qty_samples < 3) {
  # Cria um PDF com a mensagem
  #pdf("Dendrogram_of_Samples_by_SNP_Profile.pdf", width = 4, height = 3)
  jpeg("Dendrogram_of_Samples_by_SNP_Profile.jpg", width = 800, height = 600, units = "px", quality = 100)
  
  # Adiciona a mensagem ao PDF
  plot.new()  # Cria uma nova tela de plot
  text(0.5, 0.5, "Not enough samples for clustering", cex = 2, col = "black", font = 2)
  
  # Fecha o dispositivo PDF
  dev.off()
} else{

  min_value <- min(qty_samples, 10) # Use first 10 PCs for clustering
  pc_matrix <- pca_result$x[, 1:min_value]  

  # Compute the distance matrix on the PCs
  dist_matrix <- dist(pc_matrix)  

  # Perform hierarchical clustering
  hclust_result <- hclust(dist_matrix, method = "ward.D2")

 # Get sample names

  # 4. Plot the dendrogram with Individual_ID labels
  #pdf("Dendrogram_of_Samples_by_SNP_Profile.pdf", width = 16, height = 10)
  jpeg("Dendrogram_of_Samples_by_SNP_Profile.jpg", width = 1600, height = 1000, units = "px", quality = 100)


  # Plot the dendrogram directly into the PDF device
  plot(hclust_result, main = "Dendrogram of Samples by SNP Profile", 
     xlab = "Height", ylab = "Samples", cex.lab = 1.2, cex.axis = 1, 
     cex.main = 1.2)

  # Close the PDF device
  dev.off()
}

