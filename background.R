#load regions in which the background has been computed

counts_list<-readRDS("cfChromatin_in_MS/H3K27ac_ChIP_counts.rds")



estimateBackground <- function(gr) {
  # Extract background regions
  background_regions <- gr[gr$type == "background"]
  
  # Calculate the length of each region
  region_lengths <- width(background_regions)
  
  # Filter regions that are longer than 4kb
  long_background_regions <- background_regions[region_lengths >= 4000]
  
  # Extract counts from the filtered regions
  X <- mcols(long_background_regions)$counts
  
  # Find the 95th quantile of X
  T <- quantile(X, 0.95)
  
  # Restrict ourselves to values below T
  X <- X[X <= T]
  
  # Maximum likelihood of truncated poisson
  lambda_hat <- which.max(dpois(X, lambda = 1:length(X)))
  
  # Convert to reads/Kb (the median length of background windows is 5KB)
  lambda_hat <- lambda_hat / 5
  
  # Add the background estimate as a new column to the GRanges object
  gr$background <- lambda_hat
  
  return(gr)
}

# Apply the estimateBackground function to each GRanges object in the list
# Set mc.cores to the number of cores you want to use
granges_list_with_background <- mclapply(counts_list, estimateBackground, mc.cores = detectCores())


saveRDS(granges_list_with_background,"cfChromatin_in_MS/H3K27ac_ChIP_counts.rds")
