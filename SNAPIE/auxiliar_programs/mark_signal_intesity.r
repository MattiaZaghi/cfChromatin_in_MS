#!/usr/bin/env Rscript

#===============================================================================
# cfChIP Signal Intensity Analysis Pipeline
#
# Authors: Ze Zhang & Paulo Cordeiro
# Date: August 20, 2025
# 
# Description: Complete pipeline for cfChIP (cell-free Chromatin Immunoprecipitation)
#              signal intensity analysis. Processes BED files containing fragment
#              coordinates, applies blacklist filtering, calculates signal at
#              sites of interest with normalization using housekeeping regions,
#              and generates normalized signal intensity matrix.
#
# Usage: Rscript mark_signal_intesity.r <input.bed> <blacklist.bed> <sites_of_interest.bed> <housekeeping.bed> <mark>
#
# Output: <SampleName>_<MARK>_Signal_Intensity_Matrix.csv
#===============================================================================

# Load required libraries
library(stringr)
library(GenomicRanges)
library(rtracklayer)
library(dplyr)
library(tidyr)
library(ggplot2)
library(data.table)

# Function to check command line arguments
check_args <- function() {
  args <- commandArgs(trailingOnly = TRUE)
  if (length(args) != 5) {
    cat("Usage: Rscript mark_signal_intesity.r <input.bed> <blacklist.bed> <sites_of_interest.bed> <housekeeping.bed> <mark>\n")
    cat("Example: Rscript mark_signal_intesity.r sample_fragments.bed hg19-blacklist.v2.bed PRAD_vs_WBC_up.bed housekeeping_k4me3.bed H3K4me3\n")
    cat("Output: <SampleName>_<MARK>_Signal_Intensity_Matrix.csv\n")
    cat("\nAuthors: Ze Zhang & Paulo Cordeiro (2025)\n")
    quit(status = 1)
  }
  return(args)
}

# Function to create tiles across sites of interest
tile_sites <- function(sites_file_list, exclude_regions = NULL, bin_width = 40, flank_bp = 3000) {
  return_sites_list <- list()
  
  for (sites_file in sites_file_list) {
    if (!file.exists(sites_file)) {
      warning(paste("File not found:", sites_file))
      next
    }
    
    sites <- import(sites_file)
    
    # Exclude blacklisted regions if provided
    if (!is.null(exclude_regions)) {
      sites <- sites[!(sites %over% exclude_regions)]
    }
    
    # Standardize size and expand to include flanks
    sites <- resize(sites, width = 1, fix = 'center')
    sites <- sites + flank_bp
    sites <- resize(sites, flank_bp * 2, fix = 'start')
    
    # Remove overlapping sites
    sites <- GenomicRanges::reduce(sites)
    sites <- subset(sites, width == flank_bp * 2)
    
    # Create tiles
    tiles <- tile(sites, width = bin_width) %>% unlist()
    
    # Map each tile back to its original site
    ol <- findOverlaps(tiles, sites)
    agg <- tiles[from(ol)]
    tmp <- sites[to(ol)]
    agg$site_start <- start(tmp)
    agg$bin <- start(agg) - agg$site_start - flank_bp
    
    return_sites_list <- c(return_sites_list, list(agg))
  }
  
  return(return_sites_list)
}

# Function to calculate signal at sites
signal_at_sites <- function(frags, sites_list, sites_names, remove_top = 0.05, remove_peaks_wider_than = 4000) {
  return_df <- NULL
  count <- 1
  
  for (sites in sites_list) {
    if (!is.na(remove_peaks_wider_than)) {
      sites <- sites[width(sites) <= remove_peaks_wider_than]
    }
    
    sites$counts <- suppressWarnings(countOverlaps(sites, frags))
    
    # Remove top quantile of sites by fragment count
    tmp <- sites %>% 
      as.data.frame() %>% 
      group_by(site_start) %>% 
      summarize(tot = sum(counts))
    
    quantile_cutoff <- quantile(tmp$tot, 1 - remove_top)
    sites_to_remove <- tmp$site_start[tmp$tot > quantile_cutoff]
    sites <- subset(sites, !(site_start %in% sites_to_remove))
    
    counts <- sites %>%
      as.data.frame() %>%
      group_by(bin) %>%
      summarize(
        sites = sites_names[count],
        frag_counts = sum(counts)
      )
    
    counts$total_counts <- length(frags)
    return_df <- rbind(return_df, counts)
    count <- count + 1
  }
  
  return(return_df)
}

# Main pipeline function
cfchip_pipeline <- function(input_bed, blacklist_bed, sites_interest_bed, housekeeping_bed, mark) {
  
  # Generate output file name with sample name and mark
  sample_name <- tools::file_path_sans_ext(basename(input_bed))
  output_csv <- paste0(mark, "_Signal_Intensity_Matrix","_",sample_name,".csv")
  
  # Check if input files exist
  files_to_check <- c(input_bed, blacklist_bed, sites_interest_bed, housekeeping_bed)
  for (file in files_to_check) {
    if (!file.exists(file)) {
      stop(paste("File not found:", file))
    }
  }
  
  cat("Processing BED file:", input_bed, "\n")
  
  # Load input BED file
  cat("Loading fragments from BED file...\n")
  frags <- import(input_bed)
  
  # Convert to fragment centers if necessary
  if (all(width(frags) > 1)) {
    frags$frag_width <- width(frags)
    frags <- resize(frags, width = 1, fix = 'center')
  }
  
  cat("Total fragments loaded:", length(frags), "\n")
  
  # Load blacklist regions
  cat("Loading blacklisted regions...\n")
  exclude_regions <- import(blacklist_bed, format = 'bed')
  
  # Prepare sites of interest list
  sites_file_list <- c(sites_interest_bed, housekeeping_bed)
  sites_names <- c("SITES_OF_INTEREST", "MARK_HOUSEKEEPING_POSITION")
  
  cat("Processing sites of interest...\n")
  sites_list <- tile_sites(sites_file_list, exclude_regions = exclude_regions)
  
  # Calculate signal at sites
  cat("Calculating signal at sites of interest...\n")
  toplot <- signal_at_sites(
    frags = frags,
    sites_list = sites_list,
    sites_names = sites_names,
    remove_peaks_wider_than = 4000
  )
  
  # Add sample information
  sample_name <- tools::file_path_sans_ext(basename(input_bed))
  toplot$study_name <- sample_name
  
  cat("Starting normalization and signal calculation...\n")
  
  # Normalization parameters
  shoulders <- 2800
  normalize_to_these_sites <- "MARK_HOUSEKEEPING_POSITION"
  
  # Shoulder normalization
  cat("Applying shoulder normalization...\n")
  shoulder_norm <- toplot %>%
    subset(bin < -shoulders | bin > shoulders) %>%
    group_by(study_name, sites) %>% 
    summarize(shoulder_norm_factor = median(frag_counts), .groups = 'drop')
  
  toplot <- merge(toplot, shoulder_norm)
  toplot$frag_counts <- toplot$frag_counts - toplot$shoulder_norm_factor
  toplot$frag_counts[toplot$frag_counts < 0] <- 0
  
  # Initial normalization
  toplot$counts_norm <- toplot$frag_counts
  
  # In-sample normalization
  cat("Applying in-sample normalization...\n")
  insample_norm <- toplot %>% 
    subset(sites == normalize_to_these_sites) %>%
    group_by(study_name) %>%
    summarize(insample_norm_factor = mean(counts_norm), .groups = 'drop')
  
  toplot <- merge(toplot, insample_norm)
  toplot$counts_norm <- toplot$counts_norm / toplot$insample_norm_factor
  
  # Calculate AUC (Area Under the Curve)
  cat("Calculating AUC...\n")
  auc_score <- toplot %>%
    group_by(study_name, sites) %>%
    summarize(auc = sum(counts_norm), .groups = 'drop')
  
  toplot <- merge(toplot, auc_score)
  
  # Filter for bin = 0 and remove duplicates
  toplot <- toplot %>%
    subset(bin == 0) %>%
    distinct()
  
  # Calculate log AUC
  toplot$logauc <- log(toplot$auc + 1, 2)
  
  # Filter only sites of interest
  toplot_final <- toplot %>% 
    filter(sites == "SITES_OF_INTEREST") %>%
    select(study_name, logauc)
  
  # Rename columns for output
  colnames(toplot_final) <- c("Sample_Name", "Normalized_Signal")
  
  # Save results
  cat("Saving results to:", output_csv, "\n")
  write.csv(toplot_final, output_csv, row.names = FALSE)
  
  cat("Pipeline completed successfully!\n")
  cat("Output file:", output_csv, "\n")
  cat("Sample processed:", sample_name, "\n")
  cat("Normalized signal:", round(toplot_final$Normalized_Signal, 4), "\n")
  
  # Statistical summary
  cat("\nSummary of intermediate calculations:\n")
  summary_stats <- toplot %>%
    group_by(sites) %>%
    summarize(
      total_bins = n(),
      mean_auc = mean(auc),
      mean_logauc = mean(logauc),
      .groups = 'drop'
    )
  print(summary_stats)
  
  return(output_csv)
}

# Main function for command line usage
main <- function() {
  # Get command line arguments
  args <- check_args()
  input_bed <- args[1]
  blacklist_bed <- args[2]
  sites_interest_bed <- args[3]
  housekeeping_bed <- args[4]
  mark <- args[5]
  
  # Run pipeline
  tryCatch({
    cfchip_pipeline(input_bed, blacklist_bed, sites_interest_bed, housekeeping_bed, mark)
  }, error = function(e) {
    cat("Error:", e$message, "\n")
    quit(status = 1)
  })
}

# Run main function if script is called directly
if (!interactive()) {
  main()
}