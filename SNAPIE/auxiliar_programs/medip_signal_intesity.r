#!/usr/bin/env Rscript

#===============================================================================
# MeDIP Signal Intensity Analysis Pipeline
#
# Authors: Ze Zhang & Paulo Cordeiro
# Date: August 20, 2025
# 
# Description: Complete pipeline for MeDIP (Methylated DNA Immunoprecipitation)
#              signal intensity analysis. Converts BED files to RDS format,
#              processes with housekeeping genes for normalization, and 
#              calculates normalized signal intensity against sites of interest.
#
# Usage: Rscript medip_signal_intesity.r <sample.bed> <housekeeping_genes.bed> <sites_of_interest.bed>
#
# Output: MEDIP_Signal_Intensity_Matrix_<SampleName>.csv
#===============================================================================

# Load required libraries
library(stringr)
library(rtracklayer)
library(GenomicRanges)
library(dplyr)
library(tidyr)
library(ggplot2)
library(data.table)
library(readxl)
library(GGally)
library(readr)
library(gridExtra)
library(biomaRt)

# Complete MeDIP pipeline
medip_pipeline <- function(bed_file, housekeeping_bed, sites_bed) {
  
  # Step 1: Convert BED to RDS
  message("Step 1: Converting BED to RDS...")
  
  # Check if BED file exists
  if (!file.exists(bed_file)) {
    stop("BED file not found: ", bed_file)
  }
  
  # Generate RDS filename and sample name
  rds_file <- str_replace(bed_file, "\\.bed(\\.gz)?$", ".rds")
  sample_name <- tools::file_path_sans_ext(basename(bed_file))
  
  message('Making granges object for ', bed_file)
  
  # Import BED file
  tmp <- import(bed_file)
  tmp$frag_width <- width(tmp)
  tmp <- resize(tmp, width = 1, fix = 'center')
  
  # Save as RDS
  saveRDS(tmp, file = rds_file)
  message('RDS file saved as: ', rds_file)
  
  # Step 2: Process with housekeeping genes
  message("Step 2: Processing with housekeeping genes...")
  
  # Check if housekeeping BED file exists
  if (!file.exists(housekeeping_bed)) {
    stop("Housekeeping BED file not found: ", housekeeping_bed)
  }
  
  # Read housekeeping BED file
  hskp_5mc <- read.delim(housekeeping_bed, header = FALSE)
  colnames(hskp_5mc) <- c("chr", "start", "end")
  
  # Create GRanges object for housekeeping genes
  hskp_5mc <- GRanges(seqnames = hskp_5mc$chr,
                      ranges = IRanges(start = hskp_5mc$start, end = hskp_5mc$end))
  
  # Read the RDS file
  sample_data <- readRDS(rds_file)
  
  # Initialize data frame with location information
  All_count <- data.frame(location = paste0(as.character(seqnames(hskp_5mc)), ":", 
                                            start(hskp_5mc), "-", end(hskp_5mc)))
  
  # Count overlaps with housekeeping genes
  overlap_counts <- countOverlaps(hskp_5mc, sample_data)
  All_count[[sample_name]] <- overlap_counts
  
  # Set row names and keep data frame structure
  rownames(All_count) <- make.unique(All_count$location)
  All_count <- All_count[, -1, drop = FALSE]
  
  # Step 3: Process with sites of interest
  message("Step 3: Processing with sites of interest...")
  
  # Check if sites BED file exists
  if (!file.exists(sites_bed)) {
    stop("Sites BED file not found: ", sites_bed)
  }
  
  # Read sites of interest BED file
  Sites <- read.delim(sites_bed, header = FALSE)
  Sites <- GRanges(seqnames = Sites[,1],
                   ranges = IRanges(start = Sites[,2], end = Sites[,3]))
  
  # Count overlaps with Sites
  sites_overlap_counts <- countOverlaps(Sites, sample_data)
  Sites_Count <- sum(sites_overlap_counts)
  
  # Get housekeeping count for this sample
  HSK_total <- sum(All_count[, 1])
  
  # Calculate normalized signal
  Normalized_Signal <- Sites_Count / HSK_total
  
  # Create result data frame
  meta_sub <- data.frame(
    Sample_Name = sample_name,
    "Normalized_Signal" = Normalized_Signal,
    check.names = FALSE
  )
  
  # Generate output filename
  output_file <- paste0("medip_Signal_Intensity_Matrix_", sample_name, ".csv")
  
  message("Step 4: Saving final results...")
  write.csv(meta_sub, output_file, row.names = FALSE)
  
  message("Pipeline completed successfully!")
  message("Sample: ", sample_name)
  message("Sites Count: ", Sites_Count)
  message("Housekeeping Count: ", HSK_total)
  message("Normalized Signal: ", round(Normalized_Signal, 6))
  message("Output file: ", output_file)
  
  # Clean up temporary RDS file (optional)
  # file.remove(rds_file)
  
  return(output_file)
}

# Main function for command line usage
main <- function() {
  # Get command line arguments
  args <- commandArgs(trailingOnly = TRUE)
  
  if (length(args) != 3) {
    cat("Usage: Rscript medip_signal_intesity.r <sample.bed> <housekeeping_genes.bed> <sites_of_interest.bed>\n")
    cat("Example: Rscript medip_signal_intesity.r MdHP10028203.bed housekeeping_5mC.bed HeartvsWBC_Hyper.bed\n")
    cat("Output: MEDIP_Signal_Intensity_Matrix_<SampleName>.csv\n")
    cat("\nAuthors: Ze Zhang & Paulo Cordeiro (2025)\n")
    quit(status = 1)
  }
  
  bed_file <- args[1]
  housekeeping_bed <- args[2]
  sites_bed <- args[3]
  
  # Run pipeline
  tryCatch({
    medip_pipeline(bed_file, housekeeping_bed, sites_bed)
  }, error = function(e) {
    cat("Error:", e$message, "\n")
    quit(status = 1)
  })
}

# Run main function if script is called directly
if (!interactive()) {
  main()
}