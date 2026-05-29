#!/usr/bin/env Rscript

# Load the rtracklayer package
library(rtracklayer)
library(parallel)
library(easylift)
library(plyranges)
# Get the command line arguments
args <- commandArgs(trailingOnly = TRUE)

# Number of cores
num_cores <- as.integer(args[1])

# Directory containing the tagAlign files
tagalign_dir <-  args[2]

# Directory to store the lifted over tagAlign files
output_dir <- args[3]

# Specify the path to the chain file
chain_file <- "/date/gcb/gcb_MZ/hg19ToHg38.over.chain"



# Create the output directory if it doesn't exist
if (!dir.exists(output_dir)) {
  dir.create(output_dir)
}

# Get the list of tagAlign files
tagalign_files <- list.files(tagalign_dir, pattern = "*.tagAlign.gz", full.names = TRUE)

# Select the first 'num_files' files
tagalign_files <- tagalign_files[1:length(tagalign_files)]

# Define the function to lift over a single file
liftOverFile <- function(tagalign_file) {
  # Print the name of the file being processed
  print(paste("Processing", tagalign_file))
  
  # Import the tagAlign file
  tagalign <- import.bed(tagalign_file)
  genome( tagalign) <- "hg19"
  # Perform the liftover
  lifted_tagalign <-easylift(tagalign, to="hg38", chain_file)
  
  # Get the base name of the file
  base_name <- basename(tagalign_file)
  
  # Get the base name of the file
  base_name <- basename(tagalign_file)
  
  # Remove the .gz extension from the base name
  base_name <- sub("\\.gz$", "", base_name)
  
  # Define the output file name
  output_file <- file.path(output_dir, paste0(base_name))
  
  # Export the lifted over tagAlign file
  write_bed(lifted_tagalign, output_file)
  
  # Print the name of the file when done
  print(paste(tagalign_file, "done"))
}

# Use mclapply to apply the function in parallel
results <- mclapply(tagalign_files, liftOverFile, mc.cores = num_cores)
