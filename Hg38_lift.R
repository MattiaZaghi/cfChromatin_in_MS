# Load the rtracklayer package
library(rtracklayer)
library(parallel)
# Specify the path to the chain file
chain_file <- "/date/gcb/gcb_MZ/hg19ToHg38.over.chain"

# Import the chain file
chain <- import.chain(chain_file)

# Specify the directory containing the tagAlign files
tagalign_dir <- "/date/gcb/gcb_MZ/roadmap_epigenomics/Tag_align/"

# Specify the directory to store the lifted over tagAlign files
output_dir <- "/date/gcb/gcb_MZ/roadmap_epigenomics/Tag_align/Hg38_lift/"

# Create the output directory if it doesn't exist
if (!dir.exists(output_dir)) {
  dir.create(output_dir)
}

# Get the list of tagAlign files
tagalign_files <- list.files(tagalign_dir, pattern = "*.tagAlign.gz", full.names = TRUE)

# Define the number of files to process
num_files <- 2

# Select the first 'num_files' files
tagalign_files <- tagalign_files[1:num_files]

# Define the number of cores to use for parallel processing
num_cores <- detectCores() 

# Define the function to lift over a single file
liftOverFile <- function(tagalign_file) {
  # Import the tagAlign file
  tagalign <- import.bed(tagalign_file)
  
  # Perform the liftover
  lifted_tagalign <- liftOver(tagalign, chain)
  
  # Get the base name of the file
  base_name <- basename(tagalign_file)
  
  # Define the output file name
  output_file <- file.path(output_dir, paste0(base_name, ".hg38.tagAlign.gz"))
  
  # Export the lifted over tagAlign file
  export.bed(lifted_tagalign, output_file)
  
  # Generate the .tbi file
  system(paste("tabix -p vcf", output_file))
}

# Use mclapply to apply the function in parallel
results <- mclapply(tagalign_files, liftOverFile, mc.cores = num_cores)