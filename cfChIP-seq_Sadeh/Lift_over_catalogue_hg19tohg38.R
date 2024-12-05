library(rtracklayer)
library(easylift)
library(plyranges)

# Define directories
input_dir <- "/date/gcb/gcb_MZ/Analysis/BED/H3K27ac/"
output_dir <- "/date/gcb/gcb_MZ/Analysis/BED/H3K27ac_hg38/"
chain_file <- "/date/gcb/gcb_MZ/hg19ToHg38.over.chain"

# List all BED files in the input directory that start with "GSM"
bed_files <- list.files(input_dir, pattern = "^GSM.*\\.bed$", full.names = TRUE)

# Function to perform lift-over
lift_over_bed <- function(bed_file, chain_file, output_dir) {
  # Read the BED file
  bed_data <- read_bed(bed_file)
  genome(bed_data) <- "hg19"
  
  # Perform lift-over
  lifted_data <- easylift(bed_data, to = "hg38", chain = chain_file)
  
  # Create output file path
  output_file <- file.path(output_dir, basename(bed_file))
  
  # Write the lifted data to the new BED file
  write_bed(lifted_data, output_file)
}

# Loop through each BED file and perform lift-over
for (bed_file in bed_files) {
  lift_over_bed(bed_file, chain_file, output_dir)
}