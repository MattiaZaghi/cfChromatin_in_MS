library(rtracklayer)
library(easylift)
library(plyranges)

# Define directories
input_dir <- "/proj/user/mattia/Analysis/BED/H3K27ac_baca/"
output_dir <- "/proj/user/mattia/Analysis/BED/H3K27ac_hg38/"
chain_file <- "/date/gcb/gcb_MZ/hg19ToHg38.over.chain"


# List all BED files in the input directory that start with "H00"
bed_files <- list.files(input_dir, pattern = "*\\.bed$", full.names = TRUE)

# Select the first 15 files
bed_files <- head(bed_files, 15)


# Function to perform lift-over
lift_over_bed <- function(bed_file, chain_file, output_dir) {
  # Read the BED file
  bed_data <- read_narrowpeaks(bed_file)
  genome(bed_data) <- "hg19"
  
  # Perform lift-over
  lifted_data <- easylift(bed_data, to = "hg38", chain = chain_file)
  
  # Create output file path
  output_file <- file.path(output_dir, basename(bed_file))
  
  # Write the lifted data to the new BED file
  write_bed(lifted_data, output_file)
}

lift_over_narrowPeaks<- function(bed_file, chain_file, output_dir) {
  # Read the narrowPeak file (all 10 columns)
  bed_data <- read.table(bed_file, header = FALSE, sep = "\t",
                         col.names = c("chr", "start", "end", "name", "score",
                                       "strand", "signalValue", "pValue", "qValue", "peak"))
  
  # Convert "." strand to "*" for GRanges compatibility
  bed_data$strand[bed_data$strand == "."] <- "*"
  
  # Convert to GRanges, keeping metadata
  gr <- GRanges(seqnames = bed_data$chr,
                ranges = IRanges(start = bed_data$start + 1, end = bed_data$end),  # BED is 0-based
                strand = bed_data$strand,
                name = bed_data$name,
                score = bed_data$score,
                signalValue = bed_data$signalValue,
                pValue = bed_data$pValue,
                qValue = bed_data$qValue,
                peak = bed_data$peak)
  genome(gr) <- "hg19"
  
  # Perform lift-over
  lifted_data <- easylift(gr, to = "hg38", chain = chain_file)
  
  # Convert back to data frame with all columns
  # Convert "*" back to "." for standard narrowPeak format
  lifted_strand <- as.character(strand(lifted_data))
  lifted_strand[lifted_strand == "*"] <- "."
  
  lifted_df <- data.frame(
    chr = as.character(seqnames(lifted_data)),
    start = start(lifted_data) - 1,  # Convert back to 0-based
    end = end(lifted_data),
    name = mcols(lifted_data)$name,
    score = mcols(lifted_data)$score,
    strand = lifted_strand,
    signalValue = mcols(lifted_data)$signalValue,
    pValue = mcols(lifted_data)$pValue,
    qValue = mcols(lifted_data)$qValue,
    peak = mcols(lifted_data)$peak
  )
  
  # Create output file path
  output_file <- file.path(output_dir, basename(bed_file))
  
  # Write all 10 columns
  write.table(lifted_df, output_file, 
              sep = "\t", quote = FALSE, 
              row.names = FALSE, col.names = FALSE)
  
  cat("Lifted", nrow(lifted_df), "peaks to", output_file, "\n")
  
  return(lifted_df)
}



# Loop through each BED file and perform lift-over
for (bed_file in bed_files) {
  lift_over_bed(bed_file, chain_file, output_dir)
}




