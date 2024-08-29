# Load required libraries
library(rtracklayer)
library(GenomicRanges)
library(plyranges)
# Get list of BED files in the directory
bed_files <- list.files(path = "/date/gcb/gcb_MZ/Analysis/BED/H3K27ac_single_cell/", pattern = "\\.bed$", full.names = TRUE)

# Initialize an empty list to store the GRanges objects
gr_list <- list()

# Loop over BED files
for (i in seq_along(bed_files)) {
  # Import BED file as GRanges
  gr <- import.bed(bed_files[i])
  
  
  # Filter for chromosomes 1:22, X, Y
  gr <- gr[seqnames(gr) %in% c(paste0("chr", 1:22), "chrX", "chrY")]
  
  # Change the levels of the seqnames factor
  seqlevels(gr, pruning.mode = "coarse") <- c(paste0("chr", 1:22), "chrX", "chrY")
  
  # Change seqinfo to hg38
  seqinfo(gr) <- Seqinfo(genome = "hg38")
  
  # Add to list
  gr_list[[i]] <- gr
  
  # Export GRanges object back to BED file
  export.bed(gr, bed_files[i])
}

# Filter for chromosomes 1:22, X, Y
TSS.windows <- TSS.windows[seqnames(TSS.windows) %in% c(paste0("chr", 1:22), "chrX", "chrY")]
# Change the levels of the seqnames factor
saveRDS(TSS.windows,"/date/gcb/gcb_MZ/Analysis/cfChIP-seq/SetupFiles/H3K27ac_single_cell/Windows.rds")
write_bed(TSS.windows,"/date/gcb/gcb_MZ/Analysis/cfChIP-seq/SetupFiles/H3K27ac_single_cell/Windows.bed")

df_windows <- data.frame(
  seqnames=seqnames(TSS.windows),
  starts=start(TSS.windows),
  ends=end(TSS.windows),
  type=TSS.windows$type,
  name=TSS.windows$name,
  tissue=TSS.windows$tissue
)

# Create a new data frame with a comma as the first column and rownumber as the second column

write.csv(df_windows,"/date/gcb/gcb_MZ/Analysis/cfChIP-seq/SetupFiles/H3K27ac_single_cell/Windows.csv",row.names = TRUE)

