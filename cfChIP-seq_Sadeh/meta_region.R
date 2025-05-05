library(GenomicRanges)
library(dplyr)
library(readr)
library(plyranges)
options(scipen = 999)

# Load the GRanges object
TSS.windows <- readRDS("/proj/user/mattia/Analysis/cfChIP-seq/SetupFiles/H3K27ac_hg38/Windows.rds")

# Assuming TSS.windows is your GRanges object
# Filter rows where type is "Enhancer"
meta_enhancers <- TSS.windows[TSS.windows$type == "Enhancer"]

# Set the name column based on the tissue column
meta_enhancers$name <- sapply(meta_enhancers$tissue, function(tissue) {
  if (is.na(tissue)) {
    return("z.none")
  } else {
    tissues <- unlist(strsplit(tissue, ";"))
    if (length(tissues) > 4) {
      return("a.Many")
    } else {
      return("b.Middle")
    }
  }
})

# Merge adjacent regions
merged_enhancers <- IRanges::reduce(meta_enhancers)

# Sum the number of tissues for merged regions
merged_enhancers$tissue_count <- sapply(seq_along(merged_enhancers), function(i) {
  sum(sapply(meta_enhancers[subjectHits(findOverlaps(merged_enhancers[i], meta_enhancers))]$tissue, function(tissue) {
    if (is.na(tissue)) {
      return(0)
    } else {
      return(length(unlist(strsplit(tissue, ";"))))
    }
  }))
})

# Set the name column based on the tissue count
merged_enhancers$name <- sapply(merged_enhancers$tissue_count, function(count) {
  if (count == 0) {
    return("z.none")
  } else if (count > 4) {
    return("a.Many")
  } else {
    return("b.Middle")
  }
})

# Convert to data frame and remove unnecessary columns
merged_enhancers <- as.data.frame(merged_enhancers)

# Prepare the data frame for BED format
meta_enhancers_df <- data.frame(
  chr = merged_enhancers$seqnames,
  start = merged_enhancers$start,  # BED format uses 0-based start
  end = merged_enhancers$end,
  name = merged_enhancers$name,
  score = 0,
  strand = merged_enhancers$strand
) 


# Save as BED file
write_tsv(meta_enhancers_df, file = "/date/gcb/gcb_MZ/Analysis/cfChIP-seq/SetupFiles/H3K27ac_hg38/Meta-enhancers.bed", col_names = FALSE)

# Load GeneDescription and rename the first column to 'name'
GeneDescription <- read_csv("/date/gcb/gcb_MZ/Analysis/cfChIP-seq/SetupFiles/H3K4me3/GeneDescription.csv") %>% 
  dplyr::rename(name = 1)

#filter a bed file of only merged enhancers 

meta_enhacers_many <- as.data.frame(meta_enhancers) %>% dplyr::filter(name=="a.Many")%>% 
makeGRangesFromDataFrame()%>% reduce()

write_bed(meta_enhacers_many,"/proj/user/mattia/Analysis/cfChIP-seq/SetupFiles/H3K27ac_hg38/enhancers_many.bed")

 Save as BED file
write_tsv(meta_enhancers_df, file = "/date/gcb/gcb_MZ/Analysis/cfChIP-seq/SetupFiles/H3K27ac_hg38/Meta-enhancers.bed", col_names = FALSE)

library(GenomicRanges)
library(dplyr)
library(readr)

# Load the GRanges object
TSS.windows <- readRDS("/date/gcb/gcb_MZ/Analysis/cfChIP-seq/SetupFiles/H3K27ac_hg38/Windows.rds")

# Convert to data frame for filtering
TSS.windows_df <- as.data.frame(TSS.windows)

# Filter rows where type is "TSS"
meta_genes <- TSS.windows_df %>% filter(type == "TSS")

# Load GeneDescription and rename the first column to 'name'
GeneDescription <- read_csv("/date/gcb/gcb_MZ/Analysis/cfChIP-seq/SetupFiles/H3K4me3/GeneDescription.csv") %>% 
  dplyr::rename(name = 1)

# Merge with gene description data
meta_genes <- merge(meta_genes, GeneDescription, by.x = "name", by.y = "name", all.x = TRUE)

# Set the feature column based on expression and CpG status
meta_genes$feature <- apply(meta_genes, 1, function(row) {
  if (is.na(row["GTEX.Expressed"])) {
    return("e.NotExpressed")
  } else if (!is.na(row["isCPG"]) && row["isCPG"] == TRUE) {
    if (!is.na(row["GTEX.HighExpressed"])) {
      return("a.CpG.High")
    } else {
      return("b.CpG.Low")
    }
  } else if (!is.na(row["isCPG"]) && row["isCPG"] == FALSE) {
    if (!is.na(row["GTEX.HighExpressed"])) {
      return("c.NonCpG.High")
    } else {
      return("d.NonCpG.Low")
    }
  } else {
    return("d.NonCpG.Low")  # Default case if isCPG is NA
  }
})

# Convert back to GRanges
meta_genes_gr <- makeGRangesFromDataFrame(meta_genes, keep.extra.columns = TRUE)

# Merge adjacent TSS regions
merged_genes <- IRanges::reduce(meta_genes_gr)

# Preserve CpG information during the merge
merged_genes$CpG_info <- sapply(seq_along(merged_genes), function(i) {
  overlapping_genes <- meta_genes_gr[subjectHits(findOverlaps(merged_genes[i], meta_genes_gr))]
  if (any(overlapping_genes$feature %in% c("a.CpG.High", "b.CpG.Low"))) {
    if (any(overlapping_genes$feature == "a.CpG.High")) {
      return("a.CpG.High")
    } else {
      return("b.CpG.Low")
    }
  } else {
    if (any(overlapping_genes$feature == "c.NonCpG.High")) {
      return("c.NonCpG.High")
    } else {
      return("d.NonCpG.Low")
    }
  }
})

# Convert to data frame and remove unnecessary columns
merged_genes_df <- as.data.frame(merged_genes)

# Prepare the data frame for BED format
meta_genes_df <- data.frame(
  chr = merged_genes_df$seqnames,
  start = merged_genes_df$start - 1,  # BED format uses 0-based start
  end = merged_genes_df$end,
  name = merged_genes_df$CpG_info,
  score = 0,
  strand = merged_genes_df$strand
)

# Save as BED file
write_tsv(meta_genes_df, file = "/date/gcb/gcb_MZ/Analysis/cfChIP-seq/SetupFiles/H3K27ac_hg38/Meta-genes.bed", col_names = FALSE)

# Print the meta_genes DataFrame
print("meta_genes DataFrame:")
print(meta_genes_df)
