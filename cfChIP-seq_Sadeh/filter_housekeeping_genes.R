library(GenomicRanges)
library(dplyr)
library(readr)
library(biomaRt)
library(plyranges)



# Load GeneDescription and rename the first column to 'name'
GeneDescription <- read_csv("/date/gcb/gcb_MZ/Analysis/cfChIP-seq/SetupFiles/H3K4me3/GeneDescription.csv") %>% 
  dplyr::rename(name = 1)

# Merge with gene description data
House_keeping <- GeneDescription %>% 
  dplyr::filter(GTEX.HouseKeeping==TRUE)

gene_list<-House_keeping$name
mart <- useEnsembl(biomart = "genes", 
                   dataset = "hsapiens_gene_ensembl") # or remove version for latest

gene_coords <- getBM(
  attributes = c("hgnc_symbol", "chromosome_name", "start_position", "end_position", "strand"),
  filters = "hgnc_symbol",
  values = gene_list,
  mart = mart
)

# Filter rows where chromosome_name is only digits
gene_coords_clean <- gene_coords[grepl("^\\d+$", gene_coords$chromosome_name), ]

# Add prefix 'vhr' to chromosome_name
gene_coords_clean$chromosome_name <- paste0("chr", gene_coords_clean$chromosome_name)

# Create a new column for strand symbol
gene_coords_clean$strand <- ifelse(gene_coords_clean$strand == 1, "+",
                                 ifelse(gene_coords_clean$strand == -1, "-", NA))

gene_coords_clean<-gene_coords_clean %>% dplyr::rename(start=3,end=4) %>%
makeGRangesFromDataFrame(keep.extra.columns = T)


write_bed(gene_coords_clean,"/proj/user/mattia/Analysis/cfChIP-seq/SetupFiles/H3K27ac_hg38/House_keeping.bed")
# View result
head(df_clean)


