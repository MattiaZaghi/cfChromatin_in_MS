library(plyranges)
library(readr)
library(stringr)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(EnsDb.Hsapiens.v86)
library(org.Hs.eg.db)
library(parallel)
library(rtracklayer)
library(preprocessCore)
library(GenomeInfoDb)


sessionInfo()

#############################################################################################################################
###                                                          IMPORT                                                       ###
#############################################################################################################################

#working_path <- "/media/sf_ubuntu_music/cfChromatin/"
working_path <- "/date/gcb/gcb_MZ/"
annotation_dir_path <- paste0(working_path, "Chrom_HMM_hg38_annotation/")
glossary_path <- paste0(annotation_dir_path, "Glossary.txt")

core <- 20

setwd(dir = working_path)
source("/home/mattia/cfChromatin_in_MS/cfChromatin_func_new.R")

#############################################################################################################################
###                                                    CATALOGUES CREATION                                                ###
#############################################################################################################################

###### TSS
# Catalog part 1
windows_tss <- load_catalog(annotation_dir = annotation_dir_path, glossary = glossary_path, glossary_group = "ANATOMY", states = c("1_TssA", "2_TssAFlnk"))
# Catalog part 2 & 3
windows_tss <- annotate_catalog(windows_tss, feature_tag = "TSS", db = list('UCSC'=TxDb.Hsapiens.UCSC.hg38.knownGene, 'Ensembl'=EnsDb.Hsapiens.v86), db_gap = 2500, db_spe = TRUE, db_tag = "EXTRA_TSS", db_spe_center = 3000)
windows_enh<- load_catalog(annotation_dir = annotation_dir_path,glossary = glossary_path, glossary_group = "ANATOMY", states = c("6_EnhG","7_Enh"))
windows_enh$type<-"Enhancer"
windows<-c(windows_tss,windows_enh)
# Catalog part 4
windows <- add_flank(windows, feature_tag = "FLANKING", flanking_size = 1000)
windows <- add_background(windows, feature_tag = "BACKGROUND", bin_size = 5000)

# Statistics
show_statistics(windows, db_tag = "EXTRA_TSS")

df_windows<- data.frame(
  seqnames=seqnames(windows),
  start=start(windows)-1,
  end=end(windows),
  name=windows$feature,
  tissue=windows$tissue,
  type=windows$type)


df_windows <- df_windows %>%
  mutate(tissue = ifelse(type == "EXTRA_TSS", "other", tissue))



df_windows <- df_windows %>%
  mutate(name = ifelse(is.na(name) & type %in% c("TSS", "Enhancer", "EXTRA_TSS"), ".", name))


TSS_Enhancers_windows <- makeGRangesFromDataFrame(df_windows,keep.extra.columns = T)


# Create a new Seqinfo object with hg38 genome information
hg38_seqinfo <- Seqinfo(genome = "hg38")

# Assuming your_granges_object is your GRanges object
seqinfo(TSS_Enhancers_windows) <- hg38_seqinfo

# Check the updated Seqinfo object
seqinfo(TSS_Enhancers_windows)

# Save files
saveRDS(TSS_Enhancers_windows,"/date/gcb/gcb_MZ/Analysis/cfChIP-seq/SetupFiles/H3K27ac_hg38/Windows.rds")
write.table(df_windows, file="/date/gcb/gcb_MZ/Analysis/cfChIP-seq/SetupFiles/H3K27ac_hg38/Windows.csv", quote=F, sep=",", row.names=F, col.names=F)

write_bed(TSS_Enhancers_windows,"/date/gcb/gcb_MZ/Analysis/cfChIP-seq/SetupFiles/H3K27ac_hg38/Windows.bed")
