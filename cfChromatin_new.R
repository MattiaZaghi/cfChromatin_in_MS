library(plyranges)
library(readr)
library(stringr)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(EnsDb.Hsapiens.v86)
library(org.Hs.eg.db)
library(parallel)
library(rtracklayer)
library(preprocessCore)

sessionInfo()

#############################################################################################################################
###                                                          IMPORT                                                       ###
#############################################################################################################################

#working_path <- "/media/sf_ubuntu_music/cfChromatin/"
working_path <- "/date/gcb/gcb_MZ/"
annotation_dir_path <- paste0(working_path, "Chrom_HMM_hg38_annotation/")
glossary_path <- paste0(annotation_dir_path, "Glossary.txt")


core <- 30

setwd(dir = working_path)
source("/home/mattia/cfChromatin_in_MS/cfChromatin_func_new.R")

#############################################################################################################################
###                                                    CATALOGUES CREATION                                                ###
#############################################################################################################################

###### TSS
# Catalog part 1
windows_TSS <- load_catalog(annotation_dir = annotation_dir_path, glossary = glossary_path, glossary_group = "ANATOMY", states = c("1_TssA", "2_TssAFlnk"))
# Catalog part 2 & 3
windows_TSS <- annotate_catalog(windows, feature_tag = "TSS", db = list('UCSC'=TxDb.Hsapiens.UCSC.hg38.knownGene, 'Ensembl'=EnsDb.Hsapiens.v86), db_gap = 2500, db_spe = TRUE, db_tag = "EXTRA_TSS", db_spe_center = 3000)
# Catalog enhacers
windows_enh<-load_catalog(annotation_dir = annotation_dir_path, glossary = glossary_path, glossary_group = "ANATOMY",  states = c("6_EnhG","7_Enh"))
# Unify enhancers & TSS
windows<-c(windows_TSS,windows_enh)
# Catalog part 4
windows <- add_flank(windows, feature_tag = "FLANKING", flanking_size = 1000)
windows <- add_background(windows, feature_tag = "BACKGROUND", bin_size = 5000)

# Statistics
show_statistics(windows, db_tag = "EXTRA_TSS")

# Save files
saveRDS(windows,paste0(working_path,"Analysis/cfChIP-seq/SetupFiles/H3K27ac_hg38/Windows.rds"))

df_windows <- data.frame(
  seqnames=seqnames(windows),
  starts=start(windows)-1,
  ends=end(windows),
  tissue=windows$tissue,
  feature=windows$feature,
  type=windows$type)

write.table(df_windows, file=paste0(working_path,"Analysis/cfChIP-seq/SetupFiles/H3K27ac_hg38/Windows.bed"), quote=F, sep=",", row.names=F, col.names=F)
