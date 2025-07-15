library(plyranges)
library(readr)
library(stringr)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
library(EnsDb.Hsapiens.v75)
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
annotation_dir_path <- paste0(working_path, "roadmap_epigenomics_annotation/hg19/")
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
windows_tss<- dropSeqlevels(windows_tss, "chrM", pruning.mode = "coarse")
# Catalog part 2 & 3
windows_tss <- annotate_catalog(windows_tss, feature_tag = "TSS", db = list('UCSC'=TxDb.Hsapiens.UCSC.hg19.knownGene, 'Ensembl'=EnsDb.Hsapiens.v75), db_gap = 2500, db_spe = TRUE, db_tag = "EXTRA_TSS", db_spe_center = 3000)
windows_enh<- load_catalog(annotation_dir = annotation_dir_path,glossary = glossary_path, glossary_group = "ANATOMY", states = c("6_EnhG","7_Enh"))
windows_enh$type<-"Enhancer"
windows_enh<- dropSeqlevels(windows_enh, "chrM", pruning.mode = "coarse")
windows<-c(windows_tss,windows_enh)


# Catalog part 4
windows <- add_flank(windows, feature_tag = "FLANKING", flanking_size = 1000)
windows <- add_background(windows, feature_tag = "BACKGROUND", bin_size = 5000)

df_windows <- data.frame(
  seqnames=seqnames(windows),
  start=start(windows),
  end=end(windows),
  name=windows$feature,
  tissue=windows$tissue,
  type=windows$type)

df_windows <- df_windows %>%
  mutate(tissue = ifelse(type == "EXTRA_TSS", "other", tissue))



df_windows <- df_windows %>%
  mutate(name = ifelse(is.na(name) & type %in% c("TSS", "Enhancer", "EXTRA_TSS"), ".", name))

TSS_Enhancers_windows <- makeGRangesFromDataFrame(df_windows,keep.extra.columns = T)


# Remove the specific chromosome
TSS_Enhancers_windows <- dropSeqlevels(TSS_Enhancers_windows, "chrM", pruning.mode = "coarse")

# Get the current type values
current_types <- TSS_Enhancers_windows$type

# Convert uppercase types to lowercase
TSS_Enhancers_windows$type <- ifelse(current_types %in% c("BACKGROUND", "FLANKING"),
                           tolower(current_types),
                           current_types)

library(BSgenome.Hsapiens.UCSC.hg19)
hg19.si <- seqinfo(BSgenome.Hsapiens.UCSC.hg19)              # all offline

# gr is your existing GRanges
keep <- intersect(seqlevels(TSS_Enhancers_windows), seqlevels(hg19.si))
hg19.si <- hg19.si[keep]          # drop sequences you do not use

seqinfo(TSS_Enhancers_windows) <- hg19.si            # fills seqlengths and genome field
genome(TSS_Enhancers_windows)                        # check it now reads "hg19"

TSS_Enhancers_windows<-trim(TSS_Enhancers_windows)
# Save files
saveRDS(TSS_Enhancers_windows,"Analysis/cfChIP-seq/SetupFiles/H3K4me3_rodmap_hg19/Windows.rds")
write.table(df_windows, file="Analysis/cfChIP-seq/SetupFiles/H3K4me3_rodmap_hg19/Windows.csv", quote=F, sep=",", row.names=F, col.names=F)

write_bed(TSS_Enhancers_windows,"Analysis/cfChIP-seq/SetupFiles/H3K4me3_rodmap_hg19/Windows.bed")
