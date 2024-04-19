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
bed_dir_path <- "/date/gcb/gcb_MZ/Round1/bed/bedtools/normal/"

core <- 40

setwd(dir = working_path)
source("/home/mattia/cfChromatin_in_MS/cfChromatin_fun.R")

#############################################################################################################################
###                                                    CATALOGUES CREATION                                                ###
#############################################################################################################################

###### TSS
# Catalog part 1
windows <- load_catalog(annotation_dir = annotation_dir_path,glossary = glossary_path, glossary_group = "ANATOMY", states = c("1_TssA", "2_TssAFlnk"))
# Catalog part 2 & 3
windows <- annotate_catalog(windows, feature_tag = "TSS", db = list('UCSC'=TxDb.Hsapiens.UCSC.hg38.knownGene, 'Ensembl'=EnsDb.Hsapiens.v86), db_gap = 2500, db_spe = TRUE, db_tag = "EXTRA_TSS", db_spe_center = 3000)
# Catalog part 4
windows <- add_flank(windows, feature_tag = "FLANKING", flanking_size = 1000)
windows <- add_background(windows, feature_tag = "BACKGROUND", bin_size = 5000)

# Statistics
show_statistics(windows, db_tag = "EXTRA_TSS")

# Save files
saveRDS(windows,paste0(working_path,"TSS_windows.rds"))

df_windows <- data.frame(
  seqnames=seqnames(windows),
  starts=start(windows)-1,
  ends=end(windows),
  name=windows$name,
  tissue=windows$tissue,
  type=windows$type)

write.table(df_windows, file=paste0(working_path,"TSS.bed"), quote=F, sep=",", row.names=F, col.names=F)

#############################################################################################################################
###                                                 PROCESSING SEQUENCING FILES                                           ###
#############################################################################################################################

### Process fastq files separetly to get BED files

histone_mark <- "H3K27ac"

# Get a list of all BED files in the directory that end with "H3K27ac_ChIP.bed"
bed_files <- list.files(path = bed_dir_path, pattern = paste0(histone_mark,"_ChIP\\.bed$"), full.names = TRUE)
# Define regions
windows <- readRDS(paste0(working_path,"TSS_windows.rds"))
# Use mclapply() to apply the function in parallel to the list of BED files
# Set mc.cores to the number of cores you want to use
results <- mclapply(bed_files, process_bed_file, regions = windows, mc.cores = core)
# Now 'results' is a list of GRanges objects
# You can convert it to a named list of GRanges objects with the file names as the names of the list elements
counts_list <- setNames(results, gsub("\\.bed$", "", basename(bed_files)))

saveRDS(counts_list, paste0(working_path,paste0(histone_mark,"_counts.rds")))

#############################################################################################################################
###                                                 ESTIMATING BACKGROUNG SIGNAL                                          ###
#############################################################################################################################

#counts_list <- readRDS(paste0(working_path,paste0(histone_mark,"_counts.rds")))

# Apply the estimateBackground function to each GRanges object in the list
# Set mc.cores to the number of cores you want to use
granges_list_with_background <- mclapply(counts_list, estimateBackground, tag = "BACKGROUND", background_cutOff = 4000, quantile_cutOff = 0.95, mc.cores = core)

saveRDS(granges_list_with_background, paste0(working_path,paste0(histone_mark,"_counts_bg.rds")))

#############################################################################################################################
###                                                     FEATURE LEVEL SIGNAL                                              ###
#############################################################################################################################

# Remove the samples from the GRanges list
exclude_samples <- c("CD19-GSE172116-1_H3K27ac_ChIP", "CD19-GSE172116-2_H3K27ac_ChIP")
granges_list_with_background_rm <- granges_list_with_background[!names(granges_list_with_background) %in% exclude_samples]

# Apply the calculationSignal function to each GRanges object in the list
# Set mc.cores to the number of cores you want to use
granges_list_geneSignal <- mclapply(granges_list_with_background_rm, calculationSignal, tag_list = c("TSS","EXTRA_TSS"), rm_feature_list=c("UNKNOWN"), mc.cores = core)
df_geneSignal <- t(as.data.frame(do.call(rbind, granges_list_geneSignal)))

saveRDS(granges_list_geneSignal, paste0(working_path,paste0(histone_mark,"_signal.rds")))
saveRDS(df_geneSignal, paste0(working_path,paste0(histone_mark,"_df_signal.rds")))

#############################################################################################################################
###                                                    CALCULATE SIZE FACTOR                                           ###
#############################################################################################################################

#df_geneSignal <- as.data.frame(readRDS(paste0(histone_mark,"_df_signal.rds")))

# Create different dataframe for each type of sample
plasma_sample <- c("12097-P_H3K27ac_ChIP", "14223-P_H3K27ac_ChIP", "16057-P-100_H3K27ac_ChIP", "16057-P-200_H3K27ac_ChIP", "16170-P_H3K27ac_ChIP", "19019-P_H3K27ac_ChIP")
csf_sample <- c("04061-C-100_H3K27ac_ChIP", "04061-C-200_H3K27ac_ChIP", "12097-C_H3K27ac_ChIP", "14223-C_H3K27ac_ChIP", "16170-C_H3K27ac_ChIP", "19019-C_H3K27ac_ChIP")
df_plasma <- df_geneSignal[,plasma_sample]
df_csf <- df_geneSignal[,csf_sample]

# Calculate the size factor
sf_plasma <- calculation_sf(df_plasma, ctrl=c("19019-P_H3K27ac_ChIP"), housekeeping=1000)
sf_csf <- calculation_sf(df_csf, ctrl=c("19019-C_H3K27ac_ChIP"), housekeeping=1000)

# Apply size factor to genes
df_plasma_geneSignal <- t(t(df_geneSignal[,names(sf_plasma)])*sf_plasma)
df_csf_geneSignal <- t(t(df_geneSignal[,names(sf_csf)])*sf_csf)

# Apply size factor to windows
sf_joined <- c(sf_plasma,sf_csf)
for(x in names(sf_joined)){
  no_bg <- granges_list_with_background[[x]]$counts-granges_list_with_background[[x]]$background
  # Remove background signal from each window and apply the size factor
  granges_list_with_background[[x]]$norm <- ifelse(no_bg > 0,no_bg*sf_joined[x],0)
}

saveRDS(df_plasma_geneSignal, paste0(working_path,paste0(histone_mark,"_plasma_df_signal_norm.rds")))
saveRDS(df_csf_geneSignal, paste0(working_path,paste0(histone_mark,"_csf_df_signal_norm.rds")))
saveRDS(granges_list_with_background, paste0(working_path,paste0(histone_mark,"_counts_norm.rds")))
























# Define a function to create a BigWig file from a GRanges object
createBigWig <- function(gr, filename) {
  # Create a GRanges object with the normalized signal as the score
  gr_with_score <- GRanges(seqnames = seqnames(gr),
                           ranges = ranges(gr),
                           strand = strand(gr),
                           score = mcols(gr)$normalizedSignal)
  
  # Export the GRanges object to a BigWig file
  export.bw(gr_with_score, filename)
}

# Apply the createBigWig function to each GRanges object in the list
filenames <- paste0("/date/gcb/gcb_MZ/Round1/R_bigiwig/", seq_along(normalizedSignal_counts), ".bw")  # replace with your actual filenames
mapply(createBigWig, normalizedSignal_counts, filenames)