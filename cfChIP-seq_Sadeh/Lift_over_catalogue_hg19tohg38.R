library(rtracklayer)
library(easylift)
TSS.windows_hg19 <- readRDS("/date/gcb/gcb_MZ/Analysis/cfChIP-seq/SetupFiles/H3K4me3/Windows.rds")
genome(TSS.windows_hg19) <- "hg19"

# replace "hg38" with the target genome assembly
# replace "path/to/your/hg19ToHg38.over.chain.gz" with the path to your chain file
TSS.windows_hg38 <- easylift(TSS.windows_hg19, to = "hg38", chain = "/date/gcb/gcb_MZ/hg19ToHg38.over.chain.gz")
# replace "hg19" with the source genome assembly

saveRDS(TSS.windows_hg38,"/date/gcb/gcb_MZ/Analysis/cfChIP-seq/SetupFiles/H3K27ac_hg38/Windows.rds")
