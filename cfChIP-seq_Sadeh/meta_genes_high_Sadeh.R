library(plyranges)


meta_genes<-read_bed("/date/gcb/gcb_MZ/Analysis/cfChIP-seq/SetupFiles/H3K27ac/Meta-genes.bed")

# Filter rows where 'name' contains 'high'
meta_genes_high <- meta_genes[grepl("high", mcols(meta_genes)$name, ignore.case = TRUE)]

write_bed(meta_genes_high,"/date/gcb/gcb_MZ/Analysis/cfChIP-seq/SetupFiles/H3K27ac/meta_genes_high.bed")
