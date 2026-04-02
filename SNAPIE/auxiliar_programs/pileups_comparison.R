# Compare read pileups between all samples
# Description: Make one figure for each specified region that compares all samples
# Output: New folder in specified directory with a plot for each important region

#setwd('/Users/prc992/Desktop/DFCI/1-cfpipeline/2-Pileup/3-R-SetReadPileUps')

library(Gviz)
library(rtracklayer)
library(org.Hs.eg.db)

# Function to create tuple with files and name of the samples
create_list_namesamples_files <- function(strDir) {
  list_files <- list.files(strDir,pattern="*treat_pileup.bdg.bw")
  list_NameSamplesFiles <- lapply(list_files, function(str) {
    charDot <- regexpr("\\_", str)
    if (charDot == -1) {
      # If there is no point in the string consider it complete
      substring_ate_ponto <- str
    } else {
      # Otherwise get the first characters until the firs dot
      substring_ate_ponto <- substr(str, 1, charDot - 1)
    }
    lista <- list(substring_ate_ponto,str)
    return(lista)
  })
  return(list_NameSamplesFiles)
}

########### USER INPUT (only change variables in this section) ##############

#samples_dir <- "read_pileups/"
#bed_path <- "test_housekeeping.bed" # bed file of important regions to plot
#chrom_sizes_path <- "hg19_clean.chrom.sizes" # chromosome sizes for this genome model
#save_path <- "" #path to save the plots (if saving in working directory just set to "")
#buffer <- 1000 # default = 500, can change with user input
#gen <- "hg19" # genome model

args <- commandArgs( trailingOnly = TRUE )
samples_dir = args[1] # path to sample bigwig
bed_path = args[2]  # bed file of important regions to plot
chrom_sizes_path = args[3] # chromosome sizes for this genome model
gen = args[4] # genome model

buffer <- 1000 # default = 500, can change with user input
#############################################################################

# read in files/folders
samples_list <- list.files(samples_dir,pattern="*.bw")
sample_tuples <- create_list_namesamples_files(samples_dir)
bed_con <- file(bed_path, open="r") # bed file of specific regions to snapshot
chrom_sizes <- read.table(chrom_sizes_path, header = FALSE, sep="\t")

# make folder to save plots in
#system(paste("mkdir ", save_path, "pileup_comparison_plots", sep=""))

gtrack <- GenomeAxisTrack() 

# go through each region in bed file
while (TRUE) {
  # get region info
  line = readLines(bed_con, n=1)
  if (identical(line, character(0))) { # no more lines
    break
  }
  tokens <- strsplit(line, "\t")
  chr <- tokens[[1]][1]
  start <- strtoi(tokens[[1]][2])
  end <- strtoi(tokens[[1]][3])
  
  # ensure buffer within chromosome bounds
  start_buffer <- start - buffer
  end_buffer <- end + buffer
  end_chrom <- chrom_sizes[which(chrom_sizes$V1 == chr), "V2"]  
  if (start_buffer < 0) {
    start_buffer <- 0
  } else if (end_buffer > end_chrom) {
    end_buffer <- end_chrom
  }
  
  itrack <- IdeogramTrack(genome = gen, chromosome = chr)
  
  # set up sample datatracks list
  dt_list = list()
  max_height <- 0
  for (sample_file in samples_list) {
    # get sample name
    tup_index <- which(sapply(sample_tuples, function(x) sample_file %in% x))
    sample_name <- paste("", sample_tuples[[tup_index]][1], sep="")

    # make datatrack from sample
    sample_path <- paste(samples_dir, sample_file, sep = "")
    sample_gr <- import.bw(sample_path,as="GRanges") # bigwig file for sample
    sample_dt <- DataTrack(sample_gr, genome = gen, chromosome = chr, name = sample_name)
    dt_list <- append(dt_list, sample_dt) #add sample_dt to list of sample datatracks
    
    # get highest peak value
    #isolate this region of bigwig file
    interval_ir <- IRanges(start=start_buffer, width=end_buffer-start_buffer+1)
    chr_gr <- sample_gr[seqnames(sample_gr) == chr]
    chr_ir <- ranges(chr_gr)
    overlaps <- findOverlaps(chr_ir, interval_ir)
    begin_interval <- queryHits(overlaps)[1]
    end_interval <- queryHits(overlaps)[length(overlaps)]
    region_gr <- chr_gr[begin_interval:end_interval]
    
    #get highest peak (score) in region
    score_list <- mcols(region_gr)$score
    high_score <- max(score_list)
    #check if higher than max value of all samples
    if (high_score > max_height) {
      max_height <- high_score
    }
  }
  
  # NCBI RefSeq track
  ucscGenes <- UcscTrack(genome=gen, chromosome = chr, table="ncbiRefSeq", 
                         track = 'NCBI RefSeq', trackType="GeneRegionTrack", 
                         rstarts = "exonStarts", rends = "exonEnds", gene = "name", 
                         symbol = 'name', transcript = "name", strand = "strand", 
                         name = "RefSeq", stacking = 'pack', showID = T, 
                         geneSymbol = T, transcriptAnnotation="symbol")
  z <- ranges(ucscGenes)
  mcols(z)$symbol <- mapIds(org.Hs.eg.db, gsub("\\.[1-9]$", "", mcols(z)$symbol), "SYMBOL","REFSEQ")
  ucscGenes_named <- ucscGenes
  ranges(ucscGenes_named) <- z
  
  # highlight track for region without buffer
  htrack <- HighlightTrack(trackList = c(list(gtrack, ucscGenes_named), dt_list), 
                           start = start, width = end-start, chromosome = chr, 
                           col = "azure3", fill = "azure2")
  
  # plot region and save
  plot_name <- paste("comparison_",chr,"_",start,"_",end,".pdf", sep="")
  num_samples = length(samples_list)
  pdf(plot_name, height=1+num_samples*1.5)
  size_list <- c(list(0.5,1,1), as.list(rep(2, num_samples)))
  plotTracks(list(itrack,htrack), from=start_buffer,to=end_buffer, 
             chromosome=chr, type="hist", ylim = c(0, max_height + (max_height * 0.2)), 
             title.width=1.5, sizes=size_list, col.title="gray23", col.axis="gray23")
  dev.off()
}
close(bed_con)
