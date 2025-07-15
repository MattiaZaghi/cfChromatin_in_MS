#!/usr/bin/env Rscript --vanilla

# The shebang line above allows this script to be executed directly from the command line on Unix/Linux systems.
# The '--vanilla' option starts R without saved or site profiles, ensuring a clean session.

# Retrieve command-line arguments. 'trailingOnly = FALSE' includes script-related arguments.
initial.options <- commandArgs(trailingOnly = FALSE)

# Check if the script is run in an interactive mode (e.g., from RStudio) or not.
if(!any(grepl("--interactive", initial.options))) {
  # If not in interactive mode, extract the script name and determine the source directory.
  file.arg.name <- "--file="
  script.name <- sub(file.arg.name, "", initial.options[grep(file.arg.name, initial.options)])
  SourceDIR <- paste0(dirname(script.name),"/")
  DataDir = paste0(getwd(), "/")  # Use the current working directory as the data directory.
} else {
  # If in interactive mode, set the source and data directories explicitly for development purposes.
  SourceDIR = "/date/gcb/gcb_MZ/Analysis/cfChIP-seq/"
  Mod="H3K4me3_rodmap_hg19"
  DataDir = paste0(SourceDIR, "SetupFiles/", Mod,"/")
  TargetMod = "H3K4me3_rodmap_hg19"  # Specify the target modification for analysis.
}

# Load necessary libraries, suppressing warnings that might arise during loading.
suppressWarnings(library(ggplot2))
suppressWarnings(library(preprocessCore))
suppressWarnings(library(Matrix))
suppressWarnings(library(reshape2))
suppressWarnings(library(rtracklayer))

# Set the annotation directory to the data directory.

ANNOTDIR = DataDir

# Load TSS window data from an RDS file located in the annotation directory.
TSS.windows.filename = paste0(ANNOTDIR,"Windows.rds")
TSS.windows = readRDS(TSS.windows.filename)
# Retrieve genome sequence information from the TSS windows.
genome.seqinfo = seqinfo(TSS.windows)
# Define a list of chromosomes to be analyzed.
ChrList = paste0("chr", c(1:22,"X", "Y"))
# If the target modification is H3K27ac, adjust the chromosome list based on the hg38 genome.
if(TargetMod == "H3K4me3_sCER") {
  hg19.seqinfo = Seqinfo(genome="hg19")
  ChrList = (seqnames(hg19.seqinfo))
}

# Print a message indicating the start of the background model building process.
BackgroundModel.filename = paste0(DataDir,"BackgroundModel.rds")
print("Building Background model")

# Define a function to sort genomic ranges by chromosome and start position.
sortGR = function(GR) {
  GR.order = lapply(unique(chrom(GR)), function(chr) {
    W = which(as.vector(chrom(GR)) == chr)
    W[order(start(GR[W]))]
  })
  GR.order = do.call("c", GR.order)
  GR[GR.order]  
}

# Define a function to build overlapping tiles for genome coverage analysis.
# 'width' specifies the width of each tile, and 'jump' specifies the step size for overlapping tiles.
buildOverlapingTiles = function(width, jump = width/4) {
  # Generate tiles for the entire genome based on the specified width.
  T = tileGenome(hg19.seqinfo, tilewidth = width, cut.last.tile.in.chrom = TRUE)
  # Filter tiles to include only those in the specified chromosome list.
  T = T[chrom(T) %in% ChrList]
  T0 = T
  mWidth = min(width(T0))
  # Create additional tiles based on the jump parameter to ensure overlapping coverage.
  for(o in seq(jump, width-1, by=jump)) {
    T = suppressWarnings(c(T, trim(shift(T0,o))))
  }
  # Filter tiles to ensure they meet the minimum width requirement.
  T = T[width(T) >= mWidth]
  # Sort the tiles by chromosome and start position.
  T = sortGR(T)
  T
}

# Define the width of background regions to be 10 million base pairs.
Background.Regions.width =  10**7
# Generate overlapping tiles for background regions using the specified width.
Background.Regions = buildOverlapingTiles(Background.Regions.width)
# Count the number of background regions generated.
Background.Regions.num = length(Background.Regions)
# Disjoin the background regions to ensure they are non-overlapping and unique.
Background.Uniq.Regions = disjoin(Background.Regions, with.revmap=TRUE)

# Define the width of individual background intervals to be 5 million base pairs.
Background.Inds.width = 5*10**6
# Generate overlapping tiles for individual background intervals using the specified width.
Background.Inds = buildOverlapingTiles(Background.Inds.width)
# Count the number of individual background intervals generated.
Background.Inds.num = length(Background.Inds)
# Disjoin the individual background intervals to ensure they are non-overlapping and unique.
Background.Uniq.Inds = disjoin(Background.Inds, with.revmap=TRUE)

# Extract the chromosome names for each of the background regions.
Background.RegionChr = as.character(chrom(Background.Regions))
# Find the nearest individual background interval for each background region, centering the region.
Background.IndRegion = nearest(Background.Inds, subject = resize(Background.Regions, width = 10, fix = "center"))

# Extract the chromosome names for each TSS window.
WinChr = as.character(chrom(TSS.windows))

# For each background region, find overlapping TSS windows.
Background.RegionWindows =  lapply(1:Background.Regions.num, function(i) subjectHits(findOverlaps(Background.Regions[i], TSS.windows)))
# For each individual background interval, find overlapping TSS windows.
Background.IndWindows = lapply(1:Background.Inds.num, function(i) subjectHits(findOverlaps(Background.Inds[i], TSS.windows)))

# For each TSS window, find the nearest unique individual background interval.
Background.WindowInd = Rle(nearest(TSS.windows, Background.Uniq.Inds))
# For each TSS window, find the nearest unique background region.
Background.WindowRegion = Rle(nearest(TSS.windows, Background.Uniq.Regions))

# Save the background model and related information to an RDS file.
saveRDS(list(Background.Regions.width = Background.Regions.width,
             Background.Regions = Background.Regions,
             Background.Regions.num = Background.Regions.num,
             Background.Uniq.Regions = Background.Uniq.Regions,
             Background.Inds.width = Background.Inds.width,
             Background.Inds = Background.Inds,
             Background.Inds.num = Background.Inds.num,
             Background.Uniq.Inds = Background.Uniq.Inds,
             Background.RegionChr = Background.RegionChr,
             Background.IndRegion = Background.IndRegion,
             WinChr = WinChr,
             Background.RegionWindows = Background.RegionWindows,
             Background.IndWindows = Background.IndWindows,
             Background.WindowInd = Background.WindowInd,
             Background.WindowRegion = Background.WindowRegion),
        BackgroundModel.filename)

# Compute the mapping from windows to genes and save the information.
Win2Gene.Matrix.filename = paste0(DataDir,"Win2Gene.rds")
print("Computing Window to Gene mapping")

# Identify TSS windows associated with genes.
GeneWindows = which(TSS.windows$name != "." & !is.na(TSS.windows$name))
# Extract unique gene names from the TSS windows.
GeneLists = as.character(unique(TSS.windows[GeneWindows]$name))
GeneLists = GeneLists[GeneLists != "."]
# Split gene names by ";" and flatten the list to get a unique list of genes.
Genes = unique(do.call("c", strsplit(GeneLists, ";")))
Genes = Genes[!is.na(Genes)]
# Split gene names associated with each window.
GL = strsplit(as.character(TSS.windows[GeneWindows]$name), ";")
names(GL) = GeneWindows

# Compute the length of gene lists for each window.
Ls = lapply(GL, length)
# Assign IDs to each unique gene.
GenesId = 1:length(Genes)
names(GenesId) = Genes
# Find the maximum length of gene lists.
max.L =  max(do.call("c", Ls))
# Initialize vectors to store row and column indices for the sparse matrix.
rs = c()
cs = c()
# Populate the row and column indices based on gene list lengths.
for(i in 1:max.L) {
  rs = c(rs, rep(which(Ls == i), each = i))
  cs = c(cs, GenesId[do.call("c", GL[which(Ls == i)])])
}

# Create a sparse matrix representing the window-to-gene mapping.
Win2Gene.matrix = sparseMatrix(i=as.integer(rs), j=cs, x = 1, dims = c(length(GeneWindows), length(Genes)), dimnames = list(as.character(GeneWindows), Genes))

# Identify genes associated with multiple promoters.
MultiPromoterGenes = sapply(Genes, function(g) {
  I = which(Win2Gene.matrix[,g] > 0)
  if(length(I) > 0)
    max(I) - min(I) > length(I) - 1
  else
    NA
})

# Save the window-to-gene mapping and multi-promoter genes information to an RDS file.
saveRDS(list(Matrix = Win2Gene.matrix, 
             Multi = MultiPromoterGenes), Win2Gene.Matrix.filename)