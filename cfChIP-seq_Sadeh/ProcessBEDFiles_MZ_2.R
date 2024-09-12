#!/usr/bin/env Rscript
args <- commandArgs(trailingOnly = TRUE)
bed_file <- args[1]  # Get the bed file from the command-line arguments

# process_bed_files.R
catn = function(...) { cat(...,"\n") }

catn("Initializing")
suppressMessages(library(rtracklayer))
suppressMessages(library(tools))
suppressMessages(library(ggplot2))

SourceDIR <-"/date/gcb/gcb_MZ/Analysis/cfChIP-seq/"
SetupDIR<-"/date/gcb/gcb_MZ/Analysis/cfChIP-seq/SetupFiles/H3K27ac_hg38/"
RootDir =  "/date/gcb/gcb_MZ/Analysis/"
TargetMod="H3K27ac_hg38"
DataDir = paste0(RootDir, "Samples/", TargetMod, "/")
BedDir = paste0(RootDir, "BED/", TargetMod, "/")

TSS.windows = readRDS(paste0(SetupDIR,"Windows.rds"))
genome.seqinfo = seqinfo(TSS.windows)
ChrList = paste0("chr", c(1:22,"X", "Y","M"))

BaseFileName <- function( fname, extList = c(".gz$", ".bed$",".rdata$", ".bw$", ".tagAlign$", "-H3K4me3") ) {
  x = fname
  for( ext in extList )
    x = sub(ext, "", x)
  
  y = strsplit(x,"/")[[1]]
  n = length(y)
  z = y[n]
  return(z)
}

source("/date/gcb/gcb_MZ/Analysis/cfChIP-seq/cfChIP-Functions_MZ.R")

# Process the bed file
bed_file_path <- file.path(BedDir, bed_file)  # Add this line
dat <- cfChIP.ProcessFile(filename = bed_file_path, param = cfChIP.Params())  # Modify this line
