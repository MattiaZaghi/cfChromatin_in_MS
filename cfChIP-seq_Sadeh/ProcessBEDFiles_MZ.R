#!/usr/local/bin/Rscript --vanilla
#

catn = function(...) { cat(...,"\n") }


initial.options <- commandArgs(trailingOnly = FALSE)
trailing.options = commandArgs(trailingOnly = TRUE)
#Find where are the scripts at

if( !any(grepl("--interactive", initial.options)) ) {
  file.arg.name <- "--file="
  script.name <- sub(file.arg.name, "", initial.options[grep(file.arg.name, initial.options)])
  SourceDIR <- paste0(dirname(script.name),"/")
  DataDir = paste0(getwd(), "/")
  ANNOTDIR = DataDir
  DevelopmentMode = FALSE
} else {
  # development mode (inside RStudio) we can set options using these variables
  #
  # Assumes:
  #   ~/BloodChIP/Src points to source directory
  #   ~/BloodChIP/Data points to root above Analysis directories 
  #   
  # It is highly recommended to reset the environment ("Clear Objects" in RStudio) prior to running
  #
  
  
  DevelopmentMode = TRUE
  Mod = "H3K4me3/"
  SourceDIR = "~/BloodChIP/Src/"
  ANNOTDIR = paste0(SourceDIR, "SetupFiles/", Mod)
  trailing.options = paste(
    "-r ~/BloodChIP/Data/Analysis-Paper",
    "-p ~/BloodChIP/Data/Analysis-Paper",
    "-o",  "~/BloodChIP/Data/Analysis-Paper/Figures/Liver/",
    "--plotprograms=IHEC-sig",
    paste0("--geneprograms=",
           "~/BloodChIP/Data/Atlases/Programs/MSigDB-curated.csv",
           ",",
           "~/BloodChIP/Data/Analysis-Paper/Figures/BLUEprint/BLUEprint-RNA-sig.csv"
    ),
    "--select=~/BloodChIP/Data/Analysis-Paper/Figures/Liver/Select-example.txt",
    "L001.1 M002.1 L010.1 L011.1 H001.1 H002.1 H003.1 H004.1",
    
    #    "--writeforce",
    #    "--plotprograms=Fig5E-prog",
    #    paste0("--select=",
    #           "~/BloodChIP/Data/Analysis-Paper/Figures/Pathways/Select-programs.csv"),
    #    "-o", "~/BloodChIP/Data/Analysis-Paper/Figures/Pathways",
    #    paste0("--geneprograms=", 
    #           "~/BloodChIP/Data/Atlases/Programs/MSigDB-curated-small.csv"
    #"~/BloodChIP/Data/Analysis-Paper/Figures/BLUEprint/BLUEprint-RNA-sig.csv",
    #",",
    #"~/BloodChIP/Data/Analysis-Paper/Figures/CRC/TCGA_CMS-gene-sig.csv",
    #",",
    #"~/BloodChIP/Data/Analysis-Paper/Figures/CRC/CancerGroups-Sig.csv",
    #",",
    #"~/BloodChIP/Data/Analysis-Paper/Figures/CRC/ColonClassifier-crca786.csv",
    #",",
    #"~/BloodChIP/Data/Analysis-Paper/Figures/CRC/ColonClassifier-CRIS.csv",
    #  ""           
    #           ),
    #    paste0("--geneconsensus=",
    #           "~/BloodChIP/Data/Analysis-Paper/Figures/BLUEprint/BLUEprint-RNA-sig.rds"),
    #    "~/BloodChIP/Data/Analysis-Paper/Samples-Healthy-reference-cohort.txt",
    #    "H001.1", "M002.1", "M001.1", "C001.2746", "C001.2752", "C040.3606",
    #      "~/BloodChIP/Data/Analysis-Paper/Figures/Samples-Fig5E.txt",
    "" 
  )
  trailing.options = strsplit(trailing.options, " ")[[1]]
  catn(trailing.options)
}


OrigDir = getwd()

extendDir <- function(x) {
  a = substr(x,1,1)
  if( a != "/" && a != "~")
    x = paste0(OrigDir,"/",x)
  a = substr(x, nchar(x), nchar(x))
  if( a != "/")
    x = paste0(x,"/")
  x
}


outputFile = function(x) {
  a = substr(x,1,1)
  if( a != "/" && a != "~")
    x = paste0(extendDir(OutputDir), x)
  x
}

suppressMessages(library("optparse"))
option_list = list(
  # directory 
  make_option(c("-d", "--datadir"), type = "character", default = NULL,
              help="data directory"),
  make_option(c("-b", "--BEDdir"), type="character", default=NULL,
              help="location of BED files" ),
  make_option(c("-a", "--annotationdir"), type="character", default=NULL,
              help="location of annotation files" ),
  make_option(c("-t", "--trackdir"), type="character", default=NULL,
              help="location of genome browser track files" ),
  make_option(c("-o", "--outputdir"), type="character", default=NULL,
              help="location of output files" ),
  make_option(c("-m", "--mod"), type="character", default="H3K4me3",
              help="name of modification, used in combination with -r" ),
  make_option(c("-r", "--root"), type="character", default=NULL,
              help="name of root directory, input is assumed to be in Root/XXX/mod" ),
  make_option(c("-p", "--project"), type="character", default=NULL,
              help="name of project root directory, output is assumed to be in Project/XXX/mod" ),
  
  # basic processing   
  make_option(c("-B", "--background"),  action = "store_true", type="logical", default=FALSE, 
              help="Compute background model"),
  make_option(c("-C", "--count"),  action = "store_true", type="logical", default=FALSE, 
              help = "Count gene coverage"),
  make_option(c("-N", "--normalize"),  action = "store_true", type="logical", default=FALSE, 
              help = "Normalize gene counts"),
  make_option(c("-O", "--overexpressed"),  action = "store_true", type="logical", default=FALSE, 
              help = "Evaluate overexpressed genes"),
  make_option(c("-A", "--all"),  action = "store_true", type="logical", default=FALSE, 
              help="Compute all relevant data for each file (equivalent to -BCN)"),
  make_option(c("-F", "--force"),  action = "store_true", type="logical", default=FALSE, 
              help="Force recomputing"),
  make_option("--hardforce",  action = "store_true", type="logical", default=FALSE, 
              help="Force recomputing from BED file"),
  make_option(c("-X", "--extra"), action="store_true",type="logical", default=FALSE, 
              help = "Run all typical outputs per file"),
  make_option("--writeforce",  action = "store_true", type="logical", default=FALSE, 
              help="Force writing of output files such as tracks and plots")
); 

opt_parser = OptionParser(option_list=option_list);
if( exists("trailing.options") ) {
  opt = parse_args(opt_parser, args = trailing.options, positional_arguments = TRUE)
} else
  opt = parse_args(opt_parser, positional_arguments = TRUE)

Files = opt$args

if(!DevelopmentMode && length(Files) == 0)
  parse_args(opt_parser,args = "-h", print_help_and_exit = TRUE)

catn("Initializing")
suppressMessages(library(rtracklayer))
suppressMessages(library(tools))
suppressMessages(library(ggplot2))

TargetMod = opt$options$mod
if( !is.null(opt$options$root)) {
  RootDir =  extendDir(opt$options$root)
  DataDir = paste0(RootDir, "Samples/", TargetMod, "/")
  BedDir = paste0(RootDir, "BED/", TargetMod, "/")
  #  TracksDIR = paste0(RootDir, "Tracks/", TargetMod, "/")
  #  OutputDir = paste0(RootDir, "Output/", TargetMod, "/")
} 
if( !is.null(opt$options$project)) {
  ProjectDir =  extendDir(opt$options$project)
  TracksDIR = paste0(ProjectDir, "Tracks/", TargetMod, "/")
  OutputDir = paste0(ProjectDir, "Output/", TargetMod, "/")
} 
ANNOTDIR = paste0(SourceDIR, "SetupFiles/", TargetMod, "/")


if( !is.null(opt$options$datadir) )
  DataDir = extendDir(opt$options$datadir)

DIR = DataDir

if( !is.null(opt$options$annotationdir) )
  ANNOTDIR = extendDir(opt$options$annotationdir)

if( !exists("BedDir") ) 
  BedDir = DataDir


if( !is.null(opt$options$BEDdir) )
  BedDir = extendDir(opt$options$BEDdir)

if( !exists("TracksDIR") ) 
  TracksDIR = DataDir
if( !is.null(opt$options$trackdir) )
  TracksDIR = extendDir(opt$options$trackdir)

if( !is.null(opt$options$outputdir)) {
  OutputDir = extendDir(opt$options$outputdir)
} else
  if( !exists("OutputDir") ) {
    if( !DevelopmentMode  ) {
      OutputDir = OrigDir
    } else
      OutputDir = paste0("~/Data/BloodChIP/Output/", Mod)
  }

MetaDIR = OutputDir

if( !exists("SetupDIR") )
  SetupDIR = ANNOTDIR
# catn("setupDir is: ", SetupDIR)


#initialize variables
{
  TSS.windows = readRDS(paste0(SetupDIR,"Windows.rds"))
  genome.seqinfo = seqinfo(TSS.windows)
  ChrList = paste0("chr", c(1:22,"X", "Y"))
  if (TargetMod == "H3K4me3-scer") {
    #    sacCer3.seqinfo = Seqinfo(genome="sacCer3")
    ChrList = (seqnames(genome.seqinfo))
  }
  
  
  source(paste0(SourceDIR,"Background.R"))

  source(paste0(SourceDIR, "cfChIP-Functions_MZ.R"))  
  source(paste0(SourceDIR, "YieldEstimation-Functions.R"))

  

  

  WinDescription = mcols(TSS.windows)
  colnames(WinDescription) = sapply(colnames(WinDescription), function(x) paste0(toupper(substr(x,1,1)), substr(x,2,width(x))))
  colnames(WinDescription) = sub("Name", "Gene", colnames(WinDescription))
  rownames(WinDescription) = 1:length(TSS.windows)

  
}


params = cfChIP.Params()
params$DataDir = DataDir

RetainInMemory = FALSE

if( !is.null(opt$options$config) ) {
  config.filename = opt$options$config
} else
  config.filename = paste0(SetupDIR, "config.csv")

if( file.exists(config.filename) ) {
  catn("Reading configuration options from ",config.filename)
  config.data = read.csv(config.filename,as.is = TRUE, header = FALSE)
  for( i in 1:nrow(config.data) ) {
    params[[config.data[i,1]]] = config.data[i,2]
  }
}

params$Background = opt$options$background
params$GeneCounts = opt$options$count
params$Normalize = opt$options$normalize
params$OverExpressedGenes = opt$options$overexpressed
doOverExpressedGenes = opt$options$overexpressed
if( opt$options$extra ) {
  opt$options$all = TRUE
  opt$options$tracks = TRUE
  opt$options$signatures = TRUE
  opt$options$signaturesvshealthy  = TRUE
  opt$options$programs = TRUE
  opt$options$backgroundplot = TRUE
  opt$options$fraglenplot = TRUE
  #  opt$options$meta = TRUE
  #  opt$options$enrichR = TRUE
}

if( opt$options$all ) {
  params$Background = TRUE
  params$GeneBackground = TRUE
  params$GeneCounts = TRUE
  params$Normalize = TRUE
  params$OverExpressedGenes = TRUE
}

doHardForce = opt$options$hardforce
doForce = opt$options$force || doHardForce
doWriteForce = opt$options$force || opt$options$writeforce


BaseFileName <- function( fname, 
                          extList = c(".gz$", ".bed$",".rdata$", ".bw$", ".tagAlign$", "-H3K4me3") ) {
  #  x = file_path_sans_ext(fname)
  x = fname
  for( ext in extList )
    x = sub(ext, "", x)
  
  y = strsplit(x,"/")[[1]]
  n = length(y)
  z = y[n]
  return(z)
}

ProcessBEDFile  = function( BFile ) {
  if( grepl("^XXX", BFile))
    return(NULL)
  
  #catn(BFile)
  if( !file.exists(BFile) &&
      !file.exists(paste0(BFile,".bed")) &&
      !file.exists(paste0(BFile,".tagAlign")) && 
      !file.exists(paste0(BFile,".bw")) &&
      !file.exists(paste0(BFile,".bed.gz")) &&
      !file.exists(paste0(BFile,".tagAlign.gz")) )
    BFile = paste0(BedDir, BFile)
  
  catn(BFile)
  
  dat = cfChIP.ProcessFile(filename = BFile, param = params, Force = doForce, HardForce = doHardForce)
  
  
  if( RetainInMemory ) {
    # get read of the coverage since it takes too much memory....
    if(!DevelopmentMode) {
      dat$Cov = NULL
      dat$BED = NULL
      dat$BW = NULL
    }
    return(dat) 
  } else
    rm(dat)
}


ProcessBEDFileList = function( BFlist ) {
  L = lapply(BFlist, ProcessBEDFile)
  names(L) = sapply(L, function(l) l$Name)
  
  L  
}

expandFiles = function(f) {
  if( grepl(".txt$", f) | grepl(".csv$", f) ) {
    if(!file.exists(f)) {
      catn("Missing list file", f)
      return(NULL)
    }
    catn("Reading sample list from",f)
    File.list = read.table(f, as.is = TRUE)
    File.list = File.list[,1]
    names(File.list) = NULL
    as.list(File.list)
  } else
    f
}

FullFiles = unique(unlist(lapply(Files,expandFiles)))

LL = ProcessBEDFileList(FullFiles)
LL = LL[!sapply(LL, is.null)]


