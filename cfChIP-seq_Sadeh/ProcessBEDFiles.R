#!/usr/local/bin/Rscript --vanilla
#if executing in development 
setwd("/date/gcb/gcb_MZ/Analysis/")
# Define a custom function 'catn' that behaves like 'cat' but appends a newline at the end of the output.
catn = function(...) { cat(...,"\n") }

# Retrieve the command-line arguments passed to the R script.
# 'trailingOnly = FALSE' retrieves all arguments, including script name and R-specific args.
initial.options <- commandArgs(trailingOnly = FALSE)
# 'trailingOnly = TRUE' retrieves only the arguments after the script name.
trailing.options = commandArgs(trailingOnly = TRUE)

# Determine the mode of operation based on the presence of '--interactive' in the initial options.
if( !any(grepl("--interactive", initial.options)) ) {
  # If not in interactive mode, extract the script name and set directories based on its location.
  file.arg.name <- "--file="
  script.name <- sub(file.arg.name, "", initial.options[grep(file.arg.name, initial.options)])
  SourceDIR <- paste0(dirname(script.name),"/")
  DataDir = paste0(getwd(), "/")
  ANNOTDIR = DataDir
  DevelopmentMode = FALSE
} else {
  # In development mode (e.g., running inside RStudio), set directories and options manually.
  DevelopmentMode = TRUE
  Mod = "H3K27ac/"
  SourceDIR = "/date/gcb/gcb_MZ/Analysis/cfChIP-seq/"
  ANNOTDIR = paste0(SourceDIR, "SetupFiles/", Mod)
  # Simulate command-line arguments for development mode.
  trailing.options = paste(
    "-r /date/gcb/gcb_MZ/Analysis/",
    "-p /date/gcb/gcb_MZ/Analysis/",
    "-o",  "/date/gcb/gcb_MZ/Analysis/Output",
    "--plotprograms=IHEC-sig",
    paste0("--geneprograms=",
           "~/BloodChIP/Data/Atlases/Programs/MSigDB-curated.csv",
           ",",
           "~/BloodChIP/Data/Analysis-Paper/Figures/BLUEprint/BLUEprint-RNA-sig.csv"
    ),
    "--select=~/BloodChIP/Data/Analysis-Paper/Figures/Liver/Select-example.txt",
    "L001.1 M002.1 L010.1 L011.1 H001.1 H002.1 H003.1 H004.1",
    ""
  )
  # Split the simulated command-line arguments into a list.
  trailing.options = strsplit(trailing.options, " ")[[1]]
  # Print the trailing options to the console for verification.
  catn(trailing.options)
}

# Save the current working directory to revert back if needed.
OrigDir = getwd()

# Function to ensure directory paths are absolute and end with a slash.
extendDir <- function(x) {
  a = substr(x,1,1)
  if( a != "/" && a != "~")
    x = paste0(OrigDir,"/",x)
  a = substr(x, nchar(x), nchar(x))
  if( a != "/")
    x = paste0(x,"/")
  x
}

# Function to construct the full path for output files, ensuring they are placed in the correct directory.
outputFile = function(x) {
  a = substr(x,1,1)
  if( a != "/" && a != "~")
    x = paste0(extendDir(OutputDir), x)
  x
}

# Load the 'optparse' library for parsing command-line options.
suppressMessages(library("optparse"))

# Define the list of command-line options available for the script.
option_list = list(
  # Directory options
  make_option(c("-d", "--datadir"), type = "character", default = NULL, help="data directory"),
  make_option(c("-b", "--BEDdir"), type="character", default=NULL, help="location of BED files"),
  make_option(c("-a", "--annotationdir"), type="character", default=NULL, help="location of annotation files"),
  make_option(c("-t", "--trackdir"), type="character", default=NULL, help="location of genome browser track files"),
  make_option(c("-o", "--outputdir"), type="character", default=NULL, help="location of output files"),
  make_option(c("-m", "--mod"), type="character", default="H3K27ac", help="name of modification, used in combination with -r"),
  make_option(c("-r", "--root"), type="character", default=NULL, help="name of root directory, input is assumed to be in Root/XXX/mod"),
  make_option(c("-p", "--project"), type="character", default=NULL, help="name of project root directory, output is assumed to be in Project/XXX/mod"),
  
  # Basic processing options
  make_option(c("-B", "--background"),  action = "store_true", type="logical", default=FALSE, help="Compute background model"),
  make_option(c("-C", "--count"),  action = "store_true", type="logical", default=FALSE, help = "Count gene coverage"),
  make_option(c("-N", "--normalize"),  action = "store_true", type="logical", default=FALSE, help = "Normalize gene counts"),
  make_option(c("-O", "--overexpressed"),  action = "store_true", type="logical", default=FALSE, help = "Evaluate overexpressed genes"),
  make_option(c("-A", "--all"),  action = "store_true", type="logical", default=FALSE, help="Compute all relevant data for each file (equivalent to -BCN)"),
  make_option(c("-F", "--force"),  action = "store_true", type="logical", default=FALSE, help="Force recomputing"),
  make_option("--hardforce",  action = "store_true", type="logical", default=FALSE, help="Force recomputing from BED file"),
  
  # action per file 
  make_option(c("-X", "--extra"), action="store_true",type="logical", default=FALSE, 
              help = "Run all typical outputs per file"),
  make_option(c("-T", "--tracks"),  action = "store_true", type="logical", default=FALSE, 
              help = "Create normalized genome browser tracks"),
  make_option(c("-M", "--meta"),  action = "store_true", type="logical", default=FALSE, 
              help = "Generate meta plots"),
  make_option(c("-E", "--enrichR"),  action = "store_true", type="logical", default=FALSE, 
              help = "Evaluate enrichments of overexpressed genes"),
  make_option(c("-S", "--signatures"),  action = "store_true", type="logical", default=FALSE, 
              help = "Evaluate signatures of cell types"),
  make_option("--signaturesvshealthy",  action = "store_true", type="logical", default=FALSE, 
              help = "Evaluate signatures of cell types vs healthy samples"),
  
  make_option(c("-P", "--programs"),  action = "store_true", type="logical", default=FALSE, 
              help = "Evaluate expression programs"),
  make_option("--ctDNA",  action = "store_true", type="logical", default=FALSE, 
              help = "Estimate ctDNA fraction"),
  make_option("--backgroundplot",type="logical", default=FALSE, action="store_true", 
              help = "Plot background coverage in output directory"),
  make_option("--fraglenplot",type="logical", default=FALSE, action="store_true", 
              help = "Plot fragment length distribution in output directory"),
  make_option("--writeforce",  action = "store_true", type="logical", default=FALSE, 
              help="Force writing of output files such as tracks and plots"),
  
  # global output
  make_option(c("--consensus"),  type="character", default=NULL, 
              help = "Compute consensus (avg after re-normalization) and write to specified file"),
  make_option("--winconsensus",  type="character", default=NULL, 
              help = "Compute consensus for window signatures and write to specified file"),
  make_option("--geneconsensus",  type="character", default=NULL, 
              help = "Compute consensus for gene programs and write to specified file"),
  make_option("--genecounts",  type="character", default=NULL, 
              help = "output table of gene counts to specified file"),
  make_option("--genebackground",  type="character", default=NULL, 
              help = "output table of gene background to specified file"),
  make_option("--wincounts",  type="character", default=NULL, 
              help = "output table of window counts to specified file"),
  make_option("--normcounts",  type="character", default=NULL, 
              help = "output table of normalized gene counts to specified file"),
  make_option("--plotsignatures",  type="character", default=NULL, 
              help = "output plot of signature enrichment to specified file"),
  make_option("--plotprograms",  type="character", default=NULL, 
              help = "output plot of program enrichment to specified file"),
  make_option("--plotenrichments",  type="character", default=NULL, 
              help = "output plot of program enrichment (HyperGeometric) to specified file"),
  make_option("--select",  type="character", default=NULL, 
              help = "Select subset of programs/signatures for plots"),
  make_option("--plotdetails", type = "logical", action = "store_true", default = FALSE,
              help = "Generete detailed plots for each programs/signatures"),
  make_option("--export",  type="character", default=NULL, 
              help = "export data to RDS file"),
  make_option("--ExportNormCounts", type="character", default=NULL, 
              help = "export normalized count data to RDS file"),
  make_option("--QC", type = "character", default = NULL,
              help = "prepare QC report to file"),
  make_option("--cluster", type = "character", default = NULL,
              help = "prepare file for Cluster program"),
  
  
  # auxilary files (instead of default)
  make_option("--windowsignatures",  type="character", default=NULL, 
              help = "Provide window siganture file (instead of Win-sig.csv)"),
  make_option("--geneprograms",  type="character", default=NULL, 
              help = "Provide gene gene programs file (instead of Gene-sig.csv)"),
  make_option("--metagenes",  type="character", default=NULL, 
              help = "Provide gene annotation file for meta plots (instead of meta-gene.bed)"),
  make_option("--metaenhancers",  type="character", default=NULL, 
              help = "Provide enhancer annotation file for meta plots (instead of meta-enhancer.bed)"),
  make_option("--config", type="character", default = NULL,
              help = "Specify an alternative configuration file")
);
# Initialize the option parser with the defined list of options.
opt_parser = OptionParser(option_list=option_list);

# Parse the command-line arguments. If 'trailing.options' exists, use it; otherwise, use the default positional arguments.
if( exists("trailing.options") ) {
  opt = parse_args(opt_parser, args = trailing.options, positional_arguments = TRUE)
} else {
  opt = parse_args(opt_parser, positional_arguments = TRUE)
}

# Extract the files to be processed from the parsed options.
Files = opt$args

# If not in development mode and no files are provided, print the help message and exit.
if(!DevelopmentMode && length(Files) == 0)
  parse_args(opt_parser,args = "-h", print_help_and_exit = TRUE)

# Print an initialization message to indicate the start of the script's execution.
catn("Initializing")

# Load required libraries quietly to avoid cluttering the console output.
suppressMessages(library(rtracklayer))
suppressMessages(library(tools))
suppressMessages(library(ggplot2))

# Set the target modification based on the command-line option provided.
TargetMod = opt$options$mod

# If a root directory is specified, set up the directories for data, BED files, and output based on the root directory.
if( !is.null(opt$options$root)) {
  RootDir =  extendDir(opt$options$root)
  DataDir = paste0(RootDir, "Samples/", TargetMod, "/")
  BedDir = paste0(RootDir, "BED/", TargetMod, "/")
} 

# If a project directory is specified, set up the directories for tracks and output based on the project directory.
if( !is.null(opt$options$project)) {
  ProjectDir =  extendDir(opt$options$project)
  TracksDIR = paste0(ProjectDir, "Tracks/", TargetMod, "/")
  OutputDir = paste0(ProjectDir, "Output/", TargetMod, "/")
} 

# Set the annotation directory based on the source directory and target modification.
ANNOTDIR = paste0(SourceDIR, "SetupFiles/", TargetMod, "/")

# Update the DataDir, BedDir, TracksDIR, and OutputDir based on the options provided.
if( !is.null(opt$options$datadir) )
  DataDir = extendDir(opt$options$datadir)
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
} else if( !exists("OutputDir") ) {
  OutputDir = DevelopmentMode ? paste0("~/Data/BloodChIP/Output/", Mod) : OrigDir
}

# Set the MetaDIR to the OutputDir.
MetaDIR = OutputDir

# If the SetupDIR does not exist, set it to the ANNOTDIR.
if( !exists("SetupDIR") )
  SetupDIR = ANNOTDIR
# catn("setupDir is: ", SetupDIR)


# Define the ReadSignatureList function which is responsible for reading signature or program files,
# processing them, and returning a structured list containing the signatures, their references, and names.
ReadSignatureList = function(Sig.filename, type = "signature") {
  # Split the Sig.filename string by commas, allowing for multiple files to be specified.
  Sig.filename = unlist(strsplit(Sig.filename, ","))
  
  # Iterate over each file name in the Sig.filename list.
  for( f in Sig.filename ) {
    # Check if the file exists. If not, print an error message and return NULL.
    if( !file.exists(f) ) {
      catn("Cannot find", type,"file", f)
      return(NULL)
    }
  }
  
  # Initialize an empty list to store the signatures.
  Sig.list = list()
  
  # Iterate over each file name again, this time to read and process the data.
  for( f in Sig.filename ) {
    # Print a message indicating the file being read.
    catn("Reading", type, "from", f)
    # Read the CSV file into a data frame.
    ZZ = read.csv(f, as.is = TRUE)
    # Split the second column of the data frame by the first column, creating a list of signatures.
    Sig = split(ZZ[,2], ZZ[,1])
    # If the type is "programs", filter the signatures to only include those not excluded.
    if( type == "programs")
      Sig = sapply(Sig, function(s) intersect(s, Genes.notexcluded))
    # Attempt to find a reference file by replacing the csv extension with rds.
    f.Ref = sub("csv$", "rds", f)
    # Initialize Sig.Ref as NULL.
    Sig.Ref = NULL
    # If the reference file exists, read it and check for consistency with the signatures.
    if( file.exists(f.Ref) ) {
      catn("Reading", type, "reference from", f.Ref)
      Sig.Ref = readRDS(f.Ref)
      if( length(Sig) != length(Sig.Ref$avg) || any(names(Sig) != names(Sig.Ref$avg)) ) {
        catn("Mismatch in consensus of signature", f)
        Sig.Ref = NULL
      }
    }
    # If Sig.Ref is still NULL, attempt to compute averages and variances based on available data.
    if( is.null(Sig.Ref) ) {
      if( type == "programs") {
        if( exists("Healthy.GeneCount") )
          Sig.Ref = list( avg = sapply(Sig, function(p) sum(Healthy.GeneCount[p])),
                          var = sapply(Sig, function(p) sum(Healthy.GeneCount.var[p])))
      } else
        if( exists("Healthy.WinCount") )
          Sig.Ref = list( avg = sapply(Sig, function(p) sum(Healthy.WinCount[p])),
                          var = sapply(Sig, function(p) sum(Healthy.WinCount.var[p])))
    }        
    
    # Add the processed signatures and their references to the Sig.list.
    if( !is.null(Sig.Ref)) {
      Sig.list[[f]] = list(Sig = Sig, avg = Sig.Ref$avg, var = Sig.Ref$var )
    } else 
      Sig.list[[f]] = list(Sig = Sig, avg = c(), var = c() )
  } 
  # Remove the names from Sig.list to simplify further processing.
  names(Sig.list) = NULL
  # Combine all signatures into a single list.
  Sig = do.call("c", lapply(Sig.list, function(x) x$Sig))
  # Combine all references into a single list.
  Sig.Ref = list( avg = do.call("c", lapply(Sig.list, function(x) x$avg)),
                  var = do.call("c", lapply(Sig.list, function(x) x$var))
  )
  # Return a list containing the combined signatures, their references, and a list of names.
  list(Sig = Sig, Ref = Sig.Ref, List = lapply(Sig.list, function(x) names(x$Sig)))
}


{
  # Read the TSS windows data from an RDS file located in the SetupDIR directory.
  TSS.windows = readRDS(paste0(SetupDIR,"Windows.rds"))
  # Retrieve sequence information from the TSS windows data.
  genome.seqinfo = seqinfo(TSS.windows)
  # Define a list of chromosome names typically found in human genome datasets.
  ChrList = paste0("chr", c(1:22,"X", "Y"))
  # If analyzing a specific modification in Saccharomyces cerevisiae, adjust the chromosome list accordingly.
  if (TargetMod == "H3K4me3-scer") {
    ChrList = (seqnames(genome.seqinfo))
  }
  
  # Attempt to load window signatures from a CSV file located in the SetupDIR directory.
  Win.Sig.filename = paste0(SetupDIR, "Win-sig.csv")
  Win.Sig = NULL
  Win.Sig.Ref = NULL
  # If an alternative window signatures file is specified, use that instead.
  if( !is.null(opt$options$windowsignatures))
    Win.Sig.filename = opt$options$windowsignatures
  # Read the window signatures file and store the data.
  ll = ReadSignatureList(Win.Sig.filename)
  if( !is.null(ll) ) {
    Win.Sig = ll$Sig
    Win.Sig.Ref = ll$Ref
  } 
  
  # Source additional R scripts containing functions required for the analysis.
  source(paste0(SourceDIR,"Background.R"))
  source(paste0(SourceDIR,"ComputeGeneCounts.R"))
  source(paste0(SourceDIR,"CommonGenes.R"))
  # Exclude certain genes from the analysis.
  Genes.notexcluded = Genes[!Genes.excluded]
  source(paste0(SourceDIR, "cfChIP-Functions.R"))  
  source(paste0(SourceDIR, "YieldEstimation-Functions.R"))
  
  # Attempt to load gene programs from a CSV file located in the SetupDIR directory.
  GenePrograms.filename =  paste0(SetupDIR, "Gene-sig.csv")
  Gene.Programs = NULL
  Gene.Programs.Ref = NULL
  Gene.Programs.Partition = NULL
  # If an alternative gene programs file is specified, use that instead.
  if( !is.null(opt$options$geneprograms))
    GenePrograms.filename = opt$options$geneprograms
  # Read the gene programs file and store the data.
  ll = ReadSignatureList(GenePrograms.filename,"programs")
  if(!is.null(ll)) {
    Gene.Programs = ll$Sig
    Gene.Programs.Ref = ll$Ref
    Gene.Programs.Partition = ll$List
  }
  
  # Attempt to load meta gene information from a BED file located in the SetupDIR directory.
  MetaGene.filename = paste0(SetupDIR, "Meta-genes.bed")
  MetaGene = NULL
  # If an alternative meta genes file is specified, use that instead.
  if( !is.null(opt$options$metagenes))  
    MetaGene.filename = opt$options$metagenes
  # Read the meta genes file and store the data.
  if( file.exists( MetaGene.filename)){
    catn("Reading meta gene info from", MetaGene.filename )
    MetaGene = import(MetaGene.filename)
  }
  
  # Attempt to load meta enhancer information from a BED file located in the SetupDIR directory.
  MetaEnhancer.filename = paste0(SetupDIR, "Meta-enhancers.bed")
  MetaEnhancer = NULL
  # If an alternative meta enhancers file is specified, use that instead.
  if( !is.null(opt$options$metaenhancer))  {
    MetaEnhancer.filename = opt$options$metaenhancer
  }
  # Read the meta enhancers file and store the data.
  if( file.exists( MetaEnhancer.filename)) {
    catn("Reading meta enhancer info from", MetaEnhancer.filename )
    MetaEnhancer = import(MetaEnhancer.filename)
  }
  
  # Attempt to load gene description information from a CSV file located in the SetupDIR directory.
  GeneDescription.filename = paste0(SetupDIR,"GeneDescription.csv")
  GeneDescription = NULL
  # Read the gene description file and store the data.
  if( file.exists(GeneDescription.filename) ) {
    GeneDescription = read.csv(GeneDescription.filename, as.is = TRUE)
    rownames(GeneDescription) = GeneDescription[,1]
    GeneDescription=GeneDescription[,-1]
    GeneDescription$Description = gsub(",",";",GeneDescription$Description)
  }
  # Prepare window description data for analysis.
  WinDescription = mcols(TSS.windows)
  colnames(WinDescription) = sapply(colnames(WinDescription), function(x) paste0(toupper(substr(x,1,1)), substr(x,2,width(x))))
  colnames(WinDescription) = sub("Name", "Gene", colnames(WinDescription))
  rownames(WinDescription) = 1:length(TSS.windows)
  
  # Attempt to load QC information from a BED file located in the SetupDIR directory.
  QC.filename = paste0(SetupDIR,"QC.bed")
  QC.bed = TSS.windows
  # Read the QC file and store the data.
  if( file.exists(QC.filename) ) {
    catn("Reading QC bed file",QC.filename)
    QC.bed = import(QC.filename)
    QC.bed$type = QC.bed$name
  }
}



# Initialize the cfChIP parameters by calling the cfChIP.Params function.
params = cfChIP.Params()
# Set the DataDir parameter within the params object to the previously defined DataDir variable.
params$DataDir = DataDir

# Initialize a variable to control whether data should be retained in memory after processing.
RetainInMemory = FALSE

# Check if a custom configuration file has been specified through command-line options.
if( !is.null(opt$options$config) ) {
  # If specified, use the provided configuration file.
  config.filename = opt$options$config
} else {
  # If not specified, use the default configuration file located in the SetupDIR directory.
  config.filename = paste0(SetupDIR, "config.csv")
}

# If the configuration file exists, read the configuration options from the file.
if( file.exists(config.filename) ) {
  catn("Reading configuration options from ",config.filename)
  # Read the configuration data from the CSV file.
  config.data = read.csv(config.filename,as.is = TRUE, header = FALSE)
  # Iterate through each row in the configuration data.
  for( i in 1:nrow(config.data) ) {
    # Update the params object with the configuration options.
    params[[config.data[i,1]]] = config.data[i,2]
  }
}

# Update the params object with options specified through command-line arguments.
params$Background = opt$options$background
params$GeneCounts = opt$options$count
params$Normalize = opt$options$normalize
params$OverExpressedGenes = opt$options$overexpressed
# Duplicate the overexpressed genes option for use within the script.
doOverExpressedGenes = opt$options$overexpressed

# Check if the extra option is specified, which triggers multiple processing options.
if( opt$options$extra ) {
  # Enable all relevant processing options.
  opt$options$all = TRUE
  opt$options$tracks = TRUE
  opt$options$signatures = TRUE
  opt$options$signaturesvshealthy  = TRUE
  opt$options$programs = TRUE
  opt$options$backgroundplot = TRUE
  opt$options$fraglenplot = TRUE
}

# If the all option is specified, enable comprehensive data processing.
if( opt$options$all ) {
  params$Background = TRUE
  params$GeneBackground = TRUE
  params$GeneCounts = TRUE
  params$Normalize = TRUE
  params$OverExpressedGenes = TRUE
}

# Determine whether to force re-computation based on command-line options.
doHardForce = opt$options$hardforce
doForce = opt$options$force || doHardForce
doWriteForce = opt$options$force || opt$options$writeforce

# Initialize variables to control which analyses to perform based on command-line options.
doMeta = opt$options$meta
doTracks = opt$options$tracks
doEnrichR = opt$options$enrichR
doSignatures = opt$options$signatures
doSignaturesVsHealthy = opt$options$signaturesvshealthy
doPrograms = opt$options$programs
doBackgroundPlot = opt$options$backgroundplot
doFragmentPlot = opt$options$fraglenplot
doctDNA = opt$options$ctDNA

# Check for required data before proceeding with specific analyses.
if( doPrograms &&  is.null(Gene.Programs) ) {
  stop("Cannot evaluate gene programs -- missing programs!")
}

if( (doSignatures || doSignaturesVsHealthy) &&  is.null(Win.Sig) ) {
  stop("Cannot evaluate signatures -- missing signatures!")
}

# If the script is running in DevelopmentMode, certain parameters can be set to TRUE for testing or development purposes.
# Here, we ensure that data is retained in memory for further inspection or debugging.
if( DevelopmentMode ) {
  RetainInMemory = TRUE
}

# If any of the following analyses (meta analysis, tracks creation, signature evaluation, or signature evaluation vs healthy samples) are to be performed,
# normalization of data is required and thus the Normalize parameter is set to TRUE.
if( doMeta || doTracks || doSignatures || doSignaturesVsHealthy )
  params$Normalize = TRUE

# If enrichR analysis or gene programs evaluation are to be performed, evaluating overexpressed genes is required.
# Thus, the OverExpressedGenes parameter is set to TRUE.
if( doEnrichR  || doPrograms )
  params$OverExpressedGenes = TRUE

# If tracks are to be created and the directory for tracks does not exist, it is created.
if( doTracks && !dir.exists(TracksDIR))
  dir.create(TracksDIR)

# If overexpressed genes analysis or meta analysis are to be performed and the output directory does not exist, it is created.
if( doOverExpressedGenes || doMeta )
  if( !dir.exists(OutputDir)) 
    dir.create(OutputDir)

# If meta analysis is to be performed and the MetaDIR does not exist, it is created.
if( doMeta && !dir.exists(MetaDIR))
  dir.create(MetaDIR)

# If consensus analysis is to be performed, the output file for consensus is set, and necessary parameters are updated.
doConsensus = !is.null(opt$options$consensus)
if(doConsensus ) {
  outputConsensus = outputFile(opt$options$consensus)
  params$GeneBackground = TRUE
  params$GeneCounts = TRUE
  RetainInMemory = TRUE
}

# If window signature consensus analysis is to be performed, the output file for window signature consensus is set.
# Data is retained in memory for further processing, and a check is performed to ensure window signatures data is available.
doWinSignatureConsensus =  !is.null(opt$options$winconsensus)
if( doWinSignatureConsensus ) {
  outputWinSigConsensus = outputFile(opt$options$winconsensus)
  RetainInMemory = TRUE
  if( is.null(Win.Sig) ) {
    stop("Cannot evaluate window signatures -- missing signatures!")
  }
}

# If gene signature consensus analysis is to be performed, the output file for gene signature consensus is set,
# and necessary parameters are updated. A check is performed to ensure gene programs data is available.
doGeneSignatureConsensus =  !is.null(opt$options$geneconsensus)
if( doGeneSignatureConsensus ) {
  outputGeneSigConsensus = outputFile(opt$options$geneconsensus)
  params$GeneBackground = TRUE
  params$GeneCounts = TRUE
  RetainInMemory = TRUE
  if( is.null(Gene.Programs) ) {
    stop("Cannot evaluate gene programs -- missing programs!")
  }
}

# If gene counts table creation is requested, the output file for gene counts table is set,
# and the GeneCounts parameter is updated. Data is retained in memory for further processing.
doGeneCountsTable = !is.null(opt$options$genecounts)
if(doGeneCountsTable ) {
  outputGeneCountsTable = outputFile(opt$options$genecounts)
  params$GeneCounts = TRUE
  RetainInMemory = TRUE
}

# If gene background table creation is requested, the output file for gene background table is set,
# and the GeneBackground parameter is updated. Data is retained in memory for further processing.
doGeneBackgroundTable = !is.null(opt$options$genebackground)
if(doGeneBackgroundTable ) {
  outputGeneBackgrouodTable = outputFile(opt$options$genebackground)
  params$GeneBackground = TRUE
  RetainInMemory = TRUE
}

# If window counts table creation is requested, the output file for window counts table is set,
# and the WinCounts parameter is updated. Data is retained in memory for further processing.
doWinCountsTable = !is.null(opt$options$wincounts)
if(doWinCountsTable ) {
  outputWinCountsTable = outputFile(opt$options$wincounts)
  params$WinCounts = TRUE
  RetainInMemory = TRUE
}

# If normalized counts table creation is requested, the output file for normalized counts table is set,
# and the Normalize parameter is updated. Data is retained in memory for further processing.
doNormCountsTable = !is.null(opt$options$normcounts)
if( doNormCountsTable ) {
  outputNormCountsTable = outputFile(opt$options$normcounts)
  params$Normalize = TRUE
  RetainInMemory = TRUE
}




# Check if the option to plot signatures is specified and set the necessary parameters.
doPlotSignatures = !is.null(opt$options$plotsignatures)
if(doPlotSignatures) {
  # If the Signatures parameter is not set, stop the script and display an error message.
  if(is.null(params$Signatures))
    stop("Error: missing signatures for --plotsignatures")
  # Set the output file path for the plot signatures output.
  outputPlotSignatures = outputFile(opt$options$plotsignatures)
  # Ensure that data normalization and overexpressed genes analysis are enabled.
  params$Normalize = TRUE
  params$OverExpressedGenes = TRUE
  # Keep the processed data in memory after the script finishes.
  RetainInMemory = TRUE
}

# Check if the option to plot gene programs is specified and set the necessary parameters.
doPlotPrograms = !is.null(opt$options$plotprograms)
if(doPlotPrograms) {
  # If the Programs parameter is not set, stop the script and display an error message.
  if(is.null(params$Programs))
    stop("Error: missing gene programs for --plotprograms")
  # Set the output file path for the plot programs output.
  outputPlotPrograms = outputFile(opt$options$plotprograms)
  # Ensure that data normalization and overexpressed genes analysis are enabled.
  params$Normalize = TRUE
  params$OverExpressedGenes = TRUE
  # Keep the processed data in memory after the script finishes.
  RetainInMemory = TRUE
}

# Check if the option to plot enrichments is specified and set the necessary parameters.
doPlotEnrichments = !is.null(opt$options$plotenrichments)
if(doPlotEnrichments) {
  # If the Programs parameter is not set, display an error message but do not stop the script.
  if(is.null(params$Programs))
    catn("Error: missing gene programs for --plotenrichments")
  # Set the output file path for the plot enrichments output.
  outputPlotEnrichments = outputFile(opt$options$plotenrichments)
  # Ensure that data normalization and overexpressed genes analysis are enabled.
  params$Normalize = TRUE
  params$OverExpressedGenes = TRUE
  # Keep the processed data in memory after the script finishes.
  RetainInMemory = TRUE
}

# Check if detailed plots for each program/signature are requested.
if(opt$options$plotdetails) {
  # If neither plot programs nor plot signatures options are enabled, display a warning message.
  if(!(doPlotPrograms || doPlotSignatures))
    catn("Warning: --plotdetails requires either --plotprograms or --plotsignatures")
  # Enable individual heatmap, bar chart, and CSV output for signatures and programs.
  params$PlotSignatures.IndividualHeatmap = TRUE
  params$PlotPrograms.IndividualHeatmap = TRUE
  params$PlotPrograms.IndividualBarChart = TRUE
  params$PlotPrograms.IndividualCSV = TRUE
}

# Check if a selection of specific programs/signatures is provided for plotting.
if(!is.null(opt$options$select)) {
  # If neither plot programs nor plot signatures options are enabled, display a warning message.
  if(!(doPlotPrograms || doPlotSignatures))
    catn("Warning: --select requires either --plotprograms or --plotsignatures")
  # Read the list of selected programs/signatures from the specified file.
  catn("Reading list of selected programs/signatures from", opt$options$select)
  SelectNames = read.csv(opt$options$select, as.is = TRUE)
  SelectNames = unlist(SelectNames)
  # If plotting programs or enrichments, filter the programs based on the selection.
  if(doPlotPrograms || doPlotEnrichments) {
    XX = intersect(names(params$Programs), SelectNames)
    params$Programs = params$Programs[XX]
    params$Programs.Ref$avg = params$Programs.Ref$avg[XX]
    params$Programs.Ref$var = params$Programs.Ref$var[XX]
    params$Programs.Partition = lapply(params$Programs.Partition, function(x) intersect(x, SelectNames))
  }
  # If plotting signatures, filter the signatures based on the selection.
  if(doPlotSignatures) {
    XX = intersect(names(params$Signatures), SelectNames)
    params$Signatures = params$Signatures[XX]
    params$Signatures.Ref$avg = params$Signatures.Ref$avg[XX]
    params$Signatures.Ref$var = params$Signatures.Ref$var[XX]
  }
}

# Check if data export to RDS file is requested and set the necessary parameters.
doExport = !is.null(opt$options$export)
if(doExport) {
  outputExport = outputFile(opt$options$export)
  RetainInMemory = TRUE
}

# Check if export of normalized count data to RDS file is requested and set the necessary parameters.

doExportNormCounts = !is.null(opt$options$ExportNormCounts)
if(doExportNormCounts ) {
  outputExport = outputFile(opt$options$ExportNormCounts)
  RetainInMemory = TRUE
  params$Normalize = TRUE
}

doQC = !is.null(opt$options$QC)
if( doQC ) {
  outputQC = outputFile(opt$options$QC)
  QC = list()
  params$Background = TRUE
}

doCluster = !is.null(opt$options$cluster)
if( doCluster ) {
  outputCluster = outputFile(opt$options$cluster)
  Cluster = list()
  params$Background = TRUE
  params$Normalize = TRUE
}

if( doMeta && is.null(params$MetaGene) )
  catn("Missing meta gene information")


# Define a function to process BED files for cfChIP analysis.
ProcessBEDFile = function(BFile) {
  # If the filename starts with "XXX", it's likely a placeholder or invalid, so return NULL to skip processing.
  if(grepl("^XXX", BFile))
    return(NULL)
  
  # Check if the file exists in various formats by appending common extensions. If not found, prepend the BedDir path.
  if(!file.exists(BFile) &&
     !file.exists(paste0(BFile,".bed")) &&
     !file.exists(paste0(BFile,".tagAlign")) && 
     !file.exists(paste0(BFile,".bw")) &&
     !file.exists(paste0(BFile,".bed.gz")) &&
     !file.exists(paste0(BFile,".tagAlign.gz")))
    BFile = paste0(BedDir, BFile)
  
  # Print the final file path/name being processed.
  catn(BFile)
  
  # Process the BED file using the cfChIP.ProcessFile function with parameters for force and hard force re-computation.
  dat = cfChIP.ProcessFile(filename = BFile, param = params, Force = doForce, HardForce = doHardForce)
  
  # If quality control (QC) analysis is enabled, perform QC using the cfChIP.countQC function and store the results.
  if(doQC)
    QC[[dat$Name]] <<- cfChIP.countQC(dat, BFile, GR = QC.bed, param= params)
  
  # If clustering analysis is enabled, store the normalized gene counts for clustering.
  if(doCluster)
    Cluster[[dat$Name]] <<- dat$GeneCounts.QQnorm
  
  # If background plotting is enabled, generate and save the background plot.
  if(doBackgroundPlot) {
    if(params$Verbose) catn(dat$Name, ": Plotting background estimate")
    cfChIP.BackgroundPlot(dat, Dir = OutputDir, Force = doWriteForce)
  }
  
  # If fragment length plotting is enabled, generate and save the fragment length distribution plot.
  if(doFragmentPlot) {
    if(params$Verbose) catn(dat$Name, ": Plotting fragment lengths")
    cfChIP.FragmentLenPlot(dat, Dir = OutputDir, Force = doWriteForce) 
  }
  
  # If circulating tumor DNA (ctDNA) estimation is enabled, estimate and report the ctDNA contribution.
  if(doctDNA) {
    if(params$Verbose) catn(dat$Name, ": Estimating ctDNA contribution")
    cfChIP.EstimatectDNA(dat, Dir = OutputDir, Force = doWriteForce) 
  }
  
  # If analysis of overexpressed genes is enabled, identify and report overexpressed genes and windows.
  if(doOverExpressedGenes) {
    if(params$Verbose) catn(dat$Name, ": Printing overexpressed genes")
    cfChIP.OverExpressedGenes(dat, Dir = OutputDir, param = params, Force = doWriteForce)
    if(params$Verbose) catn(dat$Name, ": Printing overexpressed windows")
    cfChIP.OverExpressedWins(dat, Dir = OutputDir, param = params, Force = doWriteForce)
  }
  
  # If meta analysis is enabled and meta gene information is available, generate and save the meta coverage plot.
  if(doMeta && !is.null(params$MetaGene)) {
    if(params$Verbose) catn(dat$Name, ": Plotting meta coverage")
    cfChIP.MetaPlot(dat, BFile, Dir = MetaDIR, param = params, Force = doWriteForce)
  }
  
  # If track creation is enabled, generate and save normalized genome browser tracks.
  if(doTracks) {
    if(params$Verbose) catn(dat$Name, ": Saving normalized tracks")
    cfChIP.WriteTrack(dat, BFile, params, Dir = TracksDIR, Force = doWriteForce)
  }
  
  # If signature evaluation is enabled, perform the analysis and save the results.
  if(doSignatures) {
    if(params$Verbose) catn(dat$Name, ": Evaluating signatures")
    cfChIP.EvaluateSignatures(dat, Dir = OutputDir, param = params, Write = TRUE, Force = doWriteForce)
  }
  # If signature evaluation against healthy samples is enabled, perform the analysis and save the results.
  if(doSignaturesVsHealthy) {
    if(params$Verbose) catn(dat$Name, ": Evaluating signatures vs Healthy")
    cfChIP.EvaluateSignatures(dat, Dir = OutputDir, param = params, Write = TRUE, Force = doWriteForce, WithReference = TRUE)
  }
  
  # If gene expression program evaluation is enabled, perform the analysis and save the results.
  if(doPrograms) {
    if(params$Verbose) catn(dat$Name, ": Evaluating gene expression programs")
    cfChIP.EvaluatePrograms(dat, Dir = OutputDir, param = params, Write = TRUE, WithReference = TRUE, Force = doWriteForce)
    cfChIP.EvaluateProgramsHypG(dat, Dir = OutputDir, param = params, Write = TRUE, Force = doWriteForce)
  }
  # If gene expression program evaluation without healthy samples is enabled, perform the analysis and save the results.
  if(doSignaturesVsHealthy) {
    if(params$Verbose) catn(dat$Name, ": Evaluating gene expression programs w/o healthy")
    cfChIP.EvaluatePrograms(dat, Dir = OutputDir, param = params, Write = TRUE, Force = doWriteForce, WithReference = FALSE)
    cfChIP.EvaluateProgramsHypG(dat, Dir = OutputDir, param = params, Write = TRUE, Force = doWriteForce)
  }
  
  # If EnrichR analysis is enabled, perform the analysis and save the results.
  if(doEnrichR) {
    if(params$Verbose) catn(dat$Name, ": Checking EnrichR")
    cfChIP.WriteEnrichR(dat, Dir = OutputDir, Force = doWriteForce)
  }
  
  
  # Check if data should be retained in memory after processing.
  if(RetainInMemory) {
    # If not in development mode, clear the coverage, BED, and BW data from the 'dat' object to save memory.
    if(!DevelopmentMode) {
      dat$Cov = NULL
      dat$BED = NULL
      dat$BW = NULL
    }
    # If exporting normalized counts, clear the Counts and GeneCounts from the 'dat' object to save memory.
    if(doExportNormCounts) {
      dat$Counts = NULL
      dat$GeneCounts = NULL
    }
    # Return the processed data object.
    return(dat) 
  } else {
    # If data should not be retained, remove the 'dat' object from memory.
    rm(dat)
  }
}
  # Define a function to process a list of BED files.
  ProcessBEDFileList = function(BFlist) {
    # Apply the ProcessBEDFile function to each file in the list and store the results.
    L = lapply(BFlist, ProcessBEDFile)
    # Set the names of the list elements to the 'Name' attribute of the processed data objects.
    names(L) = sapply(L, function(l) l$Name)
    # Return the list of processed data objects.
    L  
  }
  
  # Define a function to expand file names, which can include reading from a text or CSV file containing a list of files.
  expandFiles = function(f) {
    # Check if the file name indicates a text or CSV file.
    if(grepl(".txt$", f) | grepl(".csv$", f)) {
      # If the file does not exist, print a message and return NULL.
      if(!file.exists(f)) {
        catn("Missing list file", f)
        return(NULL)
      }
      # Read the list of files from the text or CSV file.
      catn("Reading sample list from", f)
      File.list = read.table(f, as.is = TRUE)
      File.list = File.list[,1]
      names(File.list) = NULL
      # Return the list of files as a list object.
      as.list(File.list)
    } else {
      # If the file name does not indicate a list file, return the file name as is.
      f
    }
  }
  
  # Create a unique list of full file paths by expanding any list files and combining the results.
  FullFiles = unique(unlist(lapply(Files, expandFiles)))
  
  # Process the list of full file paths using the ProcessBEDFileList function.
  LL = ProcessBEDFileList(FullFiles)
  # Filter out any NULL results from the list.
  LL = LL[!sapply(LL, is.null)]
  
  # If quality control (QC) analysis is enabled, write the QC results to the specified output file.
  if(doQC) {
    catn("Writing QC numbers to ", outputQC)
    # Combine the QC data into a data frame.
    df = data.frame(do.call(rbind, QC))
    # Write the QC data to a CSV file.
    write.csv(df, file = outputQC, quote = FALSE)
  }
  
  # If clustering analysis is enabled, save the clustering data to the specified output file.
  if(doCluster) {
    catn("Saving file for clustering", outputCluster)
    # Combine the clustering data into a matrix.
    X = do.call(cbind, Cluster)
    colnames(X) = names(Cluster)
    rownames(X) = Genes
    # Create a matrix with an additional header row and column for use with clustering software.
    M = matrix(nr = nrow(X)+1, nc = ncol(X)+2)
    M[1,] = c("UID", "Name", colnames(X))
    M[2:(nrow(X)+1),1] = rownames(X)
    M[2:(nrow(X)+1),2] = rownames(X)
    M[2:(nrow(X)+1),3:(ncol(X)+2)] = formatC(log2(1+X))
    M[is.na(M)] = ""
    # Write the clustering matrix to a tab-separated file.
    write.table(M, file = outputCluster, sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)
  }
  
  # Check if data should be retained in memory after processing.
  if(RetainInMemory) {
    # Assign sample names to the processed list for easy reference.
    SampleNames = sapply(LL, function(l) l$Name)
    names(LL) = SampleNames
    SampleOrder = names(LL)  # Store the order of samples for potential future use.
    
    # If exporting the processed data is requested.
    if(doExport) {
      catn("Exporting data in RDS format")
      dat = LL[[1]]  # Take the first element as a template for the structure.
      
      # Clear out large data structures to save memory before exporting.
      dat$Cov = NULL
      if(!is.null(dat$Counts))
        dat$Counts = sapply(LL, function(l) l$Counts)  # Aggregate counts data from all samples.
      if(!is.null(dat$Background))
        dat$Background = lapply(LL, function(l) l$Background)  # Aggregate background data from all samples.
      if(!is.null(dat$GeneCounts))
        dat$GeneCounts = sapply(LL, function(l) l$GeneCounts)  # Aggregate gene counts data from all samples.
      
      catn("QQnorm")
      if(!is.null(dat$QQNorm))
        dat$QQNorm = sapply(LL, function(l) l$QQNorm)  # Aggregate normalization factors from all samples.
      if(!is.null(dat$OverExpressedGenes))
        dat$OverExpressedGenes = lapply(LL, function(l) l$OverExpressedGenes)  # Aggregate overexpressed genes data from all samples.
      if(!is.null(dat$GeneBackground))
        dat$GeneBackground = sapply(LL, function(l) l$GeneBackground)  # Aggregate gene background data from all samples.
      
      catn("GenesCounts.QQnorm")
      if(!is.null(dat$GenesCounts.QQnorm)) {
        dat$GeneCounts.QQnorm = sapply(LL, function(l) l$GeneCounts.QQnorm[,1])  # Aggregate normalized gene counts data from all samples.
      }
      if(!is.null(dat$Counts.QQnorm))
        dat$Counts.QQnorm = sapply(LL, function(l) l$Counts.QQnorm)  # Aggregate normalized counts data from all samples.
      
      saveRDS(dat, outputExport)  # Save the aggregated data to an RDS file.
    }
    
    # If exporting normalized counts is requested.
    if(doExportNormCounts) {
      catn("Exporting norm counts in RDS format")
      dat = LL[[1]]  # Use the first element as a template for the structure.
      
      # Aggregate and save normalized gene counts data.
      catn("GenesCounts.QQnorm")
      if(!is.null(dat$GenesCounts.QQnorm)) {
        dat$GeneCounts.QQnorm = sapply(LL, function(l) l$GeneCounts.QQnorm[,1])
      }
      catn("Counts.QQnorm")
      if(!is.null(dat$Counts.QQnorm))
        dat$Counts.QQnorm = sapply(LL, function(l) l$Counts.QQnorm)
      
      saveRDS(dat, outputExport)  # Save the aggregated normalized counts data to an RDS file.
    }
  }
    # Additional logic for exporting gene counts table, gene background table, window counts table, and normalized gene counts table.
    # These blocks follow a similar pattern: aggregate the specific data type from all samples, then write it to a CSV file.
    
    # The `doGeneCountsTable`, `doGeneBackgroundTable`, `doWinCountsTable`, and `doNormCountsTable` flags control whether these exports are performed.
    # For each, the script aggregates the relevant data across all samples, formats it as needed, and writes it to the specified output file.
    
    # This approach allows for a comprehensive export of processed data, ensuring that all relevant information is captured and saved for further analysis or reporting.
    
    # Check if consensus analysis is requested.
    if(doConsensus) {
      # Print a message if verbose mode is enabled.
      if(params$Verbose) catn("Writing normalized consensus to ", outputConsensus)
      
      # Aggregate window counts and background data from all samples.
      WinCounts = do.call("cbind", lapply(LL, function(l) l$Counts))
      WinBackground = do.call("cbind", lapply(LL, function(l) l$WinBackground))
      # Aggregate normalization factors from all samples.
      QQnorm = sapply(LL, function(l) l$QQNorm)
      names(QQnorm) = colnames(WinCounts)
      
      # Estimate mean and variance for window data.
      win.est = cfChIP.EstimateMeanVarianceBasis(WinCounts, WinBackground, QQnorm)
      
      # Aggregate gene counts and background data from all samples.
      GeneCounts = do.call("cbind", lapply(LL, function(l) l$GeneCounts))
      GeneBackground = do.call("cbind", lapply(LL, function(l) l$GeneBackground))
      
      # Estimate mean and variance for gene data.
      gene.est = cfChIP.EstimateMeanVarianceBasis(GeneCounts, GeneBackground, QQnorm)
      
      # Compile consensus data into a list and save it as an RDS file.
      consensus = list(Win.avg = win.est$avg, Win.var = win.est$var, Gene.avg = gene.est$avg, Gene.var = gene.est$var)
      saveRDS(consensus, outputConsensus)
    }
    
    # Check if window signature consensus analysis is requested.
    if(doWinSignatureConsensus) {
      # Print a message if verbose mode is enabled.
      if(params$Verbose) catn("Writing normalized consensus of window signature ", outputWinSigConsensus)
      
      # Estimate mean and variance for window signatures.
      win.est = cfChIP.EstimateMeanVarianceWinSig(LL, Win.Sig)
      # Save the estimated data as an RDS file.
      saveRDS(win.est, outputWinSigConsensus)
    }
    
    # Check if gene signature consensus analysis is requested.
    if(doGeneSignatureConsensus) {
      # Print a message if verbose mode is enabled.
      if(params$Verbose) catn("Writing normalized consensus of gene signature ", outputGeneSigConsensus)
      
      # Estimate mean and variance for gene signatures.
      gene.est = cfChIP.EstimateMeanVarianceGeneSig(LL, Gene.Programs)
      # Save the estimated data as an RDS file.
      saveRDS(gene.est, outputGeneSigConsensus)
    }
    
    # Check if plotting of signatures is requested.
    if(doPlotSignatures) {
      # Plot signatures using the provided data and parameters.
      cfChIP.plotSignatures(LL, outputPlotSignatures, params)
    }
    
    # Check if plotting of gene programs is requested.
    if(doPlotPrograms) {
      # Plot gene programs using the provided data and parameters.
      cfChIP.plotPrograms(LL, outputPlotPrograms, params)
    }
    
    # Check if plotting of enrichments is requested.
    if(doPlotEnrichments) {
      # Plot enrichments using the provided data and parameters.
      cfChIP.plotEnrichments(LL, outputPlotEnrichments, params)
  }
    
