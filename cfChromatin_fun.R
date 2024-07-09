# Creation of a feature catalog based on a ChromHMM annotation
#	@annotation_dir 	: directory of the ChromHMM files
#	@glossary 			: glossary metadata file for the ChromHMM files
#	@glossary_group 	: group to group the EID sample by
#	@states				: vector of states to keep for the analysis
#	@verbose			: Display verbose
load_catalog <- function(
    annotation_dir		= NULL,
    glossary			= NULL,
    glossary_group		= NULL,
    states				= NULL,
    verbose				= TRUE
){
  # Required library
  library(parallel)
  
  # Get all files from annotation directory
  if (verbose){print("Load ChromHMM files...")}
  files <- list.files(annotation_dir, pattern = "\\.bed", full.names = T)
  
  # Load each file in granges and name them
  encode <- mclapply(files, read_bed, mc.cores = detectCores())
  encode <- lapply(encode, function(gr) {
    gr[seqnames(gr) != "chrM"]})
  names(encode) <- sapply(files, function(x) {strsplit(basename(x), "_")[[1]][1]})
  
  # Load tissue metadata file
  glossary <- as.data.frame(read_delim(glossary_path, delim = "\t", escape_double = FALSE, trim_ws = TRUE))
  # Remove cell line files
  glossary <- glossary[glossary$GROUP != "ENCODE2012",]
  
  # Create a grange for each tissue and reduce it
  if (verbose){print("Reducing ChromHMM windows per group...")}
  windows_list <- mclapply(unique(glossary[,glossary_group]), function(tissue) {
    eid_tissue <- glossary[glossary[,glossary_group] == tissue,]$EID
    tmp <- bind_ranges(encode[eid_tissue])
    tmp <- tmp[tmp$name %in% states, ]
    tmp <- reduce(tmp, ignore.strand=TRUE)
    tmp$tissue <- tissue
    return(tmp)
  }, mc.cores = detectCores())
  
  
  # Aggregate the granges
  if (verbose){print("Aggregate ChromHMM windows...")}
  windows_aggr <- sort(sortSeqlevels(bind_ranges(windows_list)), ignore.strand=TRUE)
  # Split windows_aggr into chunks
  chunks <- split(windows_aggr, ceiling(seq_along(windows_aggr) / length(windows_aggr) * detectCores()))
  
  # Apply disjoin to each chunk in parallel
  result_list <- mclapply(chunks, function(chunk) {
    windows <- disjoin(chunk, with.revmap=TRUE, ignore.strand=TRUE)
    windows$tissue <- unlist(lapply(windows$revmap,function(i){paste(collapse=';',chunk$tissue[i])}))
    windows$revmap <- NULL
    return(windows)
  }, mc.cores = detectCores())
  
  # Combine the results
  windows_list <- do.call(c, result_list)
  windows <- bind_ranges(windows_list)
}


# Annotation of the catalog based on databases
#	@windows			: Windows to annotate
#	@feature_tag		: Tag name of the features of interest
#	@db					: List of database annotation to retrieve information from, the first takes priority on the second etc...
#	@db_gap				: Tolerated gap to annotate the windows
#	@db_spe				: Boolean to include or not annotation specific features
#	@db_tag				: Tag name of the annotation specific features
#	@db_spe_center		: Size of the windows centered on the feature
#	@verbose			: Display verbose
annotate_catalog <- function(
    windows				= NULL,
    feature_tag			= NULL,
    db					= NULL,
    db_gap				= NULL,
    db_spe				= NULL,
    db_tag				= NULL,
    db_spe_center		= NULL,
    verbose				= TRUE
){
  #Annotation of the windows
  if (verbose){print("Windows annotation from the databases...")}
  
  db_list <- list()
  for (db_name in names(db)){
    # Gather features annotation from UCSC
    if(toupper(db_name)=="UCSC"){
      ucscdb <- promoters(db[[db_name]], upstream=0, downstream=1, columns=c("tx_name", "gene_id"))
      ucscdb <- sort(sortSeqlevels(ucscdb), ignore.strand=TRUE)
      mcols(ucscdb)$gene_id <- as.character(mcols(ucscdb)$gene_id)
      entrezid2symbol <- select(org.Hs.eg.db, mcols(ucscdb)$gene_id, c("ENTREZID", "SYMBOL"))
      stopifnot(identical(entrezid2symbol$ENTREZID, mcols(ucscdb)$gene_id))
      mcols(ucscdb)$SYMBOL <- entrezid2symbol$SYMBOL
      ucscdb <- ucscdb[!is.na(ucscdb$SYMBOL),]
      ucscdb$gene_id <- NULL
      db_list[[db_name]] <- ucscdb
    }
    # Gather features annotation from Ensembl
    else if(toupper(db_name)=="ENSEMBL"){
      ensdb <- promoters(db[[db_name]], upstream=0, downstream=1, columns="tx_id")
      ensdb <- sort(sortSeqlevels(ensdb), ignore.strand=TRUE)
      mcols(ensdb)$tx_id <- as.character(mcols(ensdb)$tx_id)
      ensdb <- ensdb[!is.na(ensdb$tx_id),]
      chr_name <- str_replace(string=paste("chr",seqlevels(ensdb),sep=""), pattern="chrMT", replacement="chrM")
      seqlevels(ensdb) <- chr_name
      db_list[[db_name]] <- ensdb
    }
  }
  
  # First database takes priority on the second one...
  windows$name <- "UNKNOWN"
  for (db_name in names(db)){
    
    if (verbose){print(paste0("Annotation from ", db_name,"..."))}
    
    windows <- join_overlap_left(windows, db_list[[db_name]], maxgap = db_gap)
    
    #Annotation UCSC
    if(toupper(db_name)=="UCSC"){
      
      windows[windows$name == "UNKNOWN",]$name <- windows[windows$name == "UNKNOWN",]$SYMBOL
      windows[is.na(windows$name)]$name <- "UNKNOWN"
      windows$tx_name <- NULL
      windows$SYMBOL <- NULL
    }
    #Annotation Ensembl
    else if(toupper(db_name)=="ENSEMBL"){
      windows[windows$name == "UNKNOWN",]$name <- windows[windows$name == "UNKNOWN",]$tx_id
      windows[is.na(windows$name)]$name <- "UNKNOWN"
      windows$tx_id <- NULL
    }
  }
  
  # Group the ranges for the same gene names
  if (verbose){print("Group ranges...")}
  windows <- windows %>% 
    group_by(seqnames, start, end, tissue) %>%
    summarise(name = paste(unique(name), collapse = ";")) %>%
    makeGRangesFromDataFrame(keep.extra.columns = T)
  
  windows$type <- feature_tag
  
  if (db_spe){
    
    #############################################################################################################################
    ###                                       3.GET REGIONS IN DATABASES NOT IN WINDOWS                                       ###
    #############################################################################################################################
    
    if (verbose){print("Include annotation missing from catalogue")}
    
    for (db_name in names(db)){
      
      if (verbose){print(paste0("Extra annotation from ", db_name,"..."))}
      
      if (verbose){print(paste0("Find new regions..."))}
      
      # Reduce our windows to check which regions are already covered
      windows_covered <- reduce(windows)
      
      tmp_db <- db_list[[db_name]]
      
      # Extended positions covered by the database
      tmp_db_extend <- resize(tmp_db, width = width(tmp_db)+db_spe_center, fix = "center")
      tmp_db_extend <- tmp_db_extend[seqnames(tmp_db_extend) %in% seqlevels(windows)]
      seqlevels(tmp_db_extend) <- seqlevels(windows)
      tmp_db_extend <- resize_correct(tmp_db_extend)
      
      tmp_db_extend_red <- reduce(tmp_db_extend)
      tmp_db_extend_red <- sort(sortSeqlevels(tmp_db_extend_red), ignore.strand=TRUE)
      
      # Genome are not exactly the same but positions are identicial except some on chrM
      genome(tmp_db_extend_red) <- NA
      
      #Find databases specific location
      db_specific <- disjoin_specific(windows_covered,tmp_db_extend_red,db_tag)
      
      if (verbose){print(paste0("Annotate new regions..."))}
      
      # Reannotate those extra ranges from database
      db_overlaps <- findOverlaps(db_specific, tmp_db_extend)
      db_hits <- queryHits(db_overlaps)
      names(db_hits) <- subjectHits(db_overlaps)
      
      #Annotation UCSC
      if(toupper(db_name)=="UCSC"){
        db_specific$name <- unlist(lapply(seq(1:length(db_specific)),function(i){
          ids_hit <- as.integer(names(db_hits[db_hits==i]))
          paste(collapse=';',unique(tmp_db_extend$SYMBOL[ids_hit]))
        }))
      }
      #Annotation Ensembl
      else if(toupper(db_name)=="ENSEMBL"){
        db_specific$name <- unlist(lapply(seq(1:length(db_specific)),function(i){
          ids_hit <- as.integer(names(db_hits[db_hits==i]))
          paste(collapse=';',unique(tmp_db_extend$tx_id[ids_hit]))
        }))
      }
      
      # Merge the windows and new ranges
      windows <- sort(sortSeqlevels(bind_ranges(list(windows,db_specific))), ignore.strand=TRUE)
    }
  }
  
  return(windows)
}

# Add flanking regions to the catalog
#	@windows			: Windows to annotate
#	@name_tag		: Tag name of the features of interest
#	@flanking_size		: Size of the flanking regions around the features
#	@verbose			: Display verbose
add_flank <- function(
    windows				= NULL,
    feature_tag			= NULL,
    flanking_size		= NULL,
    verbose				= TRUE
){
  #############################################################################################################################
  ###                                                  4.SET UP FLANKING REGIONS                                            ###
  #############################################################################################################################
  
  if (verbose){print(paste0("Add flanking regions to features..."))}
  
  # Reduce our windows to check which regions are already covered
  windows_red <- sort(sortSeqlevels(reduce(windows)), ignore.strand=TRUE)
  # Extending both side of all features
  windows_red_flank <- resize(windows_red, width = width(windows_red)+(flanking_size*2), fix = "center")
  windows_red_flank <- resize_correct(windows_red_flank)
  #Find flank specific location
  flank_specific <- disjoin_specific(windows_red, windows_red_flank,feature_tag)
  
  windows <- sort(sortSeqlevels(bind_ranges(list(windows,flank_specific))), ignore.strand=TRUE)
  
  return(windows)
}

# Add background to the catalog
#	@windows			: Windows to annotate
#	@feature_tag		: Tag name of the features of interest
#	@bin_size			: Size of the tile of the background
#	@verbose			: Display verbose
add_background <- function(
    windows				= NULL,
    feature_tag			= NULL,
    bin_size			= NULL,
    verbose				= TRUE
){
  #############################################################################################################################
  ###                                                 4.SET UP BACKGROUND REGIONS                                           ###
  #############################################################################################################################
  
  if (verbose){print(paste0("Add background to the catalog..."))}
  
  # Tile the whole genome
  bg_genome <- tileGenome(seqlengths(windows), tilewidth=bin_size, cut.last.tile.in.chrom=TRUE)
  
  # Reduce our windows to check which regions are already covered
  windows_red <- sort(sortSeqlevels(reduce(windows)), ignore.strand=TRUE)
  
  #Find background specific location
  bg_specific <- disjoin_specific(windows_red,bg_genome,feature_tag)
  
  windows <- sort(sortSeqlevels(bind_ranges(list(windows,bg_specific))), ignore.strand=TRUE)
  
  return(windows)
}

# Statistics about the window
#	@windows			: Windows to annotate
#	@db_tag				: Tag name of the annotation specific features
show_statistics <- function(
    windows				= NULL,
    db_tag				= NULL
){
  print(table(windows$type))
  
  tmp <- windows[windows$type == db_tag,]
  nb_extra_ucsc <- length(tmp[substring(tmp$feature,1,4) != "ENST",])
  nb_extra_ensembl <- length(tmp[substring(tmp$feature,1,4) == "ENST",])
  print(paste0("UCSC : ",nb_extra_ucsc))
  print(paste0("Ensembl : ",nb_extra_ensembl))
  
  size_genome <- sum(seqlengths(windows))
  print("Genome coverage")
  for (type_reg in unique(windows$type)){
    tmp <- windows[windows$type==type_reg,]
    print(paste0(type_reg, " : ", round(sum(width(tmp))*100/size_genome, 2)))
  }
  
  print("Quantile values")
  for (type_reg in unique(windows$type)){
    tmp <- windows[windows$type==type_reg,]
    print(paste0(type_reg, " : ", length(tmp)))
    print(paste0("    Min : ", min(width(tmp))))
    print(paste0("    Quantile 25% : ", quantile(width(tmp), 0.25)))
    print(paste0("    Quantile 50% : ", quantile(width(tmp), 0.50)))
    print(paste0("    Quantile 75% : ", quantile(width(tmp), 0.75)))
    print(paste0("    Max : ", max(width(tmp))))
  }
}

# Sub function the correct ranges resizing
#	@gr					: Grange to correct
resize_correct <- function(
    gr					= NULL
){
  # Resize can go negative
  start(gr) <- ifelse(start(gr)<0, 0, start(gr))
  # Resize can go over chromosome length
  for (chr in names(seqlengths(gr))){
    end(gr[seqnames(gr)==chr,]) <- ifelse(end(gr[seqnames(gr)==chr,]) > seqlengths(gr)[chr], seqlengths(gr)[chr], end(gr[seqnames(gr)==chr,]))
  }
  return(gr)
}

# Sub function to extract ranges specific tothe second grange entry
#	@gr1				: Grange1 
#	@gr1				: Grange2
#	@tag				: Tag to use to reannotate the type
disjoin_specific <- function(
    gr1					= NULL,
    gr2					= NULL,
    tag					= NULL
){
  
  gr1$region <- "GR1"
  gr2$region <- "GR2"
  
  # Substract features from gr1 to gr2
  aggr <- sort(sortSeqlevels(bind_ranges(list(gr1,gr2))), ignore.strand=TRUE)
  disjoin <- disjoin(aggr, with.revmap=TRUE, ignore.strand=TRUE)
  disjoin <- disjoin[lengths(disjoin$revmap)==1,]
  disjoin$revmap <- unlist(disjoin$revmap)
  gr2_specific <- disjoin[disjoin$revmap %in% which(aggr$region == "GR2")]
  
  #Reannotate type with the tag
  gr2_specific$revmap <- NULL
  gr2_specific$type <- tag
  
  return(gr2_specific)
}

#############################################################################################################################
###                                                 5.GET COVERAGE READ from BED FILE                                     ###
#############################################################################################################################
Max.Matrix.Mult = 20000

ComputeGeneCounts = function (A, W2G = Win2Gene.matrix) {
  if(!is.list(A)) 
    A = list(A)
  breaks = c(seq(0,dim(W2G)[1], Max.Matrix.Mult), dim(W2G)[1])
  L = lapply(1:(length(breaks)-1), function(i) {
    Is = (breaks[i]+1):breaks[i+1]
    M = W2G[Is,] 
    do.call("cbind", lapply(A, function(v) t(v[Is] %*% M)))
  })
  G = Reduce("+", L)
  matrix(G, nr = dim(G)[1], nc = dim(G)[2], dimnames = dimnames(G))
}

MaxGeneCounts = function(A,W2G = Win2Gene.matrix) {
  W2G.triplet = as(W2G, "TsparseMatrix")
  W2G.lists = split(W2G.triplet@i+1, W2G.triplet@j+1)
  X = sapply(W2G.lists, function(l) max(A[l]))
  names(X) = colnames(W2G)
  X
}

if( !exists("Win2Gene.matrix") ) {
  Win2Gene.Matrix.filename = paste0(SetupDIR,"Win2Gene.rds")
  print("Loading Window to Gene mapping")
  L = readRDS(Win2Gene.Matrix.filename)
  Win2Gene.matrix = L$Matrix
  MultiPromoterGenes = L$Multi
  rm(L)
  #  colnames(Win2Gene.matrix)[is.na(colnames(Win2Gene.matrix))] = ""
  Genes = colnames(Win2Gene.matrix)
  GeneWindows = as.integer(rownames(Win2Gene.matrix))
  
  W = width(TSS.windows[GeneWindows])/1000
  GeneLength = ComputeGeneCounts(W)
  names(GeneLength) = Genes
  rm(W)
  
  Win2Gene.triplet = as(Win2Gene.matrix, "TsparseMatrix")
  Win2Gene.rows = Win2Gene.triplet@i +1
  Win2Gene.cols = Win2Gene.triplet@j +1
  Win2Gene.lists = split(Win2Gene.rows, Win2Gene.cols)
  names(Win2Gene.lists) = Genes
  rm(Win2Gene.cols)
  rm(Win2Gene.rows)
  rm(Win2Gene.triplet)
}


cfChIP.Params <- function() {
  list( 
    Save = TRUE,
    DataDir =NULL,
    Background = FALSE,
    GeneCounts = FALSE,
    GeneBackground = FALSE,
    Normalize = FALSE,
    OverExpressedGenes = FALSE,
    reuseSavedData = TRUE,
    Verbose = TRUE,
    TSS.windows = TSS.windows,
    Win2Gene = Win2Gene.matrix,
    #GeneWindows = GeneWindows,
    #CommonGenes = CommonGenes,
    #NormRef = Healthy.GeneCount,
    #NormRef.var = Healthy.GeneCount.var,
    #WinNormRef = Healthy.WinCount,
    #WinNormRef.var = Healthy.WinCount.var,
    #Signatures = Win.Sig,
    #Signatures.Ref = Win.Sig.Ref,
    MinFragLen = 50,
    MaxFragLen = 800,
    Heatmap.Type = "dot",
    Heatmap.MaxCount = 15,
    Heatmap.Prune.Cols = FALSE,
    Heatmap.Prune.Rows = FALSE,
    Heatmap.Cluster.Cols = FALSE,
    Heatmap.Cluster.Rows = FALSE,
    Heatmap.Qvalue.threshold = 0.0001,
    Heatmap.Filter.threshold = 6,
    PlotPrograms.GlobalHeatmap = TRUE,
    Programs.Eval.Reference = TRUE,
    Signatures.Eval.Reference = FALSE,
    Heatmap.Detailed.Type = "plain",
    Heatmap.Detailed.Prune.Cols = FALSE,
    Heatmap.Detailed.Prune.Rows = FALSE,
    Scatter.NormalizeByLength = FALSE,
    Scatter.RefSeqOnly = FALSE,
    Scatter.MarkGenePrograms = NULL,
    Scatter.UpperQuantile = 0.995,
    Scatter.LimitUnits = 2,
    #Programs = Gene.Programs,
    #Programs.Ref = Gene.Programs.Ref,
    #Programs.Partition = Gene.Programs.Partition,
    PlotEnrichments.MaxCount = 50,
    PlotSignatures.IndividualHeatmap = FALSE,
    PlotSignatures.MaxCount = 10,
    PlotSignatures.MaxZscore = 20,
    PlotSignatures.Height = 8,
    PlotSignatures.Width = 10,
    PlotPrograms.MaxCount = 10,
    PlotPrograms.IndividualHeatmap = FALSE,
    PlotPrograms.IndividualBarChart = FALSE,
    PlotPrograms.IndividualCSV = FALSE,
    PlotPrograms.IndividualCDT = FALSE,
    PlotPrograms.MaxZscore = 20,
    PlotPrograms.Height = 8,
    PlotPrograms.Width = 10,
    PlotPrograms.RowOrder = NA,
    #MetaGene = MetaGene,
    MetaGene.Offset = 5000,
    MetaGene.Width = 25000,
    MetaGene.Tick = 5000,
    MetaGene.Label = "TSS",
    MetaGene.Max = -1,
    MetaGene.Color = "red",
    MetaGene.BGColor = "black",
    #MetaEnhancer = MetaEnhancer,
    MetaEnhancer.Offset = 25000,
    MetaEnhancer.Width = 50000,
    MetaEnhancer.Tick = 10000,
    MetaEnhancer.Label = "Enhancer",
    MetaEnhancer.Max = -1,
    MetaEnhancer.Color = "red",
    MetaEnhancer.BGColor = "black",
    #    BackgroundModel = NULL
    QC.PositiveType = "TSS",
    QC.PositiveNucs = 192015, # H3K4me3 in humans
    QC.GenomeLength = 3.3e9,
    QC.SampleVolume = 2, # 2ml per sample
    QC.QuantileCutoff = 0.75, # cutoff of high peak coverage (see below)
    QC.GenomeWeight = 6.6e-12, # molecular weight of 1 genome # CHECK. see http://www.bio.net/bionet/mm/methods/1999-December/080037.html
    QC.NucleosomeLength = 200,
    QC.cfDNA.ng = 10,
    QCFields = paste("Total","Total.uniq","Total.uniq.est","TSS",
                     "%Signal Total","%Background Total",
                     "%Signal TSS","%Background TSS",
                     "lambda.mix.poiss","Mito","%Mito","Global.signal.yield",
                     "Local.signal.yield", "BG.yield","Global.SNR", "Local.SNR","Seq.factor", sep=";")
  )
}

BaseFileName <- function( fname, extList = c(".bed$", ".rdata$", ".bw$", ".tagAlign$" #, "-H3K4me3", "-H3K4me2", "-H3K4me1", "-H3K36me3"
) ) {
  #  x = file_path_sans_ext(fname)
  x = fname
  x = sub(".gz$", "",x)
  for( ext in extList )
    x = sub(ext, "", x)
  
  y = strsplit(x,"/")[[1]]
  n = length(y)
  z = y[n]
  return(z)
}

# The cfChIP.BuildFN function constructs a file name/path for saving data, based on the provided parameters.
cfChIP.BuildFN <- function(Name, param, suff = ".rdata") {
  # Check if the DataDir parameter is specified in the 'param' object.
  if(is.null(param$DataDir))
    f = ""  # If not specified, set 'f' to an empty string.
  else {
    f = DataDir  # Use the global DataDir variable as the default directory path.
    # Ensure the directory path ends with a forward slash. If not, append it.
    if(!grepl("/$", DataDir))
      f = paste0(param$DataDir, "/")
  }
  # Construct the full file name/path by combining the directory path, file name, and file suffix.
  paste0(f, Name, suff)
}

# The fixCoverage function is designed to ensure that the coverage data for each chromosome
# is of the correct length, matching the expected length of the chromosome.
fixCoverage = function(Cov) {
  # Iterate over each chromosome listed in ChrList.
  for(c in ChrList) {
    # Determine the current length of the coverage data for the chromosome.
    l = length(Cov[[c]])
    # Retrieve the expected length of the chromosome from the genome sequence information.
    m = seqlengths(genome.seqinfo[c])
    # If the current length of the coverage data is less than the expected length,
    # append zeros to the coverage data to match the expected length.
    if(l < m)
      Cov[[c]] = append(Cov[[c]], rep(0, m - l))
  }
  # Return the adjusted coverage data.
  Cov
}

cfChIP.BED.suffixes = c(".bed", ".bed.gz", ".tagAlign", ".tagAlign.gz")
cfChIP.BW.suffixes = c(".bw", ".bigWig", ".bw.gz", ".bigWig.gz")
cfChIP.File.suffixes = c(cfChIP.BED.suffixes, cfChIP.BW.suffixes)

cfChIP.FindFile <- function( filename, param = cfChIP.Params() ) {
  FileType  = NA
  BED.suffixes = cfChIP.BED.suffixes
  BW.suffixes = cfChIP.BW.suffixes
  if( any(sapply(BED.suffixes, function(s) grepl(paste0(s,"$"), filename))))
    FileType = "BED"
  if(  any(sapply(BW.suffixes, function(s) grepl(paste0(s,"$"), filename))))
    FileType = "BW"
  
  if( is.na(FileType) ) {
    for( s in BED.suffixes )
      if( file.exists(paste0(filename, s))) {
        filename = paste0(filename, s)
        FileType = "BED"
      } 
    if( is.na(FileType) ) 
      for( s in BW.suffixes )
        if( file.exists(paste0(filename, s))) {
          filename = paste0(filename, s)
          FileType = "BW"
        } 
  }
  
  if( is.na(FileType ) ) {
    catn(filename, ": Error, cannot determine file type of ",filename)
    return(NULL)
  }
  
  return(list(filename = filename, FileType = FileType))
}

cfChIP.RawData.Cache = list()
cfChIP.RawData.CacheMaxSize = 100

cfChIP.GetRawData = function(filename, param = cfChIP.Params) {
  ll = cfChIP.FindFile(filename, param)
  filename = ll$filename
  FileType = ll$FileType
  
  if( filename %in% cfChIP.RawData.Cache )  
    return(cfChIP.RawData.Cache[[filename]])
  
  if( length(cfChIP.RawData.Cache) >= cfChIP.RawData.CacheMaxSize )
    cfChIP.RawData.Cache <<- list()
  
  dat = list()
  if( FileType == "BED") {     
    if( param$Verbose ) catn(filename, ": Reading BED file")
    
    dat$RawBED = import(filename, format = "BED")
    # Annotate the genome (assuming hg19 for this example)
    seqlevelsStyle(dat$RawBED) <- "UCSC"
    genome(dat$RawBED) <- "hg38"
    
    # Filter to keep only standard chromosomes
    standard_chr <- c(paste0("chr", 1:22), "chrX", "chrY")
    dat$RawBED <- dat$RawBED[seqnames(dat$RawBED) %in% standard_chr]
    
    # remove long/short fragments and non-unique copies
    
    # check for single end reads
    if( max(width(dat$RawBED)) <  param$MinFragLen) {
      dat$RawBED = resize(dat$RawBED, width = 166)
    } else 
      dat$RawBED = dat$RawBED[width(dat$RawBED) <= param$MaxFragLen & width(dat$RawBED) > param$MinFragLen]
    
    dat$BED = unique(dat$RawBED)
    dat$Cov = coverage(dat$BED)
  } 
  cfChIP.RawData.Cache[[filename]] <<- dat
  if( FileType == "BW" ) {     
    if( param$Verbose ) catn(filename, ": Reading BigWig file", filename)
    dat$BW = import(filename)
    dat$Cov = coverage(dat$BW, weight="score")
  } 
  return(dat)
}

cfChIP.GetCoverage  = function(filename, param=cfChIP.Params()) {
  dat = cfChIP.GetRawData(filename, param)
  return(fixCoverage(dat$Cov))
}

# Define the cfChIP.ProcessFile function which processes a given file (either a new file or precomputed data) for cfChIP analysis.
cfChIP.ProcessFile <- function(filename = NULL,
                               dat = NULL,
                               param = cfChIP.Params(),
                               Force = FALSE,
                               HardForce = FALSE,
                               Change = FALSE) {
  # Check if both 'filename' and 'dat' are NULL. If so, print an error message and return NULL.
  if(is.null(filename) && is.null(dat)) {
    catn("Need one of filename or dat be assigned!")
    return(NULL)
  }
  
  # If 'dat' is NULL, extract the base name of the file without its directory path or extension.
  # This name is used for further processing and saving the data.
  if(is.null(dat)) {
    Name = BaseFileName(filename)
  } else {
    # If 'dat' is not NULL, it means precomputed data is being processed. Use the name from 'dat'.
    Name = dat$Name
    if(param$Verbose) catn(Name, ": Processing precomputed data")
  }
  
  # Construct the filename for saving or reading the processed data using the cfChIP.BuildFN function.
  fn = cfChIP.BuildFN(Name, param)
  
  # If 'dat' is NULL, check if the data should be reused (if it exists) or if new data should be processed.
  if(is.null(dat)) {
    if(param$reuseSavedData && file.exists(fn)) {
      # If reusing saved data and the file exists, read the precomputed data.
      if(param$Verbose) catn(Name, ": Reading precomputed data", fn)
      dat <- readRDS(fn)
    } else {
      # If not reusing or the file doesn't exist, initialize 'dat' with default values.
      dat = list(Name = Name, 
                 Cov = NULL,
                 Counts = NULL,
                 Heights = NULL,
                 Background = NULL,
                 GeneCounts = NULL,
                 GeneHeights = NULL,
                 GeneBackground = NULL, 
                 Counts.QQnorm = NULL,
                 GeneCounts.QQnorm = NULL, 
                 QQNorm = 1,
                 OverExpressedGenes = NULL)
      Change = TRUE  # Indicate that the data has changed (or is new).
    }
  }
  
  # The rest of the function would typically include further processing of 'dat',
  # but this snippet ends here. The function's purpose is to manage the loading and initialization
  # of data for analysis, handling both new and precomputed datasets.

  
  # Enforce dependencies between different parameters to ensure logical consistency in the analysis.
  # If normalization or overexpressed genes analysis is requested, ensure related parameters are set accordingly.
  param$Normalize = param$Normalize || param$OverExpressedGenes
  param$GeneCounts = param$GeneCounts || param$Normalize
  param$GeneBackground = param$GeneBackground || param$Normalize
  param$Background = param$Background || param$GeneBackground
  
  # If HardForce is TRUE, reset the Counts data and enforce reprocessing.
  if(HardForce) {
    dat$Counts = NULL
    Force = TRUE
  }
  
  # If Force is TRUE, clear various data fields to ensure data is reprocessed rather than reused.
  if(Force) {
    dat$Cov = NULL
    dat$Background = NULL
    dat$GeneCounts = NULL
    dat$GeneBackground = NULL
    dat$Counts.QQnorm = NULL
    dat$GeneCounts.QQnorm = NULL
    dat$QQNorm = 1
  }
  
  # Process BED file to obtain genomic data if Counts data is not already available.
  if(is.null(dat$Counts)) {
    # Retrieve raw data from the file, which includes BED, BW, and coverage data.
    dd = cfChIP.GetRawData(filename, param)
    dat$BED = dd$BED
    dat$BW = dd$BW
    dat$Cov = dd$Cov
    
    # If raw BED data is available, calculate duplicate counts.
    if(!is.null(dd$RawBED)) {
      RawBED.dups <- countOverlaps(dat$BED, dd$RawBED, type = "equal")
      dat$DupCount = as.data.frame(table(RawBED.dups))
    }
    
    # Calculate fragment counts if BED data is available.
    if(!is.null(dat$BED))
      dat$FragCount = as.data.frame(table(width(dat$BED)), stringsAsFactors = FALSE)
    
    Change = TRUE
    
    # Count fragment overlaps with TSS windows if BED data is available.
    if(!is.null(dat$BED)) {
      if(param$Verbose) catn(Name, ": counting fragment overlap")
      dat$Counts = countOverlaps(query = param$TSS.windows, 
                                 subject = resize(dat$BED, width=1, fix="center"))
    } else if(!is.null(dat$BW)) {
      # If BigWig data is available, calculate coverage and counts.
      if(param$Verbose) catn(Name, ": counting BigWig overlap")
      dat$Cov = coverage(dat$BW, weight="score")
      dat$Counts = rep(0, length(param$TSS.windows))
      
      # Ensure coverage for chromosome Y is accounted for.
      if(!("chrY" %in% names(dat$Cov)))
        dat$Cov[["chrY"]] = Rle(0, seqlengths(param$TSS.windows)["chrY"])
      
      # Aggregate coverage data across TSS windows.
      ChrRle = Rle(chrom(params$TSS.windows))
      ChrStarts = start(ChrRle)
      ChrEnds = end(ChrRle)
      ChrName = as.character(runValue(ChrRle))
      for(i in 1:nrun(ChrRle)) {
        catn(ChrName[i])
        ws = ChrStarts[i]:ChrEnds[i]
        dat$Counts[ws] = aggregate(dat$Cov[[ChrName[i]]], 
                                   ranges(params$TSS.windows)[ws], 
                                   sum)
      }
      # Normalize counts assuming a typical read length of 200bp.
      dat$Counts = dat$Counts / 200
      
      # Calculate maximum coverage heights across TSS windows.
      dat$Heights = max(dat$Cov[params$TSS.windows])
    } else {
      # If neither BED nor BigWig data could be processed, log an error.
      catn(Name, ": error! cannot compute counts")
      return(dat)
    }
    
    # Normalize coverage data across all chromosomes.
    dat$Heights = rep(0, length(params$TSS.windows))
    dat$Cov = fixCoverage(dat$Cov)
    for(chr in unique(chrom(params$TSS.windows))) {
      ww = which(chrom(params$TSS.windows) == chr)
      dat$Heights[ww] = max(dat$Cov[params$TSS.windows[ww]])
    }
  }
  # Remove large data objects from the 'dat' object to save memory.
  if( !is.null(dat$BED))  
    dat$BED = NULL  # Remove BED data.
  if( !is.null(dat$BW) )
    dat$BW = NULL  # Remove BigWig data.
  if( !is.null(dat$Cov))
    dat$Cov = NULL  # Remove coverage data.
  
  # Compute the background model if it's requested and not already computed.
  if( param$Background && is.null(dat$Background) ) {
    if( param$Verbose ) catn(Name, ": Computing background model")
    # Use the buildBackground function to compute the background model based on counts and TSS windows.
    dat$Background = buildBackground(Y = dat$Counts, TWin = params$TSS.windows)
    # Reset related fields as the background computation might affect their values.
    dat$GeneBackground = NULL
    dat$Counts.QQnorm = NULL
    dat$GeneCounts.QQnorm = NULL
    Change = TRUE  # Indicate that the data has changed.
  }
  
  # Compute gene counts if it's requested and not already computed.
  if( param$GeneCounts && is.null(dat$GeneCounts)) {
    if( param$Verbose ) catn(Name, ": Computing gene counts")
    # Compute gene counts using the ComputeGeneCounts function.
    dat$GeneCounts = ComputeGeneCounts(dat$Counts[GeneWindows], params$Win2Gene)
    # If the result is a matrix, extract the first column (assuming it's the relevant data).
    if(is.matrix(dat$GeneCounts))
      dat$GeneCounts = dat$GeneCounts[,1]
    # Assign gene names to the gene counts.
    names(dat$GeneCounts) = Genes
    # Reset related fields as the gene counts computation might affect their values.
    dat$GeneCounts.QQnorm = NULL
    dat$GeneHeights = MaxGeneCounts(dat$Heights[GeneWindows], param$Win2Gene)
    Change = TRUE
  }
  
  # Compute gene background if it's requested and not already computed.
  if( param$GeneBackground && is.null(dat$GeneBackground)) {
    if( param$Verbose ) catn(Name, ": Computing gene background")
    # Use the background model to estimate gene background.
    mu = dat$Background
    Z = getMultiBackgroundEstimate(mu, param$GeneWindows)
    dat$GeneBackground = ComputeGeneCounts(Z, param$Win2Gene)
    # If the result is a matrix, extract the first column.
    if(is.matrix(dat$GeneBackground))
      dat$GeneBackground = dat$GeneBackground[,1]
    # Assign gene names to the gene background.
    names(dat$GeneBackground) = Genes
    # Compute gene copy number variation based on gene background and genome data.
    dat$GeneCNV = dat$GeneBackground / (GeneLength * dat$Background$genome)
    dat$GeneCounts.QQnorm = NULL
    Change = TRUE
  }
  
  
  
  # Return the processed data object.
  return(dat)
}

# Normalize data if it's requested and not already normalized.
if( param$Normalize && is.null(dat$GeneCounts.QQnorm)) {
  if( param$Verbose ) catn(Name, ": Normalize")
  # Compute the difference between gene counts and gene background, ensuring no negative values.
  GeneDiff = pmax(dat$GeneCounts - dat$GeneBackground, 0)
  # Prepare data for normalization.
  A = cbind(GeneDiff, param$NormRef)
  colnames(A) = c("Sample", "Ref")
  # Normalize gene counts using the QQNormalizeGenes function.
  Qs = QQNormalizeGenes(A, CommonG = param$CommonGenes)
  dat$QQNorm = Qs[1] / Qs[2]
  dat$GeneCounts.QQnorm = GeneDiff * dat$QQNorm
  names(dat$GeneCounts.QQnorm) = Genes
  dat$OverExpressedGenes = NULL
  # Estimate background throughout the data.
  dat$WinBackground = getMultiBackgroundEstimate(dat$Background, 1:length(param$TSS.windows))
  dat$Counts.QQnorm = pmax(dat$Counts - dat$WinBackground, 0) * dat$QQNorm
  Change = TRUE
}

# Compute overexpressed genes if it's requested and not already computed.
if(param$OverExpressedGenes && is.null(dat$OverExpressedGenes)) {
  if(param$Verbose) catn(Name, ": Computing overexpressed genes")
  # Use the ComputeOverExpressed function to identify overexpressed genes.
  dat$OverExpressedGenes = cfChIP.ComputeOverExpressed(dat$GeneCounts.QQnorm, 
                                                       dat$GeneCounts, 
                                                       dat$GeneBackground, 
                                                       param$NormRef, 
                                                       param$NormRef.var, 
                                                       dat$QQNorm)
  Change = TRUE
}

# Compute overexpressed windows if it's requested, not already computed, and window normalization reference is available.
if(param$OverExpressedGenes && is.null(dat$OverExpressedWins) && !is.null(param$WinNormRef)) {
  if(param$Verbose) catn(Name, ": Computing overexpressed windows")
  # Use the ComputeOverExpressed function to identify overexpressed windows.
  dat$OverExpressedWins = cfChIP.ComputeOverExpressed(dat$Counts.QQnorm, 
                                                      dat$Counts, 
                                                      dat$WinBackground, 
                                                      param$WinNormRef, 
                                                      param$WinNormRef.var, 
                                                      dat$QQNorm)
  rownames(dat$OverExpressedWins) = 1:nrow(dat$OverExpressedWins)
  Change = TRUE
}

# Save the processed data if changes have been made and saving is enabled.
if(Change && param$Save) {
  if(param$Verbose) catn(Name, ": Saving data")
  saveRDS(dat, fn)
}