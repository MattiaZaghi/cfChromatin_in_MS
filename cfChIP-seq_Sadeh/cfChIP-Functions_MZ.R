library(Biobase)

source(paste0(SourceDIR,"cfChIP-util.R"))
source(paste0(SourceDIR, "YieldEstimation-Functions.R"))
source(paste0(SourceDIR, "EstimateGammaPoisson.R"))
source(paste0(SourceDIR, "ComputeGeneCounts.R"))
source(paste0(SourceDIR, "Background.R"))
source(paste0(SourceDIR, "NMF-util.R"))


#source(paste0(SourceDIR, "Estimate-ctDNA.R"))

mem.maxVSize(vsize = Inf)

#not the right place, but for now
WinDescription = data.frame(Type = TSS.windows$type, 
                            Gene = TSS.windows$name, 
                            Tissue = TSS.windows$tissue, 
                            stringsAsFactors = FALSE)
rownames(WinDescription) = 1:nrow(WinDescription)


logSum = function(a, b) {
  m = max(a,b)
  m + log(exp(a-m)+exp(b-m))
}

logAvg = function(a, b) {
  m = max(a,b)
  m + log(exp(a-m)+exp(b-m))-log(2)
}

twoTailedPValue = function(a,b) {
  pmin(a,b) + log(2)
}

cfChIP.Params <- function() {
  list( 
    Save = TRUE,
    DataDir = NULL,
    Background = TRUE,
    GeneCounts = TRUE,
    GeneBackground = TRUE,
    Normalize = FALSE,
    OverExpressedGenes = FALSE,
    reuseSavedData = TRUE,
    Verbose = TRUE,
    TSS.windows = TSS.windows,
    Win2Gene = Win2Gene.matrix,
    GeneWindows = GeneWindows,
    MinFragLen = 50,
    MaxFragLen = 800
  )
}



fixCoverage = function( Cov ) {
  for( c in ChrList )
  {
    l = length(Cov[[c]])
    m = seqlengths(genome.seqinfo[c])
    if( l < m )
      Cov[[c]] = append(Cov[[c]], rep(0,m-l))
  }
  Cov
}


cfChIP.BuildFN <- function(Name, param, suff = ".rdata" ) {
  if( is.null(param$DataDir) )
    f = ""
  else {
    f = DataDir
    if( !grepl("/$", DataDir))
      f = paste0(param$DataDir,"/")
  }
  paste0(f, Name, suff)
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
    # remove long/short fragments and non-unique copies
    
    # check for single end reads
    if( max(width(dat$RawBED)) <  param$MinFragLen) {
      dat$RawBED = resize(dat$RawBED, width = 166)
    } else 
      dat$RawBED = dat$RawBED[width(dat$RawBED) <= param$MaxFragLen & width(dat$RawBED) > param$MinFragLen]
    
    dat$BED = unique(dat$RawBED)
    dat$Cov = coverage(dat$BED)
  } 
  if( FileType == "BW" ) {     
    if( param$Verbose ) catn(filename, ": Reading BigWig file", filename)
    dat$BW = import(filename)
    dat$Cov = coverage(dat$BW, weight="score")
  } 
  cfChIP.RawData.Cache[[filename]] <<- dat
  
  return(dat)
}

cfChIP.GetCoverage  = function(filename, param=cfChIP.Params()) {
  dat = cfChIP.GetRawData(filename, param)
  return(fixCoverage(dat$Cov))
}

cfChIP.ComputeOverExpressed = function(X.norm, Counts, BG, NormRef, NormRef.var, QQNorm) {
  X = X.norm
  Y = Counts
  logX = log2(X+1)
  logH = log2(NormRef+1)
  Diff.log = logX - logH
  
  Lam = NormRef / QQNorm + BG 
  Lam[is.na(Lam)] = 100
  Pv.up = -computeMultiPValue(Y, BG, NormRef, NormRef.var, 1/QQNorm)
  Pv.down = -computeMultiPValue(Y, BG, NormRef, NormRef.var, 1/QQNorm, Above = FALSE)
  Pv = pmax(Pv.up,Pv.down) - log(2)
  Qv = -log10(p.adjust(exp(-Pv),method="fdr"))
  
  if(!is.null(NormRef.var))
    Zscore =  (Y - Lam)/sqrt(NormRef.var/(QQNorm**2) + Lam)
  else
    Zscore = NA
  
  data.frame(healthy = logH, 
             sample = logX, 
             obs = Y, 
             exp = Lam, 
             qvalue = Qv, 
             pvalue = Pv/log(10), 
             X = X, 
             H = NormRef, 
             Significant.up = (Qv > 3 & ((logX-logH > 1) & (Y > Lam))),
             Significant.down = (Qv > 3 & ((logX-logH < -1) & (Y < Lam)) ),
             Zscore = Zscore,
             stringsAsFactors = FALSE)
}

cfChIP.ProcessFile <- function( filename = NULL,
                                dat = NULL,
                                param = cfChIP.Params(),
                                Force = FALSE,
                                HardForce = FALSE,
                                Change = FALSE ) 
{
  if( is.null( filename ) && is.null( dat ))
  {
    catn("Need one of filename or dat be assigned!")
    return(NULL)
  }
  if( is.null( dat ) )
  {
    Name = BaseFileName(filename)
  } else {
    Name = dat$Name
    if( param$Verbose ) catn(Name, ": Processing precomputed data")
  }
  fn = cfChIP.BuildFN(Name, param )
  
  if( is.null(dat) ) {
    if( param$reuseSavedData && file.exists(fn) ) {
      if( param$Verbose ) catn(Name, ": Reading precomputed data", fn)
      dat <- readRDS(fn)
    } else  {
      dat = list(Name = Name, 
                 Cov = NULL,
                 Counts = NULL,
                 Heights = NULL,
                 Background = NULL,
                 GeneCounts = NULL,
                 GeneHeights = NULL,
                 GeneBackground = NULL)
      Change = TRUE
    }
  }
  
  
  #enforce dependencies
  param$GeneBackground = param$GeneBackground || param$Normalize
  param$Background = param$Background || param$GeneBackground
  
  
  
  if( HardForce) {
    dat$Counts = NULL
    Force = TRUE
  }
  
  if( Force ) {
    dat$Cov = NULL
    dat$Background = NULL
    dat$GeneCounts = NULL
    dat$GeneBackground = NULL 
  }
  
  # Get BED reads into a GenomicRanges object
  if( is.null(dat$Counts) ) {
    dd = cfChIP.GetRawData(filename, param)
    dat$BED = dd$BED
    dat$BW = dd$BW
    dat$Cov  = dd$Cov
    
    if(!is.null(dd$RawBED)) {
      RawBED.dups <- countOverlaps(dat$BED, dd$RawBED, type = "equal")
      dat$DupCount = as.data.frame(table(RawBED.dups))
    }
    
    if( !is.null(dat$BED))
      dat$FragCount = as.data.frame(table(width(dat$BED)), stringsAsFactors = FALSE)
    
    Change = TRUE
    
    if(!is.null(dat$BED)) {
      if( param$Verbose ) catn(Name, ": counting fragment overlap")
      dat$Counts = countOverlaps(query = param$TSS.windows, 
                                 subject = resize(dat$BED, width=1, fix="center"))
    } else
      if( !is.null(dat$BW)) {
        if( param$Verbose ) catn(Name, ": counting BigWig overlap")
        dat$Cov = coverage(dat$BW, weight="score")
        dat$Counts = rep(0,length(param$TSS.windows))
        
        if( !( "chrY" %in% names(dat$Cov)) ) 
          dat$Cov[["chrY"]] = Rle(0,seqlengths(param$TSS.windows)["chrY"])
        
        ChrRle = Rle(chrom(param$TSS.windows))
        ChrStarts = start(ChrRle)
        ChrEnds = end(ChrRle)
        ChrName = as.character(runValue(ChrRle))
        for( i in 1:nrun(ChrRle)) {
          catn(ChrName[i])
          ws = ChrStarts[i]:ChrEnds[i]
          dat$Counts[ws] = aggregate(dat$Cov[[ChrName[i]]], 
                                     ranges(param$TSS.windows)[ws], 
                                     sum)
        }
        
        # assuming a typical read is 200bp
        dat$Counts = dat$Counts/200
        
        dat$Heights = max(dat$Cov[param$TSS.windows])
        
      } else {
        catn(Name, ": error! cannot compute counts")
        return(dat)
      }
    dat$Heights = rep(0, length(param$TSS.windows))
    dat$Cov = fixCoverage(dat$Cov)
    for( chr in unique(chrom(param$TSS.windows))) {
      ww = which(chrom(TSS.windows) == chr)
      dat$Heights[ww] = max(dat$Cov[TSS.windows[ww]])
    }
  }
  # remove BED, BW, and Cov
  if( !is.null(dat$BED))  
    dat$BED = NULL
  if( !is.null(dat$BW) )
    dat$BW = NULL
  if( !is.null(dat$Cov))
    dat$Cov = NULL
  
  # Background
  if( param$Background && is.null(dat$Background) ) {
    if( param$Verbose ) catn(Name, ": Computing background model")
    dat$Background = buildBackground(Y = dat$Counts, TWin = param$TSS.windows)
    dat$GeneBackground = NULL 
    dat$Counts.QQnorm = NULL
    dat$GeneCounts.QQnorm = NULL
    Change = TRUE
  }
  
  # Gene counts
  if( param$GeneCounts && is.null(dat$GeneCounts)) {
    if( param$Verbose ) catn(Name, ": Computing gene counts")
    dat$GeneCounts =  ComputeGeneCounts(dat$Counts[GeneWindows], param$Win2Gene)
    if(is.matrix(dat$GeneCounts))
      dat$GeneCounts = dat$GeneCounts[,1]
    names(dat$GeneCounts) = Genes
    dat$GeneCounts.QQnorm = NULL
    dat$GeneHeights = MaxGeneCounts(dat$Heights[GeneWindows], param$Win2Gene)
    
    Change = TRUE
  }
  
  # GeneBackground
  if( param$GeneBackground && is.null(dat$GeneBackground)) {
    if( param$Verbose ) catn(Name, ": Computing gene background")
    mu = dat$Background
    Z = getMultiBackgroundEstimate(mu,param$GeneWindows)
    dat$GeneBackground =  ComputeGeneCounts(Z,param$Win2Gene)
    if(is.matrix(dat$GeneBackground))
      dat$GeneBackground = dat$GeneBackground[,1]
    names(dat$GeneBackground) = Genes
    dat$GeneCNV = dat$GeneBackground/(GeneLength*dat$Background$genome)
    
    dat$GeneCounts.QQnorm = NULL
    Change = TRUE
  }
  
  # we are done!
  if( Change && param$Save ) {
    if( param$Verbose ) catn(Name, ": Saving data")
    saveRDS(dat,fn)
  }
  
  return( dat )
}

cfChIP.BuildOutputName = function( dat, Dir, Prefix = NULL, Suffix = ".pdf") {
  if( is.null(Dir) ) 
    Dir = "./"
  if( !grepl("/$", Dir))
    Dir = paste0(Dir,"/")
  
  if( is.null(Prefix) )
    Prefix = ""
  
  TargetDir = paste0(Dir, Prefix)
  TargetDir = sub("/$","", TargetDir)
  if( !dir.exists(TargetDir) )
    dir.create(TargetDir)
  
  fname = paste0(TargetDir, "/", dat$Name, Suffix)
  return(fname)
}







