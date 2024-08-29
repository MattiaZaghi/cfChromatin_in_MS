# convert list of objects to dataframe of metadata
listToMeta <- function(obj) {
  for (i in 1:length(obj)) {
    mdata <- obj[[i]][[]]
    if (i==1) {
      final <- mdata
    } else {
      final <- rbind(final,mdata)
    }
  }
  return(final)
}


# Plot UMI counts per each modality and sample
plotCounts <- function(obj,quantiles,feature,ylabel=feature) {
  dataf <- listToMeta(obj)
  # plot
  pp=ggplot(dataf, aes(x="x",y=dataf[,feature],fill=sample)) +
    theme_bw() +
    geom_half_violin(draw_quantiles = .5) +
    geom_half_point(shape=1,aes(color=sample)) +
    facet_grid(sample~modality,scales = "free_x") +
    theme(strip.text.x = element_text(size = 11, colour = "black", angle = 0, face= 'bold')) +
    theme(strip.text.y = element_text(size = 11, colour = "black", face= 'bold')) +
    # x axis
    xlab("Samples") +
    theme(axis.text.x=element_text(size = 0,angle = 0, hjust = .5)) +
    theme(axis.title.x = element_text(size=0)) +
    # y axis
    ylab(ylabel) +
    theme(axis.text.y=element_text(size = 12)) +
    theme(axis.title.y = element_text(size=14)) +
    theme(legend.position = "none") +
    # add quantiles
    stat_summary(fun = "quantile", fun.args = list(probs = quantiles), 
                 geom = "hline", aes(yintercept = ..y..), linetype = "dashed")
    return(pp)
}

# modified version of the Signac DepthCor - it just adds to the plot the modality name
DepthCorMulMod <- function(obj) {
  nm <- unique(obj[[]][,"modality"])
  p1 <- DepthCor(obj) +
    ggtitle(nm) +
    theme(plot.title = element_text(size=14,hjust=0.5,face='bold'))
  return(p1)
}


# Umap with connected modalities 
plotConnectModal <- function(seurat,group) {
  nmodalities <- length(seurat)
  
  # First, get coords of UMAP, modality name, cluster name and cell barcode for each modality
  umap_embeddings <- list()
  i <- 0
  for(mod in names(seurat)) {
    i <- i + 1
    # get modality name without barcode
    mod2 <- strsplit(mod,"_")[[1]][1]
    # get UMAP 1 and 2 coords
    umap_embeddings[[mod2]]               <- as.data.frame(seurat[[mod]]@reductions[['umap']]@cell.embeddings)
    # adjust UMAP1 for plotting pursoses
    umap_embeddings[[mod2]]$UMAP_1        <- umap_embeddings[[mod2]]$UMAP_1 + (i-1)*40 
    # add modality name
    umap_embeddings[[mod2]]$modality      <- unique(seurat[[mod]]$modality) 
    # add cluster name
    umap_embeddings[[mod2]]$cluster       <- seurat[[mod]]@meta.data[,group] 
    # add cell barcode
    umap_embeddings[[mod2]]$cell_barcode  <- rownames(umap_embeddings[[mod2]]) 
  }
  
  # convert the list to dataf
  umap.embeddings.merge <- purrr::reduce(umap_embeddings,rbind)
  # get the names of the common cells among the analysed modalities
  common.cells                        <- table(umap.embeddings.merge$cell_barcode)
  common.cells                        <- names(common.cells[common.cells==nmodalities])
  # subset the umap_embeddings, selecting only the info from the common set of cells
  umap.embeddings.merge               <- umap.embeddings.merge[umap.embeddings.merge$cell_barcode %in% common.cells,]
  
  # label modality - get coords
  coords <- aggregate(umap.embeddings.merge$UMAP_1,by=list(umap.embeddings.merge$modality),max)
  names(coords) <- c("modality","max")
  coords$min <- aggregate(umap.embeddings.merge$UMAP_1,by=list(umap.embeddings.merge$modality),min)[,"x"]
  coords$mid_point <- (coords$max + coords$min) / 2
  
  plot <- ggplot(data=umap.embeddings.merge,aes(x=UMAP_1,y=UMAP_2,col=cluster)) + 
    geom_point(size=0.2) + 
    geom_line(data=umap.embeddings.merge, aes(group=cell_barcode,col=cluster),alpha=0.2,size=0.02) + 
    theme_classic() + NoAxes() +
    guides(color = guide_legend(override.aes = list(size=5),title="")) +
    geom_text(data=coords,aes(label=modality,x=mid_point,y=max(umap.embeddings.merge$UMAP_2+.1*umap.embeddings.merge$UMAP_2)),
              colour='black', fontface = "bold",size=4)
  return(plot)
  
}


# Function to plot upset on peaks of each modality
getUpsetPeaks <- function(modalities, samples, combined_peaks_ls, input_ls) {
  
  list_up <- list()
  for (mod in modalities) {
    i <- 0
    for (smp in samples) {
      i <- i + 1
      
      # overlap
      combined_mod <- as.data.frame(combined_peaks_ls[[mod]])
      overlap <- GenomicRanges::findOverlaps( toGRanges(combined_mod), input_ls[[paste0(mod,"_",smp)]], select = "first")
      # convert NA to 0s and make it binary
      overlap[is.na(overlap)] = 0
      overlap <- ifelse(overlap>0,1,0)
      combined_mod[,smp] <- overlap
      
      # add new values
      if (i==1) {
        final <- combined_mod
      } else {
        final[,smp] <- overlap
      }
      
      # if last sample, append to list
      if (smp==samples[length(samples)]) {
        list_up[[mod]] <- final
      }
      
      
    }
  }
  
  if (length(list_up)==1) {
    pfinal=upset(list_up[[1]][,6:ncol(list_up[[1]])],colnames(list_up[[1]][,6:ncol(list_up[[1]])]),min_size=10,width_ratio=.3,
                 set_sizes=(upset_set_size() + theme(axis.text.x=element_text(angle=45,hjust=1,size=8))),intersections="all",
                 base_annotations=list('Size'=(intersection_size(counts=FALSE))))  + xlab(names(list_up)[1]) + theme(axis.title.x = element_text(size=12))
  } else if (length(list_up)==2) {
    p1=upset(list_up[[1]][,6:ncol(list_up[[1]])],colnames(list_up[[1]][,6:ncol(list_up[[1]])]),min_size=10,width_ratio=.3,
             set_sizes=(upset_set_size() + theme(axis.text.x=element_text(angle=45,hjust=1,size=8))),intersections="all",
             base_annotations=list('Size'=(intersection_size(counts=FALSE)))) + xlab(names(list_up)[1])+ theme(axis.title.x = element_text(size=12)) 
    p2=upset(list_up[[2]][,6:ncol(list_up[[1]])],colnames(list_up[[2]][,6:ncol(list_up[[1]])]),min_size=10,width_ratio=.3,
             set_sizes=(upset_set_size() + theme(axis.text.x=element_text(angle=45,hjust=1,size=8))),intersections="all",
             base_annotations=list('Size'=(intersection_size(counts=FALSE)))) + xlab(names(list_up)[2]) + theme(axis.title.x = element_text(size=12))
    pfinal=ggarrange(p1,p2,ncol=2)
  } else if (length(list_up)==3) {
    p1=upset(list_up[[1]][,6:ncol(list_up[[1]])],colnames(list_up[[1]][,6:ncol(list_up[[1]])]),min_size=10,width_ratio=.3,
             set_sizes=(upset_set_size() + theme(axis.text.x=element_text(angle=45,hjust=1,size=8))),intersections="all",
             base_annotations=list('Size'=(intersection_size(counts=FALSE))))  + xlab(names(list_up)[1]) + theme(axis.title.x = element_text(size=12))
    p2=upset(list_up[[2]][,6:ncol(list_up[[1]])],colnames(list_up[[2]][,6:ncol(list_up[[1]])]),min_size=10,width_ratio=.3,
             set_sizes=(upset_set_size() + theme(axis.text.x=element_text(angle=45,hjust=1,size=8))),intersections="all",
             base_annotations=list('Size'=(intersection_size(counts=FALSE)))) + xlab(names(list_up)[2]) + theme(axis.title.x = element_text(size=12))
    p3=upset(list_up[[3]][,6:ncol(list_up[[1]])],colnames(list_up[[3]][,6:ncol(list_up[[1]])]),min_size=10,width_ratio=.3,
             set_sizes=(upset_set_size() + theme(axis.text.x=element_text(angle=45,hjust=1,size=8))),intersections="all",
             base_annotations=list('Size'=(intersection_size(counts=FALSE)))) + xlab(names(list_up)[3]) + theme(axis.title.x = element_text(size=12))
    pfinal=ggarrange(p1,p2,p3,ncol=3)
  } else {
    warning("Upset plot function implemented in this vignette is not implemented to work on more than 3 modalities")
  }
  return(pfinal)
}


plotPassed <- function(mdata.list,xaxis_text=NULL,angle_x=NULL) {
  if (is.null(xaxis_text)) { xaxis_text=12 }
  if (is.null(angle_x)) { angle_x=30 }
  i <- 0
  for (exp in names(mdata.list)) {
    i <- i + 1
    df <- mdata.list[[exp]]
    counts <- as.data.frame(table(df$passedMB))
    counts$Perc <- counts$Freq / sum(counts$Freq) * 100
    counts$sample <- exp
    # bind
    if (i==1) {
      metadata.df <- counts
    } else {
      metadata.df <- rbind(metadata.df,counts)
    }
  }
  names(metadata.df)[1] <- "passedMB"
  # now plot
  pl=ggplot(metadata.df,aes(x=sample,y=Freq,fill=passedMB)) +
    theme_bw() +
    geom_bar(stat="identity",color='black',alpha=.7) +
    # x-axis
    theme(axis.text.x=element_text(angle=angle_x, hjust=1, size = xaxis_text)) + 
    theme(axis.title.x = element_text(size= 0)) +
    xlab("") +
    theme(axis.title.x = element_text(size= 14)) +
    # y-axis
    ylab("Number of cells") +
    theme(axis.text.y=element_text(size = 12)) +
    theme(axis.title.y = element_text(size= 14)) 
  return(pl)
}


plotPassedCells <- function(obj,sample_list,mod_list) {
  # convert obj to df
  i <- 0
  for (nm in names(obj)) {
    i <- i + 1
    # get samples and barcode name
    for (s in sample_list) {
      if (grepl(s,nm)) { 
        sample <- s
      }
    }
    for (m in mod_list) {
      if (grepl(m,nm)) { 
        mod <- m
      }
    }
    obj[[nm]]$sample <- sample
    obj[[nm]]$modality <- mod
    if (i==1) {
      df <- obj[[nm]]
    } else {
      df <- rbind(df,obj[[nm]])
    }
  }
  plot <- ggplot(df,aes(x=all_unique_MB,y=peak_ratio_MB,fill=passedMB)) + 
    theme_bw() +
    geom_point(shape=21,size=.4) +
    scale_x_log10(labels=trans_format('log10',math_format(10^.x))) +
    #coord_cartesian(ylim = c(0,1),xlim = c(10,1000000)) +
    facet_grid(modality~sample) +
    theme(strip.text.x = element_text(size = 11, colour = "black", angle = 0, face= 'bold')) +
    theme(strip.text.y = element_text(size = 11, colour = "black", face= 'bold')) + 
    scale_fill_manual(values = c("#F8766D","#00BA38")) +
    # x-axis
    ylab("Fractio reads in peaks") +
    theme(axis.text.x=element_text(angle=0, hjust=.5, size = 12)) + 
    theme(axis.title.x = element_text(size= 14)) +
    # y-axis
    xlab("UMI") +
    theme(axis.text.y=element_text(size = 12)) +
    theme(axis.title.y = element_text(size= 14))+
    guides(fill = guide_legend(override.aes = list(size=5)))
  plot <- ggarrange(plot,legend='bottom')
  return(plot)
}




commonCellHistonMarks <- function(mod1,name_mod1,mod2,name_mod2,mod3=NULL,name_mod3=NULL,mod4=NULL,name_mod4=NULL,sample) {
  
  # for 2 modalities
  x <- list(name_mod1=mod1$barcode_RNA,
            name_mod2=mod2$barcode_RNA)
  names(x) <- c(name_mod1,name_mod2)
  pp=ggvenn(x) + 
    #scale_fill_gradient(low = "white", high = "white") +
    theme(legend.position = "none") +
    ggtitle(sample) +
    theme(plot.title = element_text(size=14,hjust=0.5,face='bold'))
  
  # 3 modalities
  if (!is.null(mod3)) {
    x <- list(name_mod1=mod1$barcode_RNA,
              name_mod2=mod2$barcode_RNA,
              name_mod3=mod3$barcode_RNA)
    names(x) <- c(name_mod1,name_mod2,name_mod3)
    pp=ggvenn(x) + 
      #scale_fill_gradient(low = "white", high = "white") +
      theme(legend.position = "none") +
      ggtitle(sample) +
      theme(plot.title = element_text(size=14,hjust=0.5,face='bold'))
  }
  
  # 4 modalities
  if (!is.null(mod4)) {
    x <- list(name_mod1=mod1$barcode_RNA,
              name_mod2=mod2$barcode_RNA,
              name_mod3=mod3$barcode_RNA,
              name_mod4=mod4$barcode_RNA)
    names(x) <- c(name_mod1,name_mod2,name_mod3,name_mod4)
    pp=ggvenn(x) + 
      #scale_fill_gradient(low = "white", high = "white") +
      theme(legend.position = "none") +
      ggtitle(sample) +
      theme(plot.title = element_text(size=14,hjust=0.5,face='bold'))
  }
  return(pp)
  
}

commonCellHistonMarks_all <- function(mod1,name_mod1,mod2,name_mod2,mod3=NULL,name_mod3=NULL,mod4=NULL,name_mod4=NULL) {
  
  # for 2 modalities
  x <- list(name_mod1=colnames(mod1),
            name_mod2=colnames(mod2))
  names(x) <- c(name_mod1,name_mod2)
  pp=ggvenn(x) + 
    #scale_fill_gradient(low = "white", high = "white") +
    theme(legend.position = "none") +
    ggtitle(sample) +
    theme(plot.title = element_text(size=14,hjust=0.5,face='bold'))
  
  # 3 modalities
  if (!is.null(mod3)) {
    x <- list(name_mod1=colnames(mod1),
              name_mod2=colnames(mod2),
              name_mod3=colnames(mod3))
    names(x) <- c(name_mod1,name_mod2,name_mod3)
    pp=ggvenn(x) + 
      #scale_fill_gradient(low = "white", high = "white") +
      theme(legend.position = "none") +
      ggtitle(sample) +
      theme(plot.title = element_text(size=14,hjust=0.5,face='bold'))
  }
  
  # 4 modalities
  if (!is.null(mod4)) {
    x <- list(name_mod1=colnames(mod1),
              name_mod2=colnames(mod2),
              name_mod3=colnames(mod3),
              name_mod4=colnames(mod4))
    names(x) <- c(name_mod1,name_mod2,name_mod3,name_mod4)
    pp=ggvenn(x) + 
      #scale_fill_gradient(low = "white", high = "white") +
      theme(legend.position = "none") +
      ggtitle(sample) +
      theme(plot.title = element_text(size=14,hjust=0.5,face='bold'))
  }
  return(pp)
  
}

splitGroupFragments <- function(
    object = NULL,
    assay = "ATAC",
    genome = "mm10",
    groupBy = "Sample",
    downGroupBy = "all",
    minCells = 5,
    maxCells = 1000,
    nbFrags = 10000000,
    threads=NULL,
    test=FALSE,
    test.nbFrags=NULL,
    outdir = NULL,
    tmpdir = NULL
){
  
  dir.create(file.path(outdir), showWarnings = FALSE)
  dir.create(tmpdir, showWarnings = FALSE)
  
  DefaultAssay(object) <- assay
  genome(object) <- genome
  
  #Get fragments file path
  if (length(Fragments(object)) > 1){
    print("The object contains a list of fragments files")
    frag_file <- unique(unlist(lapply(Fragments(object), function(x){GetFragmentData(x)})))
    if (frag_file > 1){
      #Could also try to run on each specific fragments files without splitting
      print("The list of fragments files are pointing to different fragment files")
      print("Please, aggregate the fragments files before submitting")
    }else{
      bed_input_path <- frag_file
    }
  }else{
    bed_input_path <- GetFragmentData(Fragments(object))
  }
  
  if(test){
    dir.create(paste0(tmpdir, "/test"), showWarnings = FALSE)
    nbFrags <- test.nbFrags
    print(paste0("Testing on a small chunk of the fragments file, selecting : ",nbFrags, " fragments"))
    system(sprintf("zcat %s | head -%d > %s/test/test_%s_atac_fragments.tsv", bed_input_path, nbFrags, tmpdir, as.character(nbFrags)))
    system(sprintf("gzip %s/test/test_%s_atac_fragments.tsv", tmpdir, as.character(nbFrags)))
    bed_input_path <- paste0(tmpdir, "/test/test_", nbFrags, "_atac_fragments.tsv.gz")
    print(paste0("Testing file located at : ",bed_input_path))
  }
  
  annot <- object@meta.data
  
  #Column to group the cells by
  Groups <- annot[, groupBy, drop=FALSE]
  
  if(unlist(downGroupBy)[1] != 'all'){
    Groups <- Groups[which(Reduce(`&`, Map(`%in%`, Groups[intersect(names(Groups), names(downGroupBy))], downGroupBy))),]
  }
  
  #Split cell names by group
  cellGroups <- lapply(split(rownames(Groups), factor(apply(Groups, 1, paste, collapse="_"))), unique)
  
  #Sample each cell group to maxCells
  if(!is.null(maxCells)){
    cellGroups <- lapply(cellGroups, function(x){if(length(x) <= maxCells){x}else{sample(x, maxCells)}})
  }
  
  #Remove group with less than maxCells
  if(!is.null(minCells)){
    cellGroups <- Filter(function(x) length(x) >= minCells, cellGroups)
  }
  
  #Create a vector of selected cells with the group as name
  cellGroups_rev_list <- unlist(cellGroups)
  names(cellGroups_rev_list) <- rep(names(cellGroups), sapply(cellGroups, length))
  
  ##If there is 1B line in the fragment file
  ##We can split the fragments file by the number of available threads (minus 1 thread not to overload)
  
  if (threads == 1){
    avai_threads <- 1
  }else{
    avai_threads <- threads - 1
  }
  theorical_nbfile_per_thread <- nbFrags/avai_threads
  suitable_nbfile_per_thread <- ceiling(theorical_nbfile_per_thread)
  print(paste0("Each sub-fragments files will contain : ", suitable_nbfile_per_thread," fragments"))
  system(sprintf("gunzip -c %s | split --lines %d --additional-suffix=.tsv - %s/sub_atac_fragments", bed_input_path, suitable_nbfile_per_thread, tmpdir))
  
  #Set number of thread in future
  plan("multicore", workers = avai_threads)
  
  files <- list.files(tmpdir)[list.files(tmpdir) %like% "sub_atac_fragments"]    
  for(group in unique(names(cellGroups_rev_list))){
    cells <- cellGroups_rev_list[names(cellGroups_rev_list)==group]
    future_lapply(files,.groupfragments, cells, group, tmpdir, future.stdout=FALSE)
    system(sprintf("cat %s/%s* | sort -k1,1 -k2,2n -k3,3n - > %s/%s.tsv", tmpdir, group, tmpdir, group))
    system(sprintf("rm %s/%s_sub*", tmpdir, group))
    system(sprintf("gzip %s/%s.tsv", tmpdir, group))
  }
  
  #Removing sub-fragments files
  system(sprintf("rm %s/sub_atac_fragments*", tmpdir))
  
  plan("multicore", workers = 1)
  
  if(test){
    system(sprintf("rm -r %s/test", tmpdir))
  }
}

#Function to generate fragments files
.groupfragments <- function(x, cellname, cellgroup, tmpdir) {
  bed_input <- gzfile(paste0(tmpdir,"/",x), "r")
  bed_output <- file(paste0(tmpdir, "/", cellgroup, "_", strsplit(x,"\\.")[[1]][1], ".tsv"), "w")
  
  while(length(line <- readLines(bed_input, n = 1)) > 0) {
    cell_name_bed <- strsplit(line, "\t")[[1]][4]
    if (cell_name_bed %in% cellname) {
      writeLines(line, bed_output)
    }
  }
  close(bed_input)
  close(bed_output)
}

ExportGroupBW  <- function(
    object,
    assay = NULL,
    group.by = NULL,
    idents = NULL,
    normMethod = "RC",
    tileSize = 100,
    minCells = 5,
    cutoff = NULL,
    chromosome = NULL,
    outdir = NULL,
    verbose=TRUE
) {
  # Check if temporary directory exist
  if (!dir.exists(outdir)){
    dir.create(outdir)
  }
  if (!requireNamespace("rtracklayer", quietly = TRUE)) { 
    message("Please install rtracklayer. http://www.bioconductor.org/packages/rtracklayer/") 
    return(NULL) 
  }
  assay <- SetIfNull(x = assay, y = DefaultAssay(object = object))
  DefaultAssay(object = object) <- assay
  group.by <- SetIfNull(x = group.by, y = 'ident')
  Idents(object = object) <- group.by
  idents <- SetIfNull(x = idents, y = levels(x = object))
  GroupsNames <- names(x = table(object[[group.by]])[table(object[[group.by]]) > minCells])
  GroupsNames <- GroupsNames[GroupsNames %in% idents]
  # Check if output files already exist
  lapply(X = GroupsNames, FUN = function(x) {
    fn <- paste0(outdir, .Platform$file.sep, x, ".bed")
    if (file.exists(fn)) {
      message(sprintf("The group \"%s\" is already present in the destination folder and will be overwritten !",x))
      file.remove(fn)
    }
  })      
  # Splitting fragments file for each idents in group.by
  SplitFragments(
    object = object,
    assay = assay,
    group.by = group.by,
    idents = idents,
    outdir = outdir,
    file.suffix = "",
    append = TRUE,
    buffer_length = 256L,
    verbose = verbose
  )
  # Column to normalized by
  if(!is.null(x = normMethod)) {
    if (tolower(x = normMethod) %in% c('rc', 'ncells', 'none')){
      normBy <- normMethod
    } else{
      normBy <- object[[normMethod, drop = FALSE]]
    }
  }
  # Get chromosome information
  if(!is.null(x = chromosome)){
    seqlevels(object) <- chromosome
  }
  availableChr <- names(x = seqlengths(object))
  chromLengths <- seqlengths(object)
  chromSizes <- GRanges(
    seqnames = availableChr,
    ranges = IRanges(
      start = rep(1, length(x = availableChr)),
      end = as.numeric(x = chromLengths)
    )
  )
  if (verbose) {
    message("Creating tiles")
  }
  # Create tiles for each chromosome, from GenomicRanges
  tiles <- unlist(
    x = slidingWindows(x = chromSizes, width = tileSize, step = tileSize)
  )
  if (verbose) {
    message("Creating bigwig files at ", outdir)
  }
  # Run the creation of bigwig for each cellgroups
  if (nbrOfWorkers() > 1) { 
    mylapply <- future_lapply 
  } else { 
    mylapply <- ifelse(test = verbose, yes = pblapply, no = lapply) 
  }
  
  covFiles <- mylapply(
    GroupsNames,
    FUN = CreateBWGroup,
    availableChr,
    chromLengths,
    tiles,
    normBy,
    tileSize,
    normMethod,
    cutoff,
    outdir
  )
  return(covFiles)
}

CreateBWGroup <- function(
    groupNamei,
    availableChr,
    chromLengths,
    tiles,
    normBy,
    tileSize,
    normMethod,
    cutoff,
    outdir
) {
  if (!requireNamespace("rtracklayer", quietly = TRUE)) { 
    message("Please install rtracklayer. http://www.bioconductor.org/packages/rtracklayer/") 
    return(NULL) 
  }
  normMethod <- tolower(x = normMethod)
  # Read the fragments file associated to the group
  fragi <- rtracklayer::import(
    paste0(outdir, .Platform$file.sep, groupNamei, ".bed"), format = "bed"
  )
  cellGroupi <- unique(x = fragi$name)
  # Open the writing bigwig file
  covFile <- file.path(
    outdir,
    paste0(groupNamei, "-TileSize-",tileSize,"-normMethod-",normMethod,".bw")
  )
  
  covList <- lapply(X = seq_along(availableChr), FUN = function(k) {
    fragik <- fragi[seqnames(fragi) == availableChr[k],]
    tilesk <- tiles[BiocGenerics::which(S4Vectors::match(seqnames(tiles), availableChr[k], nomatch = 0) > 0)]
    if (length(x = fragik) == 0) {
      tilesk$reads <- 0
      # If fragments
    } else {
      # N Tiles
      nTiles <- chromLengths[availableChr[k]] / tileSize
      # Add one tile if there is extra bases
      if (nTiles%%1 != 0) {
        nTiles <- trunc(x = nTiles) + 1
      }
      # Create Sparse Matrix
      matchID <- S4Vectors::match(mcols(fragik)$name, cellGroupi)
      
      # For each tiles of this chromosome, create start tile and end tile row,
      # set the associated counts matching with the fragments
      mat <- Matrix::sparseMatrix(
        i = c(trunc(x = start(x = fragik) / tileSize),
              trunc(x = end(x = fragik) / tileSize)) + 1,
        j = as.vector(x = c(matchID, matchID)),
        x = rep(1, 2*length(x = fragik)),
        dims = c(nTiles, length(x = cellGroupi))
      )
      
      # Max count for a cells in a tile is set to cutoff
      if (!is.null(x = cutoff)){
        mat@x[mat@x > cutoff] <- cutoff
      }
      # Sums the cells
      mat <- rowSums(x = mat)
      tilesk$reads <- mat
      # Normalization
      if (!is.null(x = normMethod)) {
        if (normMethod == "rc") {
          tilesk$reads <- tilesk$reads * 10^4 / length(fragi$name)
        } else if (normMethod == "ncells") {
          tilesk$reads <- tilesk$reads / length(cellGroupi)
        } else if (normMethod == "none") {
        } else {
          if (!is.null(x = normBy)){
            tilesk$reads <- tilesk$reads * 10^4 / sum(normBy[cellGroupi, 1])
          }
        }
      }
    }
    tilesk <- coverage(tilesk, weight = tilesk$reads)[[availableChr[k]]]
    tilesk
  })
  
  names(covList) <- availableChr
  covList <- as(object = covList, Class = "RleList")
  rtracklayer::export.bw(object = covList, con = covFile)
  return(covFile)
}

SetIfNull <- function(x, y) {
  if (is.null(x = x)) {
    return(y)
  } else {
    return(x)
  }
}