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
	#############################################################################################################################
	###                                          1.LOAD THE 111 FILES FROM CHROMHMM                                           ###
	#############################################################################################################################

	# Get all files from annotation directory
	if (verbose){print("Load ChromHMM files...")}
	files <- list.files(annotation_dir, pattern = "\\.bed", full.names = T)

	# Load each file in granges and name them
	encode <- lapply(files, read_bed)
	names(encode) <- sapply(files, function(x) {strsplit(basename(x), "_")[[1]][1]})

	# Load tissue metadata file
	glossary <- as.data.frame(read_delim(glossary_path, delim = "\t", escape_double = FALSE, trim_ws = TRUE))
	# Remove cell line files
	glossary <- glossary[glossary$GROUP != "ENCODE2012",]

	#############################################################################################################################
	###                                    2.MERGE ADJACENT RANGES BY TISSUE AND ANNOTATE                                     ###
	#############################################################################################################################

	# Create a grange for each tissue and reduce it
	if (verbose){print("Reducing ChromHMM windows per group...")}
	windows_list <- list()
	for (tissue in unique(glossary[,glossary_group])){
		eid_tissue <- glossary[glossary[,glossary_group] == tissue,]$EID
		tmp <- bind_ranges(encode[eid_tissue])
		tmp <- tmp[tmp$name %in% states, ]
		windows_list[[tissue]] <- reduce(tmp, ignore.strand=TRUE)
		windows_list[[tissue]]$tissue <- tissue
	}

	# Aggregate the granges
	if (verbose){print("Aggregate ChromHMM windows...")}
	windows_aggr <- sort(sortSeqlevels(bind_ranges(windows_list)), ignore.strand=TRUE)

	# Disjoin to get windows overlapping multiple tissues (can take a while)
	windows <- disjoin(windows_aggr, with.revmap=TRUE, ignore.strand=TRUE)
	windows$tissue <- unlist(lapply(windows$revmap,function(i){paste(collapse=';',windows_aggr$tissue[i])}))
	windows$revmap <- NULL

	return(windows)
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
	windows$feature <- "UNKNOWN"
	for (db_name in names(db)){

        if (verbose){print(paste0("Annotation from ", db_name,"..."))}

		windows <- join_overlap_left(windows, db_list[[db_name]], maxgap = db_gap)
		
		#Annotation UCSC
		if(toupper(db_name)=="UCSC"){

			windows[windows$feature == "UNKNOWN",]$feature <- windows[windows$feature == "UNKNOWN",]$SYMBOL
			windows[is.na(windows$feature)]$feature <- "UNKNOWN"
			windows$tx_name <- NULL
			windows$SYMBOL <- NULL
		}
		#Annotation Ensembl
		else if(toupper(db_name)=="ENSEMBL"){
			windows[windows$feature == "UNKNOWN",]$feature <- windows[windows$feature == "UNKNOWN",]$tx_id
			windows[is.na(windows$feature)]$feature <- "UNKNOWN"
			windows$tx_id <- NULL
		}
	}

	# Group the ranges for the same gene names
	if (verbose){print("Group ranges...")}
	windows <- windows %>% 
		group_by(seqnames, start, end, tissue) %>%
		summarise(feature = paste(unique(feature), collapse = ";")) %>%
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
				db_specific$feature <- unlist(lapply(seq(1:length(db_specific)),function(i){
					ids_hit <- as.integer(names(db_hits[db_hits==i]))
					paste(collapse=';',unique(tmp_db_extend$SYMBOL[ids_hit]))
					}))
			}
			#Annotation Ensembl
			else if(toupper(db_name)=="ENSEMBL"){
				db_specific$feature <- unlist(lapply(seq(1:length(db_specific)),function(i){
					ids_hit <- as.integer(names(db_hits[db_hits==i]))
					paste(collapse=';',unique(tmp_db_extend$tx_id[ids_hit]))
					}))
			}

			# Merge the newly annotated regions
			windows <- sort(sortSeqlevels(bind_ranges(list(windows,db_specific))), ignore.strand=TRUE)
		}
	}

	return(windows)
}

# Add flanking regions to the catalog
#	@windows			: Windows to annotate
#	@feature_tag		: Tag name of the features of interest
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
	
	# Reduce flanking regions (max flanking regions is then flanking_size*2)
	windows_red_flank <- reduce(windows_red_flank)
	
	# Find flank specific location
	flank_specific <- disjoin_specific(windows_red, windows_red_flank, feature_tag)

	# Merge the newly annotated regions
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
	#bg_genome <- tileGenome(seqlengths(windows), tilewidth=bin_size, cut.last.tile.in.chrom=TRUE)
	
	# Set whole genome as background
	bg_genome <- GRanges(names(seqlengths(windows)), IRanges(rep(1,length(seqlengths(windows))), seqlengths(windows)))

	# Reduce our windows to check which regions are already covered
	windows_red <- sort(sortSeqlevels(reduce(windows)), ignore.strand=TRUE)
	
	# Find background specific location
	bg_specific <- disjoin_specific(windows_red, bg_genome, feature_tag)
	
	# Split background regions by max=bin_size
	bg_specific_split <- split_homogeneous(bg_specific, bin_size, feature_tag)
	
	# Merge the newly annotated regions
	windows <- sort(sortSeqlevels(bind_ranges(list(windows,bg_specific_split))), ignore.strand=TRUE)

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

# Sub function to extract ranges specific to the second grange entry
#	@gr1				: Grange1 
#	@gr2				: Grange2
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

# Sub function to split ranges by nearly equal size
#	@gr				: Grange 
#	@bin_size		: Maximum size of each background tile
#	@tag    : Tag to use to reannotate the type
split_homogeneous <- function(
	gr					= NULL,
	bin_size			= NULL,
	tag					= NULL
){

	# Find how many subdivision each bin will have
	gr$nb <- ceiling(width(gr)/bin_size)
	
	
	############TRY SOMETHING BY SPLIT SUBDIVISION
	# Create a list by division
	#gr.list <- split(gr, gr$nb)
	# Subdivide Granges by each associated number of division
	#lapply(gr.list, function(x){subdivideGRanges(x, subsize=as.integer(unique(x$nb)))})
	############
	
	
	gr_1bin <- gr[gr$nb==1]
	gr_1bin$nb <- NULL
	gr <- gr[gr$nb>1]
	
	gr_split <- GRanges()
	int_split <- function(n, p) n %/% p + (sequence(p) - 1 < n %% p)
	for (i in 1:length(gr)){

		bin_sizes <- int_split(width(gr[i]),gr[i]$nb)

		incremental_ends <- c()
		for(j in 1:length(bin_sizes)){
			incremental_ends <- c(incremental_ends, sum(bin_sizes[1:j]))
		}
		new_ends <- incremental_ends+start(gr[i])-1
		new_starts <- new_ends-bin_sizes+1
		
		gr_split <- append(gr_split, GRanges(rep(as.character(seqnames(gr[i])),gr[i]$nb),IRanges(start = new_starts, end = new_ends)))
		
	}
	
	gr_split$type = tag
	
	gr <- sort(append(gr_split, gr_1bin))

	return(gr)
}

