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
	if (verbose){print("Reducing ChromHMM windows_data per group...")}
	windows_data_list <- list()
	for (tissue in unique(glossary[,glossary_group])){
		eid_tissue <- glossary[glossary[,glossary_group] == tissue,]$EID
		tmp <- bind_ranges(encode[eid_tissue])
		tmp <- tmp[tmp$name %in% states, ]
		windows_data_list[[tissue]] <- reduce(tmp, ignore.strand=TRUE)
		windows_data_list[[tissue]]$tissue <- tissue
	}

	# Aggregate the granges
	if (verbose){print("Aggregate ChromHMM windows_data...")}
	windows_data_aggr <- sort(sortSeqlevels(bind_ranges(windows_data_list)), ignore.strand=TRUE)

	# Disjoin to get windows_data overlapping multiple tissues (can take a while)
	windows_data <- disjoin(windows_data_aggr, with.revmap=TRUE, ignore.strand=TRUE)
	windows_data$tissue <- unlist(lapply(windows_data$revmap,function(i){paste(collapse=';',windows_data_aggr$tissue[i])}))
	windows_data$revmap <- NULL

	return(windows_data)
}


# Annotation of the catalog based on databases
#	@windows_data			: windows_data to annotate
#	@feature_tag		: Tag name of the features of interest
#	@db					: List of database annotation to retrieve information from, the first takes priority on the second etc...
#	@db_gap				: Tolerated gap to annotate the windows_data
#	@db_spe				: Boolean to include or not annotation specific features
#	@db_tag				: Tag name of the annotation specific features
#	@db_spe_center		: Size of the windows_data centered on the feature
#	@verbose			: Display verbose
annotate_catalog <- function(
	windows_data				= NULL,
	feature_tag			= NULL,
	db					= NULL,
	db_gap				= NULL,
	db_spe				= NULL,
	db_tag				= NULL,
	db_spe_center		= NULL,
	verbose				= TRUE
){
	#Annotation of the windows_data
	if (verbose){print("windows_data annotation from the databases...")}

	db_list <- list()
	for (db_name in names(db)){
		# Gather features annotation from UCSC
		if(toupper(db_name)=="UCSC"){
		  ucscdb <- promoters(db[[db_name]], upstream=0, downstream=1, columns=c("tx_name", "gene_id"))
		  ucscdb <- sort(sortSeqlevels(ucscdb), ignore.strand=TRUE)
		  mcols(ucscdb)$gene_id <- as.character(mcols(ucscdb)$gene_id)
		  # First, ensure gene_id is properly formatted
		  mcols(ucscdb)$gene_id <- as.character(mcols(ucscdb)$gene_id)
		  
		  # Check if any gene_ids are valid Entrez IDs
		  valid_entrez_keys <- keys(org.Hs.eg.db, keytype="ENTREZID")
		  overlap <- intersect(mcols(ucscdb)$gene_id, valid_entrez_keys)
		  print(paste("Number of matching Entrez keys:", length(overlap)))
		  
		  # If you have valid Entrez IDs, use this:
		  if(length(overlap) > 0) {
		    # Filter to keep only valid gene IDs
		    ucscdb <- ucscdb[mcols(ucscdb)$gene_id %in% valid_entrez_keys]
		    
		    # Now get the symbols
		    entrezid2symbol <- AnnotationDbi::select(org.Hs.eg.db, 
		                                             keys = mcols(ucscdb)$gene_id,
		                                             columns = c("SYMBOL"),
		                                             keytype = "ENTREZID")
		    
		    # Add symbols to the ucscdb object
		    mcols(ucscdb)$SYMBOL <- entrezid2symbol$SYMBOL[match(mcols(ucscdb)$gene_id, entrezid2symbol$ENTREZID)]
		  }
		  
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
	windows_data$feature <- "UNKNOWN"
	for (db_name in names(db)){

        if (verbose){print(paste0("Annotation from ", db_name,"..."))}

		windows_data <- join_overlap_left(windows_data, db_list[[db_name]], maxgap = db_gap)
		
		#Annotation UCSC
		if(toupper(db_name)=="UCSC"){

			windows_data[windows_data$feature == "UNKNOWN",]$feature <- windows_data[windows_data$feature == "UNKNOWN",]$SYMBOL
			windows_data[is.na(windows_data$feature)]$feature <- "UNKNOWN"
			windows_data$tx_name <- NULL
			windows_data$SYMBOL <- NULL
		}
		#Annotation Ensembl
		else if(toupper(db_name)=="ENSEMBL"){
			windows_data[windows_data$feature == "UNKNOWN",]$feature <- windows_data[windows_data$feature == "UNKNOWN",]$tx_id
			windows_data[is.na(windows_data$feature)]$feature <- "UNKNOWN"
			windows_data$tx_id <- NULL
		}
	}

	# Group the ranges for the same gene names
	if (verbose){print("Group ranges...")}
	windows_data <- windows_data %>% 
		group_by(seqnames, start, end, tissue) %>%
		summarise(feature = paste(unique(feature), collapse = ";")) %>%
		makeGRangesFromDataFrame(keep.extra.columns = T)

	windows_data$type <- feature_tag

	if (db_spe){

	#############################################################################################################################
	###                                       3.GET REGIONS IN DATABASES NOT IN windows_data                                       ###
	#############################################################################################################################

		if (verbose){print("Include annotation missing from catalogue")}

		for (db_name in names(db)){

        	if (verbose){print(paste0("Extra annotation from ", db_name,"..."))}

			if (verbose){print(paste0("Find new regions..."))}

			# Reduce our windows_data to check which regions are already covered
			windows_data_covered <- reduce(windows_data)

			tmp_db <- db_list[[db_name]]

			# Extended positions covered by the database
			tmp_db_extend <- resize(tmp_db, width = width(tmp_db)+db_spe_center, fix = "center")
			tmp_db_extend <- tmp_db_extend[seqnames(tmp_db_extend) %in% seqlevels(windows_data)]
			seqlevels(tmp_db_extend) <- seqlevels(windows_data)
			tmp_db_extend <- resize_correct(tmp_db_extend)

			tmp_db_extend_red <- reduce(tmp_db_extend)
			tmp_db_extend_red <- sort(sortSeqlevels(tmp_db_extend_red), ignore.strand=TRUE)

			# Genome are not exactly the same but positions are identicial except some on chrM
			genome(tmp_db_extend_red) <- NA

			#Find databases specific location
			db_specific <- disjoin_specific(windows_data_covered,tmp_db_extend_red,db_tag)

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
			windows_data <- sort(sortSeqlevels(bind_ranges(list(windows_data,db_specific))), ignore.strand=TRUE)
		}
	}

	return(windows_data)
}

# Add flanking regions to the catalog
#	@windows_data			: windows_data to annotate
#	@feature_tag		: Tag name of the features of interest
#	@flanking_size		: Size of the flanking regions around the features
#	@verbose			: Display verbose
add_flank <- function(
	windows_data				= NULL,
	feature_tag			= NULL,
	flanking_size		= NULL,
	verbose				= TRUE
){
	#############################################################################################################################
	###                                                  4.SET UP FLANKING REGIONS                                            ###
	#############################################################################################################################

	if (verbose){print(paste0("Add flanking regions to features..."))}

	# Reduce our windows_data to check which regions are already covered
	windows_data_red <- sort(sortSeqlevels(reduce(windows_data)), ignore.strand=TRUE)
	
	# Extending both side of all features
	windows_data_red_flank <- resize(windows_data_red, width = width(windows_data_red)+(flanking_size*2), fix = "center")
	windows_data_red_flank <- resize_correct(windows_data_red_flank)
	
	# Reduce flanking regions (max flanking regions is then flanking_size*2)
	windows_data_red_flank <- reduce(windows_data_red_flank)
	
	# Find flank specific location
	flank_specific <- disjoin_specific(windows_data_red, windows_data_red_flank, feature_tag)

	# Merge the newly annotated regions
	windows_data <- sort(sortSeqlevels(bind_ranges(list(windows_data,flank_specific))), ignore.strand=TRUE)

	return(windows_data)
}

# Add background to the catalog
#	@windows_data			: windows_data to annotate
#	@feature_tag		: Tag name of the features of interest
#	@bin_size			: Size of the tile of the background
#	@verbose			: Display verbose
add_background <- function(
	windows_data				= NULL,
	feature_tag			= NULL,
	bin_size			= NULL,
	verbose				= TRUE
){
	#############################################################################################################################
	###                                                 4.SET UP BACKGROUND REGIONS                                           ###
	#############################################################################################################################

	if (verbose){print(paste0("Add background to the catalog..."))}

	# Tile the whole genome
	#bg_genome <- tileGenome(seqlengths(windows_data), tilewidth=bin_size, cut.last.tile.in.chrom=TRUE)
	
	# Set whole genome as background
	bg_genome <- GRanges(names(seqlengths(windows_data)), IRanges(rep(1,length(seqlengths(windows_data))), seqlengths(windows_data)))

	# Reduce our windows_data to check which regions are already covered
	windows_data_red <- sort(sortSeqlevels(reduce(windows_data)), ignore.strand=TRUE)
	
	# Find background specific location
	bg_specific <- disjoin_specific(windows_data_red, bg_genome, feature_tag)
	
	# Split background regions by max=bin_size
	bg_specific_split <- split_homogeneous(bg_specific, bin_size, feature_tag)
	
	# Merge the newly annotated regions
	windows_data <- sort(sortSeqlevels(bind_ranges(list(windows_data,bg_specific_split))), ignore.strand=TRUE)

	return(windows_data)
}

# Statistics about the window
#	@windows_data			: windows_data to annotate
#	@db_tag				: Tag name of the annotation specific features
show_statistics <- function(
	windows_data				= NULL,
	db_tag				= NULL
){
	print(table(windows_data$type))

	tmp <- windows_data[windows_data$type == db_tag,]
	nb_extra_ucsc <- length(tmp[substring(tmp$feature,1,4) != "ENST",])
	nb_extra_ensembl <- length(tmp[substring(tmp$feature,1,4) == "ENST",])
	print(paste0("UCSC : ",nb_extra_ucsc))
	print(paste0("Ensembl : ",nb_extra_ensembl))

	size_genome <- sum(seqlengths(windows_data))
	print("Genome coverage")
	for (type_reg in unique(windows_data$type)){
		tmp <- windows_data[windows_data$type==type_reg,]
		print(paste0(type_reg, " : ", round(sum(width(tmp))*100/size_genome, 2)))
	}

	print("Quantile values")
	for (type_reg in unique(windows_data$type)){
		tmp <- windows_data[windows_data$type==type_reg,]
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

# Process BED file and a set of regions.
#	@bed_file			: Bed file to process 
#	@regions			: Granges of regions
#	@verbose			: Display verbose
process_bed_file <- function(
	bed_file			= NULL,
	regions				= NULL,
	verbose				= TRUE
) {

	if (verbose){print(paste0("Import bed files..."))}
	
	# Load the BED file
	bed <- import(bed_file, format="bed")

	# Resize the fragments to length 1, maintaining their center
	bed_resized <- resize(bed, width = 1, fix = "center")
	
	if (verbose){print(paste0("Count overlaps..."))}

	# Count overlaps
	counts <- countOverlaps(regions, bed_resized)

	# Add counts to the metadata columns of the regions GRanges object
	mcols(regions)$counts <- counts

	# Return the regions GRanges object
	return(regions)
}

# Estimate background for a specific histone mark
#	@gr					: Grange of counts
#	@tag				: Tag to use for background
#	@background_cutOff	: All background region bigger than x will be kept to estimate the background signal
#	@quantile_cutOff	: Quantile cut off to remove outlier
#	@verbose			: Display verbose
estimateBackground <- function(
	gr					= NULL,
	tag					= NULL,
	background_cutOff   = 4000,
	quantile_cutOff		= 0.95,
	verbose				= TRUE
) {

	if (verbose){print(paste0("Estimate background..."))}

	# Extract background regions
	background_regions <- gr[gr$type == tag]
	
	# Filter regions that are longer than background_cutOff
	long_background_regions <- background_regions[width(background_regions) >= background_cutOff]

	nb_kept <- round((length(long_background_regions)*100)/length(background_regions),2)

	if (verbose){print(paste0("Background windows_data kept ", nb_kept, "%"))}
	if (verbose){print(paste0("Median of the windows_data is ", median(width(long_background_regions)), "bp"))}
	
	# Extract counts from the filtered regions
	X <- mcols(long_background_regions)$counts
	
	# Find the 95th quantile of X
	T <- quantile(X, quantile_cutOff)
	
	# Restrict ourselves to values below T
	X <- X[X <= T]
	
	# Maximum likelihood of truncated poisson
	lambda_hat <- which.max(dpois(X, lambda = 1:length(X)))
	
	# Convert to reads/Kb (the median length of background windows_data is 5KB)
	lambda_hat <- lambda_hat / 5
	
	# Add the background estimate as a new column to the GRanges object
	gr$background <- lambda_hat
	
	return(gr)
}

# Calculate gene signal aggregating signal from windows_data while reducing background
#	@gr					: Grange of counts
#	@tag_list			: Tag to calculate feature signal
#	@rm_feature_list	: List of features to exclude from the calculation
#	@verbose			: Display verbose
calculationSignal <- function(
	gr					= NULL,
	tag_list			= NULL,
	rm_feature_list	    = NULL,
	verbose				= TRUE
) {

	if (verbose){print(paste0("Calculate signal..."))}

	# Normalize background for the width
	gr$bgNorm <- mcols(gr)$background * (width(gr)/1000)
	# Extract only tagged features
	gr <- gr[gr$type %in% tag_list,]
	# Remove some unwanted features
	gr <- gr[!gr$feature %in% rm_feature_list,]

	# Split features names by ; that we used to aggregated common windows_data
	split_features <- strsplit(gr$feature, ";")
	# Find out which rows are gonna have to be duplicated
	expand_rows <- rep(1:length(gr), lengths(split_features))

	# Dupliate the rows
	gr <- gr[expand_rows]
	# Paste back the split feature names
	gr$feature <- unlist(split_features)
	# Sort the grange
	gr <- sort(sortSeqlevels(gr), ignore.strand=TRUE)

	# Split the grange by feature (which should now by unique, without ;)
	gr_split <- split(gr, ~feature)

	# Calculate the feature signal by substracting the sum of the background to the sum of counts
	# Set signal to 0 if background > signal
	# Takes around 30min
	res <- unlist(lapply(gr_split, function(x){
		C <- sum(x$counts)
		B <- sum(x$bgNorm)
		return(ifelse(C > B, C - B, 0))
	}))

	return(res)
}

# Calculate size factor for each sample based on quantiles normalization in regard of contrl samples
#	@df					: Dataframe of sample
#	@ctrl			    : Vector of control sample names
#	@housekeeping	    : Number of top expressed genes in normal samples to be considered as housekeeping genes for the normalization
#	@verbose			: Display verbose
calculation_sf <- function(
	df					= NULL,
	ctrl			    = NULL,
	housekeeping	    = NULL,
	verbose				= TRUE
) {

	if (verbose){print(paste0("Calculate size factors..."))}

	# Create a separate df for control samples 
	df_ctrl <- df[, ctrl, drop=FALSE]

	# Fetch highly expressed genes in Ctrl sample
	housekeeping_genes <- head(sort(rowMeans(df), decreasing=TRUE), housekeeping)
	df_house <- df[housekeeping_genes,]
	df_house_norm <- df[housekeeping_genes,]
	
	# Quantiles normalization of the housekeeping genes for each sample in regard to the average of control samples
	ctrl_values <- rowMeans(df_ctrl[housekeeping_genes,, drop=FALSE])
	
	# Create a temporary dataframe to calculate the quantile normalization between control and each sample
	for (sample in colnames(df)){
		tmp_df <- df[housekeeping_genes,sample, drop=FALSE]
		tmp_df <- cbind(tmp_df, data.frame(Ctrl = ctrl_values))
		df_house_norm[,sample] <- normalize.quantiles(as.matrix(tmp_df),copy=TRUE)[,1,drop=FALSE]
	}

	# Generalize the normalization to all genes using linear regression
	size.factor_list <- c()
	for(sample in colnames(df_house)){

		if (verbose){print(paste0("Linear regression for ", sample))}

		# Create a temporary dataframe to genereate the linear model
		tmp_df <- df_house[,sample, drop=FALSE]
		colnames(tmp_df) <- "raw"
		tmp_df <- cbind(tmp_df, data.frame(norm = df_house_norm[,sample]))

		# Calculate the linear regression
		lm_tmp <- lm(raw~norm, data = tmp_df)

		# Outputs statistics
		r_squared <- summary(lm_tmp)$r.squared
		intercept <- pvalue <- summary(lm_tmp)$coefficients[,1][2]
		pvalue <- summary(lm_tmp)$coefficients[,4][2]
		# THe actual size factor is the inverse of the intercept
		size.factor <- 1/intercept

		if (verbose){print(paste0("R² : ", r_squared))}
		if (verbose){print(paste0("p-value : ", pvalue))}
		if (verbose){print(paste0("Intercept : ", intercept))}

		size.factor_list <- c(size.factor_list, size.factor)
	}

	names(size.factor_list) <- colnames(df_house)

	if (verbose){print(paste0("Size factor before scaling : "))}
	if (verbose){print(size.factor_list)}

	# Normalize the size factor to reach 1000000 counts in the control sample on average
	scaling.factor <- mean(sum(rowMeans(df_ctrl))*size.factor_list[names(size.factor_list) %in% ctrl])/1000000
	size.factor_list <- size.factor_list/scaling.factor

	if (verbose){print(paste0("Size factor after scaling : "))}
	if (verbose){print(size.factor_list)}

	return(size.factor_list)
}
