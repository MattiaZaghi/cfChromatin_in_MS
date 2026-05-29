# Load the required library
library(rtracklayer)
library(limma)
library(preprocessCore)

#load the counts list needed to obtain the normalize signal 

counts_list<-readRDS("cfChromatin_in_MS/H3K27ac_ChIP_counts.rds")

# Define the function to calculate raw signal at gene level
rawSignal <- function(gr) {
  # Extract counts and background from the GRanges object
  C <- mcols(gr)$counts
  B <- mcols(gr)$background * (width(gr)/1000)
  
  # Calculate raw signal
  S <- ifelse(C > B, C - B, 0)
  
  # Add the raw signal as a new column to the GRanges object
  gr$rawSignal <- S
  
  return(gr)
}

# Apply the rawSignal function to each GRanges object in the list
rawSignal_counts <- mclapply(counts_list, rawSignal, mc.cores = detectCores())
rm(counts_list)
#define house keeping genes
# Define the names of the samples you want to exclude
exclude_samples <- c("CD19-GSE172116-1_H3K27ac_ChIP", "CD19-GSE172116-2_H3K27ac_ChIP")

# Remove the samples from the GRanges list
rawSignal_counts_blood_only <- rawSignal_counts[!names(rawSignal_counts) %in% exclude_samples]

# Extract the 'name' and 'rawSignal' columns from the GRanges list
data <- do.call(rbind, lapply(rawSignal_counts_blood_only, function(gr) data.frame(name = gr$name, rawSignal = gr$rawSignal, stringsAsFactors = FALSE)))


# Apply the operation to all GRanges objects in the list
gr_list_modified <- lapply(rawSignal_counts_blood_only, function(gr) {
  gr$region <- paste0(seqnames(gr), ":", start(gr), "-", end(gr), "-", seq_along(gr))
  return(gr)
})

total_signal<-do.call(rbind, lapply(gr_list_modified, function(gr) data.frame(name = gr$region, rawSignal = gr$rawSignal, stringsAsFactors = FALSE)))

# Remove row names
rownames(total_signal) <- NULL


# Group by 'region' and summarize 'rawSignal'
total_signal<- total_signal %>%
  group_by(name) %>%
  summarize(rawSignal = sum(rawSignal, na.rm = TRUE))

total_signal<-as.data.frame(total_signal)


rownames(total_signal)<-total_signal$name


total_signal$name<-NULL

# Remove rows with NA in the 'name' column
data <- data[!is.na(data$name) & data$name != "", ] 

# Calculate the total signal for each gene across all samples
total_signal_gene <- tapply(data$rawSignal, data$name, sum)

# Calculate the standard deviation for each gene across all samples
std_dev <- tapply(data$rawSignal, data$name, sd)

# Calculate the median of the Poisson distribution
# Assuming 'lambda' is the mean of your data
lambda <- median(total_signal_gene)
poisson_median <- qpois(0.7, lambda)

# Define a standard deviation threshold
std_dev_threshold <- 5  # replace with your actual threshold

# Identify the housekeeping genes
# Housekeeping genes are those with a signal above the median of the Poisson distribution and a standard deviation below a certain threshold
housekeeping_genes <- names(total_signal_gene)[total_signal_gene > poisson_median & std_dev <= std_dev_threshold]

# Initialize an empty list to store the count data frames
count_dfs <- list()

# Loop over the GRanges list
for(i in seq_along(rawSignal_counts_blood_only)) {
  # Extract the 'name' and 'counts' columns from the GRanges object
  data <- data.frame(name = rawSignal_counts_blood_only[[i]]$name, rawSignal = rawSignal_counts_blood_only[[i]]$rawSignal, stringsAsFactors = FALSE)
  
  # Group by 'name' and summarize 'counts'
  Windows <- data %>%
    group_by(name) %>%
    summarize(rawSignal = mean(rawSignal))
  
  # Remove rows with NA in the 'name' column
  All_genes <- Windows[!is.na(Windows$name) & Windows$name != "", ] 
  
  All_genes<-as.data.frame(All_genes)
  
  # Let's filter only the housekeeping genes that will be used in the analysis
  housekeeping <- All_genes %>% dplyr::filter(name %in% housekeeping_genes)
  
  # Set the 'name' column as row names
  row.names(housekeeping) <- housekeeping$name
  
  # Remove the 'name' column
  housekeeping$name <- NULL
  
  # Add the data frame to the list
  count_dfs[[i]] <- housekeeping
}
# Combine all data frames into a single count matrix
count_matrix <- do.call(cbind, count_dfs)

# Set the column names of the count matrix as the names of the GRanges objects
colnames(count_matrix) <- names(rawSignal_counts_blood_only)

#filter housekeeping genes from the gene count matrix

S_housekeeping <-as.matrix(count_matrix)

normalizeSignal <- function(gr_list, housekeeping_genes, reference_set) {

# Loop over the GRanges list
for(i in seq_along(gr_list)) {
  # Apply quantile normalization to the raw signal of the housekeeping genes
  normalized_housekeeping_genes <- limma::normalizeQuantiles(S_housekeeping)
  
  # Extract the raw signal for all genes
  S_all <- as.matrix(total_signal)
  
  # Use linear regression to estimate a multiplicative normalization factor for each sample
  scaling_factor <- coef(lm(normalized_housekeeping_genes ~ S_housekeeping))[2]
  
  # Rescale the scaling factors so that the total normalized signal at the set of reference healthy samples will be a million on average
  scaling_factor <- scaling_factor * 1e6 / scaling_factor * total_signal
  
  # Calculate normalized signal
  N <- scaling_factor * S_all
  
  # Assuming 'N' is your data frame
  # Assuming df is your dataframe and you want to extract the 3rd row
  
  # Now, 'row' contains the 3rd row of df

  # Split the row names into seqnames, start, and end
  coords <- strsplit(rownames(N), "-")
  seqnames <- sapply(coords, function(x) unlist(strsplit(x[1], ":"))[1])
  starts <- sapply(coords, function(x) as.numeric(unlist(strsplit(x[1], ":"))[2]))
  ends <- sapply(coords, function(x) as.numeric(x[2]))
  
  # Remove rows with NA in 'start' or 'end'
  na_index <- is.na(starts) | is.na(ends)
  seqnames <- seqnames[!na_index]
  starts <- starts[!na_index]
  ends <- ends[!na_index]
  
  # Now create the GRanges object
  gr <-data.frame(seqnames = seqnames, start = starts, end = ends) %>% makeGRangesFromDataFrame()
  # Assuming df is your dataframe
  gr_na <- gr[apply(is.na(gr), 1, any), ]
  
  # Now, df_na contains the rows with NA values from df
  
  
  # Add the normalized signal as a metadata column
  mcols(gr)$normalizedSignal <- N$rawSignal # replace 'column_name' with your actual column name
  
  
  # Assuming 'gr' is your new GRanges object with the 'normalizedSignal' column
  # and 'gr_list' is your list of GRanges objects
  
  # Apply the operation to all GRanges objects in the list
  gr_list_modified <- lapply(gr_list, function(gr_old) {
    # Find the overlaps between the old and new GRanges objects
    olaps <- findOverlaps(gr_old, gr)
    
    # Add the 'normalizedSignal' column from the new GRanges object to the old GRanges object
    mcols(gr_old)$normalizedSignal <- mcols(gr)[subjectHits(olaps)]
    
    return(gr_old)
  })
  
    }
    
    return(gr_list)
}


# replace with your actual housekeeping genes
reference_set <- names(rawSignal_counts_blood_only)  # replace with your actual reference set

# Apply the normalizeSignal function to each GRanges object in the list
normalizedSignal_counts <-normalizeSignal(gr_list = rawSignal_counts,housekeeping_genes = housekeeping_genes,reference_set = reference_set)


# Assuming 'gr_list' is your list of GRanges objects, 'old_name' is the old column name
# and 'new_name' is the new column name you want to use

# Apply the operation to all GRanges objects in the list
normalizedSignal_counts  <- lapply(gr_list, function(gr) {
  # Get the current column names
  col_names <- colnames(mcols(gr))
  
  # Find the index of the column you want to rename
  index <- which(col_names == "old_name")
  
  # Change the column name
  col_names[index] <- "new_name"
  
  # Assign the new column names back to the GRanges object
  colnames(mcols(gr)) <- col_names
  
  return(gr)
})


# Define a function to create a BigWig file from a GRanges object
createBigWig <- function(gr, filename) {
  # Create a GRanges object with the normalized signal as the score
  gr_with_score <- GRanges(seqnames = seqnames(gr),
                           ranges = ranges(gr),
                           strand = strand(gr),
                           score = mcols(gr)$normalizedSignal)
  
  # Export the GRanges object to a BigWig file
  export.bw(gr_with_score, filename)
}

# Apply the createBigWig function to each GRanges object in the list
filenames <- paste0("/date/gcb/gcb_MZ/Round1/R_bigiwig/", seq_along(normalizedSignal_counts), ".bw")  # replace with your actual filenames
mapply(createBigWig, normalizedSignal_counts, filenames)

