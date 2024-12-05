suppressMessages(library(rtracklayer))
suppressMessages(library(tools))
suppressMessages(library(ggplot2))
suppressMessages(library(corrplot))
# Define a custom function 'catn' that behaves like 'cat' but appends a newline at the end of the output.
catn = function(...) { cat(...,"\n") }

cfChIP.BED.suffixes = c(".bed", ".bed.gz", ".tagAlign", ".tagAlign.gz")
cfChIP.BW.suffixes = c(".bw", ".bigWig", ".bw.gz", ".bigWig.gz")
cfChIP.File.suffixes = c(cfChIP.BED.suffixes, cfChIP.BW.suffixes)

cfChIP.FindFile <- function(filename) {
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

cfChIP.GetRawData = function(filename) {
  MinFragLen = 50
  MaxFragLen = 800
  Verbose = TRUE
  ll = cfChIP.FindFile(filename)
  filename = ll$filename
  FileType = ll$FileType
  
  dat = list()
  if( FileType == "BED") {     
    if(Verbose ) catn(filename, ": Reading BED file")
    
    dat$RawBED = import(filename, format = "BED")
    # remove long/short fragments and non-unique copies
    
    # check for single end reads
    if( max(width(dat$RawBED)) <  MinFragLen) {
      dat$RawBED = resize(dat$RawBED, width = 166)
    } else 
      dat$RawBED = dat$RawBED[width(dat$RawBED) <= MaxFragLen & width(dat$RawBED) > MinFragLen]
    
    dat$BED = unique(dat$RawBED)
    dat$Cov = coverage(dat$BED)
  } 
  
  return(dat)  # Return the 'dat' list
}
# Read the TSS windows data from an RDS file located in the SetupDIR directory.
TSS.windows = readRDS("/date/gcb/gcb_MZ/Analysis/cfChIP-seq/SetupFiles/H3K27ac_hg38/Windows.rds")
# Retrieve sequence information from the TSS windows data.
genome.seqinfo = seqinfo(TSS.windows)
# Define a list of chromosome names typically found in human genome datasets.
ChrList = paste0("chr", c(1:22,"X", "Y"))
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

cfChIP.GetCoverage  = function(filename) {
  dat = cfChIP.GetRawData(filename)
  return(fixCoverage(dat$Cov))
}



# Create a vector of file paths
file_paths_ref <- "/date/gcb/gcb_MZ/Analysis/Samples/H3K27ac_roadmap_hg38"
files_ref <- list.files(path = file_paths_ref, pattern = "*rdata$", full.names = TRUE)
# Remove the files that contain "_H3K4me3_ChIP"
files_ref  <- files_ref [!grepl("_H3K4me3_ChIP", files_ref)]
# Use grep to get the indices of files that start with 'H'
files_ref  <- files_ref[c(4:25,36)]
file_paths_samples <- "/date/gcb/gcb_MZ/Analysis/Samples/H3K27ac_hg38"

files_samples <- list.files(path = file_paths_samples, pattern = "*_ChIP.rdata$", full.names = TRUE)
files_samples <-files_samples[c(2,4:6,8,10,22:26,31)]
# Merge the lists
merged_list <- c(files_ref,files_samples)
# Use lapply to load all files
data_list <- lapply(merged_list ,readRDS)



# Extract the base names from the file paths
base_names <- basename(merged_list)

# Remove the pattern "__H3K27ac_ChIP-SE.rdata" from the base names
names <- sub("_H3K27ac_ChIP.rdata$", "", base_names)
# Extract the first part of each name
shortened_names <- sub("([^-_+]+).*", "\\1", names)

# Assign names to the list elements
names(data_list) <- shortened_names
# Create a vector of bed file paths
file_paths <- c("/date/gcb/gcb_MZ/Analysis/BED/H3K27ac/H1-P_H3K27ac_ChIP.bed",
                "/date/gcb/gcb_MZ/Analysis/BED/H3K27ac/H2-P_H3K27ac_ChIP.bed",
                "/date/gcb/gcb_MZ/Analysis/BED/H3K27ac/H3-P_H3K27ac_ChIP.bed",
                "/date/gcb/gcb_MZ/Analysis/BED/H3K27ac/H4-P_H3K27ac_ChIP.bed",
                "/date/gcb/gcb_MZ/Analysis/BED/H3K27ac/19019-P_H3K27ac_ChIP.bed")
# Apply the function to each file and add the result to the corresponding element in data_list
data_list <- lapply(seq_along(file_paths), function(i) {
  # Get the file path
  file_path = file_paths[i]
  
  # Calculate the new data
  new_data = cfChIP.GetCoverage(file_path)
  
  # Add the new data to the existing data
  data_list[[i]] = c(data_list[[i]], new_data)
  
  # Calculate the multiplication of QQNorm and Cov and add it to the data
  if ("QQNorm" %in% names(data_list[[i]]) && "Cov" %in% names(data_list[[i]])) {
    data_list[[i]]$QQNorm_Cov = data_list[[i]]$QQNorm * data_list[[i]]$Cov
  }
  
  # Return the updated data
  return(data_list[[i]])
})

# Assign names to the list elements
names(data_list) <- object_names

data_list[["19019"]][["Cov"]]<-tail(data_list[["19019"]][["Cov"]],2000)



# Assuming H1, H2, H3, and H4 are your objects with gene counts
gene_counts_list <- lapply(data_list, function(x) x$GeneCounts.QQnorm)

# Combine them into a data frame
gene_counts_df <- do.call(cbind, gene_counts_list)

#Create a matrix from the normalized counts

norm_counts <- as.matrix(gene_counts_df)

gene_counts_df<-gene_counts_df[rownames(gene_counts_df) != "UNKNOWN", ]
# Assuming norm_counts is your matrix of normalized counts
correlation_matrix <- cor(norm_counts)

# Create the correlation plot with squares and positive values in red
corrplot(correlation_matrix, method = "color", col = colorRampPalette(c("blue", "white", "red"))(200), tl.cex = 0.8)


# Assuming H1, H2, H3, and H4 are your objects with gene counts
gene_counts_list <- lapply(data_list, function(x) x$Counts.QQnorm)

# Combine them into a data frame
gene_counts_df <- do.call(cbind, gene_counts_list)

#Create a matrix from the normalized counts

norm_counts <- as.matrix(gene_counts_df)

# Assuming norm_counts is your matrix of normalized counts
correlation_matrix <- cor(norm_counts)

# Create the correlation plot with squares and positive values in red
corrplot(correlation_matrix, method = "color", col = colorRampPalette(c("blue", "white", "red"))(200), tl.cex = 0.8)

library(S4Vectors)

# Define a function to bin an Rle object
bin_Rle <- function(rle, bin_width) {
  # Get the start and end of each run
  starts <- start(rle)
  ends <- end(rle)
  
  # Calculate the bin number for each start and end
  start_bins <- ceiling(starts / bin_width)
  end_bins <- floor(ends / bin_width)
  
  # Initialize a list to store the binned values
  binned_values <- list()
  
  # Loop over the runs
  for (i in seq_along(rle)) {
    # If the start and end are in the same bin, add the run to that bin
    if (start_bins[i] == end_bins[i]) {
      binned_values[[start_bins[i]]] <- c(binned_values[[start_bins[i]]], rep(rle[i], ends[i] - starts[i] + 1))
    } else {
      # If the start and end are in different bins, split the run across the bins
      for (j in start_bins[i]:end_bins[i]) {
        if (j == start_bins[i]) {
          # For the start bin, add the part of the run from the start to the end of the bin
          binned_values[[j]] <- c(binned_values[[j]], rep(rle[i], j * bin_width - starts[i] + 1))
        } else if (j == end_bins[i]) {
          # For the end bin, add the part of the run from the start of the bin to the end
          binned_values[[j]] <- c(binned_values[[j]], rep(rle[i], ends[i] - (j - 1) * bin_width))
        } else {
          # For the middle bins, add the entire bin
          binned_values[[j]] <- c(binned_values[[j]], rep(rle[i], bin_width))
        }
      }
    }
  }
  
  # Convert the list of binned values to an Rle object
  binned_rle <- Rle(unlist(binned_values))
  
  return(binned_rle)
}

# Apply the function to each Rle object in your list
binned_data_list <- lapply( data_list[["19019"]][["Cov"]], function(x) bin_Rle(x, 2000))


