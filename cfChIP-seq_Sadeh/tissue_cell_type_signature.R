library(dplyr)
library(tidyverse)

# Define the directories
dir1 <- "/date/gcb/gcb_MZ/Analysis/Samples/H3K4me3_roadmap_hg19/"

# Get the list of .rdata files in each directory
files1 <- list.files(path = dir1, pattern = "\\.rdata$", full.names = TRUE)

# Combine the file lists
all_files <- c(files1)

# Initialize an empty list to store the data
data_list <- list()

# Loop over the files and load each one into the list
for (file in all_files) {
  # Load the .rdata file
  data <- readRDS(file)
  
  # Get the name of the file (without extension) to use as the list element name
  name <- tools::file_path_sans_ext(basename(file))
  
  # Add the data to the list
  data_list[[name]] <- data
}

Glossary <- read_tsv("/date/gcb/gcb_MZ/roadmap_epigenomics/Glossary.tsv")
Group <- read_tsv("/date/gcb/gcb_MZ/roadmap_epigenomics/Group.tsv")
Glossary_Group <- inner_join(Glossary, Group)

# Initialize an empty list to store the subsets
subsets <- list()

# Iterate through the dataset
for (i in 1:nrow(Glossary_Group)) {
  tissue <- Glossary_Group$EDACC_NAME[i]
  groups <- strsplit(Glossary_Group$Groups[i], ";")[[1]]
  
  for (group in groups) {
    if (!group %in% names(subsets)) {
      subsets[[group]] <- c()
    }
    subsets[[group]] <- c(subsets[[group]], tissue)
  }
}

# Load H3K27ac windows
TSS.windows<-readRDS("/date/gcb/gcb_MZ/Analysis/cfChIP-seq/SetupFiles/H3K4me3_roadmap_hg19/Windows.rds")


# Convert the GRanges object to a data frame for easier manipulation
tss_df <- as.data.frame(TSS.windows)

# Ensure the indices match between counts and tss_df
tss_indices <- which((tss_df$type == "TSS" | tss_df$type == "EXTRA_TSS") & !grepl("^ENST", tss_df$name) & tss_df$name != "UNKNOWN")
# Initialize an empty list to store the filtered data
filtered_data <- list()

# Loop over each subset
for (subset in names(subsets)) {
  tissues <- subsets[[subset]]
  
  # Loop over each tissue in the subset
  for (tissue in tissues) {
    if (!is.null(data_list[[tissue]][["Counts.QQnorm"]])) {
      # Get the normalized signal for the tissue
      signal <- data_list[[tissue]][["Counts.QQnorm"]]
      
      # Ensure the indices match between counts and tss_df
      tss_indices <- which((tss_df$type == "TSS" | tss_df$type == "EXTRA_TSS") & !grepl("^ENST", tss_df$name) & tss_df$name != "UNKNOWN")
      
      # Get the indices of the rows that are above 35 in the current tissue
      high_indices <- which(signal > 35)
      
      
      # Further subset high_indices based on TSS regions
      high_indices <- intersect(high_indices, tss_indices)
      
      # Check each of these indices against all other tissues outside the subset
      for (index in high_indices) {
        # Assume initially that the value is below 15 in all other tissues outside the subset
        below_15_in_all <- TRUE
        
        # Check the value in the same row in all other tissues outside the subset
        for (other_tissue in names(data_list)) {
          if (other_tissue != tissue && !other_tissue %in% tissues && !is.null(data_list[[other_tissue]][["Counts.QQnorm"]])) {
            if (data_list[[other_tissue]][["Counts.QQnorm"]][index] > 15) {
              # If the value is not below 15 in any other tissue outside the subset, set the flag to FALSE and break the loop
              below_15_in_all <- FALSE
              break
            }
          }
        }
        
        # If the value is below 15 in all other tissues outside the subset, add it to the filtered data for the current subset
        if (below_15_in_all) {
          if (is.null(filtered_data[[subset]])) {
            filtered_data[[subset]] <- list(values = numeric(), indices = numeric())
          }
          filtered_data[[subset]][["values"]] <- c(filtered_data[[subset]][["values"]], signal[index])
          filtered_data[[subset]][["indices"]] <- c(filtered_data[[subset]][["indices"]], index)
        }
      }
    }
  }
}


# Print the filtered data
print(filtered_data)



# Initialize an empty dataframe
df <- data.frame("signature" = character(), "window" = integer())

# Loop through each tissue in the list
for(tissue in names(filtered_data)){
  # Get the indices for the current tissue
  indices <- filtered_data[[tissue]][["indices"]]
  
  # Create a temporary dataframe with the current tissue and its indices
  temp_df <- data.frame("signature" = tissue, "window" = indices)
  
  # Append the temporary dataframe to the main dataframe
  df <- rbind(df, temp_df)
}

# Print the dataframe
print(df)
df<-unique(df)

# Sort the dataset
sorted_data <- df[order(df$signature, df$window), ]

# Print the sorted dataset
print(sorted_data)

write_csv(sorted_data,"/date/gcb/gcb_MZ/Analysis/cfChIP-seq/SetupFiles/H3K27ac/Win-sig.csv")

win_sig<- list()
# Create a vector with all values set to 0
zero_vector <- rep(0, length(unique(sorted_data$signature)))

# Get the names from the 'window' element
names_vector <- unique(sorted_data$signature)

# Create a named vector with all values set to 0
zero_vector <- setNames(rep(0, length(names_vector)), names_vector)
# Add the new 'avg' element to the 'win_sig' list
win_sig[["avg"]] <- zero_vector
win_sig[["var"]] <- zero_vector

saveRDS(win_sig,"/date/gcb/gcb_MZ/Analysis/cfChIP-seq/SetupFiles/H3K27ac/Win-sig.rds")

tss_indices <- sorted_data$window

TSS.windows_TSS<-TSS.windows[tss_indices] %>% as.data.frame()

TSS.windows_TSS<-filtered_TSS.windows_TSS %>% dplyr::select(name,tissue) %>% unique()

combined_df<-cbind(TSS.windows_TSS,sorted_data)
