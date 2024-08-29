# Define the directories
dir1 <- "/date/gcb/gcb_MZ/Analysis/Samples/H3K27ac_roadmap_hg38"


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


# Get the names of the tissues
tissues <- names(data_list)

# Initialize the data frame with the first tissue's "Counts.QQnorm" values
Windows_counts <- data.frame(data_list[[tissues[1]]][["Counts.QQnorm"]])
names(Windows_counts) <- tissues[1]

# Loop over the remaining tissues in the data list
for(tissue in tissues[-1]) {
  # Get the "Counts.QQnorm" values for the current tissue
  counts <- data_list[[tissue]][["Counts.QQnorm"]]
  
  # Add the counts as a new column in the data frame
  Windows_counts[[tissue]] <- counts
}

# Initialize an empty list to store the filtered data
filtered_data <- list()

# Loop over each column in the data frame
for(tissue in names(Windows_counts)) {
  # Get the indices of the rows that are above 35 in the current column
  high_indices <- which(Windows_counts[[tissue]] > 35)

  # Check each of these indices against all other columns
  for(index in high_indices) {
    # Assume initially that the value is below 15 in all other columns
    below_15_in_all <- TRUE
    
    # Check the value in the same row in all other columns
    for(other_tissue in names(Windows_counts)) {
      if(other_tissue != tissue) {
        if(Windows_counts[index, other_tissue] >= 15) {
          # If the value is not below 15 in any other column, set the flag to FALSE and break the loop
          below_15_in_all <- FALSE
          break
        }
      }
    }
    
    # If the value is below 15 in all other columns, add it to the filtered data for the current column
    if(below_15_in_all) {
      if(!exists(tissue, filtered_data)) {
        filtered_data[[tissue]] <- list()
      }
      filtered_data[[tissue]][["values"]] <- c(filtered_data[[tissue]][["values"]], Windows_counts[index, tissue])
      filtered_data[[tissue]][["indices"]] <- c(filtered_data[[tissue]][["indices"]], index)
    }
  }
}

# Your vector
datasets <- c("neuron","oligodendrocyte","astrocyte","oligodendrocyte_precursor_cell")



# Get the names of the tissues that are in your datasets vector
single_cells <- names(data_list)[names(data_list) %in% datasets]

# Print single_cells
print(single_cells)


# Initialize the data frame with the first tissue's "Counts.QQnorm" values
Windows_counts2 <- data.frame(data_list[[single_cells[1]]][["Counts.QQnorm"]])
names(Windows_counts2) <- single_cells[1]

# Loop over the remaining tissues in the data list
for(single_cell in single_cells[-1]) {
  # Get the "Counts.QQnorm" values for the current tissue
  counts <- data_list[[single_cell]][["Counts.QQnorm"]]
  
  # Add the counts as a new column in the data frame
  Windows_counts2[[single_cell]] <- counts
}

# Initialize an empty list to store the filtered data
filtered_data2 <- list()

# Loop over each column in the data frame
for(single_cell in names(Windows_counts2)) {
  # Get the indices of the rows that are above 35 in the current column
  high_indices <- which(Windows_counts2[[single_cell]] > 35)
  
  # Check each of these indices against all other columns
  for(index in high_indices) {
    # Assume initially that the value is below 15 in all other columns
    below_15_in_all <- TRUE
    
    # Check the value in the same row in all other columns
    for(other_single_cell in names(Windows_counts2)) {
      if(other_single_cell != single_cell) {
        if(Windows_counts2[index, other_single_cell] >= 20) {
          # If the value is not below 15 in any other column, set the flag to FALSE and break the loop
          below_15_in_all <- FALSE
          break
        }
      }
    }
    
    # If the value is below 15 in all other columns, add it to the filtered data for the current column
    if(below_15_in_all) {
      if(!exists(single_cell, filtered_data2)) {
        filtered_data2[[single_cell]] <- list()
      }
      filtered_data2[[single_cell]][["values"]] <- c(filtered_data2[[single_cell]][["values"]], Windows_counts[index, single_cell])
      filtered_data2[[single_cell]][["indices"]] <- c(filtered_data2[[single_cell]][["indices"]], index)
    }
  }
}



filtered_data_final<-c(filtered_data,filtered_data2)


# 'filtered_data' now contains the values for each column that are above 35 in the current column and below 15 in all other columns
# Find the difference
# Get the names of the elements in the lists
names1 <- names(data_list)
names2 <- names(filtered_data_final)

# Find the difference in the names
diff_names <- setdiff(names1, names2)

# Print the difference in the names
print(diff_names)


# Initialize an empty dataframe
df <- data.frame("signature" = character(), "window" = integer())

# Loop through each tissue in the list
for(tissue in names(filtered_data_final)){
  # Get the indices for the current tissue
  indices <- filtered_data_final[[tissue]][["indices"]]
  
  # Create a temporary dataframe with the current tissue and its indices
  temp_df <- data.frame("signature" = tissue, "window" = indices)
  
  # Append the temporary dataframe to the main dataframe
  df <- rbind(df, temp_df)
}

# Print the dataframe
print(df)

library(readr)
library(dplyr)
library(tidyverse)
glossary<-read_delim("/date/gcb/gcb_MZ/Chrom_HMM_hg38_annotation/Glossary.txt")


df$EDACC_NAME <- df$signature


df<-left_join(df,glossary) %>% dplyr::select(signature,window,GROUP)

# Assuming df is your dataframe and 'signature' is the column name
df$signature <- sub("Adipose_Nuclei", "Adipose", df$signature)
df$signature <- sub("Adult_Liver", "Liver", df$signature)

df$signature <- sub("CD19_Primary_Cells_Peripheral_UW", "B-Cells", df$signature)
df$signature <- sub("CD3_Primary_Cells_Peripheral_UW", "T-Cells", df$signature)
df$signature <- sub("CD3_Primary_Cells_Peripheral_UW", "T-Cells", df$signature)
df$signature <- sub("CD56_Primary_Cells", "NK-Cells", df$signature)
df$signature <- sub("Stomach_Smooth_Muscle", "Sm.Muscle", df$signature)
df$signature <- sub("Gastric", "Digestive", df$signature)
df$signature <- sub("Left_Ventricle", "Heart", df$signature)
df$signature <- sub("Fetal_Placenta", "Placenta", df$signature)
df$signature <- sub("CD4_Naive_Primary_Cells", "T-Cells", df$signature)
df$signature <- sub("CD4+_CD25+_CD127-_Treg_Primary_Cells", "T-Reg", df$signature)
df$signature <- sub("NHDF-Ad_Adult_Dermal_Fibroblasts", "Epithelial", df$signature)
df$signature <- sub("Right_Atrium", "Heart", df$signature)
df$signature <- sub("Skeletal_Muscle_Female", "Muscle", df$signature)
df$signature <- sub("CD4+_CD25-_Th_Primary_Cells", "T-Helper", df$signature)
df$signature <- sub("CD8_Naive_Primary_Cells", "T-Cells", df$signature)
df$signature <- sub("CD8_Memory_Primary_Cells", "T-Cells", df$signature)
df$signature <- sub("Placenta_Amnion", "Placenta", df$signature)
df$signature <- sub("NHEK-Epidermal_Keratinocytes", "Epithelial", df$signature)
df$signature <- sub("CD4+_CD25-_CD45RO+_Memory_Primary_Cells", "T-Cells", df$signature)
df$signature <- sub("Pancreas", "Pancreas", df$signature)
df$signature <- sub("Pancreatic_Islets", "Pancreas", df$signature)
df$signature <- sub("Psoas_Muscle", "Muscle", df$signature)
df$signature <- sub("Rectal_Smooth_Muscle", "Sm.Muscle", df$signature)
df$signature <- sub("Right_Ventricle", "Heart", df$signature)
df$signature <- sub("Rectal_Mucosa.Donor_29", "Digestive", df$signature)
df$signature <- sub("Rectal_Mucosa.Donor_31", "Digestive", df$signature)
df$signature <- sub("Colonic_Mucosa", "Digestive", df$signature)
df$signature <- sub("Esophagus", "Digestive", df$signature)
df$signature <- sub("Sigmoid_Colon", "Digestive", df$signature)
df$signature <- sub("Peripheral_Blood_Mononuclear_Primary_Cells", "PBMCs", df$signature)

# Define the entries you want to replace
entries_to_replace <- c("CD4+_CD25-_CD45RA+_Naive_Primary_Cells", 
                        "CD4+_CD25-_CD45RO+_Memory_Primary_Cells", 
                        "CD4+_CD25-_IL17-_PMA-Ionomycin_stimulated_MACS_purified_Th_Primary_Cells", 
                        "CD4+_CD25-_IL17+_PMA-Ionomcyin_stimulated_Th17_Primary_Cells", 
                        "CD4+_CD25-_Th_Primary_Cells", 
                        "CD4+_CD25+_CD127-_Treg_Primary_Cells", 
                        "CD4+_CD25int_CD127+_Tmem_Primary_Cells")

# Replace the entries in the 'signature' column
df$signature[df$signature %in% entries_to_replace] <- "T-cells"

entries_to_replace <- c("CD14_Primary_Cells","Monocytes-CD14+_RO01746")

# Replace the entries in the 'signature' column
df$signature[df$signature %in% entries_to_replace] <- "Monocytes"

print(unique(df$signature))

df<-df %>% dplyr::select(signature,window)

write_csv(df,"/date/gcb/gcb_MZ/Analysis/cfChIP-seq/SetupFiles/H3K27ac_hg38/Win-sig.csv")

win_sig<- list()
# Create a vector with all values set to 0
zero_vector <- rep(0, length(unique(df$signature)))

# Get the names from the 'window' element
names_vector <- unique(df$signature)

# Create a named vector with all values set to 0
zero_vector <- setNames(rep(0, length(names_vector)), names_vector)
# Add the new 'avg' element to the 'win_sig' list
win_sig[["avg"]] <- zero_vector
win_sig[["var"]] <- zero_vector

saveRDS(win_sig2,"/date/gcb/gcb_MZ/Analysis/cfChIP-seq/SetupFiles/H3K27ac_hg38/Win-sig.rds")

