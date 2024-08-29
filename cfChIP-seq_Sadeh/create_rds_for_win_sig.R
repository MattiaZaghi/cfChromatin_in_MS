library(tidyverse)
library(dplyr)
library(readr)


win_sig<-readRDS("/date/gcb/gcb_MZ/Analysis/cfChIP-seq/SetupFiles/H3K4me3/Win-sig.rds")
win_sig_CSV<-read_csv("/date/gcb/gcb_MZ/Analysis/cfChIP-seq/SetupFiles/H3K27ac_hg38/Win-sig.csv")



win_sig2<- list()
# Create a vector with all values set to 0
zero_vector <- rep(0, length(unique(win_sig_CSV$signature)))

# Get the names from the 'window' element
names_vector <- unique(win_sig_CSV$signature)

# Create a named vector with all values set to 0
zero_vector <- setNames(rep(0, length(names_vector)), names_vector)
# Add the new 'avg' element to the 'win_sig' list
win_sig2[["avg"]] <- zero_vector
win_sig2[["var"]] <- zero_vector


saveRDS(win_sig2,"/date/gcb/gcb_MZ/Analysis/cfChIP-seq/SetupFiles/H3K27ac_hg38/Win-sig.rds")
