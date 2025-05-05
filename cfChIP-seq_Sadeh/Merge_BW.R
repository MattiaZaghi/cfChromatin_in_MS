library(rtracklayer)

files <- list.files(path = "/proj/user/mattia/Analysis/Tracks/H3K27ac_hg38/", pattern = "^GSM.*\\.bw$", full.names = TRUE)

# Initialize an empty list to store the imported bigWig files
bw_list <- lapply(files, import)

# Merge the bigWig files by summing their scores
merged <- bw_list[[1]]
for (i in 2:length(bw_list)) {
  merged$score <- merged$score + bw_list[[i]]$score
}

export(merged, "path/to/your/directory/Baca.bw")