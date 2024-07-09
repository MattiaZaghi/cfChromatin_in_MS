# Assuming 'metadata' is your Roadmap Epigenomics metadata table
# and 'normalized_signals' is a matrix/data frame of your normalized signals N[w; s]

# Define your groups
groups <- list(
  "Lymphocytes" = c("B-Cells", "T-Cells", "NK")
  # Add other groups as needed
)

# Initialize an empty list to store the signatures
signatures <- list()

# Loop over each group
for (group_name in names(groups)) {
  group_samples <- groups[[group_name]]
  
  # Get the normalized signals for the group
  group_signals <- normalized_signals[, group_samples]
  
  # Criteria 1: The window w is on an autosomal chromosome
  # This will depend on how your data is structured. You'll need a way to identify which windows are on autosomal chromosomes.
  
  # Criteria 2: In at least one of the atlas samples in the group, N[w; s] >= 35
  criteria_2 <- rowSums(group_signals >= 35) > 0
  
  # Criteria 3: In all atlas samples outside the group, N[w; s] < 15
  non_group_signals <- normalized_signals[, !colnames(normalized_signals) %in% group_samples]
  criteria_3 <- rowSums(non_group_signals >= 15) == 0
  
  # Criteria 4: In all windows w within 1Kb of w, N[w; s] < 15
  # This will depend on how your data is structured. You'll need a way to identify which windows are within 1Kb of each other.
  
  # Apply the criteria
  specific_windows <- criteria_2 & criteria_3 # & criteria_1 & criteria_4 if you have these
  
  # Check if the group has less than 4 specific windows
  if (sum(specific_windows) < 4) {
    cat(group_name, "has no signature\n")
  } else {
    # Define the signature as the set of specific windows
    signatures[[group_name]] <- which(specific_windows)
  }
}
