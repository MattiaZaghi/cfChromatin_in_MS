library(DESeq2)
library(pwr)

# Load your normalized counts data
# Assuming your data is in a CSV file with rows as genomic bins and columns as samples
normalized_counts <- read.csv("normalized_counts.csv", row.names = 1)

# Example: Calculate the mean and standard deviation for each group
# Assuming you have two groups, Group A and Group B
group_A <- normalized_counts[, 1:5]  # First 5 samples as Group A
group_B <- normalized_counts[, 6:10] # Next 5 samples as Group B

mean_A <- rowMeans(group_A)
mean_B <- rowMeans(group_B)
sd_A <- apply(group_A, 1, sd)
sd_B <- apply(group_B, 1, sd)

# Calculate pooled standard deviation
pooled_sd <- sqrt((sd_A^2 + sd_B^2) / 2)

# Calculate effect size (Cohen's d)
effect_size <- (mean_A - mean_B) / pooled_sd

# Define parameters for power analysis
alpha <- 0.05       # Significance level
power <- 0.8        # Desired power

# Calculate required sample size
sample_size <- pwr.t.test(d = mean(effect_size), sig.level = alpha, power = power, type = "two.sample")$n

print(paste("Required sample size per group:", ceiling(sample_size)))
