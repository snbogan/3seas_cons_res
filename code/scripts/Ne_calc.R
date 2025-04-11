#!/usr/bin/env Rscript

# Load required packages
if (!require("data.table", quietly = TRUE)) {
  install.packages("data.table", repos = "https://cloud.r-project.org")
  library(data.table)
}

# Read and process LD data
ld_data <- fread("ld_output.geno.ld", header = TRUE)

# Identify the correct r² column
r2_col <- if ("R.2" %in% names(ld_data)) {
  "R.2"
} else if ("R^2" %in% names(ld_data)) {
  "R^2"
} else {
  names(ld_data)[grep("r2|R2", names(ld_data), ignore.case = TRUE)[1]]
}

if (is.na(r2_col)) {
  stop("Could not find r² column in LD output file")
}

# Read and process MAF data with proper column handling
maf_data <- fread("maf_output.frq", header = TRUE)

# Print column names and first few rows for debugging
cat("\nMAF file structure:\n")
print(head(maf_data))

# Identify MAF column - vcftools --freq typically has format:
# CHROM   POS     N_ALLELES       N_CHR   {ALLELE:FREQ}
# We need to extract the frequency of the minor allele
if (ncol(maf_data) >= 5) {
  # Extract frequency from the last column (format like "A:0.95" or "T:0.25")
  maf_data[, MAF := as.numeric(gsub("^.*:", "", maf_data[[ncol(maf_data)]]))]
} else {
  stop("MAF file doesn't have expected format. Expected at least 5 columns.")
}

# Filter data
ld_data_filtered <- ld_data[get(r2_col) > 0 & get(r2_col) < 1, ]
maf_data_filtered <- maf_data[MAF > 0 & !is.na(MAF), ]

# Calculate means
mean_r2 <- mean(ld_data_filtered[[r2_col]], na.rm = TRUE)
mean_maf <- mean(maf_data_filtered$MAF, na.rm = TRUE)

# Estimate Ne
distance_cM <- 50000 / 1e6  # 50kb ≈ 0.05 cM
Ne_estimate <- 1 / (4 * mean_r2 * distance_cM)

# Print results
cat("\nEffective population size (Ne) estimation:\n")
cat("Mean r² (excluding 0 & 1):", mean_r2, "\n")
cat("Mean MAF (excluding 0):", mean_maf, "\n")
cat("Estimated Ne (LD method):", round(Ne_estimate), "\n")

# Save results
fwrite(
  data.table(
    Mean_r2 = mean_r2,
    Mean_MAF = mean_maf,
    Ne_estimate = Ne_estimate,
    N_SNPs_LD = nrow(ld_data_filtered),
    N_SNPs_MAF = nrow(maf_data_filtered)
  ),
  "Ne_estimation_results.txt",
  sep = "\t",
  quote = FALSE
)