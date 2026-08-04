# Written by Fang Tian, Tianchu Zeng and CBIG under MIT license:
# https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md

# Script to compare fitted GAMLSS models across different distributions based on BIC

library(gamlss)

# Directory containing the GAMLSS distribution-search outputs (one
# fitted_model.rds per subdirectory). Pass it as the first argument:
#   Rscript CBIG_compare_gamlss_dists.R <base_dir>
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1) {
  stop("Usage: Rscript CBIG_compare_gamlss_dists.R <base_dir>")
}
base_dir <- args[1]
target_var <- "ei_ratio_mean"

# Find all fitted model files
model_files <- list.files(base_dir, pattern = "fitted_model.rds", recursive = TRUE, full.names = TRUE)

if (length(model_files) == 0) {
  stop(paste("No fitted models found in", base_dir))
}

results <- data.frame(
  Distribution = character(),
  AIC = numeric(),
  BIC = numeric(),
  Deviance = numeric(),
  Converged = logical(),
  stringsAsFactors = FALSE
)

cat("Comparing distributions for target:", target_var, "\n")
cat("Scanning directory:", base_dir, "\n\n")

for (f in model_files) {
  # Extract distribution name from path
  # Structure: .../gamlss_dist_search/<dist_name>/<target>/fitted_model.rds
  parts <- strsplit(f, "/")[[1]]
  
  # The distribution name is 2 levels up from the file
  dist_name <- parts[length(parts) - 2]
  
  # Check if this file belongs to the correct target variable
  target_in_path <- parts[length(parts) - 1]
  if (target_in_path != target_var) {
    next
  }
  
  tryCatch({
    obj <- readRDS(f)
    m <- obj$model
    
    results <- rbind(results, data.frame(
      Distribution = dist_name,
      AIC = AIC(m),
      BIC = BIC(m),
      Deviance = m$G.deviance,
      Converged = m$converged
    ))
  }, error = function(e) {
    cat("Error reading distribution", dist_name, ":", conditionMessage(e), "\n")
  })
}

# Sort by BIC
if (nrow(results) > 0) {
  results <- results[order(results$BIC), ]
  
  print(results)
  
  cat("\nBest distribution based on BIC:\n")
  print(results[1, ])
  
  # Save comparison to CSV
  write.csv(results, file.path(base_dir, "distribution_comparison.csv"), row.names = FALSE)
  cat("\nComparison saved to:", file.path(base_dir, "distribution_comparison.csv"), "\n")
} else {
  cat("No valid models found to compare.\n")
}
