# Written by Fang Tian, Tianchu Zeng and CBIG under MIT license:
# https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md

# Script to compare fitted GAMLSS models based on BIC

library(gamlss)

# Directory containing the GAMLSS model-selection outputs (one
# fitted_model.rds per subdirectory). Pass it as the first argument:
#   Rscript CBIG_compare_gamlss_models.R <base_dir>
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1) {
  stop("Usage: Rscript CBIG_compare_gamlss_models.R <base_dir>")
}
base_dir <- args[1]
target_var <- "ei_ratio_mean"

# Find all fitted model files
model_files <- list.files(base_dir, pattern = "fitted_model.rds", recursive = TRUE, full.names = TRUE)

results <- data.frame(
  Config = character(),
  Distribution = character(),
  AIC = numeric(),
  BIC = numeric(),
  Deviance = numeric(),
  Converged = logical(),
  stringsAsFactors = FALSE
)

cat("Comparing models for target:", target_var, "\n\n")

for (f in model_files) {
  # Extract config name from path
  # Assumes structure: .../gamlss_selection/<config_name>/<target>/fitted_model.rds
  parts <- strsplit(f, "/")[[1]]
  config_name <- parts[length(parts) - 2]
  
  tryCatch({
    obj <- readRDS(f)
    m <- obj$model
    
    results <- rbind(results, data.frame(
      Config = config_name,
      Distribution = obj$distribution,
      AIC = AIC(m),
      BIC = BIC(m),
      Deviance = m$G.deviance,
      Converged = m$converged
    ))
  }, error = function(e) {
    cat("Error reading", config_name, ":", conditionMessage(e), "\n")
  })
}

# Sort by BIC
results <- results[order(results$BIC), ]

print(results)

cat("\nBest model based on BIC:\n")
print(results[1, ])
