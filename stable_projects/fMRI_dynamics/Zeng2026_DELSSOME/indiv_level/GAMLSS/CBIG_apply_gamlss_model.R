# Written by Fang Tian, Tianchu Zeng and CBIG under MIT license:
# https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md

# Apply Fitted GAMLSS Model to New Data
# Refactored to use unified prediction logic

library(gamlss)
library(splines)
library(jsonlite)

# ==== Parse command line arguments ====
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 3) stop("Usage: Rscript CBIG_apply_gamlss_model.R <model_path> <new_data_path> <output_path> [seed]")

model_path <- args[1]
new_data_path <- args[2]
output_path <- args[3]
seed <- if (length(args) >= 4) as.integer(args[4]) else 42

set.seed(seed)

# ==== Setup logging ====
output_dir <- dirname(output_path)
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)
log_conn <- file(file.path(output_dir, "application_log.txt"), open = "wt")

log_msg <- function(...) {
  msg <- paste0(..., collapse = "")
  cat(msg, "\n"); cat(msg, "\n", file = log_conn)
}

log_msg("GAMLSS MODEL APPLICATION LOG")
log_msg("Start time: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S"))
log_msg("Seed: ", seed)

# ==== Load fitted model ====
log_msg("Loading fitted model from: ", model_path)
model_obj <- readRDS(model_path)
model <- model_obj$model
target_var <- model_obj$target_var

# Ensure model has data attached for predict() internals
if (is.null(model$data) && !is.null(model_obj$config$input_data) && file.exists(model_obj$config$input_data)) {
  training_df <- read.csv(model_obj$config$input_data)
  training_df$site_scanner <- as.factor(training_df$site_scanner)
  training_df$numeric_sex <- as.factor(training_df$numeric_sex)
  model$data <- training_df
}

# Fix call object
if (!is.null(model$call)) {
  model$call$data <- model$data
  model$call$family <- as.name(model_obj$distribution)
}

# Source the helpers from next to THIS script so it can be run from anywhere.
script_dir <- dirname(sub("^--file=", "",
                         grep("^--file=", commandArgs(FALSE), value = TRUE)[1]))
if (is.na(script_dir) || script_dir == "") script_dir <- "."
source(file.path(script_dir, "CBIG_gamlss_utils.R"))
dist_info <- get_distribution_info(model)

# ==== Load new data ====
log_msg("Loading new data from: ", new_data_path)
new_df <- read.csv(new_data_path)

if (!target_var %in% colnames(new_df)) stop(paste("Target variable", target_var, "not found."))
new_df$numeric_sex <- as.factor(new_df$numeric_sex)
new_df$site_scanner <- as.factor(new_df$site_scanner)

# ==== Run Prediction and Harmonization ====
log_msg("Running prediction and harmonization...")

app_results <- run_gamlss_application(model, new_df, target_var, dist_info)

# ==== Save Output ====
log_msg("Saving output to: ", output_path)
write.csv(app_results$result_df, output_path, row.names = FALSE)

# Summary
log_msg("Summary:")
log_msg("  Mean correction: ", round(app_results$mean_correction, 4))
log_msg("  Original Mean: ", round(app_results$original_mean, 4))
log_msg("  Corrected Mean: ", round(app_results$corrected_mean, 4))

close(log_conn)
cat("Done.\n")
