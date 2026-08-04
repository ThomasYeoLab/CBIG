# Written by Fang Tian, Tianchu Zeng and CBIG under MIT license:
# https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md

# Main GAMLSS Model Fitting Script
# This script fits a GAMLSS model with flexible specifications and saves outputs

library(gamlss)
library(splines)
library(jsonlite)

# Source utility functions from next to THIS script, so the script can be run
# from any working directory (Rscript exposes its own path via --file=).
script_dir <- dirname(sub("^--file=", "",
                         grep("^--file=", commandArgs(FALSE), value = TRUE)[1]))
if (is.na(script_dir) || script_dir == "") script_dir <- "."
source(file.path(script_dir, "CBIG_gamlss_utils.R"))

# ==== Parse command line arguments ====
args <- commandArgs(trailingOnly = TRUE)

# Defaults
config_path <- NULL
input_data_override <- NULL
target_var_override <- NULL
output_dir_override <- NULL
distribution_override <- NULL
apply_input_path <- NULL
apply_output_path <- NULL
seed <- 42
skip_fit <- FALSE

# Manual Argument Parser
i <- 1
while(i <= length(args)) {
  arg <- args[i]
  if (arg == "--skip-fit") {
    skip_fit <- TRUE
  } else if (arg == "--input_data") {
    i <- i + 1
    input_data_override <- args[i]
  } else if (arg == "--target") {
    i <- i + 1
    target_var_override <- args[i]
  } else if (arg == "--output_dir") {
    i <- i + 1
    output_dir_override <- args[i]
  } else if (arg == "--distribution") {
    i <- i + 1
    distribution_override <- args[i]
  } else if (arg == "--apply_input") {
    i <- i + 1
    apply_input_path <- args[i]
  } else if (arg == "--apply_output") {
    i <- i + 1
    apply_output_path <- args[i]
  } else if (arg == "--seed") {
    i <- i + 1
    seed <- as.integer(args[i])
  } else if (!startsWith(arg, "-")) {
    config_path <- arg
  }
  i <- i + 1
}

if (is.null(config_path)) {
  stop(paste0(
    "Usage: Rscript CBIG_fit_gamlss_model.R <config_json_path> ",
    "[--input_data <file>] [--target <var>] [--output_dir <dir>] [--distribution <dist>] ",
    "[--apply_input <file>] [--apply_output <file>] [--seed <int>] [--skip-fit]"))
}

# Set random seed
set.seed(seed)

config <- fromJSON(config_path)

# ==== Load configuration with Overrides ====
input_data <- if (!is.null(input_data_override)) input_data_override else config$input_data
if (is.null(input_data)) {
  stop("No input table: give --input_data <file> or set \"input_data\" in the config JSON.")
}
if (!file.exists(input_data)) {
  stop(paste0("Input table not found: ", input_data,
              "\n  Pass --input_data <file> to point at your own table. ",
              "See GAMLSS/README.md for the required columns."))
}
target_var <- if (!is.null(target_var_override)) target_var_override else config$target_variable
output_dir <- if (!is.null(output_dir_override)) output_dir_override else config$output_dir
distribution <- if (!is.null(distribution_override)) distribution_override else config$distribution

# Use config values for others
age_formula_mu <- config$age_formula_mu
age_formula_sigma <- config$age_formula_sigma
age_formula_nu <- config$age_formula_nu  # Optional
age_formula_tau <- config$age_formula_tau  # Optional
sample_per_year <- config$sample_per_year %||% 9
centiles_to_get <- config$centiles %||% c(2.5, 5, 25, 50, 75, 95, 97.5)
max_iterations <- config$max_iterations %||% 1000
convergence_crit <- config$convergence_criterion %||% 0.001

# Create output directory
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# ==== Setup logging ====
log_file <- file.path(output_dir, "fitting_log.txt")
log_conn <- file(log_file, open = "at")  # Append mode

# Function to log messages to both console and file
log_msg <- function(...) {
  msg <- paste0(..., collapse = "")
  cat(msg, "\n")
  cat(msg, "\n", file = log_conn)
  flush(log_conn)
}

# Start logging
log_msg(rep("=", 70))
log_msg("GAMLSS MODEL FITTING LOG")
log_msg(rep("=", 70))
log_msg("Start time: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S"))
log_msg("Random seed: ", seed)
log_msg("")

# Check if model already exists
model_path <- file.path(output_dir, "fitted_model.rds")
model_exists <- file.exists(model_path)

if (skip_fit && model_exists) {
  log_msg("--skip-fit flag detected and model exists")
  log_msg("Loading existing model from: ", model_path)
  
  saved_obj <- readRDS(model_path)
  m <- saved_obj$model
  
  log_msg("Model loaded successfully")
  
  # Skip to post-processing
  df <- read.csv(input_data)
  df$site_scanner <- as.factor(df$site_scanner)
  df$numeric_sex <- as.factor(df$numeric_sex)
  
} else {
  if (skip_fit && !model_exists) {
    log_msg("--skip-fit flag detected but model does not exist")
    log_msg("Proceeding with model fitting...")
  }
  
  # Log configuration
  log_msg("Configuration:")
  log_msg("  Config file: ", config_path)
  log_msg("  Target variable: ", target_var)
  log_msg("  Distribution: ", distribution)
  log_msg("  Output directory: ", output_dir)
  log_msg("")
  
  # ==== Read data ====
  log_msg("Loading data...")
  start_time <- Sys.time()
  
  tryCatch({
    df <- read.csv(input_data)
    log_msg("  Successfully loaded: ", input_data)
    log_msg("  Number of observations: ", nrow(df))
  }, error = function(e) {
    log_msg("ERROR loading data: ", conditionMessage(e))
    close(log_conn)
    stop(e)
  })
  
  # Make sure factors are correct
  df$site_scanner <- as.factor(df$site_scanner)
  df$numeric_sex <- as.factor(df$numeric_sex)
  
  # Check target variable exists
  if (!target_var %in% colnames(df)) {
    log_msg("ERROR: Target variable '", target_var, "' not found in data")
    close(log_conn)
    stop(paste("Target variable", target_var, "not found in data"))
  }
  
  # ==== Fit GAMLSS model ====
  log_msg(rep("-", 70))
  log_msg("Fitting GAMLSS model...")
  
  # Construct formulas
  mu_formula <- as.formula(paste0(target_var, " ~ ", age_formula_mu, " + numeric_sex + random(site_scanner)"))
  sigma_formula <- as.formula(paste0("~ ", age_formula_sigma, " + numeric_sex + random(site_scanner)"))
  
  nu_formula <- NULL
  tau_formula <- NULL
  
  if (!is.null(age_formula_nu) && age_formula_nu != "") {
    if (age_formula_nu == "1") {
      nu_formula <- ~ 1
    } else {
      nu_formula <- as.formula(paste0("~ ", age_formula_nu, " + numeric_sex + random(site_scanner)"))
    }
  }
  
  if (!is.null(age_formula_tau) && age_formula_tau != "") {
    if (age_formula_tau == "1") {
      tau_formula <- ~ 1
    } else {
      tau_formula <- as.formula(paste0("~ ", age_formula_tau, " + numeric_sex + random(site_scanner)"))
    }
  }
  
  # Try to get distribution family
  dist_family <- tryCatch({
    get(distribution)()
  }, error = function(e) {
    log_msg("ERROR: Distribution '", distribution, "' not found or invalid.")
    close(log_conn)
    stop(e)
  })
  
  fit_start <- Sys.time()
  
  # Redirect console output to file
  console_output_file <- file.path(output_dir, "gamlss_console_output.txt")
  console_conn <- file(console_output_file, open = "wt")
  sink(console_conn)
  sink(console_conn, type = "message")
  
  m <- tryCatch({
    gamlss(
      formula = mu_formula,
      sigma.formula = sigma_formula,
      nu.formula = nu_formula,
      tau.formula = tau_formula,
      family = dist_family,
      data = df,
      method = RS(),
      control = gamlss.control(n.cyc = max_iterations, c.crit = convergence_crit, trace = TRUE)
    )
  }, error = function(e) {
    sink(); sink(type = "message"); close(console_conn)
    log_msg("ERROR during model fitting: ", conditionMessage(e))
    close(log_conn)
    stop(e)
  }, finally = {
    sink(); sink(type = "message"); close(console_conn)
  })
  
  # Ensure training data and stable call components are embedded
  m$data <- df
  if (!is.null(m$call)) {
    m$call$data <- m$data
    m$call$family <- as.name(distribution)
  }

  fit_end <- Sys.time()
  fit_duration <- difftime(fit_end, fit_start, units = "secs")
  
  log_msg("Model fitting completed in: ", round(fit_duration, 2), " seconds")
  log_msg("Converged: ", m$converged)
  log_msg("Global Deviance: ", m$G.deviance)
  log_msg("AIC: ", AIC(m))
  log_msg("BIC: ", BIC(m))
  log_msg("")
  
  # ==== Save the fitted model ====
  log_msg("Saving fitted model...")
  saveRDS(list(
    model = m,
    config = config,
    target_var = target_var,
    distribution = distribution,
    fit_time = fit_duration,
    fit_date = format(Sys.time(), "%Y-%m-%d %H:%M:%S")
  ), model_path)
}

# ==== Perform site correction and generate outputs ====
log_msg(rep("-", 70))
log_msg("Generating predictions and site corrections...")

results <- tryCatch({
  process_gamlss_results(m, df, target_var, centiles_to_get, sample_per_year)
}, error = function(e) {
  log_msg("ERROR during post-processing: ", conditionMessage(e))
  close(log_conn)
  stop(e)
})

# ==== Save outputs ====
log_msg("Saving outputs...")

outputs_path <- file.path(output_dir, "gamlss_outputs.csv")
male_path <- file.path(output_dir, "male_population_centiles.csv")
female_path <- file.path(output_dir, "female_population_centiles.csv")

write.csv(results$data_with_predictions, outputs_path, row.names = FALSE)
write.csv(results$male_centiles, male_path, row.names = FALSE)
write.csv(results$female_centiles, female_path, row.names = FALSE)

# ==== Save Detailed Model Summary ====
summary_path <- file.path(output_dir, "model_summary.txt")
sink(summary_path)
cat("GAMLSS Model Summary\n")
cat("====================\n\n")
cat("Fitted:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
cat("Configuration:", config_path, "\n")
cat("Distribution:", distribution, "\n")
cat("Target Variable:", target_var, "\n\n")
print(summary(m))
cat("\n\nDistribution Information\n")
cat("========================\n")
print(results$dist_info)
cat("\n\nModel Fit Statistics\n")
cat("====================\n")
cat("Global Deviance:", m$G.deviance, "\n")
cat("AIC:", AIC(m), "\n")
cat("BIC:", BIC(m), "\n")
cat("Number of observations:", m$N, "\n")
cat("Degrees of freedom (mu):", m$mu.df, "\n")
cat("Degrees of freedom (sigma):", m$sigma.df, "\n")
if (results$dist_info$nopar >= 3) cat("Degrees of freedom (nu):", m$nu.df, "\n")
if (results$dist_info$nopar >= 4) cat("Degrees of freedom (tau):", m$tau.df, "\n")
sink()

log_msg("  Model summary saved to: ", summary_path)

# ==== Optional: Apply to new data ====
if (!is.null(apply_input_path) && !is.null(apply_output_path)) {
  log_msg(rep("-", 70))
  log_msg("Applying model to external dataset...")
  
  tryCatch({
    apply_df <- read.csv(apply_input_path)
    apply_df$numeric_sex <- as.factor(apply_df$numeric_sex)
    apply_df$site_scanner <- as.factor(apply_df$site_scanner)
    
    app_results <- run_gamlss_application(m, apply_df, target_var, results$dist_info)
    
    app_out_dir <- dirname(apply_output_path)
    if (!dir.exists(app_out_dir)) dir.create(app_out_dir, recursive = TRUE)
    write.csv(app_results$result_df, apply_output_path, row.names = FALSE)
    
    log_msg("  Application successful.")
    log_msg("  Mean correction: ", round(app_results$mean_correction, 4))
  }, error = function(e) {
    log_msg("ERROR during application: ", conditionMessage(e))
  })
}

# ==== Final summary ====
log_msg(rep("=", 70))
log_msg("FITTING COMPLETED SUCCESSFULLY")
log_msg(rep("=", 70))

total_time <- difftime(Sys.time(), start_time, units = "mins")
log_msg("Total execution time: ", round(total_time, 2), " minutes")
close(log_conn)
