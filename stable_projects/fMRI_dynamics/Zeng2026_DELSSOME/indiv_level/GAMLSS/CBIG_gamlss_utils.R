# Written by Fang Tian, Tianchu Zeng and CBIG under MIT license:
# https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md

# GAMLSS Utility Functions
# Refactored for modularity and robustness

library(gamlss)
library(splines)

# ==== Distribution Helpers ====

get_distribution_info <- function(model) {
  family_name <- as.character(model$family[1])
  family_func <- get(family_name)()
  
  list(
    family_name = family_name,
    nopar = family_func$nopar,
    mu_link = model$mu.link,
    sigma_link = model$sigma.link,
    nu_link = if (family_func$nopar >= 3) model$nu.link else NULL,
    tau_link = if (family_func$nopar >= 4) model$tau.link else NULL,
    type = family_func$type
  )
}

# CORRECTED: This now performs the INVERSE link (Link Scale -> Response Scale)
apply_inverse_link <- function(eta, link_name) {
  if (is.null(link_name)) return(NULL)
  switch(link_name,
    "identity" = eta,
    "log" = exp(eta),
    "logit" = plogis(eta), # 1 / (1 + exp(-eta))
    "probit" = pnorm(eta),
    "inverse" = 1 / eta,
    "sqrt" = eta^2,
    "1/mu^2" = 1 / sqrt(eta),
    "log(mu)" = exp(eta), # Alias for log
    "own" = eta,
    eta # Default
  )
}

# Unified function for PDF (p) and Quantile (q) calculations
get_dist_val <- function(func_type, val, mu, sigma, nu, tau, dist_info) {
  # func_type: "p" (CDF) or "q" (Quantile)
  fun_name <- paste0(func_type, dist_info$family_name)
  
  args <- list(mu = mu)
  if (func_type == "q") args$p <- val else args$q <- val
  
  if (dist_info$nopar >= 2) args$sigma <- sigma
  if (dist_info$nopar >= 3) args$nu <- nu
  if (dist_info$nopar >= 4) args$tau <- tau
  
  do.call(fun_name, args)
}

# ==== Prediction Logic ====

# Helper to create predictions without site effects (Reference Curves)
predict_without_site <- function(model, data, dist_info) {
  
  # 1. Handle New Sites: Map unknown sites to a known training site
  training_sites <- model$xlevels[["site_scanner"]]
  if (is.null(training_sites) && !is.null(model$data)) {
     training_sites <- levels(model$data$site_scanner)
  }
  
  safe_data <- data
  dummy_site <- training_sites[1]
  
  is_new_site <- !(safe_data$site_scanner %in% training_sites)
  if (any(is_new_site)) {
    safe_data$site_scanner[is_new_site] <- dummy_site
    safe_data$site_scanner <- factor(safe_data$site_scanner, levels = training_sites)
  }

  # Helper to get fixed effect part from terms
  get_fixed_pred <- function(what_par) {
    tryCatch({
      # Get terms on link scale
      terms_mat <- predict(model, newdata = safe_data, what = what_par, type = "terms")
      re_cols <- grep("random\\(", colnames(terms_mat))
      
      # Sum fixed terms + intercept to get Linear Predictor (eta)
      if (length(re_cols) > 0) {
        eta <- rowSums(terms_mat[, -re_cols, drop = FALSE]) + attr(terms_mat, "constant")
      } else {
        eta <- rowSums(terms_mat) + attr(terms_mat, "constant")
      }
      
      # Apply inverse link to get Response Scale (mu, sigma, etc.)
      link_name <- dist_info[[paste0(what_par, "_link")]]
      return(apply_inverse_link(eta, link_name))
      
    }, error = function(e) {
      # Fallback for constant parameters (e.g. nu ~ 1)
      coefs <- coef(model, what = what_par)
      if (length(coefs) == 1) {
        link_name <- dist_info[[paste0(what_par, "_link")]]
        val <- apply_inverse_link(coefs, link_name)
        return(rep(val, nrow(data))) 
      }
      return(NULL)
    })
  }

  list(
    mu = get_fixed_pred("mu"),
    sigma = get_fixed_pred("sigma"),
    nu = if (dist_info$nopar >= 3) get_fixed_pred("nu") else NULL,
    tau = if (dist_info$nopar >= 4) get_fixed_pred("tau") else NULL
  )
}

# Main wrapper: Handles Safe Data, Reference Preds, Site Preds, and Harmonization
predict_and_harmonize <- function(model, data, target_var, dist_info) {
  
  # 1. Prepare Safe Data (Map new sites to dummy to prevent crashes)
  training_sites <- levels(model$data$site_scanner)
  is_new_site <- !(data$site_scanner %in% training_sites)
  
  safe_data <- data
  # Ensure columns match training data intersection
  cols_needed <- intersect(names(safe_data), names(model$data))
  safe_data <- safe_data[, cols_needed, drop = FALSE]
  
  if (any(is_new_site)) {
    safe_data$site_scanner <- as.character(safe_data$site_scanner)
    safe_data$site_scanner[is_new_site] <- training_sites[1]
  }
  safe_data$site_scanner <- factor(safe_data$site_scanner, levels = training_sites)
  
  # 2. Get Reference Predictions (No Site Effects)
  preds_no_site <- predict_without_site(model, safe_data, dist_info)
  
  data$pred_mu_no_site <- preds_no_site$mu
  data$pred_sigma_no_site <- preds_no_site$sigma
  if (!is.null(preds_no_site$nu)) data$pred_nu_no_site <- preds_no_site$nu
  if (!is.null(preds_no_site$tau)) data$pred_tau_no_site <- preds_no_site$tau
  
  # 3. Get Site-Specific Predictions
  get_site_pred <- function(what) {
    if (what == "nu" && dist_info$nopar < 3) return(NULL)
    if (what == "tau" && dist_info$nopar < 4) return(NULL)
    # Note: predict() with random() effects might trigger warnings about re-fitting.
    # This is expected behavior in GAMLSS when predicting on new data structures.
    predict(model, newdata = safe_data, what = what, type = "response")
  }
  
  data$pred_mu_with_site <- get_site_pred("mu")
  data$pred_sigma_with_site <- get_site_pred("sigma")
  data$pred_nu_with_site <- get_site_pred("nu")
  data$pred_tau_with_site <- get_site_pred("tau")
  
  # 4. Overwrite New Sites with Population Average
  if (any(is_new_site)) {
    data$pred_mu_with_site[is_new_site] <- data$pred_mu_no_site[is_new_site]
    data$pred_sigma_with_site[is_new_site] <- data$pred_sigma_no_site[is_new_site]
    if (!is.null(data$pred_nu_with_site)) data$pred_nu_with_site[is_new_site] <- data$pred_nu_no_site[is_new_site]
    if (!is.null(data$pred_tau_with_site)) data$pred_tau_with_site[is_new_site] <- data$pred_tau_no_site[is_new_site]
  }
  
  # 5. Calculate Centiles and Harmonized Values
  data$centile_with_site <- get_dist_val("p", data[[target_var]], 
                                         data$pred_mu_with_site, data$pred_sigma_with_site, 
                                         data$pred_nu_with_site, data$pred_tau_with_site, dist_info)
  
  data$corrected_value <- get_dist_val("q", data$centile_with_site, 
                                       data$pred_mu_no_site, data$pred_sigma_no_site, 
                                       data$pred_nu_no_site, data$pred_tau_no_site, dist_info)
  
  data$residuals_raw_with_site <- data[[target_var]] - data$pred_mu_with_site
  
  return(data)
}

# ==== High-Level Processing ====

process_gamlss_results <- function(model, data, target_var, centiles_to_get, sample_per_year = 9) {
  
  dist_info <- get_distribution_info(model)
  
  # 1. Generate all predictions and harmonized values
  data <- predict_and_harmonize(model, data, target_var, dist_info)
  
  # 2. Add Z-scores (only available for training data via residuals)
  if (nrow(data) == nrow(model$data)) {
    data$residuals_z_with_site <- residuals(model, what = "z-scores")
  }
  
  # 3. Generate Population Centile Curves (Reference)
  # We use the 'no_site' predictions already in 'data' to calculate specific centile lines
  for (p in centiles_to_get) {
    cname <- paste0("C", p, "_no_site")
    data[[cname]] <- get_dist_val("q", p/100, 
                                  data$pred_mu_no_site, data$pred_sigma_no_site, 
                                  data$pred_nu_no_site, data$pred_tau_no_site, dist_info)
  }
  
  # 4. Generate Sex-Specific Growth Curves
  age_seq <- create_age_sequence(data$age_in_years, sample_per_year)
  
  male_centiles <- generate_centiles_by_sex(model, 0, "male", age_seq, centiles_to_get, dist_info)
  female_centiles <- generate_centiles_by_sex(model, 1, "female", age_seq, centiles_to_get, dist_info)
  
  list(
    data_with_predictions = data,
    male_centiles = male_centiles,
    female_centiles = female_centiles,
    dist_info = dist_info
  )
}

# Helper function to create age sequence
create_age_sequence <- function(ages, sample_per_year = 9) {
  min_age <- floor(min(ages))
  max_age <- ceiling(max(ages))
  
  # Create integer years
  integer_years <- seq(min_age, max_age, by = 1)
  
  # Create sub-year intervals
  step_size <- 1 / (sample_per_year + 1)
  
  age_seq <- c()
  for (i in 1:(length(integer_years) - 1)) {
    year_start <- integer_years[i]
    year_end <- integer_years[i + 1]
    
    # Add points between year_start and year_end (excluding year_end to avoid duplication)
    sub_points <- seq(year_start, year_end - step_size, by = step_size)
    age_seq <- c(age_seq, sub_points)
  }
  
  # Add the final integer year
  age_seq <- c(age_seq, max_age)
  
  return(age_seq)
}

generate_centiles_by_sex <- function(model, sex_value, sex_label, age_seq, centiles_to_get, dist_info) {
  
  pred_data <- data.frame(
    age_in_years = age_seq,
    numeric_sex = factor(rep(sex_value, length(age_seq)), levels = c(0, 1))
  )
  
  # Add dummy site for predict() compatibility
  if (!is.null(model$xlevels[["site_scanner"]])) {
    pred_data$site_scanner <- factor(model$xlevels[["site_scanner"]][1], levels = model$xlevels[["site_scanner"]])
  } else if (!is.null(model$data)) {
    pred_data$site_scanner <- model$data$site_scanner[1]
  }

  preds <- predict_without_site(model, pred_data, dist_info)
  
  centile_df <- data.frame(age_in_years = age_seq, sex = sex_label, mu = preds$mu, sigma = preds$sigma)
  if (!is.null(preds$nu)) centile_df$nu <- preds$nu
  if (!is.null(preds$tau)) centile_df$tau <- preds$tau
  
  for (p in centiles_to_get) {
    centile_df[[paste0("C", p)]] <- get_dist_val("q", p/100, preds$mu, preds$sigma, preds$nu, preds$tau, dist_info)
  }
  
  # Calculate growth rate (derivative of median)
  median_col <- if ("C50" %in% colnames(centile_df)) centile_df$C50 else 
                get_dist_val("q", 0.5, preds$mu, preds$sigma, preds$nu, preds$tau, dist_info)
  
  centile_df$growth_rate <- calculate_growth_rate(age_seq, median_col)
  
  return(centile_df)
}

calculate_growth_rate <- function(age_seq, centile_values) {
  n <- length(age_seq)
  growth_rate <- rep(NA_real_, n)
  
  # Central differences
  for (i in 2:(n-1)) {
    growth_rate[i] <- (centile_values[i+1] - centile_values[i-1]) / (age_seq[i+1] - age_seq[i-1])
  }
  # Edges
  if (n >= 2) {
    growth_rate[1] <- (centile_values[2] - centile_values[1]) / (age_seq[2] - age_seq[1])
    growth_rate[n] <- (centile_values[n] - centile_values[n-1]) / (age_seq[n] - age_seq[n-1])
  }
  return(growth_rate)
}

# ==== Application Wrapper ====

run_gamlss_application <- function(model, new_data, target_var, dist_info) {
  # Check for new sites
  training_sites <- levels(model$data$site_scanner)
  num_new_sites <- sum(!(new_data$site_scanner %in% training_sites))
  
  if (num_new_sites > 0) {
    message(paste("  WARNING: Found", num_new_sites, "observations from sites not in training data."))
    message("  These will be harmonized using population-average effects.")
  }
  
  # Call the unified function
  result_df <- predict_and_harmonize(model, new_data, target_var, dist_info)
  
  # Calculate correction summary
  correction <- result_df$corrected_value - result_df[[target_var]]
  
  list(
    result_df = result_df,
    mean_correction = mean(correction, na.rm=TRUE),
    original_mean = mean(result_df[[target_var]], na.rm=TRUE),
    corrected_mean = mean(result_df$corrected_value, na.rm=TRUE)
  )
}
