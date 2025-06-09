# MCMC K-Means Clustering Analysis

# --- 1. Load Libraries ---
library(dplyr)
library(readr)
# Add any other necessary libraries here, e.g., for specific distributions or optimization
# library(DistributionFitR) # May be useful for fitting distributions to CI


# --- 2. Configuration ---
N_MCMC_ITERATIONS <- 5000 # Number of MCMC iterations (Increased for robust run)
N_PRESYMP_SAMPLES <- 5000 # Number of samples for presymptomatic proportion estimation (Increased for robust run)
set.seed(123) # For reproducibility


# --- 3. Load Data ---
# Assuming the script is run from the 'archetypes' directory or the path is relative to it
params_df <- read_csv("Clustering/mcmc/Kmeans/pathogen_params_kmeans.csv")

# Optional: Quick check of the data
# print(head(params_df))
# print(str(params_df))


# --- 4. Helper Functions ---

# Function to derive parameters for rbeta from mean and 95% CI
get_beta_params_from_mean_ci <- function(mean_val, lower_ci, upper_ci, n_eff_guess = 1000) {
    # This is a simplified approach. A robust version would use optimization.
    alpha <- mean_val * n_eff_guess
    beta <- (1 - mean_val) * n_eff_guess
    
    # Ensure positive parameters, using a small positive as a fallback
    min_param_val <- 1e-6 # Small positive value
    if (is.na(alpha) || !is.finite(alpha) || alpha <= 0) alpha <- min_param_val 
    if (is.na(beta) || !is.finite(beta) || beta <= 0) beta <- min_param_val

    # Warn if inputs were problematic for a more robust method
    if (mean_val <= 0 || mean_val >= 1 || lower_ci < 0 || upper_ci > 1 || lower_ci >= upper_ci) {
        warning(paste("Input values for get_beta_params_from_mean_ci may be problematic for robust optimization. Mean:", mean_val, "CI: [", lower_ci, ",", upper_ci, "]"))
    }
    return(list(shape1 = alpha, shape2 = beta))
}

# Function to derive parameters for rgamma from mean and 95% CI
get_gamma_params_from_mean_ci <- function(mean_val, lower_ci, upper_ci) {
    # This is a simplified approach. A robust version would use optimization.
    min_param_val <- 1e-6 # Small positive value

    if (is.na(mean_val) || mean_val <= 0) {
        warning(paste("Mean value for get_gamma_params_from_mean_ci is NA or non-positive:", mean_val, ". Attempting to use CI midpoint."))
        if(!is.na(lower_ci) && !is.na(upper_ci) && lower_ci < upper_ci && (lower_ci + upper_ci)/2 > 0 ){
            mean_val <- (lower_ci + upper_ci) / 2
        } else {
            warning("Cannot robustly determine mean for Gamma. Returning minimal positive parameters.")
            return(list(shape = min_param_val, rate = min_param_val)) 
        }
    }
    
    sd_from_ci <- (upper_ci - lower_ci) / (2 * qnorm(0.975))
    if (is.na(sd_from_ci) || !is.finite(sd_from_ci) || sd_from_ci <= 0) {
       sd_from_ci <- mean_val * 0.5 # fallback SD
       if (sd_from_ci <= 0) sd_from_ci <- min_param_val # ensure positive
       warning(paste("SD from CI was NA or non-positive for Gamma params. Fallback SD used:", sd_from_ci))
    }

    shape <- (mean_val / sd_from_ci)^2
    rate <- mean_val / sd_from_ci^2
    
    if (is.na(shape) || !is.finite(shape) || shape <= 0) shape <- min_param_val
    if (is.na(rate) || !is.finite(rate) || rate <= 0) rate <- min_param_val
    
    return(list(shape = shape, rate = rate))
}


# Function to sample a value for a _Clust_ parameter (R0, SI, CFR)
sample_clust_parameter <- function(pathogen_row, param_prefix) {
  val_col <- paste0(param_prefix, "_Clust_ReportedValue")
  uncert_col <- paste0(param_prefix, "_Clust_UncertaintyType")
  low_col <- paste0(param_prefix, "_Clust_LowerBound")
  high_col <- paste0(param_prefix, "_Clust_UpperBound")
  dist_col <- paste0(param_prefix, "_Clust_SamplingDist")
  cv_col <- paste0(param_prefix, "_Clust_PointEstimateCV")

  val <- as.numeric(pathogen_row[[val_col]])
  uncert_type <- pathogen_row[[uncert_col]]
  sampling_dist <- pathogen_row[[dist_col]]
  lower_bound <- as.numeric(pathogen_row[[low_col]])
  upper_bound <- as.numeric(pathogen_row[[high_col]])
  
  if (is.na(uncert_type) || uncert_type == "") { 
      if (!is.na(val)) return(val)
      stop(paste("Missing value and uncertainty for", param_prefix, pathogen_row$Pathogen_Name))
  }

  sampled_val <- NA

  if (uncert_type == "PointEstimate") {
    cv <- as.numeric(pathogen_row[[cv_col]])
    if (is.na(cv) || cv == 0) return(val)
    sd_val <- val * cv
    if (sd_val <= 0) { # if val is 0 (e.g. CCHF R0 can be near 0), sd_val can be 0 or negative if val is negative
        warning(paste("Calculated SD for PointEstimate is non-positive for", param_prefix, pathogen_row$Pathogen_Name, "Value:", val, "CV:", cv, ". Using reported value."))
        return(val)
    }

    if (sampling_dist == "rnorm" || is.na(sampling_dist) || sampling_dist == "") {
        sampled_val <- rnorm(1, mean = val, sd = sd_val)
        if (param_prefix == "CFR") sampled_val <- pmax(0, pmin(1, sampled_val))
        else if (param_prefix %in% c("R0", "SI")) sampled_val <- pmax(0.00001, sampled_val) # Ensure positive, avoid exactly zero
    } else if (sampling_dist == "rbeta" && param_prefix == "CFR") {
        m <- val
        s <- sd_val
        if (m > 0 && m < 1 && s > 0 && (m * (1 - m) > s^2) && ((m * (1 - m) / s^2) - 1) > 0) {
            alpha <- m * (((m * (1 - m)) / s^2) - 1)
            beta <- (1 - m) * (((m * (1 - m)) / s^2) - 1)
            if (alpha > 0 && beta > 0) {
              sampled_val <- rbeta(1, shape1 = alpha, shape2 = beta)
            } else { 
              warning(paste("Derived alpha/beta for rbeta PointEstimate not positive for", param_prefix, pathogen_row$Pathogen_Name, ". Using reported value."))
              sampled_val <- val 
            }
        } else { 
          warning(paste("Mean/SD for rbeta PointEstimate not suitable for", param_prefix, pathogen_row$Pathogen_Name, "Mean:",m,", SD:",s, ". Using reported value."))
          sampled_val <- val 
        }
    } else {
        warning(paste("Unhandled sampling_dist for PointEstimate:", sampling_dist, "for", param_prefix, pathogen_row$Pathogen_Name))
        sampled_val <- val
    }
  } else if (uncert_type %in% c("95CI", "IQR", "Range", "range")) {
    if (is.na(lower_bound) || is.na(upper_bound)) {
        warning(paste("Missing bounds for", uncert_type, param_prefix, pathogen_row$Pathogen_Name, "- using reported value."))
        return(val)
    }
    if (lower_bound > upper_bound) { 
        warning(paste("Lower bound > Upper bound for", param_prefix, pathogen_row$Pathogen_Name, ". Swapping them."))
        tmp <- lower_bound; lower_bound <- upper_bound; upper_bound <- tmp;
    }

    if (sampling_dist == "rlnorm") {
      current_lower_bound <- lower_bound
      if (current_lower_bound <= 0) { 
          warning(paste("Lower bound for rlnorm is non-positive (", current_lower_bound, ") for", param_prefix, pathogen_row$Pathogen_Name, ". Adjusting to 0.00001."));
          current_lower_bound <- 0.00001;
      }
      if (upper_bound <= current_lower_bound) { # Check if upper bound is not greater than (adjusted) lower bound
           warning(paste("Upper bound (", upper_bound, ") not greater than lower bound (", current_lower_bound, ") for rlnorm for", param_prefix, pathogen_row$Pathogen_Name, ". Using reported value."));
           return(val);
      }
      log_lower <- log(current_lower_bound)
      log_upper <- log(upper_bound)
      sdlog <- (log_upper - log_lower) / (2 * qnorm(0.975))
      meanlog <- (log_upper + log_lower) / 2 
      if(is.na(sdlog) || !is.finite(sdlog) || sdlog <= 0) { 
          warning(paste("Calculated sdlog is non-positive/non-finite for rlnorm for", param_prefix, pathogen_row$Pathogen_Name, ". sdlog:", sdlog, ". Using reported value."));
          return(val);
      }
      sampled_val <- rlnorm(1, meanlog = meanlog, sdlog = sdlog)

    } else if (sampling_dist == "rgamma") {
      # Using val (ReportedValue) as the mean for get_gamma_params_from_mean_ci
      if(is.na(val)) { 
          warning(paste("ReportedValue (mean) is NA for rgamma for", param_prefix, pathogen_row$Pathogen_Name, ". Using midpoint of CI if available, else failing."))
          if(!is.na(lower_bound) && !is.na(upper_bound)) val <- (lower_bound + upper_bound) / 2 else return(NA) # or some other fallback
      }
      gamma_params <- get_gamma_params_from_mean_ci(val, lower_bound, upper_bound)
      if(is.na(gamma_params$shape) || is.na(gamma_params$rate) || gamma_params$shape <=0 || gamma_params$rate <=0 || !is.finite(gamma_params$shape) || !is.finite(gamma_params$rate)){
          warning(paste("Derived gamma params are NA/non-positive/non-finite for", param_prefix, pathogen_row$Pathogen_Name, ". Shape:", gamma_params$shape, "Rate:", gamma_params$rate, ". Using reported value."))
          return(val) # return original val, not NA, as a fallback if params are bad
      }
      sampled_val <- rgamma(1, shape = gamma_params$shape, rate = gamma_params$rate)
      if (sampled_val <= 0 && val > 0) sampled_val <- val # Fallback if sampling gives non-positive but original val was positive
      else if (sampled_val <= 0) sampled_val <- 0.00001 # Ensure positive if original was also non-positive

    } else if (sampling_dist == "rbeta") { 
      if(is.na(val)) { 
          warning(paste("ReportedValue (mean) is NA for rbeta for", param_prefix, pathogen_row$Pathogen_Name, ". Using midpoint of CI if available, else failing."))
          if(!is.na(lower_bound) && !is.na(upper_bound)) val <- (lower_bound + upper_bound) / 2 else return(NA) 
      }
      if(val < 0 || val > 1) { 
          warning(paste("ReportedValue (mean) for rbeta is outside [0,1] for", param_prefix, pathogen_row$Pathogen_Name, ":", val, ". Clamping to [0.00001, 0.99999]."))
          val <- pmax(0.00001, pmin(0.99999, val))
      }
       clamped_lower <- pmax(0, lower_bound)
       clamped_upper <- pmin(1, upper_bound)
       if (clamped_lower >= clamped_upper) {
            warning(paste("Clamped bounds for rbeta are invalid for", param_prefix, pathogen_row$Pathogen_Name, ". Lower:",clamped_lower,"Upper:",clamped_upper,". Using reported value."))
            return(val)
       }

      beta_params <- get_beta_params_from_mean_ci(val, clamped_lower, clamped_upper)
      if(is.na(beta_params$shape1) || is.na(beta_params$shape2) || beta_params$shape1 <=0 || beta_params$shape2 <=0 || !is.finite(beta_params$shape1) || !is.finite(beta_params$shape2)){
          warning(paste("Derived beta params are NA/non-positive/non-finite for", param_prefix, pathogen_row$Pathogen_Name, ". Shape1:", beta_params$shape1, "Shape2:", beta_params$shape2, ". Using reported value."))
          return(val)
      }
      sampled_val <- rbeta(1, shape1 = beta_params$shape1, shape2 = beta_params$shape2)

    } else if (sampling_dist == "rnorm") {
      sd_val <- (upper_bound - lower_bound) / (2 * qnorm(0.975))
      if(is.na(sd_val) || !is.finite(sd_val) || sd_val <= 0){
           warning(paste("Calculated SD for rnorm is NA/non-positive/non-finite for", param_prefix, pathogen_row$Pathogen_Name, ". SD:", sd_val, ". Using reported value."));
           return(val);
      }
      sampled_val <- rnorm(1, mean = val, sd = sd_val)
      if (param_prefix == "CFR") sampled_val <- pmax(0, pmin(1, sampled_val))
      else if (param_prefix %in% c("R0", "SI")) sampled_val <- pmax(0.00001, sampled_val)

    } else if (sampling_dist == "runif") {
      if (lower_bound >= upper_bound) { # runif gives error if min == max
            warning(paste("Min >= Max for runif for", param_prefix, pathogen_row$Pathogen_Name, ". Returning min value."))
            return(lower_bound)
      }
      sampled_val <- runif(1, min = lower_bound, max = upper_bound)
    } else {
      warning(paste("Unhandled sampling distribution:", sampling_dist, "for", param_prefix, pathogen_row$Pathogen_Name))
      sampled_val <- val 
    }
  } else {
    warning(paste("Unhandled uncertainty type:", uncert_type, "for", param_prefix, pathogen_row$Pathogen_Name))
    sampled_val <- val 
  }
  
  if (is.na(sampled_val)) {
    warning(paste("Sampled value is NA for", param_prefix, pathogen_row$Pathogen_Name, "(Uncertainty type:", uncert_type, ", Sampling dist:", sampling_dist, ") - returning reported value."))
    return(val)
  }
  return(sampled_val)
}


get_full_dist_params <- function(pathogen_row, type_prefix) { 
  family <- pathogen_row[[paste0(type_prefix, "_FullDist_Family")]]
  raw_params <- list()

  for (i in 1:2) { 
    param_name_col <- paste0(type_prefix, "_FullDist_Param", i, "_Name")
    param_val_col <- paste0(type_prefix, "_FullDist_Param", i, "_Value")
    param_uncert_col <- paste0(type_prefix, "_FullDist_Param", i, "_UncertaintySource")
    
    if (!is.na(pathogen_row[[param_name_col]]) && pathogen_row[[param_name_col]] != "") {
      param_name_from_csv <- pathogen_row[[param_name_col]]
      param_val <- as.numeric(pathogen_row[[param_val_col]])
      param_uncert <- pathogen_row[[param_uncert_col]]

      if (is.na(param_uncert) || param_uncert == "Fixed" || param_uncert == "") {
        raw_params[[param_name_from_csv]] <- param_val
      } else {
        if (param_uncert == "95CI") {
             param_lower <- as.numeric(pathogen_row[[paste0(type_prefix, "_FullDist_Param", i, "_LowerBound")]])
             param_upper <- as.numeric(pathogen_row[[paste0(type_prefix, "_FullDist_Param", i, "_UpperBound")]])
             param_sample_dist <- pathogen_row[[paste0(type_prefix, "_FullDist_Param", i, "_SamplingDist")]]
             
             if (!is.na(param_lower) && !is.na(param_upper) && !is.na(param_sample_dist) && param_sample_dist != "") {
                 if (param_sample_dist == "rnorm") {
                     sd_p <- (param_upper - param_lower) / (2 * qnorm(0.975))
                     if(is.na(sd_p) || sd_p <= 0){ 
                         warning(paste("Derived SD for sampling FullDist parameter", param_name_from_csv, "is NA or non-positive. Using reported value."))
                         sampled_param_val <- param_val
                     } else {
                        sampled_param_val <- rnorm(1, mean = param_val, sd = sd_p)
                     }
                     if (param_name_from_csv %in% c("SD", "sd", "shape", "rate", "scale", "sdlog") && sampled_param_val <= 0) {
                         warning(paste("Sampled FullDist parameter", param_name_from_csv, "was non-positive (", sampled_param_val, "), using reported value (", param_val, ") instead."))
                         sampled_param_val <- param_val 
                     }
                     raw_params[[param_name_from_csv]] <- sampled_param_val
                 } else {
                     raw_params[[param_name_from_csv]] <- param_val 
                     warning(paste("Sampling for uncertain FullDist parameter", param_name_from_csv, "with dist", param_sample_dist, "not yet implemented. Using mean."))
                 }
             } else {
                 raw_params[[param_name_from_csv]] <- param_val 
                 warning(paste("Missing info for uncertain FullDist parameter", param_name_from_csv, ". Using mean."))
             }
        } else {
            raw_params[[param_name_from_csv]] <- param_val 
            warning(paste("Unhandled uncertainty source '",param_uncert,"' for FullDist parameter:", param_name_from_csv, ". Using reported value."))
        }
      }
    }
  }

  final_params <- raw_params
  if (!is.na(family)) {
    if (family == "Gamma") {
      if (all(c("mean", "SD") %in% names(raw_params)) && !("shape" %in% names(raw_params))) {
        mean_val <- raw_params$mean
        sd_val <- raw_params$SD
        if (!is.na(mean_val) && !is.na(sd_val) && sd_val > 0 && mean_val > 0) {
          final_params$shape <- (mean_val / sd_val)^2
          final_params$rate <- mean_val / (sd_val^2)
          final_params$mean <- NULL; final_params$SD <- NULL
        } else {
          warning(paste("Cannot convert Gamma mean/SD to shape/rate for", type_prefix, pathogen_row$Pathogen_Name, "due to invalid values. Mean:", mean_val, "SD:", sd_val))
          return(list(family = family, params = list())) 
        }
      }
      if (!is.null(final_params$shape) && (!is.finite(final_params$shape) || final_params$shape <= 1e-9)) { final_params$shape <- 1e-9; warning("Gamma shape adjusted to be positive.") }
      if (!is.null(final_params$rate) && (!is.finite(final_params$rate) || final_params$rate <= 1e-9)) { final_params$rate <- 1e-9; warning("Gamma rate adjusted to be positive.") }
      if (!is.null(final_params$scale) && (!is.finite(final_params$scale) || final_params$scale <= 1e-9)) { final_params$scale <- 1e-9; warning("Gamma scale adjusted to be positive.") }
    
    } else if (family == "Lognormal") {
      if (all(c("mean", "SD") %in% names(raw_params)) && !("meanlog" %in% names(raw_params))) {
        mean_val <- raw_params$mean
        sd_val <- raw_params$SD
        if(!is.na(mean_val) && !is.na(sd_val) && mean_val > 0 && sd_val > 0){
          final_params$sdlog <- sqrt(log(sd_val^2 / mean_val^2 + 1))
          final_params$meanlog <- log(mean_val) - 0.5 * final_params$sdlog^2
          final_params$mean <- NULL; final_params$SD <- NULL
        } else {
           warning(paste("Cannot convert Lognormal mean/SD to meanlog/sdlog for", type_prefix, pathogen_row$Pathogen_Name, "due to invalid values. Mean:", mean_val, "SD:", sd_val))
           return(list(family = family, params = list()))
        }
      }
      if (!is.null(final_params$sdlog) && (!is.finite(final_params$sdlog) || final_params$sdlog <= 1e-9)) { final_params$sdlog <- 1e-9; warning("Lognormal sdlog adjusted to be positive.") }
    
    } else if (family == "Normal") {
      if (!is.null(final_params$SD) && (!is.finite(final_params$SD) || final_params$SD <= 1e-9)) { final_params$SD <- 1e-9; warning("Normal SD adjusted to be positive.") }
    
    } else if (family == "Weibull") {
      if (!is.null(final_params$shape) && (!is.finite(final_params$shape) || final_params$shape <= 1e-9)) { final_params$shape <- 1e-9; warning("Weibull shape adjusted to be positive.") }
      if (!is.null(final_params$scale) && (!is.finite(final_params$scale) || final_params$scale <= 1e-9)) { final_params$scale <- 1e-9; warning("Weibull scale adjusted to be positive.") }
    }
  } else {
      warning(paste("Missing distribution family for", type_prefix, "for pathogen:", pathogen_row$Pathogen_Name))
      return(list(family = NA, params = list()))
  }
  return(list(family = family, params = final_params))
}


estimate_presymptomatic_proportion <- function(si_dist_info, ip_dist_info, n_samples = N_PRESYMP_SAMPLES) {
  if (is.null(si_dist_info$family) || is.na(si_dist_info$family) || length(si_dist_info$params) == 0 ||
      is.null(ip_dist_info$family) || is.na(ip_dist_info$family) || length(ip_dist_info$params) == 0) {
    warning(paste("SI or IP distribution info incomplete for", si_dist_info$Pathogen_Name))
    return(NA)
  }

  # Helper to create a quantile function from dist_info
  get_q_func <- function(d_info) {
    params <- d_info$params
    switch(d_info$family,
      "Lognormal" = function(p) qlnorm(p, meanlog = params$meanlog, sdlog = params$sdlog),
      "Gamma"     = function(p) qgamma(p, shape = params$shape, rate = if (!is.null(params$rate)) params$rate else 1/params$scale),
      "Normal"    = function(p) qnorm(p, mean = params$mean, sd = params$SD),
      "Weibull"   = function(p) qweibull(p, shape = params$shape, scale = params$scale),
      NULL
    )
  }

  q_si <- get_q_func(si_dist_info)
  q_ip <- get_q_func(ip_dist_info)

  if (is.null(q_si) || is.null(q_ip)) {
    warning(paste("Could not get quantile function for", si_dist_info$Pathogen_Name))
    return(NA)
  }

  # --- Correlated Sampling using Quantiles ---
  # This method assumes a positive rank correlation between SI and IP,
  # addressing the unrealistic outcomes from the previous independent sampling method.
  random_quantiles <- runif(n_samples)
  
  si_samples <- q_si(random_quantiles)
  ip_samples <- q_ip(random_quantiles)
  
  # Truncate any negative samples, which can occur with Normal distributions
  if(si_dist_info$family == "Normal") si_samples[si_samples < 0] <- 0
  if(ip_dist_info$family == "Normal") ip_samples[ip_samples < 0] <- 0

  mean(si_samples < ip_samples, na.rm = TRUE)
}


# --- 5. Main MCMC Loop ---
mcmc_results <- list()

for (iter in 1:N_MCMC_ITERATIONS) {
  if (iter %% max(1, (N_MCMC_ITERATIONS/10)) == 0) { # Ensure non-zero divisor for modulo
    print(paste("MCMC Iteration:", iter, "/", N_MCMC_ITERATIONS))
  }
  
  current_iteration_params_list <- list()
  
  for (i in 1:nrow(params_df)) {
    pathogen_row <- params_df[i, ]
    pathogen_name <- pathogen_row$Pathogen_Name
    
    r0_sampled <- sample_clust_parameter(pathogen_row, "R0")
    si_clust_sampled <- sample_clust_parameter(pathogen_row, "SI")
    cfr_sampled <- sample_clust_parameter(pathogen_row, "CFR")
    
    si_full_dist_info <- get_full_dist_params(pathogen_row, "SI")
    ip_full_dist_info <- get_full_dist_params(pathogen_row, "IP")
    
    presymp_prop_sampled <- NA 
    if (!is.null(si_full_dist_info) && !is.null(ip_full_dist_info) && 
        !is.na(si_full_dist_info$family) && !is.na(ip_full_dist_info$family) &&
        length(si_full_dist_info$params) > 0 && length(ip_full_dist_info$params) > 0) {
      
      si_full_dist_info$Pathogen_Name <- pathogen_name # Pass name for warnings
      ip_full_dist_info$Pathogen_Name <- pathogen_name
      presymp_prop_sampled <- estimate_presymptomatic_proportion(si_full_dist_info, ip_full_dist_info)
    } else {
      warning(paste("Could not estimate presymptomatic proportion for", pathogen_name, "due to incomplete SI/IP full dist info."))
    }
    
    pathogen_params_df <- data.frame(
      Pathogen_Name = pathogen_name,
      R0_sampled = ifelse(is.null(r0_sampled) || is.na(r0_sampled) || !is.finite(r0_sampled), NA_real_, r0_sampled),
      SI_Clust_sampled = ifelse(is.null(si_clust_sampled) || is.na(si_clust_sampled) || !is.finite(si_clust_sampled), NA_real_, si_clust_sampled),
      CFR_sampled = ifelse(is.null(cfr_sampled) || is.na(cfr_sampled) || !is.finite(cfr_sampled), NA_real_, cfr_sampled),
      Presymp_Proportion_sampled = ifelse(is.null(presymp_prop_sampled) || is.na(presymp_prop_sampled) || !is.finite(presymp_prop_sampled), NA_real_, presymp_prop_sampled),
      Route_resp = as.integer(pathogen_row$Route_resp),
      Route_direct = as.integer(pathogen_row$Route_direct),
      Route_sexual = as.integer(pathogen_row$Route_sexual),
      Route_animal = as.integer(pathogen_row$Route_animal),
      Route_vector = as.integer(pathogen_row$Route_vector)
    )
    current_iteration_params_list[[pathogen_name]] <- pathogen_params_df
  }
  
  if(length(current_iteration_params_list) > 0){
      mcmc_results[[iter]] <- bind_rows(current_iteration_params_list)
  } else {
      warning(paste("No pathogen data generated for MCMC iteration", iter))
  }
}

if (length(mcmc_results) > 0 && !all(sapply(mcmc_results, is.null))){
    all_mcmc_samples_df <- bind_rows(mcmc_results[!sapply(mcmc_results, is.null)], .id = "MCMC_Iteration") # Filter out NULL iterations
    if(nrow(all_mcmc_samples_df) > 0){
        all_mcmc_samples_df$MCMC_Iteration <- as.integer(all_mcmc_samples_df$MCMC_Iteration)

        print("MCMC sampling complete.")
        print(paste("Generated", length(unique(all_mcmc_samples_df$MCMC_Iteration)), "valid sets of parameters for", length(unique(all_mcmc_samples_df$Pathogen_Name)), "pathogens."))
    
        # --- 6. Post-MCMC Analysis (Clustering) ---
        # Example: Save the samples for now. Clustering can be a separate script or section.

        # --- 7. Save Results ---
        write_csv(all_mcmc_samples_df, "Clustering/mcmc/Kmeans/mcmc_parameter_samples.csv")
        print("MCMC parameter samples saved to Clustering/mcmc/Kmeans/mcmc_parameter_samples.csv")
    } else {
        print("MCMC sampling completed, but no rows in the final combined data frame.")
    }
} else {
    print("MCMC sampling did not produce any results or all iterations resulted in NULL.")
}

print("Script finished.")

# TODO:
# - Refine get_beta_params_from_mean_ci and get_gamma_params_from_mean_ci with robust optimization (e.g., using optim or specific packages).
# - Further test sampling for uncertain _FullDist_Param_Values, especially for distributions other than rnorm.
# - Add more specific error handling or fallbacks in sampling functions if parameters are still problematic.
# - Strategy for NAs in Presymp_Proportion_sampled: current script passes them on. Clustering part will need to handle them (e.g. imputation, exclusion, or use of clustering algorithms robust to NAs).
# - Expand the post-MCMC clustering analysis section (currently only saves samples).
# - Consider using helper functions from 'epitools' or similar for CI to parameter conversion if appropriate for _Clust_ params.

# --- 8. Clustering Analysis (Revised for per-iteration clustering and K=6) ---
if (exists("all_mcmc_samples_df") && nrow(all_mcmc_samples_df) > 0) {
  print("--- Starting Clustering Analysis (Per MCMC Iteration, K=6) ---")
  
  CHOSEN_K <- 6 # User specified K=6
  
  # --- 8.1. Prepare Data for Clustering (from all_mcmc_samples_df) ---
  # Select features that will be used. Pathogen_Name and MCMC_Iteration are for tracking.
  clustering_data_full <- all_mcmc_samples_df %>%
    select(
      MCMC_Iteration,
      Pathogen_Name,
      R0_sampled, 
      SI_Clust_sampled, 
      CFR_sampled, 
      Presymp_Proportion_sampled,
      Route_resp,
      Route_direct,
      Route_sexual,
      Route_animal,
      Route_vector
    )

  # Handle NAs in the sampled parameters (e.g., if Presymp_Proportion_sampled was NA for some)
  # Impute with the overall median for that parameter column
  numerical_param_cols <- c("R0_sampled", "SI_Clust_sampled", "CFR_sampled", "Presymp_Proportion_sampled")
  for(col in numerical_param_cols){
    if(any(is.na(clustering_data_full[[col]]))){
      na_count <- sum(is.na(clustering_data_full[[col]]))
      warning(paste("Column", col, "in all_mcmc_samples_df has", na_count, "NA values. Imputing with overall median before scaling."))
      # Ensure the column is numeric before attempting median calculation if it became all NAs somehow
      if(is.numeric(clustering_data_full[[col]])){
        clustering_data_full[[col]][is.na(clustering_data_full[[col]])] <- median(clustering_data_full[[col]], na.rm = TRUE)
      } else {
        warning(paste("Column", col, "is not numeric or became all NAs, cannot impute with median."))
      }
    }
  }
  
  # Scale numerical features across the entire dataset
  scaled_clustering_data_full <- clustering_data_full
  # Ensure columns to scale are indeed numeric and exist
  cols_to_scale_exist <- numerical_param_cols[numerical_param_cols %in% names(scaled_clustering_data_full)]
  if(length(cols_to_scale_exist) > 0) {
      scaled_clustering_data_full[cols_to_scale_exist] <- scale(scaled_clustering_data_full[cols_to_scale_exist])
  } else {
      warning("No numerical columns found for scaling as specified in numerical_param_cols.")
  }

  print("Scaled features for per-iteration clustering (first 6 rows of full dataset):")
  print(head(scaled_clustering_data_full))

  # Lists to store results from each iteration
  all_iteration_centroids <- list() # Will store UN SCALED centroids
  all_iteration_assignments <- list()

  mcmc_iterations <- unique(scaled_clustering_data_full$MCMC_Iteration)

  # --- 8.2. K-Means Clustering for Each MCMC Iteration ---
  if (length(mcmc_iterations) > 0) {
    for (iter_val in mcmc_iterations) {
      if (as.integer(iter_val) %% max(1, round(length(mcmc_iterations)/10)) == 0) {
        print(paste("Clustering MCMC Iteration:", iter_val, "/", length(mcmc_iterations)))
      }
      
      current_iter_data <- scaled_clustering_data_full %>%
        filter(MCMC_Iteration == iter_val)
      
      # Also get the original (unscaled) data for this iteration for calculating unscaled centroids later
      current_iter_data_unscaled <- clustering_data_full %>%
        filter(MCMC_Iteration == iter_val)

      current_iter_features_scaled <- current_iter_data %>%
        select(-MCMC_Iteration, -Pathogen_Name)
          
      if (nrow(current_iter_features_scaled) < CHOSEN_K) {
          warning(paste("Skipping K-means for iteration", iter_val, "due to insufficient data points (", nrow(current_iter_features_scaled), ") for K=", CHOSEN_K))
          next
      }
      if (ncol(current_iter_features_scaled) == 0) {
          warning(paste("Skipping K-means for iteration", iter_val, "due to no features selected."))
          next
      }

      set.seed(123 + as.integer(iter_val)) 
      kmeans_iter_result <- tryCatch({
          kmeans(current_iter_features_scaled, centers = CHOSEN_K, nstart = 25)
      }, error = function(e) {
          warning(paste("Error in kmeans for iteration", iter_val, ":", e$message, "Skipping this iteration."))
          NULL
      })
      
      if (!is.null(kmeans_iter_result)) {
          assignments_df <- data.frame(
              MCMC_Iteration = iter_val,
              Pathogen_Name = current_iter_data$Pathogen_Name,
              Cluster_Assigned = kmeans_iter_result$cluster
          )
          all_iteration_assignments[[as.character(iter_val)]] <- assignments_df

          # Calculate centroids from the ORIGINAL UNSCALED data for this iteration
          # Directly use all_mcmc_samples_df to be certain of unscaled source
          original_data_this_iter_for_centroids <- all_mcmc_samples_df %>%
            filter(MCMC_Iteration == iter_val) %>%
            mutate(Cluster_Assigned_Iter = kmeans_iter_result$cluster) # Add cluster assignments
          
          iter_unscaled_centroids <- original_data_this_iter_for_centroids %>%
            group_by(Cluster_Assigned_Iter) %>%
            summarise(
              # Explicitly list the original columns from all_mcmc_samples_df to average
              R0_sampled = mean(R0_sampled, na.rm = TRUE),
              SI_Clust_sampled = mean(SI_Clust_sampled, na.rm = TRUE),
              CFR_sampled = mean(CFR_sampled, na.rm = TRUE),
              Presymp_Proportion_sampled = mean(Presymp_Proportion_sampled, na.rm = TRUE),
              Route_resp = mean(Route_resp, na.rm = TRUE),
              Route_direct = mean(Route_direct, na.rm = TRUE),
              Route_sexual = mean(Route_sexual, na.rm = TRUE),
              Route_animal = mean(Route_animal, na.rm = TRUE),
              Route_vector = mean(Route_vector, na.rm = TRUE),
              .groups = 'drop'
            ) %>%
            rename(Cluster_ID_Iter = Cluster_Assigned_Iter)
          
          iter_unscaled_centroids$MCMC_Iteration <- iter_val
          all_iteration_centroids[[as.character(iter_val)]] <- iter_unscaled_centroids
      }
    }
  } else {
    print("No MCMC iterations found in the data.")
  }

  if (length(all_iteration_centroids) > 0 && length(all_iteration_assignments) > 0) {
      combined_centroids_df <- bind_rows(all_iteration_centroids)
      combined_assignments_df <- bind_rows(all_iteration_assignments)

      # --- 8.3. Summarize Cluster Centroids (Mean and 95% CI) ---
      # These centroids are now from UNSCALED data
      feature_cols_for_centroids <- names(combined_centroids_df)[!names(combined_centroids_df) %in% c("Cluster_ID_Iter", "MCMC_Iteration")]
      
      if(length(feature_cols_for_centroids) > 0 && nrow(combined_centroids_df) > 0) {
          centroid_summary <- combined_centroids_df %>%
            group_by(Cluster_ID_Iter) %>%
            summarise(
              across(all_of(feature_cols_for_centroids), list(
                mean = ~mean(.x, na.rm = TRUE),
                lowerCI = ~quantile(.x, 0.025, na.rm = TRUE),
                upperCI = ~quantile(.x, 0.975, na.rm = TRUE)
              )),
              N_iterations_in_summary = n(), # Count how many MCMC iterations contribute to each cluster summary
              .groups = 'drop'
            )
          
          print("Summary of Cluster Centroids (Mean and 95% CI across MCMC iterations):")
          print(as.data.frame(centroid_summary))
          write_csv(centroid_summary, "Clustering/mcmc/Kmeans/cluster_centroids_summary_with_ci_k6.csv")
          print("Cluster centroids summary saved to Clustering/mcmc/Kmeans/cluster_centroids_summary_with_ci_k6.csv")
      } else {
          print("No feature columns or data available for centroid summary.")
      }

      # --- 8.4. Determine Modal Cluster Assignment for each Pathogen ---
      if (nrow(combined_assignments_df) > 0) {
          modal_assignments <- combined_assignments_df %>%
            group_by(Pathogen_Name, Cluster_Assigned) %>%
            summarise(count = n(), .groups = 'drop_last') %>%
            slice_max(order_by = count, n = 1, with_ties = FALSE) %>%
            ungroup() %>%
            select(Pathogen_Name, Modal_Cluster = Cluster_Assigned, Modal_Cluster_Count = count)
            
          print("Modal cluster assignments for each pathogen:")
          print(modal_assignments)
          write_csv(modal_assignments, "Clustering/mcmc/Kmeans/pathogen_modal_cluster_assignments_k6.csv")
      } else {
          print("No assignment data available for modal cluster calculation.")
          modal_assignments <- data.frame(Pathogen_Name = character(), Modal_Cluster = integer(), Modal_Cluster_Count = integer())
      }

      # --- 8.5. Visualization (using modal assignments and overall data) ---
      pathogen_summary_for_plot <- all_mcmc_samples_df %>%
          group_by(Pathogen_Name) %>%
          summarise(
              R0_mean_overall = mean(R0_sampled, na.rm = TRUE),
              SI_Clust_mean_overall = mean(SI_Clust_sampled, na.rm = TRUE),
              CFR_mean_overall = mean(CFR_sampled, na.rm = TRUE),
              Presymp_Proportion_mean_overall = mean(Presymp_Proportion_sampled, na.rm = TRUE),
              Route_resp_mean = mean(Route_resp, na.rm=TRUE), # Mean for binary gives proportion
              Route_direct_mean = mean(Route_direct, na.rm=TRUE),
              Route_sexual_mean = mean(Route_sexual, na.rm=TRUE),
              Route_animal_mean = mean(Route_animal, na.rm=TRUE),
              Route_vector_mean = mean(Route_vector, na.rm=TRUE),
              .groups = 'drop'
          ) %>%
          left_join(modal_assignments, by = "Pathogen_Name")

      plot_features_numerical <- c("R0_mean_overall", "SI_Clust_mean_overall", "CFR_mean_overall", "Presymp_Proportion_mean_overall")
      plot_features_routes <- c("Route_resp_mean", "Route_direct_mean", "Route_sexual_mean", "Route_animal_mean", "Route_vector_mean")

      features_for_pca_plot <- pathogen_summary_for_plot %>%
          select(all_of(plot_features_numerical), all_of(plot_features_routes))
      
      # Scale these summarized numerical features for PCA (routes are proportions 0-1, treat as somewhat scaled)
      scaled_features_for_pca_plot <- features_for_pca_plot
      scaled_features_for_pca_plot[plot_features_numerical] <- scale(scaled_features_for_pca_plot[plot_features_numerical])

      if (!requireNamespace("ggplot2", quietly = TRUE)) {
          print("Package 'ggplot2' is not installed.")
      } else if (!requireNamespace("factoextra", quietly = TRUE)){
          print("Package 'factoextra' is not installed.")
      } else if (!requireNamespace("ggrepel", quietly = TRUE)){
          print("Package 'ggrepel' is not installed. Please install it: install.packages('ggrepel')")
      } else if (nrow(scaled_features_for_pca_plot) > 1 && !any(is.na(pathogen_summary_for_plot$Modal_Cluster))) {
          library(ggplot2)
          library(factoextra)
          library(ggrepel)

          pca_result_modal <- prcomp(scaled_features_for_pca_plot, center = FALSE, scale. = FALSE) 
          pca_data_modal <- as.data.frame(pca_result_modal$x)
          pca_data_modal$Cluster <- factor(pathogen_summary_for_plot$Modal_Cluster, levels = 1:CHOSEN_K) # Ensure all K levels are present for color mapping
          pca_data_modal$Pathogen_Name <- pathogen_summary_for_plot$Pathogen_Name

          pca_plot_modal <- ggplot(pca_data_modal, aes(x = PC1, y = PC2, color = Cluster, label = Pathogen_Name)) +
            geom_point(size = 3) +
            geom_text_repel(size = 3.5, max.overlaps = Inf, fontface = "bold") +
            labs(title = paste0("Pathogen Clusters (Modal Assignment, K=", CHOSEN_K, ")"),
                 x = paste0("PC1 (", round(summary(pca_result_modal)$importance[2,1]*100, 1), "%)"),
                 y = paste0("PC2 (", round(summary(pca_result_modal)$importance[2,2]*100, 1), "%)")) +
            theme_minimal(base_size = 12) +
            scale_color_brewer(palette = "Set1", name = "Modal Cluster") # Example palette
          
          print(pca_plot_modal)
          ggsave(paste0("Clustering/mcmc/Kmeans/pca_modal_cluster_plot_k",CHOSEN_K,".png"), plot = pca_plot_modal, width=10, height=8)
          print(paste0("PCA plot with modal cluster assignments saved to Clustering/mcmc/Kmeans/pca_modal_cluster_plot_k",CHOSEN_K,".png"))
          
          # Silhouette plot for an example MCMC iteration
          if (length(mcmc_iterations) > 0) {
            example_iter_val <- mcmc_iterations[1]
            example_iter_data <- scaled_clustering_data_full %>% filter(MCMC_Iteration == example_iter_val)
            example_iter_features <- example_iter_data %>% select(-MCMC_Iteration, -Pathogen_Name)
            
            if (nrow(example_iter_features) >= CHOSEN_K) {
                set.seed(123 + as.integer(example_iter_val))
                kmeans_example_iter_result <- tryCatch({
                    kmeans(example_iter_features, centers = CHOSEN_K, nstart = 25)
                }, error = function(e) NULL)

                if(!is.null(kmeans_example_iter_result)){
                    sil_plot_example_iter <- fviz_silhouette(kmeans_example_iter_result, data = example_iter_features, ggtheme = theme_minimal()) +
                                          labs(title=paste("Silhouette Plot for K=", CHOSEN_K, "(Example MCMC Iteration:", example_iter_val, ")"))
                    print(sil_plot_example_iter)
                    ggsave(paste0("Clustering/mcmc/Kmeans/silhouette_plot_k",CHOSEN_K,"_iter",example_iter_val,".png"), plot = sil_plot_example_iter, width=8, height=6)
                    print(paste0("Silhouette plot for K=", CHOSEN_K, " saved for example iteration ", example_iter_val))
                } else {
                     print(paste("Kmeans failed for example iteration ", example_iter_val, ", skipping silhouette plot for it."))
                }
            } else {
                 print(paste("Not enough data points in example iteration ", example_iter_val, " for K=", CHOSEN_K,", skipping silhouette plot."))
            }
          } else {
             print("No MCMC iterations available to generate an example silhouette plot.")
          }
      } else {
          print("Skipping PCA plot due to insufficient data, NA in modal clusters, or missing graphics packages.")
      }
  } else {
      print("No centroids or assignments were generated from per-iteration clustering. Cannot summarize or visualize.")
  }

} else {
  print("MCMC samples data frame 'all_mcmc_samples_df' not found or empty. Skipping clustering.")
}

print("Clustering analysis section (per-iteration) finished.")


# --- 9. Ensemble Clustering / Consensus Clustering --- 
if (exists("combined_assignments_df") && nrow(combined_assignments_df) > 0 && 
    exists("all_mcmc_samples_df") && nrow(all_mcmc_samples_df) > 0) {
  
  print("--- Starting Ensemble/Consensus Clustering Analysis ---")
  
  if (!requireNamespace("stats", quietly = TRUE)) { install.packages("stats", quiet = TRUE) }
  if (!requireNamespace("dendextend", quietly = TRUE)) { install.packages("dendextend", quiet = TRUE) }
  library(stats) # For hclust, cutree
  library(dendextend) # Added dendextend

  # Define the mapping for full pathogen names
  pathogen_full_name_map <- c(
    "COVID-19_WT" = "SARS-CoV-2 (WT)",
    "COVID-19_A" = "SARS-CoV-2 (Alpha)",
    "COVID-19_D" = "SARS-CoV-2 (Delta)",
    "COVID-19_O" = "SARS-CoV-2 (Omicron)", 
    "H1N1_18" = "A/H1N1",
    "H2N2" = "A/H2N2",
    "H3N2" = "A/H3N2",
    "H1N1_09" = "A/H1N1/09",
    "H5N1" = "A/H5N1",
    "Ebola" = "EBOV",
    "Marburg" = "MARV",
    "Mpox" = "MPV", # Note: User provided MPV for Mpox
    "Lassa" = "LASV",
    "Nipah" = "NiV",
    "Zika" = "ZIKV",
    "SARS" = "SARS-CoV-1",
    "MERS" = "MERS-CoV",
    "CCHF" = "CCHFV"
  )
  
  pathogen_names <- sort(unique(combined_assignments_df$Pathogen_Name))
  num_pathogens <- length(pathogen_names)
  num_mcmc_iterations <- length(unique(combined_assignments_df$MCMC_Iteration))

  # --- 9.1. Construct Co-assignment Matrix ---
  print("Constructing co-assignment matrix...")
  coassignment_matrix <- matrix(0, nrow = num_pathogens, ncol = num_pathogens,
                                dimnames = list(pathogen_names, pathogen_names))

  # Efficiently get assignments per iteration
  assignments_by_iter <- split(combined_assignments_df[, c("Pathogen_Name", "Cluster_Assigned")], 
                               combined_assignments_df$MCMC_Iteration)

  for (iter_idx in 1:length(assignments_by_iter)) {
    if (iter_idx %% max(1, round(length(assignments_by_iter)/10)) == 0) {
      print(paste("Processing co-assignments for MCMC iteration", iter_idx, "/", length(assignments_by_iter)))
    }
    current_iter_assignments_df <- assignments_by_iter[[iter_idx]]
    # Ensure pathogen order for consistent indexing if necessary, though direct matching is used here
    # Pathogens might not all be present in every K-means iteration if K-means failed for some with few data points etc.
    # However, combined_assignments_df should only contain successful assignments.

    for (i in 1:(num_pathogens - 1)) {
      for (j in (i + 1):num_pathogens) {
        pathogen1_name <- pathogen_names[i]
        pathogen2_name <- pathogen_names[j]
        
        p1_assignment_row <- current_iter_assignments_df[current_iter_assignments_df$Pathogen_Name == pathogen1_name, ]
        p2_assignment_row <- current_iter_assignments_df[current_iter_assignments_df$Pathogen_Name == pathogen2_name, ]

        # Only increment if both pathogens have assignments in this iteration and are in the same cluster
        if (nrow(p1_assignment_row) > 0 && nrow(p2_assignment_row) > 0) {
          if (p1_assignment_row$Cluster_Assigned == p2_assignment_row$Cluster_Assigned) {
            coassignment_matrix[pathogen1_name, pathogen2_name] <- coassignment_matrix[pathogen1_name, pathogen2_name] + 1
            coassignment_matrix[pathogen2_name, pathogen1_name] <- coassignment_matrix[pathogen1_name, pathogen2_name] # Keep symmetric
          }
        }
      }
    }
  }
  # Diagonal should be total number of iterations a pathogen was part of clustering
  # For simplicity here, we assume all pathogens are in all iterations in combined_assignments_df.
  # If not, a pathogen-specific count would be more precise for the diagonal.
  diag(coassignment_matrix) <- num_mcmc_iterations 
  
  print("Co-assignment matrix (first 6x6):")
  print(head(coassignment_matrix[,1:min(6, ncol(coassignment_matrix))]))
  write.csv(coassignment_matrix, "Clustering/mcmc/Kmeans/coassignment_matrix.csv")

  # --- 9.2. Calculate Dissimilarity Matrix ---
  print("Calculating dissimilarity matrix...")
  # N_MCMC_ITERATIONS is defined at the top of the script
  dissimilarity_matrix <- 1 - (coassignment_matrix / N_MCMC_ITERATIONS) 
  # Ensure diagonal is 0 after division, if any NA from 0/0 (though N_MCMC_ITERATIONS is likely >0)
  diag(dissimilarity_matrix) <- 0 

  print("Dissimilarity matrix (first 6x6):")
  print(head(dissimilarity_matrix[,1:min(6, ncol(dissimilarity_matrix))]))
  write.csv(dissimilarity_matrix, "Clustering/mcmc/Kmeans/dissimilarity_matrix.csv")

  # --- 9.3. Perform Hierarchical Clustering ---
  print("Performing hierarchical clustering...")
  # Convert to dist object for hclust
  pathogen_dist <- as.dist(dissimilarity_matrix)
  hierarchical_clust <- hclust(pathogen_dist, method = "average") # Changed back to average

  # --- 9.4. Plot Dendrogram (Revised with dendextend) ---
  print("Plotting enhanced dendrogram with dendextend...")
  
  # Define the number of consensus clusters (K) and a color palette before plotting.
  # This value is also used in section 9.5 to cut the tree.
  K_CONSENSUS <- 6
  if (!requireNamespace("RColorBrewer", quietly = TRUE)) { install.packages("RColorBrewer", quiet = TRUE) }
  library(RColorBrewer)
  
  cluster_colors <- brewer.pal(n = max(3, K_CONSENSUS), name = "Set2") # Ensure at least 3 colors
  if (K_CONSENSUS > length(cluster_colors)) { # If K is larger than palette size
      cluster_colors <- rep(cluster_colors, length.out = K_CONSENSUS) # Repeat colors
  }

  dend <- as.dendrogram(hierarchical_clust)
  
  # Apply full names to labels
  original_labels <- labels(dend)
  new_labels <- pathogen_full_name_map[original_labels]
  # Handle any names not in the map by keeping their original short name
  new_labels[is.na(new_labels)] <- original_labels[is.na(new_labels)]
  labels(dend) <- new_labels
  
  dend_colored <- color_branches(dend, k = K_CONSENSUS, col = cluster_colors[1:K_CONSENSUS])
  dend_colored <- set(dend_colored, "labels_cex", 0.7) # Adjust label size
  dend_colored <- set(dend_colored, "branches_lwd", 3) # Adjust branch line width
  
  png_filename <- paste0("Clustering/mcmc/Kmeans/consensus_dendrogram_colored_k", K_CONSENSUS, ".png")
  png(png_filename, width=1200, height=800, units="px", res=100)
  par(mar = c(5, 4, 4, 10)) # Adjust right margin for labels
  plot(dend_colored, horiz = TRUE, 
       main = paste("Consensus Clustering Dendrogram (K=", K_CONSENSUS, ", Average Linkage)"), 
       xlab = "Dissimilarity (1 - Proportion Co-assigned)")

  dev.off()
  print(paste("Colored consensus dendrogram saved to", png_filename))

  # --- 9.5. Extract Consensus Cluster Assignments ---
  # User should inspect the dendrogram to choose K_consensus
  # K_CONSENSUS <- 6 # This is now defined in section 9.4 before the dendrogram plot.
  print(paste("Cutting tree to get", K_CONSENSUS, "consensus clusters..."))
  consensus_clusters <- cutree(hierarchical_clust, k = K_CONSENSUS)
  consensus_assignments_df <- data.frame(Pathogen_Name = names(consensus_clusters),
                                         Consensus_Cluster = consensus_clusters)
  print("Consensus cluster assignments:")
  print(consensus_assignments_df)
  write.csv(consensus_assignments_df, paste0("Clustering/mcmc/Kmeans/consensus_cluster_assignments_k", K_CONSENSUS, ".csv"))

  # --- 9.6. Characterize Consensus Clusters ---
  print("Characterizing consensus clusters (means and 95% CIs from all_mcmc_samples_df)...")
  
  # Select the numerical and route columns to summarize from all_mcmc_samples_df
  # These are the original sampled values, not the means from previous steps.
  params_to_summarize <- c("R0_sampled", "SI_Clust_sampled", "CFR_sampled", "Presymp_Proportion_sampled",
                           "Route_resp", "Route_direct", "Route_sexual", "Route_animal", "Route_vector")

  # Join consensus assignments with the full MCMC samples
  pathogen_mcmc_samples_with_consensus_clusters <- all_mcmc_samples_df %>% 
    filter(Pathogen_Name %in% consensus_assignments_df$Pathogen_Name) %>% # Ensure we only use pathogens that got consensus assignment
    left_join(consensus_assignments_df, by = "Pathogen_Name")

  if(nrow(pathogen_mcmc_samples_with_consensus_clusters) > 0 && "Consensus_Cluster" %in% names(pathogen_mcmc_samples_with_consensus_clusters)){
    consensus_cluster_summary_list <- list()
    for (param in params_to_summarize) {
      if (param %in% names(pathogen_mcmc_samples_with_consensus_clusters)) {
        summary_for_param <- pathogen_mcmc_samples_with_consensus_clusters %>%
          group_by(Consensus_Cluster) %>% # Group by the new stable consensus cluster ID
          summarise(
            mean_val = mean(get(param), na.rm = TRUE),
            lower_ci = quantile(get(param), 0.025, na.rm = TRUE),
            upper_ci = quantile(get(param), 0.975, na.rm = TRUE),
            median_val = median(get(param), na.rm = TRUE),
            sd_val = sd(get(param), na.rm=TRUE),
            n_mcmc_samples_in_calc = sum(!is.na(get(param))), # count non-NA samples for this param in this cluster
            .groups = 'drop'
          ) %>% 
          rename_with(~paste0(param, "_", .), .cols = c(mean_val, lower_ci, upper_ci, median_val, sd_val, n_mcmc_samples_in_calc))
        consensus_cluster_summary_list[[param]] <- summary_for_param
      } else {
        warning(paste("Parameter", param, "not found in pathogen_mcmc_samples_with_consensus_clusters."))
      }
    }
    
    if(length(consensus_cluster_summary_list) > 0){
      # Merge all parameter summaries by Consensus_Cluster
      final_consensus_summary_df <- Reduce(function(x, y) full_join(x, y, by = "Consensus_Cluster"), consensus_cluster_summary_list)
      print("Summary of Consensus Cluster Characteristics (Mean, 95% CI, Median, SD):")
      print(as.data.frame(final_consensus_summary_df))
      write_csv(final_consensus_summary_df, paste0("Clustering/mcmc/Kmeans/consensus_clusters_summary_k", K_CONSENSUS, ".csv"))
      print(paste0("Consensus clusters summary saved to Clustering/mcmc/Kmeans/consensus_clusters_summary_k", K_CONSENSUS, ".csv"))
    } else {
      print("No parameters were summarized for consensus clusters.")
    }
  } else {
    print("Could not prepare data for characterizing consensus clusters (pathogen_mcmc_samples_with_consensus_clusters is empty or missing Consensus_Cluster column).")
  }
  
  print("Ensemble/Consensus Clustering Analysis finished.")

} else {
  print("Skipping Ensemble/Consensus Clustering: combined_assignments_df or all_mcmc_samples_df not found or empty.")
}