# MCMC K-Means Clustering Analysis

# --- 1. Load Libraries ---
library(dplyr)
library(readr)
# Add any other necessary libraries here, e.g., for specific distributions or optimization
# library(DistributionFitR) # May be useful for fitting distributions to CI


# --- 2. Configuration ---
N_MCMC_ITERATIONS <- 100 # Number of MCMC iterations (Reduced for testing)
N_PRESYMP_SAMPLES <- 1000 # Number of samples for presymptomatic proportion estimation (Reduced for testing)
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
          # mean = exp(meanlog + sdlog^2/2)
          # sd^2 = (exp(sdlog^2) - 1) * exp(2*meanlog + sdlog^2)
          # sd^2/mean^2 = exp(sdlog^2) - 1
          # sdlog^2 = log(sd^2/mean^2 + 1)
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
    warning(paste("SI or IP distribution info incomplete. SI Family:", si_dist_info$family, "IP Family:", ip_dist_info$family,
                  "SI Params:", paste(names(si_dist_info$params), collapse=", "), 
                  "IP Params:", paste(names(ip_dist_info$params), collapse=", ")))
    return(NA)
  }

  si_samples <- NULL; ip_samples <- NULL

  # Generate samples from SI distribution
  if(si_dist_info$family == "Lognormal"){
    if(all(c("meanlog", "sdlog") %in% names(si_dist_info$params)) && si_dist_info$params$sdlog > 0){
      si_samples <- rlnorm(n_samples, meanlog = si_dist_info$params$meanlog, sdlog = si_dist_info$params$sdlog)
    } else { warning(paste("Missing/invalid params for SI Lognormal:", si_dist_info$Pathogen_Name, names(si_dist_info$params))); return(NA) }
  } else if (si_dist_info$family == "Gamma"){
    if(all(c("shape") %in% names(si_dist_info$params)) && si_dist_info$params$shape > 0 && 
       (( "rate" %in% names(si_dist_info$params) && si_dist_info$params$rate > 0) || 
        ( "scale" %in% names(si_dist_info$params) && si_dist_info$params$scale > 0 ))){
      si_samples <- rgamma(n_samples, shape = si_dist_info$params$shape, 
                                  rate = if (!is.null(si_dist_info$params$rate)) si_dist_info$params$rate else 1/si_dist_info$params$scale)
    } else { warning(paste("Missing/invalid params for SI Gamma:", si_dist_info$Pathogen_Name, names(si_dist_info$params))); return(NA) }
  } else if (si_dist_info$family == "Normal"){
     if(all(c("mean", "SD") %in% names(si_dist_info$params)) && si_dist_info$params$SD > 0){
      si_samples <- rnorm(n_samples, mean = si_dist_info$params$mean, sd = si_dist_info$params$SD)
    } else { warning(paste("Missing/invalid params for SI Normal:", si_dist_info$Pathogen_Name, names(si_dist_info$params))); return(NA) }
  } else if (si_dist_info$family == "Weibull"){
    if(all(c("shape", "scale") %in% names(si_dist_info$params)) && si_dist_info$params$shape > 0 && si_dist_info$params$scale > 0){
      si_samples <- rweibull(n_samples, shape = si_dist_info$params$shape, scale = si_dist_info$params$scale)
    } else { warning(paste("Missing/invalid params for SI Weibull:", si_dist_info$Pathogen_Name, names(si_dist_info$params))); return(NA) }
  } else {
    warning(paste("Unsupported SI distribution family for sampling:", si_dist_info$Pathogen_Name, si_dist_info$family))
    return(NA)
  }

  # Generate samples from IP distribution
  if(ip_dist_info$family == "Lognormal"){
    if(all(c("meanlog", "sdlog") %in% names(ip_dist_info$params)) && ip_dist_info$params$sdlog > 0){
      ip_samples <- rlnorm(n_samples, meanlog = ip_dist_info$params$meanlog, sdlog = ip_dist_info$params$sdlog)
    } else { warning(paste("Missing/invalid params for IP Lognormal:", ip_dist_info$Pathogen_Name, names(ip_dist_info$params))); return(NA) }
  } else if (ip_dist_info$family == "Gamma"){
    if(all(c("shape") %in% names(ip_dist_info$params)) && ip_dist_info$params$shape > 0 &&
       (( "rate" %in% names(ip_dist_info$params) && ip_dist_info$params$rate > 0) || 
        ( "scale" %in% names(ip_dist_info$params) && ip_dist_info$params$scale > 0 ))){
      ip_samples <- rgamma(n_samples, shape = ip_dist_info$params$shape, 
                                  rate = if (!is.null(ip_dist_info$params$rate)) ip_dist_info$params$rate else 1/ip_dist_info$params$scale)
    } else { warning(paste("Missing/invalid params for IP Gamma:", ip_dist_info$Pathogen_Name, names(ip_dist_info$params))); return(NA) }
  } else if (ip_dist_info$family == "Normal"){
    if(all(c("mean", "SD") %in% names(ip_dist_info$params)) && ip_dist_info$params$SD > 0){
      ip_samples <- rnorm(n_samples, mean = ip_dist_info$params$mean, sd = ip_dist_info$params$SD)
    } else { warning(paste("Missing/invalid params for IP Normal:", ip_dist_info$Pathogen_Name, names(ip_dist_info$params))); return(NA) }
  } else if (ip_dist_info$family == "Weibull"){
    if(all(c("shape", "scale") %in% names(ip_dist_info$params)) && ip_dist_info$params$shape > 0 && ip_dist_info$params$scale > 0){
      ip_samples <- rweibull(n_samples, shape = ip_dist_info$params$shape, scale = ip_dist_info$params$scale)
    } else { warning(paste("Missing/invalid params for IP Weibull:", ip_dist_info$Pathogen_Name, names(ip_dist_info$params))); return(NA) }
  } else {
    warning(paste("Unsupported IP distribution family for sampling:", ip_dist_info$Pathogen_Name, ip_dist_info$family))
    return(NA)
  }
  
  si_samples[si_samples < 0] <- 0
  ip_samples[ip_samples < 0] <- 0

  presymp_proportion <- mean(si_samples < ip_samples, na.rm = TRUE)
  return(presymp_proportion)
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
        !is.null(si_full_dist_info$family) && !is.null(ip_full_dist_info$family) &&
        !is.na(si_full_dist_info$family) && !is.na(ip_full_dist_info$family) &&
        length(si_full_dist_info$params) > 0 && length(ip_full_dist_info$params) > 0) {
      # Add Pathogen_Name to dist_info for more informative warnings inside estimate_presymptomatic_proportion
      si_full_dist_info$Pathogen_Name <- pathogen_name 
      ip_full_dist_info$Pathogen_Name <- pathogen_name
      presymp_prop_sampled <- estimate_presymptomatic_proportion(si_full_dist_info, ip_full_dist_info)
    } else {
      warning(paste("Could not estimate presymptomatic proportion for", pathogen_name, "due to incomplete SI/IP full dist info from get_full_dist_params."))
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