# MCMC K-Means Clustering Analysis

# --- 1. Load Libraries & Setup ---
# Helper function to install and load packages
install_and_load <- function(packages) {
  for (pkg in packages) {
    if (!require(pkg, character.only = TRUE)) {
      install.packages(pkg, dependencies = TRUE)
      library(pkg, character.only = TRUE)
    }
  }
}

# List of required packages
required_packages <- c(
  "dplyr", "readr", "ggplot2", "factoextra", "ggrepel", 
  "cluster", "dendextend", "RColorBrewer", "future", "furrr",
  "patchwork", "ggplotify"
)
install_and_load(required_packages)

# Setup for parallel processing
plan(multisession)


# --- 2. Configuration ---
N_MCMC_ITERATIONS <- 5000 # Number of MCMC iterations 
MCMC_SEED <- 123 # Seed for reproducibility of the MCMC parameter sampling
CLUSTERING_SEED <- 456 # Seed for reproducibility of the K-means clustering
set.seed(MCMC_SEED) # Set seed for the MCMC part


# --- 3. Load Data ---
# Load the parameter data and transmission routes
params_long_df <- read_csv("Clustering/mcmc/Kmeans/data/pathogen_params.csv")
transmission_df <- read_csv("Clustering/mcmc/Kmeans/data/transmission_route.csv")

# Filter for the specific pathogens of interest 
pathogens_of_interest <- c(
  "COVID-19_WT", "COVID-19_A", "COVID-19_D", "COVID-19_O", "Ebola", "Marburg", "Lassa", "CCHF", "Nipah", "Zika", "Mpox", 
  "H1N1_18", "H2N2", "H3N2", "H1N1_09", "H5N1", "SARS", "MERS", "HIV", "Cholera", "Measles", "Norovirus", "SFTS", 
  "EVA71", "hMPV", "Plague", "Smallpox", "Chikungunya", "RVFV", "Anthrax"
)
params_long_df <- params_long_df %>% filter(Pathogen_Name %in% pathogens_of_interest)
transmission_df <- transmission_df %>% filter(Pathogen_Name %in% pathogens_of_interest)



# print(head(params_long_df))
# print(str(params_long_df))
# print(head(transmission_df))


# --- 4. Helper Functions ---

# ---  Quantile-Matching Helper Functions for Beta and Gamma Distributions ---
# These functions use optimization to find distribution parameters that best match a given 95% CI.


# Function to derive parameters for rbeta by matching quantiles of the 95% CI
get_beta_params_from_ci <- function(lower_ci, upper_ci, mean_val) {
    min_param_val <- 1e-6

    # Objective function to minimize: squared error between target CIs and beta quantiles
    objective_beta <- function(params, lower_ci, upper_ci) {
        # Use exp() to ensure parameters are positive, a common practice in optimization
        shape1 <- exp(params[1])
        shape2 <- exp(params[2])
        if (is.na(shape1) || is.na(shape2) || !is.finite(shape1) || !is.finite(shape2)) return(1e10)
        
        # Calculate error
        q_lower <- qbeta(0.025, shape1, shape2)
        q_upper <- qbeta(0.975, shape1, shape2)
        err <- (q_lower - lower_ci)^2 + (q_upper - upper_ci)^2
        return(err)
    }

    # Use mean to generate a reasonable starting guess for optimization
    if (!is.na(mean_val) && mean_val > 0 && mean_val < 1) {
        # Crude guess for variance, assuming symmetric CI initially
        sd_guess <- (upper_ci - lower_ci) / (2 * qnorm(0.975))
        var_guess <- sd_guess^2
        # Method of moments estimates for alpha and beta
        alpha_guess <- ((1 - mean_val) / var_guess - 1 / mean_val) * mean_val^2
        beta_guess <- alpha_guess * (1 / mean_val - 1)
        
        if (is.na(alpha_guess) || alpha_guess <= 0) alpha_guess <- 1
        if (is.na(beta_guess) || beta_guess <= 0) beta_guess <- 1
        start_params <- log(c(alpha_guess, beta_guess))
    } else {
        start_params <- log(c(1, 1)) # Generic starting point
    }

    # Run optimization
    opt_result <- tryCatch({
        optim(start_params, objective_beta, lower_ci = lower_ci, upper_ci = upper_ci, control = list(maxit = 1000, reltol = 1e-8))
    }, error = function(e) NULL)

    # If optimization is successful, return the parameters
    if (!is.null(opt_result) && opt_result$convergence == 0) {
        return(list(shape1 = exp(opt_result$par[1]), shape2 = exp(opt_result$par[2])))
    } else {
        warning(paste("Beta CI-matching optimization failed for CI [", lower_ci, ",", upper_ci, "]. Falling back to simpler method."))
        # Fallback to the simpler, mean-based method
        return(get_beta_params_from_mean_ci_fallback(mean_val, lower_ci, upper_ci))
    }
}

# Function to derive parameters for rgamma by matching quantiles of the 95% CI
get_gamma_params_from_ci <- function(lower_ci, upper_ci, mean_val) {
    min_param_val <- 1e-6

    # Objective function to minimize for gamma
    objective_gamma <- function(params, lower_ci, upper_ci) {
        shape <- exp(params[1])
        rate <- exp(params[2])
        if (is.na(shape) || is.na(rate) || !is.finite(shape) || !is.finite(rate)) return(1e10)

        q_lower <- qgamma(0.025, shape = shape, rate = rate)
        q_upper <- qgamma(0.975, shape = shape, rate = rate)
        err <- (q_lower - lower_ci)^2 + (q_upper - upper_ci)^2
        return(err)
    }

    # Use mean and CI to generate a reasonable starting guess
    if(!is.na(mean_val) && mean_val > 0){
        sd_guess <- (upper_ci - lower_ci) / (2 * qnorm(0.975))
        if(is.na(sd_guess) || sd_guess <= 0) sd_guess <- mean_val # fallback sd guess
        shape_guess <- (mean_val / sd_guess)^2
        rate_guess <- mean_val / sd_guess^2
        if (is.na(shape_guess) || shape_guess <= 0) shape_guess <- 1
        if (is.na(rate_guess) || rate_guess <= 0) rate_guess <- 1
        start_params <- log(c(shape_guess, rate_guess))
    } else {
        start_params <- log(c(1,1)) # Generic starting point
    }

    # Run optimization
    opt_result <- tryCatch({
        optim(start_params, objective_gamma, lower_ci = lower_ci, upper_ci = upper_ci, control = list(maxit = 1000, reltol = 1e-8))
    }, error = function(e) NULL)

    # If successful, return params
    if (!is.null(opt_result) && opt_result$convergence == 0) {
        return(list(shape = exp(opt_result$par[1]), rate = exp(opt_result$par[2])))
    } else {
        warning(paste("Gamma CI-matching optimization failed for CI [", lower_ci, ",", upper_ci, "]. Falling back to simpler method."))
        # Fallback to the simpler, mean-based method
        return(get_gamma_params_from_mean_ci_fallback(mean_val, lower_ci, upper_ci))
    }
}


# --- FALLBACK Helper Functions  ---

# Function to derive parameters for rbeta from mean and 95% CI
get_beta_params_from_mean_ci_fallback <- function(mean_val, lower_ci, upper_ci, n_eff_guess = 1000) {
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
get_gamma_params_from_mean_ci_fallback <- function(mean_val, lower_ci, upper_ci) {
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

# --- BOOTSTRAP AGGREGATION SAMPLING FUNCTION ---
# This function uses bootstrap aggregation to create a robust parameter estimate for each MCMC iteration.
sample_parameter_bootstrap_aggregation <- function(param_name, pathogen_name, data_df) {
  # 1. Filter for all studies for the given pathogen and parameter
  studies <- data_df %>% filter(Pathogen_Name == pathogen_name, Parameter == param_name)
  if (nrow(studies) == 0) return(NA)

  # 2. Bootstrap: Resample studies with replacement. If only 1 study, it will be used.
  if (nrow(studies) > 1) {
      bootstrapped_indices <- sample(seq_len(nrow(studies)), size = nrow(studies), replace = TRUE)
      bootstrapped_studies <- studies[bootstrapped_indices, ]
  } else {
      bootstrapped_studies <- studies
  }

  # 3. Generate one sample from each study in the bootstrapped set
  samples <- apply(bootstrapped_studies, 1, function(study_row) {
    mean_val <- as.numeric(study_row["ReportedValue"])
    lower_ci <- as.numeric(study_row["LowerBound"])
    upper_ci <- as.numeric(study_row["UpperBound"])
    
    # If CI is missing or invalid, fall back to the point estimate
    if (is.na(lower_ci) || is.na(upper_ci) || lower_ci >= upper_ci) {
        return(mean_val) # Return point estimate directly for this study
    }
    
    sampled_val_internal <- NA
    
    if (param_name == "CFR") {
        if (!is.na(mean_val) && mean_val > 0 && mean_val < 1 && !is.na(lower_ci) && !is.na(upper_ci) && lower_ci < upper_ci) {
            beta_params <- get_beta_params_from_ci(lower_ci, upper_ci, mean_val)
            if (!is.na(beta_params$shape1) && !is.na(beta_params$shape2)) {
                sampled_val_internal <- rbeta(1, shape1 = beta_params$shape1, shape2 = beta_params$shape2)
            }
        }
    } else if (param_name %in% c("R0", "SI", "k", "IP", "LP", "InfP")) {
        if (!is.na(lower_ci) && !is.na(upper_ci) && lower_ci < upper_ci) {
            gamma_params <- get_gamma_params_from_ci(lower_ci, upper_ci, mean_val)
            if (!is.na(gamma_params$shape) && !is.na(gamma_params$rate)) {
                sampled_val_internal <- rgamma(1, shape = gamma_params$shape, rate = gamma_params$rate)
            }
        }
    } else {
        sd_from_ci <- (upper_ci - lower_ci) / (2 * qnorm(0.975))
        if (!is.na(mean_val) && !is.na(sd_from_ci) && sd_from_ci > 0) {
            sampled_val_internal <- rnorm(1, mean = mean_val, sd = sd_from_ci)
        }
    }
    
    # Fallback to mean if sampling failed
    if (is.na(sampled_val_internal) || !is.finite(sampled_val_internal)) {
        return(mean_val)
    }
    return(sampled_val_internal)
  })

  # 4. Aggregate the samples from the bootstrapped studies by taking the mean
  final_value <- mean(samples, na.rm = TRUE)
  
  # 5. Apply final constraints to the aggregated value to ensure plausibility
  if (is.na(final_value)) return(NA)
  
  if (param_name == "CFR") {
    final_value <- pmax(0, pmin(1, final_value))
  } else if (param_name %in% c("R0", "SI", "k", "IP", "LP", "InfP")) {
    final_value <- pmax(0.00001, final_value)
  }
  
  return(final_value)
}


get_full_dist_params_from_presym <- function(study_row) {
  if (is.null(study_row) || nrow(study_row) == 0) {
    return(list(family = NA, params = list()))
  }

  family <- study_row$FullDist_Family
  raw_params <- list()

  for (i in 1:2) {
    param_name_col <- paste0("FullDist_Param", i, "_Name")
    param_val_col <- paste0("FullDist_Param", i, "_Value")
    
    if (!is.na(study_row[[param_name_col]]) && study_row[[param_name_col]] != "") {
      param_name <- study_row[[param_name_col]]
      param_val <- as.numeric(study_row[[param_val_col]])
      raw_params[[param_name]] <- param_val
    }
  }

  # Convert mean/sd to shape/rate for Gamma, etc. if needed
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
        }
      }
    } else if (family == "Lognormal") {
      if (all(c("mean", "SD") %in% names(raw_params)) && !("meanlog" %in% names(raw_params))) {
        mean_val <- raw_params$mean
        sd_val <- raw_params$SD
        if(!is.na(mean_val) && !is.na(sd_val) && mean_val > 0 && sd_val > 0){
          final_params$sdlog <- sqrt(log(sd_val^2 / mean_val^2 + 1))
          final_params$meanlog <- log(mean_val) - 0.5 * final_params$sdlog^2
          final_params$mean <- NULL; final_params$SD <- NULL
        }
      }
    }
  }

  return(list(family = family, params = final_params))
}


get_dist_info_from_study_row <- function(study_row) {
    if (is.null(study_row) || nrow(study_row) == 0) return(list(family = NA, params = list()))

    family_map <- c("rgamma" = "Gamma", "rlnorm" = "Lognormal", "rnorm" = "Normal")
    family <- family_map[study_row$SamplingDist]
    
    if (is.na(family)) {
        warning(paste("Unsupported distribution family:", study_row$SamplingDist, "for", study_row$Pathogen_Name))
        return(list(family = NA, params = list()))
    }
    
    mean_val <- as.numeric(study_row$ReportedValue)
    lower_ci <- as.numeric(study_row$LowerBound)
    upper_ci <- as.numeric(study_row$UpperBound)
    
    params <- list()
    if (family == "Gamma") {
        if (!is.na(mean_val) && !is.na(lower_ci) && !is.na(upper_ci)) {
            gamma_params <- get_gamma_params_from_ci(lower_ci, upper_ci, mean_val)
            params$shape <- gamma_params$shape
            params$rate <- gamma_params$rate
        } else {
            warning(paste("Incomplete info for Gamma dist for", study_row$Pathogen_Name))
            return(list(family = family, params = list()))
        }
    } else if (family == "Lognormal") {
        if (!is.na(lower_ci) && !is.na(upper_ci) && lower_ci > 0 && upper_ci > 0) {
            log_lower <- log(lower_ci)
            log_upper <- log(upper_ci)
            sdlog <- (log_upper - log_lower) / (2 * qnorm(0.975))
            meanlog <- (log_upper + log_lower) / 2
            if(is.na(sdlog) || !is.finite(sdlog) || sdlog <= 0){
                warning(paste("Could not derive valid sdlog for Lognormal for", study_row$Pathogen_Name))
                return(list(family=family, params=list()))
            }
            params$meanlog <- meanlog
            params$sdlog <- sdlog
        } else {
            warning(paste("Incomplete/invalid info for Lognormal dist for", study_row$Pathogen_Name))
            return(list(family = family, params = list()))
        }
    } else if (family == "Normal") {
        if (!is.na(mean_val) && !is.na(lower_ci) && !is.na(upper_ci)) {
            sd_val <- (upper_ci - lower_ci) / (2 * qnorm(0.975))
            if(is.na(sd_val) || !is.finite(sd_val) || sd_val <= 0){
                warning(paste("Could not derive valid SD for Normal dist for", study_row$Pathogen_Name))
                return(list(family=family, params=list()))
            }
            params$mean <- mean_val
            params$SD <- sd_val
        } else {
             warning(paste("Incomplete info for Normal dist for", study_row$Pathogen_Name))
             return(list(family = family, params = list()))
        }
    }
    
    return(list(family = family, params = params))
}


# --- 5. Main MCMC Loop ---
print("Starting parallel MCMC sampling...")

# Helper function to run a single MCMC iteration for all pathogens
run_single_mcmc_iteration <- function(iter_num, .progress = FALSE) {
  
  # The progress bar can be enabled for interactive sessions
  if (.progress) {
    if (iter_num %% max(1, (N_MCMC_ITERATIONS/10)) == 0) {
      print(paste("MCMC Iteration:", iter_num, "/", N_MCMC_ITERATIONS))
    }
  }
  
  unique_pathogens <- unique(params_long_df$Pathogen_Name)
  flu_group <- c("H1N1_09", "H1N1_18", "H2N2", "H3N2")
  
  # Use map_dfr to iterate over pathogens and bind results into a dataframe for this iteration
  iteration_results <- purrr::map_dfr(unique_pathogens, function(pathogen_name) {
    
    # --- Sample parameters using bootstrap aggregation ---
    r0_sampled <- sample_parameter_bootstrap_aggregation("R0", pathogen_name, params_long_df)
    si_clust_sampled <- sample_parameter_bootstrap_aggregation("SI", pathogen_name, params_long_df)
    cfr_sampled <- sample_parameter_bootstrap_aggregation("CFR", pathogen_name, params_long_df)
    k_sampled <- sample_parameter_bootstrap_aggregation("k", pathogen_name, params_long_df)
    ip_sampled <- sample_parameter_bootstrap_aggregation("IP", pathogen_name, params_long_df)
    lp_sampled <- sample_parameter_bootstrap_aggregation("LP", pathogen_name, params_long_df)
    infp_sampled <- sample_parameter_bootstrap_aggregation("InfP", pathogen_name, params_long_df)
    
    # --- Get transmission routes ---
    pathogen_routes <- transmission_df %>% filter(Pathogen_Name == pathogen_name)
    if(nrow(pathogen_routes) == 0){
        warning(paste("No transmission route data for pathogen:", pathogen_name))
        pathogen_routes <- data.frame(Route_resp=NA, Route_direct=NA, Route_sexual=NA, Route_animal=NA, Route_vector=NA)
    }

    # Return a single row data.frame (tibble)
    tibble(
      Pathogen_Name = pathogen_name,
      R0_sampled = ifelse(is.null(r0_sampled) || !is.finite(r0_sampled), NA_real_, r0_sampled),
      SI_Clust_sampled = ifelse(is.null(si_clust_sampled) || !is.finite(si_clust_sampled), NA_real_, si_clust_sampled),
      CFR_sampled = ifelse(is.null(cfr_sampled) || !is.finite(cfr_sampled), NA_real_, cfr_sampled),
      k_sampled = ifelse(is.null(k_sampled) || !is.finite(k_sampled), NA_real_, k_sampled),
      IP_sampled = ifelse(is.null(ip_sampled) || !is.finite(ip_sampled), NA_real_, ip_sampled),
      LP_sampled = ifelse(is.null(lp_sampled) || !is.finite(lp_sampled), NA_real_, lp_sampled),
      InfP_sampled = ifelse(is.null(infp_sampled) || !is.finite(infp_sampled), NA_real_, infp_sampled),
      Route_resp = as.integer(pathogen_routes$Route_resp),
      Route_direct = as.integer(pathogen_routes$Route_direct),
      Route_sexual = as.integer(pathogen_routes$Route_sexual),
      Route_animal = as.integer(pathogen_routes$Route_animal),
      Route_vector = as.integer(pathogen_routes$Route_vector)
    )
  })
  
  return(iteration_results)
}

# Use future_map_dfr for parallel execution. 
# It runs the function for each iteration in parallel and combines the results into a single data frame.
# .id = "MCMC_Iteration" automatically creates a column with the iteration number.
all_mcmc_samples_df <- future_map_dfr(
  1:N_MCMC_ITERATIONS, 
  run_single_mcmc_iteration, 
  .options = furrr_options(seed = MCMC_SEED), # Ensures reproducibility
  .progress = TRUE, # Shows a progress bar
  .id = "MCMC_Iteration"
)

if (nrow(all_mcmc_samples_df) > 0) {
    all_mcmc_samples_df$MCMC_Iteration <- as.integer(all_mcmc_samples_df$MCMC_Iteration)

    print("MCMC sampling complete.")
    print(paste("Generated", length(unique(all_mcmc_samples_df$MCMC_Iteration)), "valid sets of parameters for", length(unique(all_mcmc_samples_df$Pathogen_Name)), "pathogens."))

    # --- 6. Post-MCMC Analysis (Clustering) ---
    # The results are already in a tidy data frame, ready for clustering.
    
    # --- 7. Save Results ---
    write_csv(all_mcmc_samples_df, "Clustering/mcmc/Kmeans/S11_outputs/mcmc_parameter_samples.csv")
    print("MCMC parameter samples saved to Clustering/mcmc/Kmeans/S11_outputs/mcmc_parameter_samples.csv")
} else {
    print("MCMC sampling did not produce any results.")
}

print("Script finished.")


# --- 8. Clustering Analysis ---
if (exists("all_mcmc_samples_df") && nrow(all_mcmc_samples_df) > 0) {
  
  CHOSEN_K <- 4 #specified K=4
  
  print(paste0("--- Starting Clustering Analysis (Per MCMC Iteration, K=", CHOSEN_K, ") ---"))
  
  # --- 8.1. Prepare Data for Clustering (from all_mcmc_samples_df) ---
  # Select features that will be used. Pathogen_Name and MCMC_Iteration are for tracking.
  clustering_data_full <- all_mcmc_samples_df %>%
    select(
      MCMC_Iteration,
      Pathogen_Name,
      R0_sampled, 
      CFR_sampled, 
      IP_sampled,
      Route_resp,
      Route_direct,
      Route_sexual,
      Route_animal,
      Route_vector
    )

  # Handle NAs in the sampled parameters. Impute with the pathogen-specific median first,
  # then use the overall median as a fallback for any remaining NAs.
  numerical_param_cols <- c("R0_sampled", "CFR_sampled", "IP_sampled")
  
  # Log initial NA counts
  for(col in numerical_param_cols){
    if(any(is.na(clustering_data_full[[col]]))){
      na_count <- sum(is.na(clustering_data_full[[col]]))
      warning(paste("Column", col, "has", na_count, "NA values. Imputing with pathogen-specific median..."))
    }
  }

  # Step 1: Impute with pathogen-specific median
  clustering_data_full <- clustering_data_full %>%
    group_by(Pathogen_Name) %>%
    mutate(across(all_of(numerical_param_cols), 
                  ~ifelse(is.na(.), median(., na.rm = TRUE), .))) %>%
    ungroup()

  # After pathogen-specific imputation, remove any rows that still contain NA values.
  # This occurs if a pathogen had all NA values for a parameter across all MCMC iterations,
  # making imputation impossible without using a global metric.
  rows_before_na_removal <- nrow(clustering_data_full)
  clustering_data_full <- clustering_data_full %>%
      filter(if_all(all_of(numerical_param_cols), ~!is.na(.)))
  rows_after_na_removal <- nrow(clustering_data_full)

  if (rows_before_na_removal > rows_after_na_removal) {
      removed_count <- rows_before_na_removal - rows_after_na_removal
      warning(paste(removed_count, "rows were removed from clustering data due to remaining NA values after pathogen-specific imputation."))
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

  # --- 8.2. K-Means Clustering for Each MCMC Iteration (Parallelized) ---
  
  # This function encapsulates the logic for clustering a single iteration's data
  perform_kmeans_for_iteration <- function(iter_val, scaled_data, original_data) {
    
    current_iter_data_scaled <- scaled_data %>% filter(MCMC_Iteration == iter_val)
    current_iter_features_scaled <- current_iter_data_scaled %>% select(-MCMC_Iteration, -Pathogen_Name)
    
    if (nrow(current_iter_features_scaled) < CHOSEN_K) {
      warning(paste("Skipping K-means for iter", iter_val, "due to insufficient data points."))
      return(NULL) # Return NULL for failed iterations
    }

    set.seed(CLUSTERING_SEED + as.integer(iter_val))
    kmeans_result <- tryCatch({
        kmeans(current_iter_features_scaled, centers = CHOSEN_K, nstart = 25)
    }, error = function(e) NULL)
    
    if (is.null(kmeans_result)) return(NULL)
    
    # Get assignments
    assignments_df <- tibble(
      MCMC_Iteration = iter_val,
      Pathogen_Name = current_iter_data_scaled$Pathogen_Name,
      Cluster_Assigned = kmeans_result$cluster
    )
    
    # Calculate unscaled centroids
    unscaled_centroids_df <- original_data %>%
      filter(MCMC_Iteration == iter_val) %>%
      mutate(Cluster_Assigned_Iter = kmeans_result$cluster) %>%
      group_by(Cluster_Assigned_Iter) %>%
      summarise(
        across(all_of(numerical_param_cols), ~mean(.x, na.rm = TRUE)),
        across(starts_with("Route_"), ~mean(.x, na.rm = TRUE)),
        .groups = 'drop'
      ) %>%
      rename(Cluster_ID_Iter = Cluster_Assigned_Iter) %>%
      mutate(MCMC_Iteration = iter_val)
      
    return(list(assignments = assignments_df, centroids = unscaled_centroids_df))
  }
  
  mcmc_iterations <- unique(scaled_clustering_data_full$MCMC_Iteration)
  
  # Use future_map for parallel execution
  print("Starting parallel per-iteration K-means clustering...")
  clustering_results_list <- future_map(
    mcmc_iterations,
    ~perform_kmeans_for_iteration(.x, scaled_clustering_data_full, clustering_data_full),
    .options = furrr_options(seed = CLUSTERING_SEED),
    .progress = TRUE
  )
  
  # --- Post-processing results from parallel execution ---
  # Filter out any NULL results from failed iterations and separate assignments/centroids
  successful_results <- purrr::compact(clustering_results_list)
  
  if (length(successful_results) > 0) {
    all_iteration_assignments <- purrr::map_dfr(successful_results, "assignments")
    all_iteration_centroids <- purrr::map_dfr(successful_results, "centroids")

    # --- 8.3. Summarize Cluster Centroids (Mean and 95% CI) ---
    feature_cols_for_centroids <- names(all_iteration_centroids)[!names(all_iteration_centroids) %in% c("Cluster_ID_Iter", "MCMC_Iteration")]
    
    if (length(feature_cols_for_centroids) > 0 && nrow(all_iteration_centroids) > 0) {
        centroid_summary <- all_iteration_centroids %>%
          group_by(Cluster_ID_Iter) %>%
          summarise(
            across(all_of(feature_cols_for_centroids), list(
              mean = ~mean(.x, na.rm = TRUE),
              lowerCI = ~quantile(.x, 0.025, na.rm = TRUE),
              upperCI = ~quantile(.x, 0.975, na.rm = TRUE)
            )),
            N_iterations_in_summary = n(),
            .groups = 'drop'
          )
        
        print("Summary of Cluster Centroids (Mean and 95% CI across MCMC iterations):")
        print(as.data.frame(centroid_summary))
        summary_filename <- paste0("Clustering/mcmc/Kmeans/S11_outputs/cluster_centroids_summary_with_ci_k", CHOSEN_K, ".csv")
        write_csv(centroid_summary, summary_filename)
        print(paste("Cluster centroids summary saved to", summary_filename))
    } else {
        print("No feature columns or data available for centroid summary.")
    }

    # --- 8.4. Determine Modal Cluster Assignment for each Pathogen ---
    if (nrow(all_iteration_assignments) > 0) {
        modal_assignments <- all_iteration_assignments %>%
          group_by(Pathogen_Name, Cluster_Assigned) %>%
          summarise(count = n(), .groups = 'drop_last') %>%
          slice_max(order_by = count, n = 1, with_ties = FALSE) %>%
          ungroup() %>%
          select(Pathogen_Name, Modal_Cluster = Cluster_Assigned, Modal_Cluster_Count = count)
          
        print("Modal cluster assignments for each pathogen:")
        print(modal_assignments)
        modal_filename <- paste0("Clustering/mcmc/Kmeans/S11_outputs/pathogen_modal_cluster_assignments_k", CHOSEN_K, ".csv")
        write_csv(modal_assignments, modal_filename)
        print(paste("Modal cluster assignments saved to", modal_filename))
    } else {
        print("No assignment data available for modal cluster calculation.")
        modal_assignments <- data.frame(Pathogen_Name = character(), Modal_Cluster = integer(), Modal_Cluster_Count = integer())
    }

      # --- 8.5. Visualization (using modal assignments and overall data) ---
      pathogen_summary_for_plot <- all_mcmc_samples_df %>%
          group_by(Pathogen_Name) %>%
          summarise(
              R0_mean_overall = mean(R0_sampled, na.rm = TRUE),
              CFR_mean_overall = mean(CFR_sampled, na.rm = TRUE),
              IP_mean_overall = mean(IP_sampled, na.rm = TRUE),
              Route_resp_mean = mean(Route_resp, na.rm=TRUE), # Mean for binary gives proportion
              Route_direct_mean = mean(Route_direct, na.rm=TRUE),
              Route_sexual_mean = mean(Route_sexual, na.rm=TRUE),
              Route_animal_mean = mean(Route_animal, na.rm=TRUE),
              Route_vector_mean = mean(Route_vector, na.rm=TRUE),
              .groups = 'drop'
          ) %>%
          left_join(modal_assignments, by = "Pathogen_Name")

      plot_features_numerical <- c("R0_mean_overall", "CFR_mean_overall", "IP_mean_overall")
      plot_features_routes <- c("Route_resp_mean", "Route_direct_mean", "Route_sexual_mean", "Route_animal_mean", "Route_vector_mean")

      features_for_pca_plot <- pathogen_summary_for_plot %>%
          select(all_of(plot_features_numerical), all_of(plot_features_routes))
      
      # Scale these summarized numerical features for PCA (routes are proportions 0-1)
      scaled_features_for_pca_plot <- features_for_pca_plot
      scaled_features_for_pca_plot[plot_features_numerical] <- scale(scaled_features_for_pca_plot[plot_features_numerical])

      # --- Check for infinite or NA values before running PCA ---
      if (any(!is.finite(as.matrix(scaled_features_for_pca_plot)))) {
          warning("Skipping PCA plot due to non-finite values (NA, NaN, Inf) in the scaled data for plotting. This can happen if imputation or scaling fails due to excessive NAs.")
          scaled_features_for_pca_plot <- features_for_pca_plot # Reset to unscaled for safety in subsequent checks
      }


      if (!requireNamespace("ggplot2", quietly = TRUE)) {
          print("Package 'ggplot2' is not installed.")
      } else if (!requireNamespace("factoextra", quietly = TRUE)){
          print("Package 'factoextra' is not installed.")
      } else if (!requireNamespace("ggrepel", quietly = TRUE)){
          print("Package 'ggrepel' is not installed. Please install it: install.packages('ggrepel')")
      } else if (nrow(scaled_features_for_pca_plot) > 1 && !any(is.na(pathogen_summary_for_plot$Modal_Cluster)) && all(is.finite(as.matrix(scaled_features_for_pca_plot)))) {
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
          ggsave(paste0("Clustering/mcmc/Kmeans/S11_outputs/pca_modal_cluster_plot_k",CHOSEN_K,".png"), plot = pca_plot_modal, width=10, height=8)
          print(paste0("PCA plot with modal cluster assignments saved to Clustering/mcmc/Kmeans/S11_outputs/pca_modal_cluster_plot_k",CHOSEN_K,".png"))
          
          # Silhouette plot for an example MCMC iteration
          if (length(mcmc_iterations) > 0) {
            example_iter_val <- mcmc_iterations[1]
            example_iter_data <- scaled_clustering_data_full %>% filter(MCMC_Iteration == example_iter_val)
            example_iter_features <- example_iter_data %>% select(-MCMC_Iteration, -Pathogen_Name)
            
            if (nrow(example_iter_features) >= CHOSEN_K) {
                set.seed(CLUSTERING_SEED + as.integer(example_iter_val))
                kmeans_example_iter_result <- tryCatch({
                    kmeans(example_iter_features, centers = CHOSEN_K, nstart = 25)
                }, error = function(e) NULL)

                if(!is.null(kmeans_example_iter_result)){
                    if (!requireNamespace("cluster", quietly = TRUE)) {
                        print("Package 'cluster' is not installed, cannot generate silhouette plot.")
                    } else {
                        library(cluster)
                        sil <- silhouette(kmeans_example_iter_result$cluster, dist(example_iter_features))
                        sil_plot_example_iter <- fviz_silhouette(sil, ggtheme = theme_minimal()) +
                                              labs(title=paste("Silhouette Plot for K=", CHOSEN_K, "(Example MCMC Iteration:", example_iter_val, ")"))
                        print(sil_plot_example_iter)
                        ggsave(paste0("Clustering/mcmc/Kmeans/S11_outputs/silhouette_plot_k",CHOSEN_K,"_iter",example_iter_val,".png"), plot = sil_plot_example_iter, width=8, height=6)
                        print(paste0("Silhouette plot for K=", CHOSEN_K, " saved for example iteration ", example_iter_val))
                    }
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
if (exists("all_iteration_assignments") && nrow(all_iteration_assignments) > 0 && 
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
    "Mpox" = "MPV", 
    "Lassa" = "LASV",
    "Nipah" = "NiV",
    "Zika" = "ZIKV",
    "SARS" = "SARS-CoV-1",
    "MERS" = "MERS-CoV",
    "CCHF" = "CCHFV",
    "HIV" = "HIV",
    "Cholera" = "V. cholerae",
    "Measles" = "Measles virus",
    "Norovirus" = "Norovirus",
    "SFTS" = "SFTS virus",
    "EVA71" = "Enterovirus A71",
    "hMPV" = "hMPV",
    "Plague" = "Y. pestis",
    "Smallpox" = "Variola virus",
    "Chikungunya" = "Chikungunya virus",
    "RVFV" = "RVFV",
    "Anthrax" = "B. anthracis"
  )
  
  pathogen_names <- sort(unique(all_iteration_assignments$Pathogen_Name))
  num_pathogens <- length(pathogen_names)
  num_mcmc_iterations <- length(unique(all_iteration_assignments$MCMC_Iteration))

  # --- 9.1. Construct Co-assignment Matrix ---
  print("Constructing co-assignment matrix...")
  coassignment_matrix <- matrix(0, nrow = num_pathogens, ncol = num_pathogens,
                                dimnames = list(pathogen_names, pathogen_names))

  # Efficiently get assignments per iteration
  assignments_by_iter <- split(all_iteration_assignments[, c("Pathogen_Name", "Cluster_Assigned")], 
                               all_iteration_assignments$MCMC_Iteration)

  for (iter_idx in seq_along(assignments_by_iter)) {
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
  print(head(coassignment_matrix[,seq_len(min(6, ncol(coassignment_matrix)))]))
  write.csv(coassignment_matrix, "Clustering/mcmc/Kmeans/S11_outputs/coassignment_matrix.csv")

  # --- 9.2. Calculate Dissimilarity Matrix ---
  print("Calculating dissimilarity matrix...")
  # N_MCMC_ITERATIONS is defined at the top of the script
  dissimilarity_matrix <- 1 - (coassignment_matrix / N_MCMC_ITERATIONS) 
  # Ensure diagonal is 0 after division, if any NA from 0/0 (though N_MCMC_ITERATIONS is likely >0)
  diag(dissimilarity_matrix) <- 0 

  print("Dissimilarity matrix (first 6x6):")
  print(head(dissimilarity_matrix[,seq_len(min(6, ncol(dissimilarity_matrix)))]))
  write.csv(dissimilarity_matrix, "Clustering/mcmc/Kmeans/S11_outputs/dissimilarity_matrix.csv")

  # --- 9.3. Perform Hierarchical Clustering ---
  print("Performing hierarchical clustering...")
  # Convert to dist object for hclust
  pathogen_dist <- as.dist(dissimilarity_matrix)
  hierarchical_clust <- hclust(pathogen_dist, method = "average") # Changed back to average

  # --- 9.3a. Find Optimal K for Consensus Clusters (Silhouette Method) ---
  print("--- Starting Optimal K Analysis (Silhouette Method) ---")
  if (!requireNamespace("cluster", quietly = TRUE)) { install.packages("cluster", quiet = TRUE) }
  library(cluster)

  k_range <- 2:10 # Range of K to test
  avg_silhouette_widths <- numeric(length(k_range))
  names(avg_silhouette_widths) <- k_range
  
  print(paste("Calculating average silhouette width for K from", min(k_range), "to", max(k_range), "..."))
  
  for (k in k_range) {
    cluster_assignments <- cutree(hierarchical_clust, k = k)
    # The silhouette function can take the original distance matrix directly.
    silhouette_info <- silhouette(cluster_assignments, dmatrix = dissimilarity_matrix)
    
    if (is.null(silhouette_info) || !is.matrix(silhouette_info) || nrow(silhouette_info) == 0) {
        avg_width <- NA
        warning(paste("Could not compute silhouette info for k =", k))
    } else {
        avg_width <- summary(silhouette_info)$avg.width
    }
    avg_silhouette_widths[as.character(k)] <- avg_width
  }
  
  silhouette_results_df <- data.frame(
    K = k_range,
    Average_Silhouette_Width = avg_silhouette_widths
  )
  
  print("Silhouette analysis complete. Results:")
  print(silhouette_results_df)
  
  # Find and report the optimal K
  if(any(!is.na(silhouette_results_df$Average_Silhouette_Width))){
      optimal_k_value <- silhouette_results_df$K[which.max(silhouette_results_df$Average_Silhouette_Width)]
      print(paste("Optimal K based on silhouette analysis:", optimal_k_value))
  } else {
      optimal_k_value <- NA
      print("Could not determine optimal K from silhouette analysis.")
  }

  # Plot the results
  optimal_k_plot <- ggplot(silhouette_results_df, aes(x = K, y = Average_Silhouette_Width)) +
    geom_line(color = "steelblue", size = 1) +
    geom_point(color = "steelblue", size = 3) +
    geom_point(data = silhouette_results_df[which.max(silhouette_results_df$Average_Silhouette_Width), ],
               aes(x = K, y = Average_Silhouette_Width), color = "red", size = 5, shape = 8) +
    scale_x_continuous(breaks = k_range) +
    labs(
      title = "",
      subtitle = "",
      x = "Number of Clusters (K)",
      y = "Average Silhouette Width",
      tag = "A"
    ) +
    theme_minimal(base_size = 14) +
    theme(plot.title = element_text(face = "bold"))
  
  
  # Save the results data
  results_filename_optimal_k <- "Clustering/mcmc/Kmeans/S11_outputs/optimal_k_silhouette_results.csv"
  write.csv(silhouette_results_df, results_filename_optimal_k, row.names = FALSE)
  print(paste("Optimal K analysis results saved to", results_filename_optimal_k))
  print("--- Finished Optimal K Analysis ---")

  # --- 9.4. Plot Dendrogram (Revised with dendextend) ---
  print("Plotting enhanced dendrogram with dendextend...")
  
  # Define the number of consensus clusters (K) and a color palette before plotting.
  # This value is also used in section 9.5 to cut the tree.
  if (!is.na(optimal_k_value)) {
    K_CONSENSUS <- optimal_k_value
    print(paste("Using optimal K=", K_CONSENSUS, "for consensus clustering."))
  } else {
    K_CONSENSUS <- 4 # Fallback value
    warning("Optimal K could not be determined. Falling back to K=4 for consensus clustering.")
  }
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
  
  # Convert the base R dendrogram plot to a ggplot object
  if (!requireNamespace("ggplotify", quietly = TRUE)) { install.packages("ggplotify", quiet = TRUE) }
  library(ggplotify)
  dendro_grob <- as.grob(~{
    par(mar = c(5, 4, 4, 10)) # Adjust right margin for labels
    plot(dend_colored, horiz = TRUE, 
         main = "", 
         xlab = "Dissimilarity (1 - Proportion Co-assigned)")
  })
  dendro_ggplot <- as.ggplot(dendro_grob) + labs(tag = "B")

  # --- 9.4b. Combine plots with patchwork and save ---
  if (!requireNamespace("patchwork", quietly = TRUE)) { install.packages("patchwork", quiet = TRUE) }
  library(patchwork)
  
  combined_plot <- optimal_k_plot / dendro_ggplot + plot_layout(heights = c(1, 1.5))
  
  combined_filename <- "Clustering/mcmc/Kmeans/figures/figure.S11.png"
  ggsave(combined_filename, plot = combined_plot, width = 8, height = 10, bg = "white")
  print(paste("Combined plot saved to", combined_filename))


  # --- 9.5. Extract Consensus Cluster Assignments ---
  # User should inspect the dendrogram to choose K_consensus
  # K_CONSENSUS <- 6 # This is now defined in section 9.4 before the dendrogram plot.
  print(paste("Cutting tree to get", K_CONSENSUS, "consensus clusters..."))
  consensus_clusters <- cutree(hierarchical_clust, k = K_CONSENSUS)
  consensus_assignments_df <- data.frame(Pathogen_Name = names(consensus_clusters),
                                         Consensus_Cluster = consensus_clusters)
  print("Consensus cluster assignments:")
  print(consensus_assignments_df)
  write.csv(consensus_assignments_df, paste0("Clustering/mcmc/Kmeans/S11_outputs/consensus_cluster_assignments_k", K_CONSENSUS, ".csv"))

  # --- 9.6. Characterize Consensus Clusters ---
  print("Characterizing consensus clusters (means and 95% CIs from all_mcmc_samples_df)...")
  
  # Select the numerical and route columns to summarize from all_mcmc_samples_df
  # These are the original sampled values, not the means from previous steps.
  params_to_summarize <- c("R0_sampled", "CFR_sampled", "IP_sampled",
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
      write_csv(final_consensus_summary_df, paste0("Clustering/mcmc/Kmeans/S11_outputs/consensus_clusters_summary_k", K_CONSENSUS, ".csv"))
      print(paste0("Consensus clusters summary saved to Clustering/mcmc/Kmeans/S11_outputs/consensus_clusters_summary_k", K_CONSENSUS, ".csv"))
    } else {
      print("No parameters were summarized for consensus clusters.")
    }
  } else {
    print("Could not prepare data for characterizing consensus clusters (pathogen_mcmc_samples_with_consensus_clusters is empty or missing Consensus_Cluster column).")
  }
  
  print("Ensemble/Consensus Clustering Analysis finished.")

} else {
  print("Skipping Ensemble/Consensus Clustering: 'all_iteration_assignments' or 'all_mcmc_samples_df' not found or empty.")
}
