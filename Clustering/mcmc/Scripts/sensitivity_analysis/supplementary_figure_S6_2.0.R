# MCMC K-Means Clustering Analysis for Supplementary Figure S6 2.0

# --- 1. Load Libraries ---
library(dplyr)
library(readr)
# Add any other necessary libraries here, e.g., for specific distributions or optimization
# library(DistributionFitR) # May be useful for fitting distributions to CI


# --- 2. Configuration ---
N_MCMC_ITERATIONS <- 500 # Number of MCMC iterations 
N_PRESYMP_SAMPLES <- 5000 # Number of samples for presymptomatic proportion estimation 
set.seed(123) # For reproducibility

# Define output directory for this specific analysis
output_dir <- "Clustering/mcmc/Kmeans/S6_2.0_outputs"
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)


# --- 3. Load Data ---
# Load the new long-format parameter data and transmission routes
params_long_df <- read_csv("Clustering/mcmc/Kmeans/data/pathogen_params.csv")
transmission_df <- read_csv("Clustering/mcmc/Kmeans/data/transmission_route.csv")
presym_dist_df <- read_csv("Clustering/mcmc/Kmeans/data/pathogen_presym.csv")

# Filter for the specific pathogens of interest for this supplementary figure
pathogens_of_interest <- c(
  "COVID-19_WT", "COVID-19_D", "COVID-19_O", "Ebola", 
  "Marburg", "Mpox", "H1N1_18", "H1N1_09", "SARS", "MERS"
)
params_long_df <- params_long_df %>% filter(Pathogen_Name %in% pathogens_of_interest)
transmission_df <- transmission_df %>% filter(Pathogen_Name %in% pathogens_of_interest)
presym_dist_df <- presym_dist_df %>% filter(Pathogen_Name %in% pathogens_of_interest)


# Optional: Quick check of the data
# print(head(params_long_df))
# print(str(params_long_df))
# print(head(transmission_df))


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

# --- NEW BOOTSTRAP AGGREGATION SAMPLING FUNCTION ---
# This function uses bootstrap aggregation to create a robust parameter estimate for each MCMC iteration.
# It addresses concerns about overemphasizing single studies (whether outliers or precise-but-biased)
# by synthesizing evidence from all available studies in each sampling step.
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
        if (!is.na(mean_val) && mean_val > 0 && mean_val < 1) {
            beta_params <- get_beta_params_from_mean_ci(mean_val, lower_ci, upper_ci)
            if (!is.na(beta_params$shape1) && !is.na(beta_params$shape2)) {
                sampled_val_internal <- rbeta(1, shape1 = beta_params$shape1, shape2 = beta_params$shape2)
            }
        }
    } else if (param_name %in% c("R0", "SI")) {
        gamma_params <- get_gamma_params_from_mean_ci(mean_val, lower_ci, upper_ci)
        if (!is.na(gamma_params$shape) && !is.na(gamma_params$rate)) {
            sampled_val_internal <- rgamma(1, shape = gamma_params$shape, rate = gamma_params$rate)
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
  } else if (param_name %in% c("R0", "SI")) {
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

  random_quantiles <- runif(n_samples)
  
  si_samples <- q_si(random_quantiles)
  ip_samples <- q_ip(random_quantiles)
  
  if(si_dist_info$family == "Normal") si_samples[si_samples < 0] <- 0
  if(ip_dist_info$family == "Normal") ip_samples[ip_samples < 0] <- 0

  mean(si_samples < ip_samples, na.rm = TRUE)
}


# --- 5. Main MCMC Loop ---
mcmc_results <- list()
unique_pathogens <- unique(params_long_df$Pathogen_Name)
flu_group <- c("H1N1_09", "H1N1_18", "H2N2", "H3N2")

for (iter in 1:N_MCMC_ITERATIONS) {
  if (iter %% max(1, (N_MCMC_ITERATIONS/10)) == 0) {
    print(paste("MCMC Iteration:", iter, "/", N_MCMC_ITERATIONS))
  }
  
  current_iteration_params_list <- list()
  
  for (pathogen_name in unique_pathogens) {
    
    # --- Sample parameters ---
    r0_sampled <- sample_parameter_bootstrap_aggregation("R0", pathogen_name, params_long_df)
    si_clust_sampled <- sample_parameter_bootstrap_aggregation("SI", pathogen_name, params_long_df)
    cfr_sampled <- sample_parameter_bootstrap_aggregation("CFR", pathogen_name, params_long_df)
    ip_sampled <- sample_parameter_bootstrap_aggregation("IP", pathogen_name, params_long_df)
    lp_sampled <- sample_parameter_bootstrap_aggregation("LP", pathogen_name, params_long_df)
    infp_sampled <- sample_parameter_bootstrap_aggregation("InfP", pathogen_name, params_long_df)
    k_sampled <- sample_parameter_bootstrap_aggregation("k", pathogen_name, params_long_df)
    
    # --- Estimate presymptomatic proportion ---
    presymp_prop_sampled <- NA
    
    if (pathogen_name %in% flu_group) {
      si_studies <- presym_dist_df %>% filter(Pathogen_Name %in% flu_group, Parameter == "SI")
      ip_studies <- presym_dist_df %>% filter(Pathogen_Name %in% flu_group, Parameter == "IP")
    } else {
      si_studies <- presym_dist_df %>% filter(Pathogen_Name == pathogen_name, Parameter == "SI")
      ip_studies <- presym_dist_df %>% filter(Pathogen_Name == pathogen_name, Parameter == "IP")
    }
    
    if (nrow(si_studies) > 0 && nrow(ip_studies) > 0) {
      si_study_sampled <- si_studies[sample(seq_len(nrow(si_studies)), 1), ]
      ip_study_sampled <- ip_studies[sample(seq_len(nrow(ip_studies)), 1), ]
      
      si_full_dist_info <- get_full_dist_params_from_presym(si_study_sampled)
      ip_full_dist_info <- get_full_dist_params_from_presym(ip_study_sampled)
      
      if (!is.null(si_full_dist_info) && !is.null(ip_full_dist_info) && 
          !is.na(si_full_dist_info$family) && !is.na(ip_full_dist_info$family) &&
          length(si_full_dist_info$params) > 0 && length(ip_full_dist_info$params) > 0) {
        
        si_full_dist_info$Pathogen_Name <- pathogen_name
        ip_full_dist_info$Pathogen_Name <- pathogen_name
        presymp_prop_sampled <- estimate_presymptomatic_proportion(si_full_dist_info, ip_full_dist_info)
      } else {
        warning(paste("Could not estimate presymptomatic proportion for", pathogen_name, "due to incomplete SI/IP dist info."))
      }
    } else {
        warning(paste("Not enough SI or IP studies to estimate presymptomatic proportion for", pathogen_name))
    }
    
    # --- Get transmission routes ---
    pathogen_routes <- transmission_df %>% filter(Pathogen_Name == pathogen_name)
    if(nrow(pathogen_routes) == 0){
        warning(paste("No transmission route data for pathogen:", pathogen_name))
        pathogen_routes <- data.frame(Route_resp=NA, Route_direct=NA, Route_sexual=NA, Route_animal=NA, Route_vector=NA)
    }

    pathogen_params_df <- data.frame(
      Pathogen_Name = pathogen_name,
      R0_sampled = ifelse(is.null(r0_sampled) || is.na(r0_sampled) || !is.finite(r0_sampled), NA_real_, r0_sampled),
      SI_Clust_sampled = ifelse(is.null(si_clust_sampled) || is.na(si_clust_sampled) || !is.finite(si_clust_sampled), NA_real_, si_clust_sampled),
      CFR_sampled = ifelse(is.null(cfr_sampled) || is.na(cfr_sampled) || !is.finite(cfr_sampled), NA_real_, cfr_sampled),
      Presymp_Proportion_sampled = ifelse(is.null(presymp_prop_sampled) || is.na(presymp_prop_sampled) || !is.finite(presymp_prop_sampled), NA_real_, presymp_prop_sampled),
      IP_sampled = ifelse(is.null(ip_sampled) || is.na(ip_sampled) || !is.finite(ip_sampled), NA_real_, ip_sampled),
      LP_sampled = ifelse(is.null(lp_sampled) || is.na(lp_sampled) || !is.finite(lp_sampled), NA_real_, lp_sampled),
      InfP_sampled = ifelse(is.null(infp_sampled) || is.na(infp_sampled) || !is.finite(infp_sampled), NA_real_, infp_sampled),
      k_sampled = ifelse(is.null(k_sampled) || is.na(k_sampled) || !is.finite(k_sampled), NA_real_, k_sampled),
      Route_resp = as.integer(pathogen_routes$Route_resp),
      Route_direct = as.integer(pathogen_routes$Route_direct),
      Route_sexual = as.integer(pathogen_routes$Route_sexual),
      Route_animal = as.integer(pathogen_routes$Route_animal),
      Route_vector = as.integer(pathogen_routes$Route_vector)
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
    all_mcmc_samples_df <- bind_rows(mcmc_results[!sapply(mcmc_results, is.null)], .id = "MCMC_Iteration")
    if(nrow(all_mcmc_samples_df) > 0){
        all_mcmc_samples_df$MCMC_Iteration <- as.integer(all_mcmc_samples_df$MCMC_Iteration)

        print("MCMC sampling complete.")
        print(paste("Generated", length(unique(all_mcmc_samples_df$MCMC_Iteration)), "valid sets of parameters for", length(unique(all_mcmc_samples_df$Pathogen_Name)), "pathogens."))
    
        # --- Save Results ---
        write_csv(all_mcmc_samples_df, file.path(output_dir, "mcmc_parameter_samples.csv"))
        print(paste("MCMC parameter samples saved to", file.path(output_dir, "mcmc_parameter_samples.csv")))
    } else {
        print("MCMC sampling completed, but no rows in the final combined data frame.")
    }
} else {
    print("MCMC sampling did not produce any results or all iterations resulted in NULL.")
}

print("Script finished.")


# --- 8. Clustering Analysis ---
if (exists("all_mcmc_samples_df") && nrow(all_mcmc_samples_df) > 0) {
  
  CHOSEN_K <- 4 #specified K=4
  
  print(paste0("--- Starting Clustering Analysis (Per MCMC Iteration, K=", CHOSEN_K, ") ---"))
  
  # --- 8.1. Prepare Data for Clustering ---
  clustering_data_full <- all_mcmc_samples_df %>%
    select(
      MCMC_Iteration,
      Pathogen_Name,
      R0_sampled, 
      SI_Clust_sampled, 
      CFR_sampled, 
      Presymp_Proportion_sampled,
      IP_sampled,
      LP_sampled,
      InfP_sampled,
      k_sampled,
      Route_resp,
      Route_direct,
      Route_sexual,
      Route_animal,
      Route_vector
    )

  numerical_param_cols <- c("R0_sampled", "SI_Clust_sampled", "CFR_sampled", "Presymp_Proportion_sampled", "IP_sampled", "LP_sampled", "InfP_sampled", "k_sampled")
  
  # --- IMPORTANT: Do NOT impute NA values. Leave them as NA. ---
  
  # Scale numerical features across the entire dataset
  scaled_clustering_data_full <- clustering_data_full
  cols_to_scale_exist <- numerical_param_cols[numerical_param_cols %in% names(scaled_clustering_data_full)]
  if(length(cols_to_scale_exist) > 0) {
      scaled_clustering_data_full[cols_to_scale_exist] <- scale(scaled_clustering_data_full[cols_to_scale_exist])
  } else {
      warning("No numerical columns found for scaling.")
  }

  print("Scaled features for per-iteration clustering (first 6 rows of full dataset):")
  print(head(scaled_clustering_data_full))

  # Lists to store results from each iteration
  all_iteration_centroids <- list()
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
      
      current_iter_features_scaled <- current_iter_data %>%
        select(-MCMC_Iteration, -Pathogen_Name) %>%
        na.omit() # Remove rows with NA values for this iteration's clustering
          
      if (nrow(current_iter_features_scaled) < CHOSEN_K) {
          warning(paste("Skipping K-means for iteration", iter_val, "due to insufficient data points after NA removal."))
          next
      }

      set.seed(123 + as.integer(iter_val)) 
      kmeans_iter_result <- tryCatch({
          kmeans(current_iter_features_scaled, centers = CHOSEN_K, nstart = 25)
      }, error = function(e) {
          warning(paste("Error in kmeans for iteration", iter_val, ":", e$message))
          NULL
      })
      
      if (!is.null(kmeans_iter_result)) {
          # Get pathogen names corresponding to the rows that were NOT omitted
          pathogen_names_for_iter <- current_iter_data %>% na.omit() %>% pull(Pathogen_Name)
          
          assignments_df <- data.frame(
              MCMC_Iteration = iter_val,
              Pathogen_Name = pathogen_names_for_iter,
              Cluster_Assigned = kmeans_iter_result$cluster
          )
          all_iteration_assignments[[as.character(iter_val)]] <- assignments_df

          original_data_this_iter_for_centroids <- all_mcmc_samples_df %>%
            filter(MCMC_Iteration == iter_val, Pathogen_Name %in% pathogen_names_for_iter) %>%
            mutate(Cluster_Assigned_Iter = kmeans_iter_result$cluster)
          
          iter_unscaled_centroids <- original_data_this_iter_for_centroids %>%
            group_by(Cluster_Assigned_Iter) %>%
            summarise(
              across(all_of(c(numerical_param_cols, "Route_resp", "Route_direct", "Route_sexual", "Route_animal", "Route_vector")), ~mean(.x, na.rm = TRUE)),
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
              N_iterations_in_summary = n(),
              .groups = 'drop'
            )
          
          print("Summary of Cluster Centroids (Mean and 95% CI across MCMC iterations):")
          print(as.data.frame(centroid_summary))
          summary_filename <- file.path(output_dir, paste0("cluster_centroids_summary_with_ci_k", CHOSEN_K, ".csv"))
          write_csv(centroid_summary, summary_filename)
          print(paste("Cluster centroids summary saved to", summary_filename))
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
          modal_filename <- file.path(output_dir, paste0("pathogen_modal_cluster_assignments_k", CHOSEN_K, ".csv"))
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
              across(all_of(c(numerical_param_cols, "Route_resp", "Route_direct", "Route_sexual", "Route_animal", "Route_vector")),
                     ~mean(.x, na.rm = TRUE),
                     .names = "{.col}_mean_overall"),
              .groups = 'drop'
          ) %>%
          left_join(modal_assignments, by = "Pathogen_Name")

      plot_features_numerical <- paste0(numerical_param_cols, "_mean_overall")
      plot_features_routes <- c("Route_resp_mean_overall", "Route_direct_mean_overall", "Route_sexual_mean_overall", "Route_animal_mean_overall", "Route_vector_mean_overall")

      features_for_pca_plot <- pathogen_summary_for_plot %>%
          select(all_of(plot_features_numerical), all_of(plot_features_routes))
      
      scaled_features_for_pca_plot <- features_for_pca_plot
      scaled_features_for_pca_plot[plot_features_numerical] <- scale(scaled_features_for_pca_plot[plot_features_numerical])

      if (any(!is.finite(as.matrix(scaled_features_for_pca_plot)))) {
          warning("Skipping PCA plot due to non-finite values.")
      } else if (requireNamespace("ggplot2", quietly = TRUE) && requireNamespace("factoextra", quietly = TRUE) && requireNamespace("ggrepel", quietly = TRUE)) {
          library(ggplot2)
          library(factoextra)
          library(ggrepel)

          pca_result_modal <- prcomp(scaled_features_for_pca_plot, center = FALSE, scale. = FALSE) 
          pca_data_modal <- as.data.frame(pca_result_modal$x)
          pca_data_modal$Cluster <- factor(pathogen_summary_for_plot$Modal_Cluster, levels = 1:CHOSEN_K)
          pca_data_modal$Pathogen_Name <- pathogen_summary_for_plot$Pathogen_Name

          pca_plot_modal <- ggplot(pca_data_modal, aes(x = PC1, y = PC2, color = Cluster, label = Pathogen_Name)) +
            geom_point(size = 3) +
            geom_text_repel(size = 3.5, max.overlaps = Inf, fontface = "bold") +
            labs(title = paste0("Pathogen Clusters (Modal Assignment, K=", CHOSEN_K, ")"),
                 x = paste0("PC1 (", round(summary(pca_result_modal)$importance[2,1]*100, 1), "%)"),
                 y = paste0("PC2 (", round(summary(pca_result_modal)$importance[2,2]*100, 1), "%)")) +
            theme_minimal(base_size = 12) +
            scale_color_brewer(palette = "Set1", name = "Modal Cluster")
          
          print(pca_plot_modal)
          ggsave(file.path(output_dir, paste0("pca_modal_cluster_plot_k", CHOSEN_K, ".png")), plot = pca_plot_modal, width=10, height=8)
      }

      # --- 9. Ensemble Clustering / Consensus Clustering --- 
      if (!requireNamespace("stats", quietly = TRUE)) { install.packages("stats", quiet = TRUE) }
      if (!requireNamespace("dendextend", quietly = TRUE)) { install.packages("dendextend", quiet = TRUE) }
      library(stats)
      library(dendextend)

      pathogen_names <- sort(unique(combined_assignments_df$Pathogen_Name))
      num_pathogens <- length(pathogen_names)
      num_mcmc_iterations <- length(unique(combined_assignments_df$MCMC_Iteration))

      coassignment_matrix <- matrix(0, nrow = num_pathogens, ncol = num_pathogens,
                                    dimnames = list(pathogen_names, pathogen_names))

      assignments_by_iter <- split(combined_assignments_df[, c("Pathogen_Name", "Cluster_Assigned")], 
                                   combined_assignments_df$MCMC_Iteration)

      for (iter_assignments in assignments_by_iter) {
        for (i in 1:(num_pathogens - 1)) {
          for (j in (i + 1):num_pathogens) {
            p1_name <- pathogen_names[i]
            p2_name <- pathogen_names[j]
            
            p1_assign <- iter_assignments$Cluster_Assigned[iter_assignments$Pathogen_Name == p1_name]
            p2_assign <- iter_assignments$Cluster_Assigned[iter_assignments$Pathogen_Name == p2_name]

            if (length(p1_assign) > 0 && length(p2_assign) > 0 && p1_assign == p2_assign) {
              coassignment_matrix[p1_name, p2_name] <- coassignment_matrix[p1_name, p2_name] + 1
            }
          }
        }
      }
      coassignment_matrix <- coassignment_matrix + t(coassignment_matrix)
      diag(coassignment_matrix) <- num_mcmc_iterations

      write.csv(coassignment_matrix, file.path(output_dir, "coassignment_matrix.csv"))

      dissimilarity_matrix <- 1 - (coassignment_matrix / num_mcmc_iterations)
      diag(dissimilarity_matrix) <- 0
      write.csv(dissimilarity_matrix, file.path(output_dir, "dissimilarity_matrix.csv"))

      pathogen_dist <- as.dist(dissimilarity_matrix)
      hierarchical_clust <- hclust(pathogen_dist, method = "average")

      K_CONSENSUS <- 4
      if (requireNamespace("RColorBrewer", quietly = TRUE)) {
          library(RColorBrewer)
          cluster_colors <- brewer.pal(n = max(3, K_CONSENSUS), name = "Set2")

          dend <- as.dendrogram(hierarchical_clust)
          
          pathogen_full_name_map <- c("COVID-19_WT" = "SARS-CoV-2 (WT)", "COVID-19_D" = "SARS-CoV-2 (Delta)", "COVID-19_O" = "SARS-CoV-2 (Omicron)", "Ebola" = "EBOV", "Marburg" = "MARV", "Mpox" = "MPV", "H1N1_18" = "A/H1N1", "H1N1_09" = "A/H1N1/09", "SARS" = "SARS-CoV-1", "MERS" = "MERS-CoV")
          labels(dend) <- pathogen_full_name_map[labels(dend)]
          
          dend_colored <- color_branches(dend, k = K_CONSENSUS, col = cluster_colors)
          
          png(file.path(output_dir, "figure.S6_2.0.png"), width=1200, height=800, units="px", res=100)
          par(mar = c(5, 4, 4, 10))
          plot(dend_colored, horiz = TRUE, main = "", xlab = "Dissimilarity (1 - Proportion Co-assigned)")
          dev.off()
      }

      consensus_clusters <- cutree(hierarchical_clust, k = K_CONSENSUS)
      consensus_assignments_df <- data.frame(Pathogen_Name = names(consensus_clusters), Consensus_Cluster = consensus_clusters)
      write.csv(consensus_assignments_df, file.path(output_dir, paste0("consensus_cluster_assignments_k", K_CONSENSUS, ".csv")))
  }
}
# The rest of the script (consensus clustering, etc.) should be pasted here,
# ensuring all file output paths use the `output_dir` variable.
