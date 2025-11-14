# MCMC K-Means Clustering Analysis - Bootstrap First Approach
# This script implements an alternative methodology where bootstrapping of studies
# is performed *before* the MCMC-style parameter sampling begins. This allows for
# assessing the stability of the entire clustering pipeline to study selection.

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
  "cluster", "dendextend", "RColorBrewer", "future", "furrr"
)
install_and_load(required_packages)

# Setup for parallel processing
plan(multisession)


# --- 2. Configuration ---
# Reduced numbers for feasible test execution
N_BOOTSTRAP_REPLICATES <- 10 # Number of bootstrap replicates
N_MCMC_ITERATIONS <- 100    # Number of MCMC iterations per bootstrap
N_PRESYMP_SAMPLES <- 1000   # Number of samples for presymptomatic proportion estimation
MCMC_SEED <- 123            # Base seed for reproducibility
CLUSTERING_SEED <- 456      # Seed for reproducibility of K-means

set.seed(MCMC_SEED)


# --- 3. Load Data ---
# --- NEW: 3. Output Directory Setup ---
main_output_dir <- "Clustering/mcmc/Kmeans/bootstrap_first_outputs"
if (!dir.exists(main_output_dir)) {
  dir.create(main_output_dir, recursive = TRUE)
  print(paste("Created output directory:", main_output_dir))
}
params_long_df <- read_csv("Clustering/mcmc/Kmeans/data/pathogen_params.csv")
transmission_df <- read_csv("Clustering/mcmc/Kmeans/data/transmission_route.csv")
presym_dist_df <- read_csv("Clustering/mcmc/Kmeans/data/pathogen_presym.csv")

# Filter for the specific pathogens of interest
pathogens_of_interest <- c(
  "COVID-19_WT", "COVID-19_D", "COVID-19_O", "Ebola", "Marburg", "Mpox", 
  "H1N1_18", "H1N1_09", "SARS", "MERS"
)
params_long_df <- params_long_df %>% filter(Pathogen_Name %in% pathogens_of_interest)
transmission_df <- transmission_df %>% filter(Pathogen_Name %in% pathogens_of_interest)
presym_dist_df <- presym_dist_df %>% filter(Pathogen_Name %in% pathogens_of_interest)


# --- 4. Helper Functions ---

# --- Quantile-Matching Helper Functions (get_beta_params_from_ci, get_gamma_params_from_ci, etc.) ---
# These functions are unchanged from the original script.
# (Code for these functions is omitted here for brevity, but would be included in the actual script)
# Function to derive parameters for rbeta by matching quantiles of the 95% CI
get_beta_params_from_ci <- function(lower_ci, upper_ci, mean_val) {
    min_param_val <- 1e-6
    objective_beta <- function(params, lower_ci, upper_ci) {
        shape1 <- exp(params[1]); shape2 <- exp(params[2])
        if (is.na(shape1) || is.na(shape2) || !is.finite(shape1) || !is.finite(shape2)) return(1e10)
        q_lower <- qbeta(0.025, shape1, shape2); q_upper <- qbeta(0.975, shape1, shape2)
        return((q_lower - lower_ci)^2 + (q_upper - upper_ci)^2)
    }
    if (!is.na(mean_val) && mean_val > 0 && mean_val < 1) {
        sd_guess <- (upper_ci - lower_ci) / (2 * qnorm(0.975)); var_guess <- sd_guess^2
        alpha_guess <- ((1 - mean_val) / var_guess - 1 / mean_val) * mean_val^2
        beta_guess <- alpha_guess * (1 / mean_val - 1)
        if (is.na(alpha_guess) || alpha_guess <= 0) alpha_guess <- 1
        if (is.na(beta_guess) || beta_guess <= 0) beta_guess <- 1
        start_params <- log(c(alpha_guess, beta_guess))
    } else { start_params <- log(c(1, 1)) }
    opt_result <- tryCatch({ optim(start_params, objective_beta, lower_ci = lower_ci, upper_ci = upper_ci) }, error = function(e) NULL)
    if (!is.null(opt_result) && opt_result$convergence == 0) {
        return(list(shape1 = exp(opt_result$par[1]), shape2 = exp(opt_result$par[2])))
    } else { return(get_beta_params_from_mean_ci_fallback(mean_val, lower_ci, upper_ci)) }
}
get_gamma_params_from_ci <- function(lower_ci, upper_ci, mean_val) {
    min_param_val <- 1e-6
    objective_gamma <- function(params, lower_ci, upper_ci) {
        shape <- exp(params[1]); rate <- exp(params[2])
        if (is.na(shape) || is.na(rate) || !is.finite(shape) || !is.finite(rate)) return(1e10)
        q_lower <- qgamma(0.025, shape = shape, rate = rate); q_upper <- qgamma(0.975, shape = shape, rate = rate)
        return((q_lower - lower_ci)^2 + (q_upper - upper_ci)^2)
    }
    if(!is.na(mean_val) && mean_val > 0){
        sd_guess <- (upper_ci - lower_ci) / (2 * qnorm(0.975))
        if(is.na(sd_guess) || sd_guess <= 0) sd_guess <- mean_val
        shape_guess <- (mean_val / sd_guess)^2; rate_guess <- mean_val / sd_guess^2
        if (is.na(shape_guess) || shape_guess <= 0) shape_guess <- 1
        if (is.na(rate_guess) || rate_guess <= 0) rate_guess <- 1
        start_params <- log(c(shape_guess, rate_guess))
    } else { start_params <- log(c(1,1)) }
    opt_result <- tryCatch({ optim(start_params, objective_gamma, lower_ci = lower_ci, upper_ci = upper_ci) }, error = function(e) NULL)
    if (!is.null(opt_result) && opt_result$convergence == 0) {
        return(list(shape = exp(opt_result$par[1]), rate = exp(opt_result$par[2])))
    } else { return(get_gamma_params_from_mean_ci_fallback(mean_val, lower_ci, upper_ci)) }
}
get_beta_params_from_mean_ci_fallback <- function(mean_val, lower_ci, upper_ci, n_eff_guess = 1000) {
    alpha <- mean_val * n_eff_guess; beta <- (1 - mean_val) * n_eff_guess
    min_param_val <- 1e-6
    if (is.na(alpha) || !is.finite(alpha) || alpha <= 0) alpha <- min_param_val 
    if (is.na(beta) || !is.finite(beta) || beta <= 0) beta <- min_param_val
    return(list(shape1 = alpha, shape2 = beta))
}
get_gamma_params_from_mean_ci_fallback <- function(mean_val, lower_ci, upper_ci) {
    min_param_val <- 1e-6
    if (is.na(mean_val) || mean_val <= 0) {
        if(!is.na(lower_ci) && !is.na(upper_ci) && lower_ci < upper_ci && (lower_ci + upper_ci)/2 > 0 ){
            mean_val <- (lower_ci + upper_ci) / 2
        } else { return(list(shape = min_param_val, rate = min_param_val)) }
    }
    sd_from_ci <- (upper_ci - lower_ci) / (2 * qnorm(0.975))
    if (is.na(sd_from_ci) || !is.finite(sd_from_ci) || sd_from_ci <= 0) { sd_from_ci <- mean_val * 0.5 }
    shape <- (mean_val / sd_from_ci)^2; rate <- mean_val / sd_from_ci^2
    if (is.na(shape) || !is.finite(shape) || shape <= 0) shape <- min_param_val
    if (is.na(rate) || !is.finite(rate) || rate <= 0) rate <- min_param_val
    return(list(shape = shape, rate = rate))
}
get_full_dist_params_from_presym <- function(study_row) {
    if (is.null(study_row) || nrow(study_row) == 0) return(list(family = NA, params = list()))
    family <- study_row$FullDist_Family; raw_params <- list()
    for (i in 1:2) {
        param_name_col <- paste0("FullDist_Param", i, "_Name"); param_val_col <- paste0("FullDist_Param", i, "_Value")
        if (!is.na(study_row[[param_name_col]]) && study_row[[param_name_col]] != "") {
            param_name <- study_row[[param_name_col]]; param_val <- as.numeric(study_row[[param_val_col]])
            raw_params[[param_name]] <- param_val
        }
    }
    final_params <- raw_params
    if (!is.na(family)) {
        if (family == "Gamma" && all(c("mean", "SD") %in% names(raw_params))) {
            mean_val <- raw_params$mean; sd_val <- raw_params$SD
            if (!is.na(mean_val) && !is.na(sd_val) && sd_val > 0 && mean_val > 0) {
                final_params$shape <- (mean_val / sd_val)^2; final_params$rate <- mean_val / (sd_val^2)
                final_params$mean <- NULL; final_params$SD <- NULL
            }
        } else if (family == "Lognormal" && all(c("mean", "SD") %in% names(raw_params))) {
            mean_val <- raw_params$mean; sd_val <- raw_params$SD
            if(!is.na(mean_val) && !is.na(sd_val) && mean_val > 0 && sd_val > 0){
                final_params$sdlog <- sqrt(log(sd_val^2 / mean_val^2 + 1))
                final_params$meanlog <- log(mean_val) - 0.5 * final_params$sdlog^2
                final_params$mean <- NULL; final_params$SD <- NULL
            }
        }
    }
    return(list(family = family, params = final_params))
}
estimate_presymptomatic_proportion <- function(si_dist_info, ip_dist_info, n_samples = N_PRESYMP_SAMPLES) {
    if (is.null(si_dist_info$family) || is.na(si_dist_info$family) || length(si_dist_info$params) == 0 ||
        is.null(ip_dist_info$family) || is.na(ip_dist_info$family) || length(ip_dist_info$params) == 0) {
        return(NA)
    }
    get_q_func <- function(d_info) {
        params <- d_info$params
        switch(d_info$family,
               "Lognormal" = function(p) qlnorm(p, meanlog = params$meanlog, sdlog = params$sdlog),
               "Gamma"     = function(p) qgamma(p, shape = params$shape, rate = if (!is.null(params$rate)) params$rate else 1/params$scale),
               "Normal"    = function(p) qnorm(p, mean = params$mean, sd = params$SD),
               "Weibull"   = function(p) qweibull(p, shape = params$shape, scale = params$scale),
               NULL)
    }
    q_si <- get_q_func(si_dist_info); q_ip <- get_q_func(ip_dist_info)
    if (is.null(q_si) || is.null(q_ip)) return(NA)
    random_quantiles <- runif(n_samples)
    si_samples <- q_si(random_quantiles); ip_samples <- q_ip(random_quantiles)
    if(si_dist_info$family == "Normal") si_samples[si_samples < 0] <- 0
    if(ip_dist_info$family == "Normal") ip_samples[ip_samples < 0] <- 0
    mean(si_samples < ip_samples, na.rm = TRUE)
}

# --- NEW Parameter Sampling Function (No Inner Bootstrap) ---
# This function generates a single aggregated parameter estimate for one MCMC iteration.
# It works on a pre-bootstrapped dataset of studies.
sample_parameter_aggregate <- function(param_name, pathogen_name, data_df) {
  # 1. Filter for all studies for the given pathogen and parameter from the provided dataframe
  studies <- data_df %>% filter(Pathogen_Name == pathogen_name, Parameter == param_name)
  if (nrow(studies) == 0) return(NA)

  # 2. Generate one sample from EACH study in the (bootstrapped) set
  samples <- apply(studies, 1, function(study_row) {
    mean_val <- as.numeric(study_row["ReportedValue"])
    lower_ci <- as.numeric(study_row["LowerBound"])
    upper_ci <- as.numeric(study_row["UpperBound"])
    
    if (is.na(lower_ci) || is.na(upper_ci) || lower_ci >= upper_ci) return(mean_val)
    
    sampled_val_internal <- NA
    if (param_name == "CFR") {
        if (!is.na(mean_val) && mean_val > 0 && mean_val < 1) {
            beta_params <- get_beta_params_from_ci(lower_ci, upper_ci, mean_val)
            if (!is.na(beta_params$shape1)) sampled_val_internal <- rbeta(1, beta_params$shape1, beta_params$shape2)
        }
    } else if (param_name %in% c("R0", "SI", "k", "IP", "LP", "InfP")) {
        gamma_params <- get_gamma_params_from_ci(lower_ci, upper_ci, mean_val)
        if (!is.na(gamma_params$shape)) sampled_val_internal <- rgamma(1, gamma_params$shape, gamma_params$rate)
    } else {
        sd_from_ci <- (upper_ci - lower_ci) / (2 * qnorm(0.975))
        if (!is.na(mean_val) && !is.na(sd_from_ci) && sd_from_ci > 0) {
            sampled_val_internal <- rnorm(1, mean = mean_val, sd = sd_from_ci)
        }
    }
    
    if (is.na(sampled_val_internal) || !is.finite(sampled_val_internal)) return(mean_val)
    return(sampled_val_internal)
  })

  # 3. Aggregate the samples from all studies by taking the mean
  final_value <- mean(samples, na.rm = TRUE)
  
  # 4. Apply final constraints
  if (is.na(final_value)) return(NA)
  if (param_name == "CFR") final_value <- pmax(0, pmin(1, final_value))
  else if (param_name %in% c("R0", "SI", "k", "IP", "LP", "InfP")) final_value <- pmax(1e-5, final_value)
  
  return(final_value)
}


# --- 5. Main Bootstrap and MCMC Loop ---

# List to store key results from each bootstrap replicate
all_bootstrap_results <- list()

for (b in 1:N_BOOTSTRAP_REPLICATES) {
  
  print(paste("--- Starting Bootstrap Replicate:", b, "/", N_BOOTSTRAP_REPLICATES, "---"))
  
  # --- 5.1. Create a bootstrapped dataset of studies for this replicate ---
  # We perform a stratified bootstrap, resampling studies within each pathogen-parameter group.
  set.seed(MCMC_SEED + b) # Ensure bootstrap is reproducible
  bootstrapped_params_df <- params_long_df %>%
    group_by(Pathogen_Name, Parameter) %>%
    sample_n(size = n(), replace = TRUE) %>%
    ungroup()

  # --- 5.2. MCMC Sampling using the bootstrapped dataset ---
  
  # Helper function to run a single MCMC iteration (now takes a data_df)
  run_single_mcmc_iteration <- function(iter_num, data_df, presym_data_df, transmission_data_df) {
    unique_pathogens <- unique(data_df$Pathogen_Name)
    flu_group <- c("H1N1_09", "H1N1_18", "H2N2", "H3N2")
    
    iteration_results <- purrr::map_dfr(unique_pathogens, function(pathogen_name) {
      r0_sampled <- sample_parameter_aggregate("R0", pathogen_name, data_df)
      si_clust_sampled <- sample_parameter_aggregate("SI", pathogen_name, data_df)
      cfr_sampled <- sample_parameter_aggregate("CFR", pathogen_name, data_df)
      k_sampled <- sample_parameter_aggregate("k", pathogen_name, data_df)
      ip_sampled <- sample_parameter_aggregate("IP", pathogen_name, data_df)
      lp_sampled <- sample_parameter_aggregate("LP", pathogen_name, data_df)
      infp_sampled <- sample_parameter_aggregate("InfP", pathogen_name, data_df)
      
      presymp_prop_sampled <- NA
      # (Presymptomatic estimation logic remains the same as it uses a different dataframe)
      if (pathogen_name %in% flu_group) {
        si_studies <- presym_data_df %>% filter(Pathogen_Name %in% flu_group, Parameter == "SI")
        ip_studies <- presym_data_df %>% filter(Pathogen_Name %in% flu_group, Parameter == "IP")
      } else {
        si_studies <- presym_data_df %>% filter(Pathogen_Name == pathogen_name, Parameter == "SI")
        ip_studies <- presym_data_df %>% filter(Pathogen_Name == pathogen_name, Parameter == "IP")
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
        }
      }
      
      pathogen_routes <- transmission_data_df %>% filter(Pathogen_Name == pathogen_name)
      if(nrow(pathogen_routes) == 0) pathogen_routes <- data.frame(Route_resp=NA, Route_direct=NA, Route_sexual=NA, Route_animal=NA, Route_vector=NA)

      tibble(
        Pathogen_Name = pathogen_name, R0_sampled = r0_sampled, SI_Clust_sampled = si_clust_sampled,
        CFR_sampled = cfr_sampled, Presymp_Proportion_sampled = presymp_prop_sampled, k_sampled = k_sampled,
        IP_sampled = ip_sampled, LP_sampled = lp_sampled, InfP_sampled = infp_sampled,
        Route_resp = as.integer(pathogen_routes$Route_resp), Route_direct = as.integer(pathogen_routes$Route_direct),
        Route_sexual = as.integer(pathogen_routes$Route_sexual), Route_animal = as.integer(pathogen_routes$Route_animal),
        Route_vector = as.integer(pathogen_routes$Route_vector)
      )
    })
    return(iteration_results)
  }

  print(paste("Starting parallel MCMC sampling for bootstrap replicate", b))
  all_mcmc_samples_df <- future_map_dfr(
    1:N_MCMC_ITERATIONS, 
    ~run_single_mcmc_iteration(.x, 
                               data_df = bootstrapped_params_df, 
                               presym_data_df = presym_dist_df,
                               transmission_data_df = transmission_df),
    .options = furrr_options(seed = MCMC_SEED + b),
    .progress = TRUE,
    .id = "MCMC_Iteration"
  )

  if (nrow(all_mcmc_samples_df) == 0) {
    warning(paste("MCMC sampling produced no results for bootstrap replicate", b, ". Skipping."))
    next # Skip to the next bootstrap replicate
  }
  
  # --- 5.3. Clustering Analysis for this Bootstrap Replicate ---
  
  # Create a dedicated output directory for this replicate's results
  output_dir <- file.path(main_output_dir, paste0("replicate_", b))
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }

  CHOSEN_K <- 4 # K for the per-iteration k-means
  
  # --- 8. Clustering Analysis (Adapted from original script) ---
  if (exists("all_mcmc_samples_df") && nrow(all_mcmc_samples_df) > 0) {
    
    # --- 8.1. Prepare Data for Clustering ---
    clustering_data_full <- all_mcmc_samples_df %>%
      select( MCMC_Iteration, Pathogen_Name, R0_sampled, SI_Clust_sampled, CFR_sampled, 
               Presymp_Proportion_sampled, k_sampled, IP_sampled, LP_sampled, InfP_sampled,
               Route_resp, Route_direct, Route_sexual, Route_animal, Route_vector )

    numerical_param_cols <- c("R0_sampled", "SI_Clust_sampled", "CFR_sampled", "Presymp_Proportion_sampled", "k_sampled", "IP_sampled", "LP_sampled", "InfP_sampled")
    
    # Impute missing values
    clustering_data_full <- clustering_data_full %>%
      group_by(Pathogen_Name) %>%
      mutate(across(all_of(numerical_param_cols), ~ifelse(is.na(.), median(., na.rm = TRUE), .))) %>%
      ungroup() %>%
      filter(if_all(all_of(numerical_param_cols), ~!is.na(.)))

    # Scale numerical features
    scaled_clustering_data_full <- clustering_data_full
    scaled_clustering_data_full[numerical_param_cols] <- scale(scaled_clustering_data_full[numerical_param_cols])

    # --- 8.2. K-Means Clustering for Each MCMC Iteration ---
    perform_kmeans_for_iteration <- function(iter_val, scaled_data, original_data) {
      current_iter_data_scaled <- scaled_data %>% filter(MCMC_Iteration == iter_val)
      current_iter_features_scaled <- current_iter_data_scaled %>% select(-MCMC_Iteration, -Pathogen_Name)
      if (nrow(current_iter_features_scaled) < CHOSEN_K) return(NULL)
      set.seed(CLUSTERING_SEED + as.integer(iter_val))
      kmeans_result <- tryCatch({ kmeans(current_iter_features_scaled, centers = CHOSEN_K, nstart = 25) }, error = function(e) NULL)
      if (is.null(kmeans_result)) return(NULL)
      assignments_df <- tibble( MCMC_Iteration = iter_val, Pathogen_Name = current_iter_data_scaled$Pathogen_Name, Cluster_Assigned = kmeans_result$cluster )
      return(list(assignments = assignments_df))
    }
    
    mcmc_iterations <- unique(scaled_clustering_data_full$MCMC_Iteration)
    clustering_results_list <- future_map( mcmc_iterations, ~perform_kmeans_for_iteration(.x, scaled_clustering_data_full, clustering_data_full), .options = furrr_options(seed = CLUSTERING_SEED), .progress = TRUE )
    
    successful_results <- purrr::compact(clustering_results_list)
    
    if (length(successful_results) > 0) {
      all_iteration_assignments <- purrr::map_dfr(successful_results, "assignments")
      
      # --- 9. Ensemble Clustering / Consensus Clustering (for this replicate) ---
      pathogen_names_replicate <- sort(unique(all_iteration_assignments$Pathogen_Name))
      num_pathogens_replicate <- length(pathogen_names_replicate)
      num_mcmc_replicate <- length(unique(all_iteration_assignments$MCMC_Iteration))

      # --- 9.1. Co-assignment Matrix ---
      coassignment_matrix <- matrix(0, nrow = num_pathogens_replicate, ncol = num_pathogens_replicate, dimnames = list(pathogen_names_replicate, pathogen_names_replicate))
      assignments_by_iter <- split(all_iteration_assignments[, c("Pathogen_Name", "Cluster_Assigned")], all_iteration_assignments$MCMC_Iteration)
      for (iter_idx in seq_along(assignments_by_iter)) {
        current_iter_assignments_df <- assignments_by_iter[[iter_idx]]
        for (i in 1:(num_pathogens_replicate - 1)) {
          for (j in (i + 1):num_pathogens_replicate) {
            pathogen1_name <- pathogen_names_replicate[i]; pathogen2_name <- pathogen_names_replicate[j]
            p1_assignment_row <- current_iter_assignments_df[current_iter_assignments_df$Pathogen_Name == pathogen1_name, ]
            p2_assignment_row <- current_iter_assignments_df[current_iter_assignments_df$Pathogen_Name == pathogen2_name, ]
            if (nrow(p1_assignment_row) > 0 && nrow(p2_assignment_row) > 0 && p1_assignment_row$Cluster_Assigned == p2_assignment_row$Cluster_Assigned) {
              coassignment_matrix[pathogen1_name, pathogen2_name] <- coassignment_matrix[pathogen1_name, pathogen2_name] + 1
              coassignment_matrix[pathogen2_name, pathogen1_name] <- coassignment_matrix[pathogen1_name, pathogen2_name]
            }
          }
        }
      }
      diag(coassignment_matrix) <- num_mcmc_replicate 
      
      # --- 9.2. Dissimilarity Matrix ---
      dissimilarity_matrix <- 1 - (coassignment_matrix / N_MCMC_ITERATIONS)
      diag(dissimilarity_matrix) <- 0

      # --- 9.3. Hierarchical Clustering ---
      pathogen_dist <- as.dist(dissimilarity_matrix)
      hierarchical_clust <- hclust(pathogen_dist, method = "average")
      
      # --- 9.4. Find Optimal K for this replicate using Silhouette Method ---
      k_range <- 2:10 # Range of K to test
      avg_silhouette_widths <- numeric(length(k_range))
      names(avg_silhouette_widths) <- k_range

      for (k in k_range) {
        cluster_assignments <- cutree(hierarchical_clust, k = k)
        silhouette_info <- silhouette(cluster_assignments, dmatrix = dissimilarity_matrix)
        if (is.null(silhouette_info) || !is.matrix(silhouette_info) || nrow(silhouette_info) == 0) {
            avg_width <- NA
        } else {
            avg_width <- summary(silhouette_info)$avg.width
        }
        avg_silhouette_widths[as.character(k)] <- avg_width
      }
      
      silhouette_results_df <- data.frame(K = k_range, Average_Silhouette_Width = avg_silhouette_widths)
      
      if(any(!is.na(silhouette_results_df$Average_Silhouette_Width))){
          optimal_k_value <- silhouette_results_df$K[which.max(silhouette_results_df$Average_Silhouette_Width)]
          print(paste("Optimal K for replicate", b, "based on silhouette analysis:", optimal_k_value))
      } else {
          optimal_k_value <- 4 # Fallback if analysis fails
          warning(paste("Could not determine optimal K for replicate", b, ". Falling back to K=4."))
      }

      K_CONSENSUS <- optimal_k_value
      
      # --- 9.5. Extract Consensus Cluster Assignments ---
      consensus_clusters <- cutree(hierarchical_clust, k = K_CONSENSUS)
      consensus_assignments_df <- data.frame(Pathogen_Name = names(consensus_clusters),
                                             Consensus_Cluster = consensus_clusters,
                                             Optimal_K_for_Replicate = K_CONSENSUS)

      # --- [End of embedded clustering analysis] ---
      print(paste("Finished clustering for bootstrap replicate", b))
    } else {
      warning(paste("Per-iteration K-means failed for all iterations in bootstrap replicate", b))
      consensus_assignments_df <- data.frame() # Empty dataframe
    }
  } else {
    warning(paste("No valid MCMC data for clustering in bootstrap replicate", b))
    consensus_assignments_df <- data.frame() # Empty dataframe
  }

  # --- 5.4. Store Results of the Replicate ---
  if (exists("consensus_assignments_df") && nrow(consensus_assignments_df) > 0) {
    # We now store the assignments which include the Optimal K for that replicate
    all_bootstrap_results[[b]] <- consensus_assignments_df 
  } else {
    warning(paste("Consensus clustering failed for bootstrap replicate", b))
  }
} # --- End of Bootstrap Loop ---


# --- 6. Analyze Results Across All Bootstrap Replicates ---
if (length(all_bootstrap_results) > 0) {
  print("--- Analyzing Stability Across Bootstrap Replicates ---")
  
  # Combine all results into a single dataframe
  final_bootstrap_summary_df <- bind_rows(all_bootstrap_results)
  
  # Save the raw bootstrap results
  write_csv(final_bootstrap_summary_df, file.path(main_output_dir, "bootstrap_all_assignments.csv"))
  
  # --- NEW: 6.0 Summarize Optimal K ---
  optimal_k_summary <- final_bootstrap_summary_df %>%
    group_by(Bootstrap_Replicate) %>%
    summarise(Optimal_K = first(Optimal_K_for_Replicate)) %>%
    ungroup()
  
  print("--- Summary of Optimal K Across Bootstrap Replicates ---")
  print(table(optimal_k_summary$Optimal_K))
  
  # Save the optimal K summary
  write_csv(optimal_k_summary, file.path(main_output_dir, "bootstrap_optimal_k_summary.csv"))

  # --- 6.1. Construct a final co-assignment matrix from bootstrap results ---
  pathogen_names <- sort(unique(final_bootstrap_summary_df$Pathogen_Name))
  num_pathogens <- length(pathogen_names)
  
  final_coassignment_matrix <- matrix(0, nrow = num_pathogens, ncol = num_pathogens,
                                      dimnames = list(pathogen_names, pathogen_names))

  assignments_by_replicate <- split(final_bootstrap_summary_df, final_bootstrap_summary_df$Bootstrap_Replicate)

  for (replicate_assignments in assignments_by_replicate) {
    for (i in 1:(num_pathogens - 1)) {
      for (j in (i + 1):num_pathogens) {
        p1_name <- pathogen_names[i]
        p2_name <- pathogen_names[j]
        
        p1_cluster <- replicate_assignments$Consensus_Cluster[replicate_assignments$Pathogen_Name == p1_name]
        p2_cluster <- replicate_assignments$Consensus_Cluster[replicate_assignments$Pathogen_Name == p2_name]
        
        if (length(p1_cluster) > 0 && length(p2_cluster) > 0 && p1_cluster == p2_cluster) {
          final_coassignment_matrix[p1_name, p2_name] <- final_coassignment_matrix[p1_name, p2_name] + 1
          final_coassignment_matrix[p2_name, p1_name] <- final_coassignment_matrix[p1_name, p2_name]
        }
      }
    }
  }
  
  # Normalize by the number of successful bootstrap replicates
  final_coassignment_matrix <- final_coassignment_matrix / length(all_bootstrap_results)
  diag(final_coassignment_matrix) <- 1
  
  print("Final Co-assignment Matrix (Proportion of times pathogens clustered together across bootstraps):")
  print(head(final_coassignment_matrix))
  write.csv(final_coassignment_matrix, file.path(main_output_dir, "bootstrap_coassignment_matrix.csv"))

  # --- 6.2. Final "most likely" clustering ---
  # Can use this final co-assignment matrix to perform one last hierarchical clustering
  # to get the most stable, final archetypes.
  final_dissimilarity <- 1 - final_coassignment_matrix
  final_hclust <- hclust(as.dist(final_dissimilarity), method = "average")
  
  # --- 6.3. Plot Final Dendrogram ---
  print("Plotting final consensus dendrogram...")

  # Determine the modal optimal K from the summary
  k_counts <- table(optimal_k_summary$Optimal_K)
  modal_k <- as.numeric(names(k_counts)[which.max(k_counts)])
  print(paste("Using modal optimal K (K=", modal_k, ") for final dendrogram and assignments.", sep=""))
  
  # Define the mapping for full pathogen names
  pathogen_full_name_map <- c(
    "COVID-19_WT" = "SARS-CoV-2 (WT)", "COVID-19_A" = "SARS-CoV-2 (Alpha)",
    "COVID-19_D" = "SARS-CoV-2 (Delta)", "COVID-19_O" = "SARS-CoV-2 (Omicron)", 
    "H1N1_18" = "A/H1N1", "H2N2" = "A/H2N2", "H3N2" = "A/H3N2",
    "H1N1_09" = "A/H1N1/09", "H5N1" = "A/H5N1", "Ebola" = "EBOV",
    "Marburg" = "MARV", "Mpox" = "MPV", "Lassa" = "LASV",
    "Nipah" = "NiV", "Zika" = "ZIKV", "SARS" = "SARS-CoV-1",
    "MERS" = "MERS-CoV", "CCHF" = "CCHFV"
  )

  dend <- as.dendrogram(final_hclust)

  # Apply full names to labels
  original_labels <- labels(dend)
  new_labels <- pathogen_full_name_map[original_labels]
  new_labels[is.na(new_labels)] <- original_labels[is.na(new_labels)]
  labels(dend) <- new_labels

  # Color the dendrogram based on the modal K
  K_FINAL_CONSENSUS <- modal_k
  if (!requireNamespace("RColorBrewer", quietly = TRUE)) { install.packages("RColorBrewer", quiet = TRUE) }
  cluster_colors <- RColorBrewer::brewer.pal(n = max(3, K_FINAL_CONSENSUS), name = "Set2")
  if (K_FINAL_CONSENSUS > length(cluster_colors)) {
      cluster_colors <- rep(cluster_colors, length.out = K_FINAL_CONSENSUS)
  }

  dend_colored <- dendextend::color_branches(dend, k = K_FINAL_CONSENSUS, col = cluster_colors[1:K_FINAL_CONSENSUS])
  dend_colored <- dendextend::set(dend_colored, "labels_cex", 0.8)
  dend_colored <- dendextend::set(dend_colored, "branches_lwd", 2)

  # Save the plot
  png_filename <- file.path(main_output_dir, paste0("final_consensus_dendrogram_k", K_FINAL_CONSENSUS, ".png"))
  png(png_filename, width=1200, height=800, units="px", res=100)
  par(mar = c(5, 4, 4, 12)) # Adjust right margin for labels
  plot(dend_colored, horiz = TRUE, 
       main = "Final Consensus Clustering from Bootstrap Analysis", 
       xlab = "Dissimilarity (1 - Proportion Co-assigned across bootstraps)")
  dev.off()
  print(paste("Final consensus dendrogram saved to", png_filename))
  
  # Cut tree at the modal K to get final assignments
  final_assignments <- cutree(final_hclust, k = modal_k)
  final_assignments_df <- data.frame(Pathogen_Name = names(final_assignments),
                                     Final_Consensus_Cluster = final_assignments)
  
  print("Final stable cluster assignments based on bootstrap analysis:")
  print(final_assignments_df)
  write.csv(final_assignments_df, file.path(main_output_dir, "bootstrap_final_assignments.csv"))
  
} else {
  print("Bootstrap analysis did not produce any results.")
}

print("Script finished.")
