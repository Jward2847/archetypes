# MCMC K-Means Clustering Analysis - TEST SCRIPT FOR NIEL'S SUGGESTION
# This script explores an alternative methodological approach suggested by Niel Hens.
# The standard approach uses a Monte Carlo loop where each iteration generates a new
# parameter set by bootstrapping the available literature (studies).
# This script reverses that order:
# 1. OUTER LOOP: Bootstrap the entire dataset of studies N times (N_BOOTSTRAP_ITERATIONS).
# 2. INNER LOOP: For each bootstrapped dataset, run the full Monte Carlo + Consensus
#    Clustering analysis.
# This allows us to assess the stability of the final consensus clusters to the
# specific set of studies available in the literature.

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
N_BOOTSTRAP_ITERATIONS <- 10 # Number of bootstrap iterations of the entire dataset
N_MCMC_ITERATIONS <- 500 # Number of MCMC iterations for each bootstrap
N_PRESYMP_SAMPLES <- 5000 # Number of samples for presymptomatic proportion estimation
MCMC_SEED <- 123 # Seed for reproducibility of the MCMC parameter sampling
CLUSTERING_SEED <- 456 # Seed for reproducibility of the K-means clustering
set.seed(MCMC_SEED) # Set seed for the MCMC part


# --- 3. Load Data ---
# Load the new long-format parameter data and transmission routes
params_long_df_orig <- read_csv("Clustering/mcmc/Kmeans/data/pathogen_params.csv")
transmission_df_orig <- read_csv("Clustering/mcmc/Kmeans/data/transmission_route.csv")
presym_dist_df_orig <- read_csv("Clustering/mcmc/Kmeans/data/pathogen_presym.csv")

# Filter for the specific pathogens of interest from the image
pathogens_of_interest <- c(
  "COVID-19_WT", "COVID-19_D", "COVID-19_O", "Ebola", "Marburg", "Mpox",
  "H1N1_18", "H1N1_09", "SARS", "MERS"
)
params_long_df_orig <- params_long_df_orig %>% filter(Pathogen_Name %in% pathogens_of_interest)
transmission_df_orig <- transmission_df_orig %>% filter(Pathogen_Name %in% pathogens_of_interest)
presym_dist_df_orig <- presym_dist_df_orig %>% filter(Pathogen_Name %in% pathogens_of_interest)


# --- 4. Helper Functions ---
# NOTE: Helper functions (get_beta_params_from_ci, get_gamma_params_from_ci, etc.) are
# assumed to be the same as in the original script and are included here for completeness.

# Function to derive parameters for rbeta by matching quantiles of the 95% CI
get_beta_params_from_ci <- function(lower_ci, upper_ci, mean_val) {
    min_param_val <- 1e-6
    objective_beta <- function(params, lower_ci, upper_ci) {
        shape1 <- exp(params[1])
        shape2 <- exp(params[2])
        if (is.na(shape1) || is.na(shape2) || !is.finite(shape1) || !is.finite(shape2)) return(1e10)
        q_lower <- qbeta(0.025, shape1, shape2)
        q_upper <- qbeta(0.975, shape1, shape2)
        err <- (q_lower - lower_ci)^2 + (q_upper - upper_ci)^2
        return(err)
    }
    if (!is.na(mean_val) && mean_val > 0 && mean_val < 1) {
        sd_guess <- (upper_ci - lower_ci) / (2 * qnorm(0.975))
        var_guess <- sd_guess^2
        alpha_guess <- ((1 - mean_val) / var_guess - 1 / mean_val) * mean_val^2
        beta_guess <- alpha_guess * (1 / mean_val - 1)
        if (is.na(alpha_guess) || alpha_guess <= 0) alpha_guess <- 1
        if (is.na(beta_guess) || beta_guess <= 0) beta_guess <- 1
        start_params <- log(c(alpha_guess, beta_guess))
    } else {
        start_params <- log(c(1, 1))
    }
    opt_result <- tryCatch({
        optim(start_params, objective_beta, lower_ci = lower_ci, upper_ci = upper_ci, control = list(maxit = 1000, reltol = 1e-8))
    }, error = function(e) NULL)
    if (!is.null(opt_result) && opt_result$convergence == 0) {
        return(list(shape1 = exp(opt_result$par[1]), shape2 = exp(opt_result$par[2])))
    } else {
        return(get_beta_params_from_mean_ci_fallback(mean_val, lower_ci, upper_ci))
    }
}

# Function to derive parameters for rgamma by matching quantiles of the 95% CI
get_gamma_params_from_ci <- function(lower_ci, upper_ci, mean_val) {
    min_param_val <- 1e-6
    objective_gamma <- function(params, lower_ci, upper_ci) {
        shape <- exp(params[1])
        rate <- exp(params[2])
        if (is.na(shape) || is.na(rate) || !is.finite(shape) || !is.finite(rate)) return(1e10)
        q_lower <- qgamma(0.025, shape = shape, rate = rate)
        q_upper <- qgamma(0.975, shape = shape, rate = rate)
        err <- (q_lower - lower_ci)^2 + (q_upper - upper_ci)^2
        return(err)
    }
    if(!is.na(mean_val) && mean_val > 0){
        sd_guess <- (upper_ci - lower_ci) / (2 * qnorm(0.975))
        if(is.na(sd_guess) || sd_guess <= 0) sd_guess <- mean_val
        shape_guess <- (mean_val / sd_guess)^2
        rate_guess <- mean_val / sd_guess^2
        if (is.na(shape_guess) || shape_guess <= 0) shape_guess <- 1
        if (is.na(rate_guess) || rate_guess <= 0) rate_guess <- 1
        start_params <- log(c(shape_guess, rate_guess))
    } else {
        start_params <- log(c(1,1))
    }
    opt_result <- tryCatch({
        optim(start_params, objective_gamma, lower_ci = lower_ci, upper_ci = upper_ci, control = list(maxit = 1000, reltol = 1e-8))
    }, error = function(e) NULL)
    if (!is.null(opt_result) && opt_result$convergence == 0) {
        return(list(shape = exp(opt_result$par[1]), rate = exp(opt_result$par[2])))
    } else {
        return(get_gamma_params_from_mean_ci_fallback(mean_val, lower_ci, upper_ci))
    }
}

# --- FALLBACK Helper Functions (Original Methods) ---
get_beta_params_from_mean_ci_fallback <- function(mean_val, lower_ci, upper_ci, n_eff_guess = 1000) {
    alpha <- mean_val * n_eff_guess
    beta <- (1 - mean_val) * n_eff_guess
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
        } else {
            return(list(shape = min_param_val, rate = min_param_val))
        }
    }
    sd_from_ci <- (upper_ci - lower_ci) / (2 * qnorm(0.975))
    if (is.na(sd_from_ci) || !is.finite(sd_from_ci) || sd_from_ci <= 0) {
       sd_from_ci <- mean_val * 0.5
       if (sd_from_ci <= 0) sd_from_ci <- min_param_val
    }
    shape <- (mean_val / sd_from_ci)^2
    rate <- mean_val / sd_from_ci^2
    if (is.na(shape) || !is.finite(shape) || shape <= 0) shape <- min_param_val
    if (is.na(rate) || !is.finite(rate) || rate <= 0) rate <- min_param_val
    return(list(shape = shape, rate = rate))
}

# The main sampling function, now simplified. It NO LONGER BOOTSTRAPS internally.
# It just samples from the provided studies.
sample_parameter_from_studies <- function(param_name, pathogen_name, data_df) {
  studies <- data_df %>% filter(Pathogen_Name == pathogen_name, Parameter == param_name)
  if (nrow(studies) == 0) return(NA)

  # Generate one sample from each available study
  samples <- apply(studies, 1, function(study_row) {
    mean_val <- as.numeric(study_row["ReportedValue"])
    lower_ci <- as.numeric(study_row["LowerBound"])
    upper_ci <- as.numeric(study_row["UpperBound"])
    if (is.na(lower_ci) || is.na(upper_ci) || lower_ci >= upper_ci) {
        return(mean_val)
    }
    sampled_val_internal <- NA
    if (param_name == "CFR") {
        if (!is.na(mean_val) && mean_val > 0 && mean_val < 1) {
            beta_params <- get_beta_params_from_ci(lower_ci, upper_ci, mean_val)
            if (!is.na(beta_params$shape1)) {
                sampled_val_internal <- rbeta(1, shape1 = beta_params$shape1, shape2 = beta_params$shape2)
            }
        }
    } else if (param_name %in% c("R0", "SI", "k", "IP", "LP", "InfP")) {
        gamma_params <- get_gamma_params_from_ci(lower_ci, upper_ci, mean_val)
        if (!is.na(gamma_params$shape)) {
            sampled_val_internal <- rgamma(1, shape = gamma_params$shape, rate = gamma_params$rate)
        }
    } else {
        sd_from_ci <- (upper_ci - lower_ci) / (2 * qnorm(0.975))
        if (!is.na(mean_val) && !is.na(sd_from_ci) && sd_from_ci > 0) {
            sampled_val_internal <- rnorm(1, mean = mean_val, sd = sd_from_ci)
        }
    }
    if (is.na(sampled_val_internal) || !is.finite(sampled_val_internal)) {
        return(mean_val)
    }
    return(sampled_val_internal)
  })
  final_value <- mean(samples, na.rm = TRUE)
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
    return(NA)
  }
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
    return(NA)
  }
  random_quantiles <- runif(n_samples)
  si_samples <- q_si(random_quantiles)
  ip_samples <- q_ip(random_quantiles)
  if(si_dist_info$family == "Normal") si_samples[si_samples < 0] <- 0
  if(ip_dist_info$family == "Normal") ip_samples[ip_samples < 0] <- 0
  mean(si_samples < ip_samples, na.rm = TRUE)
}

# --- 5. Main Analysis Function (to be called for each bootstrap) ---
# This function encapsulates the entire analysis process: MCMC sampling and clustering.
run_full_analysis_on_bootstrap <- function(bootstrap_iter_num, params_df_boot, presym_df_boot, transmission_df) {

  print(paste("--- Starting Analysis for Bootstrap Iteration:", bootstrap_iter_num, "---"))

  # --- 5.1. MCMC Sampling ---
  run_single_mcmc_iteration <- function(iter_num, .progress = FALSE) {
    unique_pathogens <- unique(params_df_boot$Pathogen_Name)
    flu_group <- c("H1N1_09", "H1N1_18", "H2N2", "H3N2")

    iteration_results <- purrr::map_dfr(unique_pathogens, function(pathogen_name) {
      r0_sampled <- sample_parameter_from_studies("R0", pathogen_name, params_df_boot)
      si_clust_sampled <- sample_parameter_from_studies("SI", pathogen_name, params_df_boot)
      cfr_sampled <- sample_parameter_from_studies("CFR", pathogen_name, params_df_boot)
      k_sampled <- sample_parameter_from_studies("k", pathogen_name, params_df_boot)
      ip_sampled <- sample_parameter_from_studies("IP", pathogen_name, params_df_boot)
      lp_sampled <- sample_parameter_from_studies("LP", pathogen_name, params_df_boot)
      infp_sampled <- sample_parameter_from_studies("InfP", pathogen_name, params_df_boot)
      presymp_prop_sampled <- NA
      if (pathogen_name %in% flu_group) {
        si_studies <- presym_df_boot %>% filter(Pathogen_Name %in% flu_group, Parameter == "SI")
        ip_studies <- presym_df_boot %>% filter(Pathogen_Name %in% flu_group, Parameter == "IP")
      } else {
        si_studies <- presym_df_boot %>% filter(Pathogen_Name == pathogen_name, Parameter == "SI")
        ip_studies <- presym_df_boot %>% filter(Pathogen_Name == pathogen_name, Parameter == "IP")
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
      pathogen_routes <- transmission_df %>% filter(Pathogen_Name == pathogen_name)
      if(nrow(pathogen_routes) == 0){
          pathogen_routes <- data.frame(Route_resp=NA, Route_direct=NA, Route_sexual=NA, Route_animal=NA, Route_vector=NA)
      }
      tibble(
        Pathogen_Name = pathogen_name,
        R0_sampled = ifelse(is.null(r0_sampled), NA_real_, r0_sampled),
        SI_Clust_sampled = ifelse(is.null(si_clust_sampled), NA_real_, si_clust_sampled),
        CFR_sampled = ifelse(is.null(cfr_sampled), NA_real_, cfr_sampled),
        Presymp_Proportion_sampled = ifelse(is.null(presymp_prop_sampled), NA_real_, presymp_prop_sampled),
        k_sampled = ifelse(is.null(k_sampled), NA_real_, k_sampled),
        IP_sampled = ifelse(is.null(ip_sampled), NA_real_, ip_sampled),
        LP_sampled = ifelse(is.null(lp_sampled), NA_real_, lp_sampled),
        InfP_sampled = ifelse(is.null(infp_sampled), NA_real_, infp_sampled),
        Route_resp = as.integer(pathogen_routes$Route_resp),
        Route_direct = as.integer(pathogen_routes$Route_direct),
        Route_sexual = as.integer(pathogen_routes$Route_sexual),
        Route_animal = as.integer(pathogen_routes$Route_animal),
        Route_vector = as.integer(pathogen_routes$Route_vector)
      )
    })
    return(iteration_results)
  }

  all_mcmc_samples_df <- future_map_dfr(
    1:N_MCMC_ITERATIONS,
    run_single_mcmc_iteration,
    .options = furrr_options(seed = MCMC_SEED + bootstrap_iter_num),
    .progress = TRUE,
    .id = "MCMC_Iteration"
  )

  # --- 5.2. Consensus Clustering ---
  # This part is simplified from the main script. It will run K-means and then consensus clustering.
  if (nrow(all_mcmc_samples_df) < 1) return(NULL)
  
  numerical_param_cols <- c("R0_sampled", "SI_Clust_sampled", "CFR_sampled", "Presymp_Proportion_sampled", "k_sampled", "IP_sampled", "LP_sampled", "InfP_sampled")
  clustering_data_full <- all_mcmc_samples_df %>%
    group_by(Pathogen_Name) %>%
    mutate(across(all_of(numerical_param_cols), ~ifelse(is.na(.), median(., na.rm = TRUE), .))) %>%
    ungroup() %>%
    filter(if_all(all_of(numerical_param_cols), ~!is.na(.)))
    
  if(nrow(clustering_data_full) == 0) return(NULL)

  scaled_clustering_data_full <- clustering_data_full
  scaled_clustering_data_full[numerical_param_cols] <- scale(clustering_data_full[numerical_param_cols])

  CHOSEN_K <- 4 # For this test, we'll keep K fixed.

  perform_kmeans_for_iteration <- function(iter_val, scaled_data) {
    current_iter_data_scaled <- scaled_data %>% filter(MCMC_Iteration == iter_val)
    current_iter_features_scaled <- current_iter_data_scaled %>% select(-MCMC_Iteration, -Pathogen_Name)
    if (nrow(current_iter_features_scaled) < CHOSEN_K) return(NULL)
    set.seed(CLUSTERING_SEED + as.integer(iter_val))
    kmeans_result <- tryCatch(kmeans(current_iter_features_scaled, centers = CHOSEN_K, nstart = 25), error = function(e) NULL)
    if (is.null(kmeans_result)) return(NULL)
    tibble(
      MCMC_Iteration = iter_val,
      Pathogen_Name = current_iter_data_scaled$Pathogen_Name,
      Cluster_Assigned = kmeans_result$cluster
    )
  }

  mcmc_iterations <- unique(scaled_clustering_data_full$MCMC_Iteration)
  all_iteration_assignments <- future_map_dfr(
    mcmc_iterations,
    ~perform_kmeans_for_iteration(.x, scaled_clustering_data_full)
  )

  if(nrow(all_iteration_assignments) == 0) return(NULL)

  pathogen_names <- sort(unique(all_iteration_assignments$Pathogen_Name))
  num_pathogens <- length(pathogen_names)
  coassignment_matrix <- matrix(0, nrow = num_pathogens, ncol = num_pathogens, dimnames = list(pathogen_names, pathogen_names))
  assignments_by_iter <- split(all_iteration_assignments, all_iteration_assignments$MCMC_Iteration)

  for (iter_assignments in assignments_by_iter) {
    for (i in 1:(num_pathogens - 1)) {
      for (j in (i + 1):num_pathogens) {
        p1_name <- pathogen_names[i]
        p2_name <- pathogen_names[j]
        p1_assign <- iter_assignments$Cluster_Assigned[iter_assignments$Pathogen_Name == p1_name]
        p2_assign <- iter_assignments$Cluster_Assigned[iter_assignments$Pathogen_Name == p2_name]
        if (length(p1_assign) > 0 && length(p2_assign) > 0 && p1_assign == p2_assign) {
          coassignment_matrix[p1_name, p2_name] <- coassignment_matrix[p1_name, p2_name] + 1
          coassignment_matrix[p2_name, p1_name] <- coassignment_matrix[p1_name, p2_name]
        }
      }
    }
  }
  diag(coassignment_matrix) <- length(assignments_by_iter)
  
  dissimilarity_matrix <- 1 - (coassignment_matrix / length(assignments_by_iter))
  diag(dissimilarity_matrix) <- 0

  hierarchical_clust <- hclust(as.dist(dissimilarity_matrix), method = "average")

  # Return the final consensus assignments for this bootstrap iteration
  k_range <- 2:min(9, num_pathogens - 1)
  avg_silhouette_widths <- sapply(k_range, function(k) {
    assignments <- cutree(hierarchical_clust, k = k)
    summary(silhouette(assignments, dmatrix = dissimilarity_matrix))$avg.width
  })
  optimal_k <- k_range[which.max(avg_silhouette_widths)]
  
  consensus_clusters <- cutree(hierarchical_clust, k = optimal_k)
  
  result <- tibble(
      Bootstrap_Iteration = bootstrap_iter_num,
      Pathogen_Name = names(consensus_clusters),
      Consensus_Cluster = consensus_clusters,
      Optimal_K = optimal_k
  )

  return(result)
}


# --- 6. Main Bootstrap Loop ---
all_bootstrap_results <- list()

for (boot_i in 1:N_BOOTSTRAP_ITERATIONS) {
  set.seed(MCMC_SEED + boot_i)
  
  # --- Bootstrap the studies ---
  # We group by pathogen and parameter, then resample the studies within each group
  params_df_bootstrapped <- params_long_df_orig %>%
    group_by(Pathogen_Name, Parameter) %>%
    sample_n(size = n(), replace = TRUE) %>%
    ungroup()

  presym_df_bootstrapped <- presym_dist_df_orig %>%
    group_by(Pathogen_Name, Parameter) %>%
    sample_n(size = n(), replace = TRUE) %>%
    ungroup()

  # Run the full analysis on this bootstrapped dataset
  # Using tryCatch to ensure the loop continues even if one bootstrap fails
  bootstrap_result <- tryCatch({
    run_full_analysis_on_bootstrap(
      bootstrap_iter_num = boot_i,
      params_df_boot = params_df_bootstrapped,
      presym_df_boot = presym_df_bootstrapped,
      transmission_df = transmission_df_orig # Transmission routes are fixed
    )
  }, error = function(e) {
    warning(paste("Analysis failed for bootstrap iteration", boot_i, ":", e$message))
    return(NULL)
  })
  
  if (!is.null(bootstrap_result)) {
    all_bootstrap_results[[boot_i]] <- bootstrap_result
  }
}

# --- 7. Analyze Bootstrap Results ---
if (length(all_bootstrap_results) > 0) {
  final_results_df <- bind_rows(all_bootstrap_results)
  
  # Save the raw results
  write_csv(final_results_df, "Clustering/mcmc/Kmeans/test_outputs/bootstrap_test_all_assignments.csv")
  print("Saved all raw bootstrap assignments to Clustering/mcmc/Kmeans/test_outputs/bootstrap_test_all_assignments.csv")

  # --- Analyze stability ---
  # 1. How often was a certain K chosen as optimal?
  optimal_k_summary <- final_results_df %>%
    distinct(Bootstrap_Iteration, Optimal_K) %>%
    count(Optimal_K) %>%
    mutate(Proportion = n / sum(n))
  
  print("Summary of Optimal K chosen across bootstrap iterations:")
  print(optimal_k_summary)
  write_csv(optimal_k_summary, "Clustering/mcmc/Kmeans/test_outputs/bootstrap_test_optimal_k_summary.csv")

  # 2. Co-assignment matrix across bootstraps
  # This tells us the proportion of bootstrap iterations where two pathogens ended up in the same cluster
  pathogen_names <- sort(unique(final_results_df$Pathogen_Name))
  num_pathogens <- length(pathogen_names)
  
  bootstrap_coassignment <- matrix(0, nrow = num_pathogens, ncol = num_pathogens, dimnames = list(pathogen_names, pathogen_names))
  
  assignments_by_bootstrap <- split(final_results_df, final_results_df$Bootstrap_Iteration)
  
  for (boot_assignments in assignments_by_bootstrap) {
    for (i in 1:(num_pathogens - 1)) {
      for (j in (i + 1):num_pathogens) {
        p1_name <- pathogen_names[i]
        p2_name <- pathogen_names[j]
        p1_cluster <- boot_assignments$Consensus_Cluster[boot_assignments$Pathogen_Name == p1_name]
        p2_cluster <- boot_assignments$Consensus_Cluster[boot_assignments$Pathogen_Name == p2_name]
        
        if (length(p1_cluster) > 0 && length(p2_cluster) > 0 && p1_cluster == p2_cluster) {
          bootstrap_coassignment[p1_name, p2_name] <- bootstrap_coassignment[p1_name, p2_name] + 1
          bootstrap_coassignment[p2_name, p1_name] <- bootstrap_coassignment[p1_name, p2_name]
        }
      }
    }
  }
  
  bootstrap_coassignment_prob <- bootstrap_coassignment / length(assignments_by_bootstrap)
  diag(bootstrap_coassignment_prob) <- 1
  
  print("Co-assignment probability matrix across all bootstraps:")
  print(bootstrap_coassignment_prob)
  write.csv(bootstrap_coassignment_prob, "Clustering/mcmc/Kmeans/test_outputs/bootstrap_test_coassignment_matrix.csv")

  # 3. Final consensus clustering based on the bootstrap co-assignment matrix
  final_dissimilarity <- 1 - bootstrap_coassignment_prob
  diag(final_dissimilarity) <- 0
  final_hclust <- hclust(as.dist(final_dissimilarity), method = "average")
  
  # Use the most frequent optimal K from the summary
  if(nrow(optimal_k_summary) > 0) {
      final_optimal_k <- optimal_k_summary$Optimal_K[which.max(optimal_k_summary$n)]
  } else {
      final_optimal_k <- 4 # Fallback
  }

  final_assignments <- cutree(final_hclust, k = final_optimal_k)
  final_assignments_df <- data.frame(Pathogen_Name = names(final_assignments), Final_Consensus_Cluster = final_assignments)
  
  print("Final consensus assignments from bootstrap analysis:")
  print(final_assignments_df)
  write.csv(final_assignments_df, "Clustering/mcmc/Kmeans/test_outputs/bootstrap_test_final_assignments.csv")
  
} else {
  print("No successful bootstrap iterations to analyze.")
}

print("Test script finished.")
