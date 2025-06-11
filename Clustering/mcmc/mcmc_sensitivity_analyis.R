# MCMC K-Means Clustering Sensitivity Analysis
#
# This script conducts two separate sensitivity analyses for the MCMC K-means clustering.
#
# 1. RVFV Sensitivity Analysis: Loads data from 'pathogen_params_RVFV.csv', which includes
#    RVFV. It then performs clustering excluding serial interval and presymptomatic transmission
#    proportion from the clustering variables.
#
# 2. HIV Sensitivity Analysis: Loads data from 'pathogen_params_HIV.csv', which includes
#    HIV. It then performs clustering excluding serial interval and case fatality rate (CFR)
#    from the clustering variables.
#
# Each analysis is self-contained, running its own MCMC sampling and saving results
# to uniquely named files.

# --- 1. Load Libraries ---
library(dplyr)
library(readr)
library(stats)
library(dendextend)
library(RColorBrewer)
library(ggplot2)
library(factoextra)
library(ggrepel)

# --- 2. Configuration ---
N_MCMC_ITERATIONS <- 5000 # Number of MCMC iterations
N_PRESYMP_SAMPLES <- 5000 # Number of samples for presymptomatic proportion estimation
set.seed(123) # For reproducibility

# --- 3. Helper Functions (Copied from main analysis script) ---

# Function to derive parameters for rbeta from mean and 95% CI
get_beta_params_from_mean_ci <- function(mean_val, lower_ci, upper_ci, n_eff_guess = 1000) {
    alpha <- mean_val * n_eff_guess
    beta <- (1 - mean_val) * n_eff_guess
    min_param_val <- 1e-6
    if (is.na(alpha) || !is.finite(alpha) || alpha <= 0) alpha <- min_param_val 
    if (is.na(beta) || !is.finite(beta) || beta <= 0) beta <- min_param_val
    if (mean_val <= 0 || mean_val >= 1 || lower_ci < 0 || upper_ci > 1 || lower_ci >= upper_ci) {
        warning(paste("Input values for get_beta_params_from_mean_ci may be problematic. Mean:", mean_val, "CI: [", lower_ci, ",", upper_ci, "]"))
    }
    return(list(shape1 = alpha, shape2 = beta))
}

# Function to derive parameters for rgamma from mean and 95% CI
get_gamma_params_from_mean_ci <- function(mean_val, lower_ci, upper_ci) {
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

# Function to sample a value for a _Clust_ parameter (R0, SI, CFR)
sample_clust_parameter <- function(pathogen_row, param_prefix) {
  val_col <- paste0(param_prefix, "_Clust_ReportedValue"); uncert_col <- paste0(param_prefix, "_Clust_UncertaintyType"); low_col <- paste0(param_prefix, "_Clust_LowerBound"); high_col <- paste0(param_prefix, "_Clust_UpperBound"); dist_col <- paste0(param_prefix, "_Clust_SamplingDist"); cv_col <- paste0(param_prefix, "_Clust_PointEstimateCV")
  val <- as.numeric(pathogen_row[[val_col]]); uncert_type <- pathogen_row[[uncert_col]]; sampling_dist <- pathogen_row[[dist_col]]; lower_bound <- as.numeric(pathogen_row[[low_col]]); upper_bound <- as.numeric(pathogen_row[[high_col]])
  if (is.na(uncert_type) || uncert_type == "") { if (!is.na(val)) return(val) else return(NA) } # Return NA if no info
  sampled_val <- NA
  if (uncert_type == "PointEstimate") {
    cv <- as.numeric(pathogen_row[[cv_col]]); if (is.na(cv) || cv == 0) return(val); sd_val <- val * cv
    if (sd_val <= 0) { return(val) }
    if (sampling_dist == "rnorm" || is.na(sampling_dist) || sampling_dist == "") {
        sampled_val <- rnorm(1, mean = val, sd = sd_val)
        if (param_prefix %in% c("CFR", "Presymp")) sampled_val <- pmax(0, pmin(1, sampled_val)) else if (param_prefix %in% c("R0", "SI")) sampled_val <- pmax(1e-5, sampled_val)
    } else if (sampling_dist == "rbeta" && param_prefix %in% c("CFR", "Presymp")) {
        m <- val; s <- sd_val
        if (m > 0 && m < 1 && s > 0 && (m * (1 - m) > s^2) && ((m * (1 - m) / s^2) - 1) > 0) {
            alpha <- m * (((m * (1 - m)) / s^2) - 1); beta <- (1 - m) * (((m * (1 - m)) / s^2) - 1)
            if (alpha > 0 && beta > 0) { sampled_val <- rbeta(1, shape1 = alpha, shape2 = beta) } else { sampled_val <- val }
        } else { sampled_val <- val }
    } else { sampled_val <- val }
  } else if (uncert_type %in% c("95CI", "IQR", "Range", "range")) {
    if (is.na(lower_bound) || is.na(upper_bound)) { return(val) }
    if (lower_bound > upper_bound) { tmp <- lower_bound; lower_bound <- upper_bound; upper_bound <- tmp; }
    if (sampling_dist == "rlnorm") {
      current_lower_bound <- lower_bound; if (current_lower_bound <= 0) { current_lower_bound <- 1e-5; }
      if (upper_bound <= current_lower_bound) { return(val); }
      log_lower <- log(current_lower_bound); log_upper <- log(upper_bound); sdlog <- (log_upper - log_lower) / (2 * qnorm(0.975)); meanlog <- (log_upper + log_lower) / 2
      if(is.na(sdlog) || !is.finite(sdlog) || sdlog <= 0) { return(val); }
      sampled_val <- rlnorm(1, meanlog = meanlog, sdlog = sdlog)
    } else if (sampling_dist == "rgamma") {
      if(is.na(val)) { if(!is.na(lower_bound) && !is.na(upper_bound)) val <- (lower_bound + upper_bound) / 2 else return(NA) }
      gamma_params <- get_gamma_params_from_mean_ci(val, lower_bound, upper_bound)
      if(any(is.na(gamma_params)) || any(unlist(gamma_params) <= 0) || any(!is.finite(unlist(gamma_params)))){ return(val) }
      sampled_val <- rgamma(1, shape = gamma_params$shape, rate = gamma_params$rate); if (sampled_val <= 0) sampled_val <- 1e-5
    } else if (sampling_dist == "rbeta") { 
      if(is.na(val)) { if(!is.na(lower_bound) && !is.na(upper_bound)) val <- (lower_bound + upper_bound) / 2 else return(NA) }
      if(val < 0 || val > 1) { val <- pmax(1e-5, pmin(1 - 1e-5, val)) }
      clamped_lower <- pmax(0, lower_bound); clamped_upper <- pmin(1, upper_bound); if (clamped_lower >= clamped_upper) { return(val) }
      beta_params <- get_beta_params_from_mean_ci(val, clamped_lower, clamped_upper)
      if(any(is.na(beta_params)) || any(unlist(beta_params) <= 0) || any(!is.finite(unlist(beta_params)))){ return(val) }
      sampled_val <- rbeta(1, shape1 = beta_params$shape1, shape2 = beta_params$shape2)
    } else if (sampling_dist == "rnorm") {
      sd_val <- (upper_bound - lower_bound) / (2 * qnorm(0.975)); if(is.na(sd_val) || !is.finite(sd_val) || sd_val <= 0){ return(val); }
      sampled_val <- rnorm(1, mean = val, sd = sd_val); if (param_prefix %in% c("CFR", "Presymp")) sampled_val <- pmax(0, pmin(1, sampled_val)) else if (param_prefix %in% c("R0", "SI")) sampled_val <- pmax(1e-5, sampled_val)
    } else if (sampling_dist == "runif") {
      if (lower_bound >= upper_bound) { return(lower_bound) }; sampled_val <- runif(1, min = lower_bound, max = upper_bound)
    } else { sampled_val <- val }
  } else { sampled_val <- val }
  if (is.na(sampled_val)) { return(val) }; return(sampled_val)
}

get_full_dist_params <- function(pathogen_row, type_prefix) { 
  family <- pathogen_row[[paste0(type_prefix, "_FullDist_Family")]]; raw_params <- list()
  for (i in 1:2) { 
    param_name_col <- paste0(type_prefix, "_FullDist_Param", i, "_Name"); param_val_col <- paste0(type_prefix, "_FullDist_Param", i, "_Value"); param_uncert_col <- paste0(type_prefix, "_FullDist_Param", i, "_UncertaintySource")
    if (!is.na(pathogen_row[[param_name_col]]) && pathogen_row[[param_name_col]] != "") {
      param_name_from_csv <- pathogen_row[[param_name_col]]; param_val <- as.numeric(pathogen_row[[param_val_col]]); param_uncert <- pathogen_row[[param_uncert_col]]
      if (is.na(param_uncert) || param_uncert %in% c("Fixed", "")) { raw_params[[param_name_from_csv]] <- param_val
      } else if (param_uncert == "95CI") {
        param_lower <- as.numeric(pathogen_row[[paste0(type_prefix, "_FullDist_Param", i, "_LowerBound")]]); param_upper <- as.numeric(pathogen_row[[paste0(type_prefix, "_FullDist_Param", i, "_UpperBound")]]); param_sample_dist <- pathogen_row[[paste0(type_prefix, "_FullDist_Param", i, "_SamplingDist")]]
        if (!any(is.na(c(param_lower, param_upper, param_sample_dist))) && param_sample_dist != "") {
          if (param_sample_dist == "rnorm") {
            sd_p <- (param_upper - param_lower) / (2 * qnorm(0.975))
            if(is.na(sd_p) || sd_p <= 0){ sampled_param_val <- param_val } else { sampled_param_val <- rnorm(1, mean = param_val, sd = sd_p) }
            if (param_name_from_csv %in% c("SD", "sd", "shape", "rate", "scale", "sdlog") && sampled_param_val <= 0) { sampled_param_val <- param_val }
            raw_params[[param_name_from_csv]] <- sampled_param_val
          } else { raw_params[[param_name_from_csv]] <- param_val }
        } else { raw_params[[param_name_from_csv]] <- param_val }
      } else { raw_params[[param_name_from_csv]] <- param_val }
    }
  }
  final_params <- raw_params
  if (!is.na(family)) {
    if (family == "Gamma") {
      if (all(c("mean", "SD") %in% names(raw_params)) && !("shape" %in% names(raw_params))) {
        mean_val <- raw_params$mean; sd_val <- raw_params$SD
        if (all(!is.na(c(mean_val, sd_val))) && sd_val > 0 && mean_val > 0) {
          final_params$shape <- (mean_val / sd_val)^2; final_params$rate <- mean_val / (sd_val^2); final_params$mean <- NULL; final_params$SD <- NULL
        } else { return(list(family = family, params = list())) }
      }
      for(p in c("shape", "rate", "scale")) { if (!is.null(final_params[[p]]) && (!is.finite(final_params[[p]]) || final_params[[p]] <= 1e-9)) { final_params[[p]] <- 1e-9 } }
    } else if (family == "Lognormal") {
      if (all(c("mean", "SD") %in% names(raw_params)) && !("meanlog" %in% names(raw_params))) {
        mean_val <- raw_params$mean; sd_val <- raw_params$SD
        if(all(!is.na(c(mean_val, sd_val))) && mean_val > 0 && sd_val > 0){
          final_params$sdlog <- sqrt(log(sd_val^2 / mean_val^2 + 1)); final_params$meanlog <- log(mean_val) - 0.5 * final_params$sdlog^2; final_params$mean <- NULL; final_params$SD <- NULL
        } else { return(list(family = family, params = list())) }
      }
      if (!is.null(final_params$sdlog) && (!is.finite(final_params$sdlog) || final_params$sdlog <= 1e-9)) { final_params$sdlog <- 1e-9 }
    } else if (family == "Normal") {
      if (!is.null(final_params$SD) && (!is.finite(final_params$SD) || final_params$SD <= 1e-9)) { final_params$SD <- 1e-9 }
    } else if (family == "Weibull") {
      for(p in c("shape", "scale")) { if (!is.null(final_params[[p]]) && (!is.finite(final_params[[p]]) || final_params[[p]] <= 1e-9)) { final_params[[p]] <- 1e-9 } }
    }
  } else { return(list(family = NA, params = list())) }
  return(list(family = family, params = final_params))
}

estimate_presymptomatic_proportion <- function(si_dist_info, ip_dist_info, n_samples = N_PRESYMP_SAMPLES) {
  if (is.null(si_dist_info$family) || is.na(si_dist_info$family) || length(si_dist_info$params) == 0 ||
      is.null(ip_dist_info$family) || is.na(ip_dist_info$family) || length(ip_dist_info$params) == 0) { return(NA) }
  get_q_func <- function(d_info) {
    params <- d_info$params
    switch(d_info$family,
      "Lognormal" = function(p) qlnorm(p, meanlog = params$meanlog, sdlog = params$sdlog),
      "Gamma"     = function(p) qgamma(p, shape = params$shape, rate = if (!is.null(params$rate)) params$rate else 1/params$scale),
      "Normal"    = function(p) qnorm(p, mean = params$mean, sd = params$SD),
      "Weibull"   = function(p) qweibull(p, shape = params$shape, scale = params$scale), NULL)
  }
  q_si <- get_q_func(si_dist_info); q_ip <- get_q_func(ip_dist_info); if (is.null(q_si) || is.null(q_ip)) { return(NA) }
  random_quantiles <- runif(n_samples); si_samples <- q_si(random_quantiles); ip_samples <- q_ip(random_quantiles)
  if(si_dist_info$family == "Normal") si_samples[si_samples < 0] <- 0
  if(ip_dist_info$family == "Normal") ip_samples[ip_samples < 0] <- 0
  mean(si_samples < ip_samples, na.rm = TRUE)
}


# --- 4. Main Analysis Function ---
run_full_analysis <- function(data_file_path, clustering_vars, k, output_suffix) {
  
  print(paste("--- Starting Full Analysis for scenario:", output_suffix, "---"))
  print(paste("--- Loading data from:", data_file_path, "---"))
  
  params_df <- read_csv(data_file_path)

  # --- MCMC Loop ---
  print("--- Running MCMC Sampling ---")
  mcmc_results <- list()
  for (iter in 1:N_MCMC_ITERATIONS) {
    if (iter %% max(1, (N_MCMC_ITERATIONS/10)) == 0) { print(paste("MCMC Iteration:", iter, "/", N_MCMC_ITERATIONS)) }
    current_iteration_params_list <- list()
    for (i in 1:nrow(params_df)) {
      pathogen_row <- params_df[i, ]; pathogen_name <- pathogen_row$Pathogen_Name
      r0_sampled <- sample_clust_parameter(pathogen_row, "R0"); si_clust_sampled <- sample_clust_parameter(pathogen_row, "SI"); cfr_sampled <- sample_clust_parameter(pathogen_row, "CFR")
      
      presymp_prop_sampled <- NA
      # First, try to calculate from full distributions
      si_full_dist_info <- get_full_dist_params(pathogen_row, "SI"); ip_full_dist_info <- get_full_dist_params(pathogen_row, "IP")
      if (all(!sapply(list(si_full_dist_info, ip_full_dist_info), function(x) is.null(x) || is.na(x$family) || length(x$params) == 0))) {
        si_full_dist_info$Pathogen_Name <- pathogen_name; ip_full_dist_info$Pathogen_Name <- pathogen_name
        presymp_prop_sampled <- estimate_presymptomatic_proportion(si_full_dist_info, ip_full_dist_info)
      } else {
        # If calculation fails, try to sample it directly as a cluster parameter
        presymp_prop_sampled <- sample_clust_parameter(pathogen_row, "Presymp")
      }

      pathogen_params_df <- data.frame(
        Pathogen_Name = pathogen_name,
        R0_sampled = ifelse(is.null(r0_sampled) || !is.finite(r0_sampled), NA_real_, r0_sampled),
        SI_Clust_sampled = ifelse(is.null(si_clust_sampled) || !is.finite(si_clust_sampled), NA_real_, si_clust_sampled),
        CFR_sampled = ifelse(is.null(cfr_sampled) || !is.finite(cfr_sampled), NA_real_, cfr_sampled),
        Presymp_Proportion_sampled = ifelse(is.null(presymp_prop_sampled) || !is.finite(presymp_prop_sampled), NA_real_, presymp_prop_sampled),
        Route_resp = as.integer(pathogen_row$Route_resp), Route_direct = as.integer(pathogen_row$Route_direct),
        Route_sexual = as.integer(pathogen_row$Route_sexual), Route_animal = as.integer(pathogen_row$Route_animal), Route_vector = as.integer(pathogen_row$Route_vector)
      )
      current_iteration_params_list[[pathogen_name]] <- pathogen_params_df
    }
    if(length(current_iteration_params_list) > 0){ mcmc_results[[iter]] <- bind_rows(current_iteration_params_list) }
  }
  
  if (length(mcmc_results) == 0 || all(sapply(mcmc_results, is.null))){
    print(paste("MCMC sampling produced no results for scenario", output_suffix, ". Aborting analysis."))
    return(NULL)
  }
  all_mcmc_samples_df <- bind_rows(mcmc_results[!sapply(mcmc_results, is.null)], .id = "MCMC_Iteration")
  all_mcmc_samples_df$MCMC_Iteration <- as.integer(all_mcmc_samples_df$MCMC_Iteration)
  print("MCMC sampling complete.")
  
  # --- Clustering Analysis (Per-iteration clustering) ---
  print(paste("--- Starting Per-iteration Clustering (K=", k, ") ---", sep=""))
  
  clustering_data_full <- all_mcmc_samples_df %>% select(MCMC_Iteration, Pathogen_Name, all_of(clustering_vars))
  numerical_param_cols <- intersect(clustering_vars, c("R0_sampled", "SI_Clust_sampled", "CFR_sampled", "Presymp_Proportion_sampled"))
  for(col in numerical_param_cols){
    if(any(is.na(clustering_data_full[[col]]))){
      if(is.numeric(clustering_data_full[[col]])){ clustering_data_full[[col]][is.na(clustering_data_full[[col]])] <- median(clustering_data_full[[col]], na.rm = TRUE) }
    }
  }
  scaled_clustering_data_full <- clustering_data_full
  if(length(numerical_param_cols) > 0) { scaled_clustering_data_full[numerical_param_cols] <- scale(scaled_clustering_data_full[numerical_param_cols]) }

  all_iteration_assignments <- list()
  mcmc_iterations <- unique(scaled_clustering_data_full$MCMC_Iteration)
  
  if (length(mcmc_iterations) > 0) {
    for (iter_val in mcmc_iterations) {
      current_iter_data <- scaled_clustering_data_full %>% filter(MCMC_Iteration == iter_val)
      current_iter_features_scaled <- current_iter_data %>% select(-MCMC_Iteration, -Pathogen_Name)
      if (nrow(current_iter_features_scaled) < k || ncol(current_iter_features_scaled) == 0) { next }
      set.seed(123 + as.integer(iter_val)) 
      kmeans_iter_result <- tryCatch({ kmeans(current_iter_features_scaled, centers = k, nstart = 25) }, error = function(e) { NULL })
      if (!is.null(kmeans_iter_result)) {
        all_iteration_assignments[[as.character(iter_val)]] <- data.frame(MCMC_Iteration = iter_val, Pathogen_Name = current_iter_data$Pathogen_Name, Cluster_Assigned = kmeans_iter_result$cluster)
      }
    }
  }
  
  if (length(all_iteration_assignments) == 0) {
    print(paste("No cluster assignments generated for scenario", output_suffix, ". Aborting analysis."))
    return(NULL)
  }
  combined_assignments_df <- bind_rows(all_iteration_assignments)

  # --- Ensemble Clustering / Consensus Clustering --- 
  print(paste("--- Starting Ensemble/Consensus Clustering ---"))
  pathogen_names <- sort(unique(combined_assignments_df$Pathogen_Name))
  num_pathogens <- length(pathogen_names)
  num_mcmc_iterations_successful <- length(unique(combined_assignments_df$MCMC_Iteration))
  
  coassignment_matrix <- matrix(0, nrow = num_pathogens, ncol = num_pathogens, dimnames = list(pathogen_names, pathogen_names))
  assignments_by_iter <- split(combined_assignments_df[, c("Pathogen_Name", "Cluster_Assigned")], combined_assignments_df$MCMC_Iteration)
  for (iter_assignments in assignments_by_iter) {
    for (i in 1:(num_pathogens - 1)) {
      for (j in (i + 1):num_pathogens) {
        p1_name <- pathogen_names[i]; p2_name <- pathogen_names[j]
        p1_assign <- iter_assignments$Cluster_Assigned[iter_assignments$Pathogen_Name == p1_name]
        p2_assign <- iter_assignments$Cluster_Assigned[iter_assignments$Pathogen_Name == p2_name]
        if (length(p1_assign) > 0 && length(p2_assign) > 0 && p1_assign == p2_assign) {
          coassignment_matrix[p1_name, p2_name] <- coassignment_matrix[p1_name, p2_name] + 1
          coassignment_matrix[p2_name, p1_name] <- coassignment_matrix[p1_name, p2_name]
        }
      }
    }
  }
  diag(coassignment_matrix) <- num_mcmc_iterations_successful
  dissimilarity_matrix <- 1 - (coassignment_matrix / num_mcmc_iterations_successful)
  diag(dissimilarity_matrix) <- 0
  
  pathogen_dist <- as.dist(dissimilarity_matrix)
  hierarchical_clust <- hclust(pathogen_dist, method = "average")
  
  # --- Plot Dendrogram ---
  pathogen_full_name_map <- c("COVID-19_WT"="SARS-CoV-2 (WT)", "COVID-19_A"="SARS-CoV-2 (Alpha)", "COVID-19_D"="SARS-CoV-2 (Delta)", "COVID-19_O"="SARS-CoV-2 (Omicron)", "H1N1_18"="A/H1N1", "H2N2"="A/H2N2", "H3N2"="A/H3N2", "H1N1_09"="A/H1N1/09", "H5N1"="A/H5N1", "Ebola"="EBOV", "Marburg"="MARV", "Mpox"="MPV", "Lassa"="LASV", "Nipah"="NiV", "Zika"="ZIKV", "SARS"="SARS-CoV-1", "MERS"="MERS-CoV", "CCHF"="CCHFV", "RVFV"="RVFV", "HIV"="HIV")
  
  K_CONSENSUS <- k
  cluster_colors <- brewer.pal(n = max(3, K_CONSENSUS), name = "Set2")
  if (K_CONSENSUS > length(cluster_colors)) { cluster_colors <- rep(cluster_colors, length.out = K_CONSENSUS) }
  dend <- as.dendrogram(hierarchical_clust)
  original_labels <- labels(dend); new_labels <- pathogen_full_name_map[original_labels]; new_labels[is.na(new_labels)] <- original_labels[is.na(new_labels)]; labels(dend) <- new_labels
  dend_colored <- color_branches(dend, k = K_CONSENSUS, col = cluster_colors[1:K_CONSENSUS]) %>% set("labels_cex", 0.7) %>% set("branches_lwd", 3)
  
  png_filename <- paste0("Clustering/mcmc/Kmeans/consensus_dendrogram_k", K_CONSENSUS, output_suffix, ".png")
  png(png_filename, width=1200, height=800, units="px", res=100)
  par(mar = c(5, 4, 4, 10)) 
  plot(dend_colored, horiz = TRUE, main = paste("Consensus Dendrogram (K=", K_CONSENSUS, ", Scenario: ", output_suffix, ")"), xlab = "Dissimilarity (1 - Proportion Co-assigned)")
  dev.off()
  print(paste("Consensus dendrogram saved to", png_filename))
  
  # --- Save Results ---
  consensus_clusters <- cutree(hierarchical_clust, k = K_CONSENSUS)
  consensus_assignments_df <- data.frame(Pathogen_Name = names(consensus_clusters), Consensus_Cluster = consensus_clusters)
  write.csv(consensus_assignments_df, paste0("Clustering/mcmc/Kmeans/consensus_assignments_k", K_CONSENSUS, output_suffix, ".csv"))
  
  print(paste("--- Finished Full Analysis for scenario:", output_suffix, "---"))
  return(TRUE)
}

# --- 5. Execute Analyses ---
CHOSEN_K <- 6 # Assuming K=6 for all sensitivity analyses.

# --- Analysis 1: Including RVFV (excluding SI and Presymptomatic Proportion) ---
run_full_analysis(
  data_file_path = "Clustering/mcmc/Kmeans/pathogen_params_RVFV.csv",
  clustering_vars = c("R0_sampled", "CFR_sampled", "Route_resp", "Route_direct", "Route_sexual", "Route_animal", "Route_vector"),
  k = CHOSEN_K,
  output_suffix = "_sens_rvfv"
)

# --- Analysis 2: Including HIV (excluding SI and CFR) ---
run_full_analysis(
  data_file_path = "Clustering/mcmc/Kmeans/pathogen_params_HIV.csv",
  clustering_vars = c("R0_sampled", "Presymp_Proportion_sampled", "Route_resp", "Route_direct", "Route_sexual", "Route_animal", "Route_vector"),
  k = CHOSEN_K,
  output_suffix = "_sens_hiv"
)

print("--- Sensitivity Analysis Script finished. ---") 