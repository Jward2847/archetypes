# MCMC K-Means Clustering Analysis for Supplementary Figure S7 2.0

# --- 1. Load Libraries ---
library(dplyr)
library(readr)
# Add any other necessary libraries here, e.g., for specific distributions or optimization
# library(DistributionFitR) # May be useful for fitting distributions to CI


# --- 2. Configuration ---
N_MCMC_ITERATIONS <- 500 # Number of MCMC iterations 
set.seed(123) # For reproducibility

# Define output directory for this specific analysis
output_dir <- "Clustering/mcmc/Kmeans/S7_2.0_outputs"
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)


# --- 3. Load Data ---
# Load the new long-format parameter data and transmission routes
params_long_df <- read_csv("Clustering/mcmc/Kmeans/data/pathogen_params.csv")
transmission_df <- read_csv("Clustering/mcmc/Kmeans/data/transmission_route.csv")

# Filter for the specific pathogens of interest for this supplementary figure
pathogens_of_interest <- c(
  "COVID-19_WT", "COVID-19_D", "COVID-19_O", "Ebola", 
  "Marburg", "Mpox", "H1N1_18", "H1N1_09", "SARS", "MERS"
)
params_long_df <- params_long_df %>% filter(Pathogen_Name %in% pathogens_of_interest)
transmission_df <- transmission_df %>% filter(Pathogen_Name %in% pathogens_of_interest)


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
    } else if (param_name %in% c("R0", "SI", "IP", "LP", "InfP", "k")) {
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
  } else if (param_name %in% c("R0", "SI", "IP", "LP", "InfP", "k")) {
    final_value <- pmax(0.00001, final_value)
  }
  
  return(final_value)
}


# --- 5. Main MCMC Loop ---
mcmc_results <- list()
unique_pathogens <- unique(params_long_df$Pathogen_Name)
flu_group <- c("H1N1_09", "H1N1_18", "H2N2", "H3N2") # Define Influenza group for pooling

for (iter in 1:N_MCMC_ITERATIONS) {
  if (iter %% max(1, (N_MCMC_ITERATIONS/10)) == 0) { # Ensure non-zero divisor for modulo
    print(paste("MCMC Iteration:", iter, "/", N_MCMC_ITERATIONS))
  }
  
  current_iteration_params_list <- list()
  
  for (pathogen_name in unique_pathogens) {
    
    # --- Sample parameters using bootstrap aggregation ---
    r0_sampled <- sample_parameter_bootstrap_aggregation("R0", pathogen_name, params_long_df)
    si_clust_sampled <- sample_parameter_bootstrap_aggregation("SI", pathogen_name, params_long_df)
    cfr_sampled <- sample_parameter_bootstrap_aggregation("CFR", pathogen_name, params_long_df)
    ip_sampled <- sample_parameter_bootstrap_aggregation("IP", pathogen_name, params_long_df)
    lp_sampled <- sample_parameter_bootstrap_aggregation("LP", pathogen_name, params_long_df)
    infp_sampled <- sample_parameter_bootstrap_aggregation("InfP", pathogen_name, params_long_df)
    k_sampled <- sample_parameter_bootstrap_aggregation("k", pathogen_name, params_long_df)
    
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
    all_mcmc_samples_df <- bind_rows(mcmc_results[!sapply(mcmc_results, is.null)], .id = "MCMC_Iteration") # Filter out NULL iterations
    if(nrow(all_mcmc_samples_df) > 0){
        all_mcmc_samples_df$MCMC_Iteration <- as.integer(all_mcmc_samples_df$MCMC_Iteration)

        print("MCMC sampling complete.")
        print(paste("Generated", length(unique(all_mcmc_samples_df$MCMC_Iteration)), "valid sets of parameters for", length(unique(all_mcmc_samples_df$Pathogen_Name)), "pathogens."))
    
        # --- 6. Post-MCMC Analysis (Clustering) ---
        # Example: Save the samples for now. Clustering can be a separate script or section.

        # --- 7. Save Results ---
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
  
  # --- 8.1. Prepare Data for Clustering (from all_mcmc_samples_df) ---
  # Select features that will be used. Pathogen_Name and MCMC_Iteration are for tracking.
  clustering_data_full <- all_mcmc_samples_df %>%
    select(
      MCMC_Iteration,
      Pathogen_Name,
      R0_sampled, 
      SI_Clust_sampled, 
      CFR_sampled, 
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

  # Handle NAs in the sampled parameters (e.g., if Presymp_Proportion_sampled was NA for some)
  # Impute with the overall median for that parameter column
  numerical_param_cols <- c("R0_sampled", "SI_Clust_sampled", "CFR_sampled", "IP_sampled", "LP_sampled", "InfP_sampled", "k_sampled")
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
              IP_sampled = mean(IP_sampled, na.rm = TRUE),
              LP_sampled = mean(LP_sampled, na.rm = TRUE),
              InfP_sampled = mean(InfP_sampled, na.rm = TRUE),
              k_sampled = mean(k_sampled, na.rm = TRUE),
              Route_resp = mean(Route_resp, na.rm=TRUE), # Mean for binary gives proportion
              Route_direct = mean(Route_direct, na.rm=TRUE),
              Route_sexual = mean(Route_sexual, na.rm=TRUE),
              Route_animal = mean(Route_animal, na.rm=TRUE),
              Route_vector = mean(Route_vector, na.rm=TRUE),
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
          summary_filename <- paste0("cluster_centroids_summary_with_ci_k", CHOSEN_K, ".csv")
          write_csv(centroid_summary, file.path(output_dir, summary_filename))
          print(paste("Cluster centroids summary saved to", file.path(output_dir, summary_filename)))
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
          modal_filename <- paste0("pathogen_modal_cluster_assignments_k", CHOSEN_K, ".csv")
          write_csv(modal_assignments, file.path(output_dir, modal_filename))
          print(paste("Modal cluster assignments saved to", file.path(output_dir, modal_filename)))
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
              IP_mean_overall = mean(IP_sampled, na.rm = TRUE),
              LP_mean_overall = mean(LP_sampled, na.rm = TRUE),
              InfP_mean_overall = mean(InfP_sampled, na.rm = TRUE),
              k_mean_overall = mean(k_sampled, na.rm = TRUE),
              Route_resp_mean = mean(Route_resp, na.rm=TRUE), # Mean for binary gives proportion
              Route_direct_mean = mean(Route_direct, na.rm=TRUE),
              Route_sexual_mean = mean(Route_sexual, na.rm=TRUE),
              Route_animal_mean = mean(Route_animal, na.rm=TRUE),
              Route_vector_mean = mean(Route_vector, na.rm=TRUE),
              .groups = 'drop'
          ) %>%
          left_join(modal_assignments, by = "Pathogen_Name")

      plot_features_numerical <- c("R0_mean_overall", "SI_Clust_mean_overall", "CFR_mean_overall", "IP_mean_overall", "LP_mean_overall", "InfP_mean_overall", "k_mean_overall")
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
          pca_filename <- paste0("pca_modal_cluster_plot_k",CHOSEN_K,".png")
          ggsave(file.path(output_dir, pca_filename), plot = pca_plot_modal, width=10, height=8)
          print(paste0("PCA plot with modal cluster assignments saved to ", file.path(output_dir, pca_filename)))
          
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
                    if (!requireNamespace("cluster", quietly = TRUE)) {
                        print("Package 'cluster' is not installed, cannot generate silhouette plot.")
                    } else {
                        library(cluster)
                        sil <- silhouette(kmeans_example_iter_result$cluster, dist(example_iter_features))
                        sil_plot_example_iter <- fviz_silhouette(sil, ggtheme = theme_minimal()) +
                                              labs(title=paste("Silhouette Plot for K=", CHOSEN_K, "(Example MCMC Iteration:", example_iter_val, ")"))
                        print(sil_plot_example_iter)
                        sil_iter_filename <- paste0("silhouette_plot_k",CHOSEN_K,"_iter",example_iter_val,".png")
                        ggsave(file.path(output_dir, sil_iter_filename), plot = sil_plot_example_iter, width=8, height=6)
                        print(paste0("Silhouette plot for K=", CHOSEN_K, " saved for example iteration ", example_iter_val, " to ", file.path(output_dir, sil_iter_filename)))
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
    "Mpox" = "MPV", 
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
  write.csv(coassignment_matrix, file.path(output_dir, "coassignment_matrix.csv"))

  # --- 9.2. Calculate Dissimilarity Matrix ---
  print("Calculating dissimilarity matrix...")
  # N_MCMC_ITERATIONS is defined at the top of the script
  dissimilarity_matrix <- 1 - (coassignment_matrix / N_MCMC_ITERATIONS) 
  # Ensure diagonal is 0 after division, if any NA from 0/0 (though N_MCMC_ITERATIONS is likely >0)
  diag(dissimilarity_matrix) <- 0 

  print("Dissimilarity matrix (first 6x6):")
  print(head(dissimilarity_matrix[,seq_len(min(6, ncol(dissimilarity_matrix)))]))
  write.csv(dissimilarity_matrix, file.path(output_dir, "dissimilarity_matrix.csv"))

  # --- 9.3. Perform Hierarchical Clustering ---
  print("Performing hierarchical clustering...")
  # Convert to dist object for hclust
  pathogen_dist <- as.dist(dissimilarity_matrix)
  hierarchical_clust <- hclust(pathogen_dist, method = "average") # Changed back to average

  # --- 9.4. Find Optimal K for Consensus Clustering ---
  print("--- Starting Optimal K Analysis for Consensus Clustering (Silhouette Method) ---")
  
  if (!requireNamespace("cluster", quietly = TRUE)) { install.packages("cluster", quiet = TRUE) }
  if (!requireNamespace("ggplot2", quietly = TRUE)) { install.packages("ggplot2", quiet = TRUE) }
  library(cluster)
  library(ggplot2)
  
  k_range <- 2:10 # Range of K to test
  avg_silhouette_widths <- numeric(length(k_range))
  names(avg_silhouette_widths) <- k_range
  
  for (k in k_range) {
    cluster_assignments <- cutree(hierarchical_clust, k = k)
    silhouette_info <- silhouette(cluster_assignments, dmatrix = as.matrix(dissimilarity_matrix))
    if (is.null(silhouette_info) || !is.matrix(silhouette_info) || nrow(silhouette_info) == 0) {
        avg_width <- NA
        warning(paste("Could not compute silhouette info for k =", k))
    } else {
        avg_width <- summary(silhouette_info)$avg.width
    }
    avg_silhouette_widths[as.character(k)] <- avg_width
  }
  
  silhouette_results_df <- data.frame(K = k_range, Average_Silhouette_Width = avg_silhouette_widths)
  
  optimal_k <- silhouette_results_df$K[which.max(silhouette_results_df$Average_Silhouette_Width)]
  if (length(optimal_k) == 0) {
      warning("Could not determine optimal K. Defaulting to 4.")
      optimal_k <- 4
  }
  print(paste("Optimal K determined as:", optimal_k))

  optimal_k_plot <- ggplot(silhouette_results_df, aes(x = K, y = Average_Silhouette_Width)) +
    geom_line(color = "steelblue", size = 1) + geom_point(color = "steelblue", size = 3) +
    geom_point(data = silhouette_results_df[which.max(silhouette_results_df$Average_Silhouette_Width), ],
               aes(x = K, y = Average_Silhouette_Width), color = "red", size = 5, shape = 8) +
    scale_x_continuous(breaks = k_range) +
    labs(title = "Optimal Number of Clusters (K)", x = "Number of Clusters (K)", y = "Average Silhouette Width") +
    theme_minimal(base_size = 14)
  
  plot_filename <- file.path(output_dir, "optimal_k_plot.png")
  ggsave(plot_filename, plot = optimal_k_plot, width = 8, height = 6)
  print(paste("Optimal K plot saved to", plot_filename))
  
  results_filename <- file.path(output_dir, "optimal_k_silhouette_results.csv")
  write.csv(silhouette_results_df, results_filename, row.names = FALSE)
  print(paste("Optimal K analysis results saved to", results_filename))


  # --- 9.5. Plot Dendrogram (Revised with dendextend) ---
  print("Plotting enhanced dendrogram with dendextend...")
  
  K_CONSENSUS <- optimal_k
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
  
  png_filename <- file.path(output_dir, "figure.S7_2.0.png")
  png(png_filename, width=1200, height=800, units="px", res=100)
  par(mar = c(5, 4, 4, 10)) # Adjust right margin for labels
  plot(dend_colored, horiz = TRUE, 
       main = paste(""), 
       xlab = "Dissimilarity (1 - Proportion Co-assigned)")

  dev.off()
  print(paste("Colored consensus dendrogram saved to", png_filename))

  # --- 9.6. Extract Consensus Cluster Assignments ---
  print(paste("Cutting tree to get", K_CONSENSUS, "consensus clusters..."))
  consensus_clusters <- cutree(hierarchical_clust, k = K_CONSENSUS)
  consensus_assignments_df <- data.frame(Pathogen_Name = names(consensus_clusters),
                                         Consensus_Cluster = consensus_clusters)
  print("Consensus cluster assignments:")
  print(consensus_assignments_df)
  write.csv(consensus_assignments_df, file.path(output_dir, paste0("consensus_cluster_assignments_k", K_CONSENSUS, ".csv")))
  
  # --- 9.7. Silhouette Plot for Consensus Clustering ---
  if (requireNamespace("cluster", quietly = TRUE) && requireNamespace("factoextra", quietly = TRUE)) {
      library(cluster)
      library(factoextra)
      
      print("Generating silhouette plot for consensus clustering...")
      
      # Use the dissimilarity matrix directly
      sil_consensus <- silhouette(consensus_clusters, dmatrix = dissimilarity_matrix)
      
      sil_plot_consensus <- fviz_silhouette(sil_consensus, ggtheme = theme_minimal()) +
          labs(title = paste("Silhouette Plot for Consensus Clustering (K=", K_CONSENSUS, ")"))
          
      print(sil_plot_consensus)
      
      sil_plot_filename <- file.path(output_dir, paste0("silhouette_plot_consensus_k", K_CONSENSUS, ".png"))
      ggsave(sil_plot_filename, plot = sil_plot_consensus, width = 8, height = 6)
      print(paste("Consensus silhouette plot saved to", sil_plot_filename))
  } else {
      print("Packages 'cluster' or 'factoextra' not installed. Skipping consensus silhouette plot.")
  }

  # --- 9.8. Characterize Consensus Clusters ---
  print("Characterizing consensus clusters (means and 95% CIs from all_mcmc_samples_df)...")
  
  # Select the numerical and route columns to summarize from all_mcmc_samples_df
  # These are the original sampled values, not the means from previous steps.
  params_to_summarize <- c("R0_sampled", "SI_Clust_sampled", "CFR_sampled",
                           "IP_sampled", "LP_sampled", "InfP_sampled", "k_sampled",
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
      summary_filename <- paste0("consensus_clusters_summary_k", K_CONSENSUS, ".csv")
      write_csv(final_consensus_summary_df, file.path(output_dir, summary_filename))
      print(paste0("Consensus clusters summary saved to ", file.path(output_dir, summary_filename)))
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
