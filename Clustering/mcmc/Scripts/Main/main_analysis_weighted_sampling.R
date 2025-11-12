# Main Analysis Script with Uniform Sampling

# --- 1. Load Libraries ---
library(dplyr)
library(readr)

# --- 2. Configuration ---
N_MCMC_ITERATIONS <- 500
N_PRESYMP_SAMPLES <- 500
set.seed(123) # For reproducibility

# --- 3. Load Data ---
params_df <- read_csv("Clustering/mcmc/Kmeans/data/pathogen_params.csv")
routes_df <- read_csv("Clustering/mcmc/Kmeans/data/transmission_route.csv")

# --- 4. Data Preprocessing ---

# Add a check for study usability before any other processing
is_study_parametrizable <- function(df_row) {
  # A study is parametrizable if it has a defined sampling distribution AND
  # sufficient uncertainty information to derive distribution parameters.
  has_dist <- !is.na(df_row$SamplingDist) && df_row$SamplingDist != ""
  has_ci_iqr_range <- !is.na(df_row$UncertaintyType) && df_row$UncertaintyType %in% c("95CI", "IQR", "Range", "range") && !is.na(df_row$LowerBound) && !is.na(df_row$UpperBound)
  has_cv <- !is.na(df_row$PointEstimateCV) && df_row$PointEstimateCV > 0
  
  return(has_dist && (has_ci_iqr_range || has_cv))
}

# Filter the dataset to only include studies from which we can parameterize a distribution for SI and IP
original_rows <- nrow(params_df %>% filter(Parameter %in% c("SI", "IP")))

parametrizable_studies <- params_df %>%
  filter(Parameter %in% c("SI", "IP")) %>%
  rowwise() %>%
  filter(is_study_parametrizable(cur_data())) %>%
  ungroup()

non_parametrizable <- params_df %>%
  filter(Parameter %in% c("SI", "IP")) %>%
  anti_join(parametrizable_studies, by = c("Pathogen_Name", "Parameter", "citation", "ReportedValue"))

if(nrow(non_parametrizable) > 0) {
  print("--- Warning: The following SI/IP studies were excluded because they have a sampling distribution but are missing necessary uncertainty data (e.g., CI, IQR, or CV) to parameterize it: ---")
  print(as.data.frame(non_parametrizable %>% select(Pathogen_Name, Parameter, citation, ReportedValue, SamplingDist, UncertaintyType)))
}

# Keep all other parameters (R0, CFR etc) and combine with the usable SI/IP studies
params_df_filtered <- params_df %>%
  filter(!Parameter %in% c("SI", "IP")) %>%
  bind_rows(parametrizable_studies)

# List of pathogens for the main analysis
main_pathogens <- c(
    "COVID-19_WT", "COVID-19_A", "COVID-19_D", "COVID-19_O", "Ebola", 
    "Marburg", "Lassa", "CCHF", "Nipah", "Zika", "Mpox", "H1N1_18", 
    "H2N2", "H3N2", "H1N1_09", "H5N1", "SARS", "MERS"
)

# Merge and filter data using the pre-filtered params_df
data <- params_df_filtered %>%
  filter(Pathogen_Name %in% main_pathogens) %>%
  left_join(routes_df, by = "Pathogen_Name")

# Define influenza group for pooled presymptomatic analysis
flu_group <- c("H1N1_18", "H2N2", "H3N2", "H1N1_09")

# --- 5. Diagnostic Pre-check for Presymptomatic Calculation Data ---
print("--- Running Diagnostic Pre-check for Presymptomatic Data ---")

# Define what makes a study usable for distribution parameterization
is_study_usable <- function(df) {
  nrow(df %>% filter(
    !is.na(ReportedValue) & 
    !is.na(SamplingDist) & SamplingDist != "" &
    (
      (!is.na(UncertaintyType) & UncertaintyType %in% c("95CI", "IQR", "Range", "range") & !is.na(LowerBound) & !is.na(UpperBound)) |
      (!is.na(PointEstimateCV) & PointEstimateCV > 0)
    )
  )) > 0
}

for (pathogen in main_pathogens) {
  pathogen_data <- data %>% filter(Pathogen_Name == pathogen)
  
  has_usable_si <- FALSE
  has_usable_ip <- FALSE

  # The check must be specific for each pathogen, even if they are in a pooled group.
  si_studies <- pathogen_data %>% filter(Parameter == "SI")
  ip_studies <- pathogen_data %>% filter(Parameter == "IP")
  has_usable_si <- is_study_usable(si_studies)
  has_usable_ip <- is_study_usable(ip_studies)

  # For pathogens in the flu group, we add a note if their data is missing 
  # but they will proceed by using pooled data from other flu strains.
  if (pathogen %in% flu_group) {
    pooled_si_studies <- data %>% filter(Parameter == "SI", Pathogen_Name %in% flu_group)
    pooled_ip_studies <- data %>% filter(Parameter == "IP", Pathogen_Name %in% flu_group)
    can_use_pooled_si <- is_study_usable(pooled_si_studies)
    can_use_pooled_ip <- is_study_usable(pooled_ip_studies)
    
    if (!has_usable_si && can_use_pooled_si) {
        print(paste0("[INFO] ", pathogen, ": Missing specific SI data, but will use pooled flu group data."))
        has_usable_si <- TRUE # It will proceed using the pool
    }
    if (!has_usable_ip && can_use_pooled_ip) {
        print(paste0("[INFO] ", pathogen, ": Missing specific IP data, but will use pooled flu group data."))
        has_usable_ip <- TRUE # It will proceed using the pool
    }
  }
  
  if (has_usable_si && has_usable_ip) {
    print(paste0("[OK] ", pathogen, ": Found usable SI and IP data. Presymptomatic proportion will be calculated."))
  } else {
    print(paste0("[WARNING] ", pathogen, ": Missing usable data. Has SI: ", has_usable_si, ", Has IP: ", has_usable_ip, ". Presymptomatic proportion will be NA."))
  }
}
print("--- Diagnostic Pre-check Finished ---")


# --- 6. Helper Functions ---

# Function to calculate variance from uncertainty
calculate_variance <- function(row) {
  uncert_type <- row["UncertaintyType"]
  lower <- as.numeric(row["LowerBound"])
  upper <- as.numeric(row["UpperBound"])
  val <- as.numeric(row["ReportedValue"])
  cv <- as.numeric(row["PointEstimateCV"])
  
  if (is.na(uncert_type) || uncert_type == "") {
    return(NA)
  }
  
  sd_val <- NA
  if (uncert_type == "95CI") {
    if (!is.na(lower) && !is.na(upper)) {
      sd_val <- (upper - lower) / (2 * qnorm(0.975)) # More precise than 1.96
    }
  } else if (uncert_type == "IQR") {
    if (!is.na(lower) && !is.na(upper)) {
      sd_val <- (upper - lower) / 1.349
    }
  } else if (uncert_type %in% c("Range", "range")) {
     if (!is.na(lower) && !is.na(upper)) {
      sd_val <- (upper - lower) / 4 # Approximation
    }
  } else if (uncert_type == "PointEstimate") {
    if (!is.na(val) && !is.na(cv) && cv > 0) {
      sd_val <- val * cv
    }
  }
  
  if (!is.na(sd_val) && sd_val > 0) {
    return(sd_val^2)
  } else {
    return(NA)
  }
}

# Apply variance calculation
data <- data %>%
  mutate(variance = apply(., 1, calculate_variance))

# Helper function to sample a single study with uniform probability
sample_study <- function(df_group) {
  # sample_n without a weight argument performs uniform sampling
  sample_n(df_group, 1, replace = TRUE)
}

# Function to get parameters for a distribution from a sampled study row
get_dist_params_from_study <- function(study_row) {
    mean_val <- as.numeric(study_row[["ReportedValue"]])
    lower_ci <- as.numeric(study_row[["LowerBound"]])
    upper_ci <- as.numeric(study_row[["UpperBound"]])
    uncert_type <- study_row[["UncertaintyType"]]
    cv <- as.numeric(study_row[["PointEstimateCV"]])
    dist_family <- study_row[["SamplingDist"]]

    if (is.na(mean_val) || is.na(dist_family) || dist_family == "") return(NULL)

    sd_val <- NA
    # Prioritize CI/Range for calculating SD
    if (!is.na(uncert_type) && uncert_type %in% c("95CI", "IQR", "Range", "range") && !is.na(lower_ci) && !is.na(upper_ci)) {
        if (uncert_type == "95CI") {
            sd_val <- (upper_ci - lower_ci) / (2 * qnorm(0.975))
        } else if (uncert_type == "IQR") {
            sd_val <- (upper_ci - lower_ci) / 1.349
        } else { # Range
            sd_val <- (upper_ci - lower_ci) / 4 # Approximation
        }
    } 
    
    # Fallback to CV if SD is still NA
    if (is.na(sd_val) && !is.na(cv) && cv > 0) {
        sd_val <- mean_val * cv
    }

    # If we still couldn't derive an SD, we can't proceed
    if (is.na(sd_val) || sd_val <= 0) return(NULL)

    params <- list()
    family_clean <- gsub("r", "", dist_family)
    
    if (family_clean == "gamma") {
        shape <- (mean_val / sd_val)^2
        rate <- mean_val / sd_val^2
        if(is.na(shape) || shape <= 0 || is.na(rate) || rate <= 0) return(NULL)
        params <- list(shape = shape, rate = rate)
    } else if (family_clean == "lnorm") {
        # Need to convert mean/sd to meanlog/sdlog for lnorm
        if (mean_val <= 0) return(NULL) # Lognormal must have a positive mean
        sdlog <- sqrt(log(sd_val^2 / mean_val^2 + 1))
        meanlog <- log(mean_val) - 0.5 * sdlog^2
        if(is.na(sdlog) || sdlog <= 0 || is.na(meanlog)) return(NULL)
        params <- list(meanlog = meanlog, sdlog = sdlog)
    } else if (family_clean == "weibull") {
        # This conversion is more complex, often requiring optimization.
        # Using a simpler approximation for now.
        # Note: This is a rough approximation.
        k <- (sd_val / mean_val)^(-1.086)
        lambda <- mean_val / (gamma(1 + 1/k))
        if(is.na(k) || k <= 0 || is.na(lambda) || lambda <= 0) return(NULL)
        params <- list(shape = k, scale = lambda)
    } else {
        # For other distributions like 'norm', 'unif', etc.
        # We might not have enough info, but can handle 'norm'
        if (family_clean == "norm"){
            params <- list(mean = mean_val, sd = sd_val)
        } else {
            return(NULL) # Can't handle other distributions for now
        }
    }
    
    return(list(family = family_clean, params = params))
}


# Presymptomatic proportion estimation
estimate_presymptomatic_proportion <- function(si_dist_info, ip_dist_info, n_samples = N_PRESYMP_SAMPLES) {
  # This function is adapted from the original script
  if (is.null(si_dist_info) || is.null(ip_dist_info) || is.null(si_dist_info$family) || is.null(ip_dist_info$family)) {
    return(NA)
  }

  get_q_func <- function(d_info) {
    params <- d_info$params
    switch(d_info$family,
      "lnorm" = function(p) qlnorm(p, meanlog = params$meanlog, sdlog = params$sdlog),
      "gamma"     = function(p) qgamma(p, shape = params$shape, rate = params$rate),
      "weibull"   = function(p) qweibull(p, shape = params$shape, scale = params$scale),
      "norm"      = function(p) qnorm(p, mean = params$mean, sd = params$sd),
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

  # Truncate negative values that might arise from Normal distributions
  if(si_dist_info$family == "norm") si_samples[si_samples < 0] <- 0
  if(ip_dist_info$family == "norm") ip_samples[ip_samples < 0] <- 0
  
  mean(si_samples < ip_samples, na.rm = TRUE)
}


# --- 6. Main MCMC Loop ---
mcmc_results <- list()

for (iter in 1:N_MCMC_ITERATIONS) {
  if (iter %% (N_MCMC_ITERATIONS / 10) == 0) {
    print(paste("MCMC Iteration:", iter, "/", N_MCMC_ITERATIONS))
  }
  
  current_iteration_params_list <- list()
  
  for (pathogen in main_pathogens) {
    pathogen_data <- data %>% filter(Pathogen_Name == pathogen)
    
    # Weighted sampling for R0, SI, CFR
    r0_study <- sample_study(pathogen_data %>% filter(Parameter == "R0"))
    si_study_for_clust <- sample_study(pathogen_data %>% filter(Parameter == "SI"))
    cfr_study <- sample_study(pathogen_data %>% filter(Parameter == "CFR"))
    
    r0_sampled <- rnorm(1, mean = r0_study$ReportedValue, sd = sqrt(r0_study$variance))
    si_clust_sampled <- rnorm(1, mean = si_study_for_clust$ReportedValue, sd = sqrt(si_study_for_clust$variance))
    cfr_sampled <- rnorm(1, mean = cfr_study$ReportedValue, sd = sqrt(cfr_study$variance))

    # Clamp values to be realistic
    r0_sampled <- max(0, r0_sampled)
    si_clust_sampled <- max(0, si_clust_sampled)
    cfr_sampled <- pmax(0, pmin(1, cfr_sampled))
    
    # Weighted sampling for SI and IP for presymptomatic estimation
    if (pathogen %in% flu_group) {
        # Pool studies for the specified flu strains
        si_studies_for_presymp <- data %>% filter(Parameter == "SI", Pathogen_Name %in% flu_group)
        ip_studies_for_presymp <- data %>% filter(Parameter == "IP", Pathogen_Name %in% flu_group)
    } else {
        # Use pathogen-specific studies for all other pathogens
        si_studies_for_presymp <- pathogen_data %>% filter(Parameter == "SI")
        ip_studies_for_presymp <- pathogen_data %>% filter(Parameter == "IP")
    }
    
    presymp_prop_sampled <- NA
    
    if(nrow(si_studies_for_presymp) > 0 && nrow(ip_studies_for_presymp) > 0) {
        si_study_presymp <- sample_study(si_studies_for_presymp)
        ip_study_presymp <- sample_study(ip_studies_for_presymp)
        
        si_dist_info <- get_dist_params_from_study(si_study_presymp)
        ip_dist_info <- get_dist_params_from_study(ip_study_presymp)
        
        if (!is.null(si_dist_info) && !is.null(ip_dist_info)) {
            presymp_prop_sampled <- estimate_presymptomatic_proportion(si_dist_info, ip_dist_info)
        }
    }
    
    pathogen_routes <- pathogen_data[1,] %>% select(starts_with("Route_"))

    pathogen_params_df <- data.frame(
      Pathogen_Name = pathogen,
      R0_sampled = r0_sampled,
      SI_Clust_sampled = si_clust_sampled,
      CFR_sampled = cfr_sampled,
      Presymp_Proportion_sampled = presymp_prop_sampled,
      Route_resp = pathogen_routes$Route_resp,
      Route_direct = pathogen_routes$Route_direct,
      Route_sexual = pathogen_routes$Route_sexual,
      Route_animal = pathogen_routes$Route_animal,
      Route_vector = pathogen_routes$Route_vector
    )
    current_iteration_params_list[[pathogen]] <- pathogen_params_df
  }
  
  mcmc_results[[iter]] <- bind_rows(current_iteration_params_list)
}

all_mcmc_samples_df <- bind_rows(mcmc_results, .id = "MCMC_Iteration")
all_mcmc_samples_df$MCMC_Iteration <- as.integer(all_mcmc_samples_df$MCMC_Iteration)

# --- 7. Save MCMC Samples ---
write_csv(all_mcmc_samples_df, "Clustering/mcmc/Kmeans/main_outputs/mcmc_parameter_samples_uniform.csv")
print("Uniformly sampled MCMC parameter samples saved.")

# --- 8. Clustering and Consensus Analysis (Adapted from original script) ---
# This part of the script would be largely the same as mcmc_Kmeans_analysis.R,
# starting from section --- 8. Clustering Analysis (K=6) ---.
# It would use `all_mcmc_samples_df` as input. For brevity, I am not repeating
# the entire clustering, visualization, and consensus analysis code here, but it
# should be copied from the original script and appended here.

print("Script finished. The subsequent clustering analysis can now be run on the weighted samples.")


# --- 8. Clustering Analysis (K=6) ---
if (exists("all_mcmc_samples_df") && nrow(all_mcmc_samples_df) > 0) {
  print("--- Starting Clustering Analysis (Per MCMC Iteration, K=6) ---")
  
  CHOSEN_K <- 6 #specified K=6
  
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
          write_csv(centroid_summary, "Clustering/mcmc/Kmeans/main_outputs/cluster_centroids_summary_with_ci_k6_uniform.csv")
          print("Cluster centroids summary saved to Clustering/mcmc/Kmeans/main_outputs/cluster_centroids_summary_with_ci_k6_uniform.csv")
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
          write_csv(modal_assignments, "Clustering/mcmc/Kmeans/main_outputs/pathogen_modal_cluster_assignments_k6_uniform.csv")
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
      
      # Scale these summarized numerical features for PCA (routes are proportions 0-1)
      scaled_features_for_pca_plot <- features_for_pca_plot
      scaled_features_for_pca_plot[plot_features_numerical] <- scale(scaled_features_for_pca_plot[plot_features_numerical])
      
      # Final check to remove any NA/NaN/Inf after scaling, replacing with 0
      scaled_features_for_pca_plot[!is.finite(as.matrix(scaled_features_for_pca_plot))] <- 0

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
          ggsave(paste0("Clustering/mcmc/Kmeans/main_outputs/pca_modal_cluster_plot_k",CHOSEN_K,"_uniform.png"), plot = pca_plot_modal, width=10, height=8)
          print(paste0("PCA plot with modal cluster assignments saved to Clustering/mcmc/Kmeans/main_outputs/pca_modal_cluster_plot_k",CHOSEN_K,"_uniform.png"))
          
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


# --- 9. Ensemble Clustering / Consensus Clustering (for K=4 through K=8) --- 
if (exists("combined_assignments_df") && nrow(combined_assignments_df) > 0 && 
    exists("all_mcmc_samples_df") && nrow(all_mcmc_samples_df) > 0) {
  
  print("--- Starting Ensemble/Consensus Clustering Analysis for K=4 to K=8 ---")
  
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

  for (iter_idx in 1:length(assignments_by_iter)) {
    if (iter_idx %% max(1, round(length(assignments_by_iter)/10)) == 0) {
      print(paste("Processing co-assignments for MCMC iteration", iter_idx, "/", length(assignments_by_iter)))
    }
    current_iter_assignments_df <- assignments_by_iter[[iter_idx]]
    
    for (i in 1:(num_pathogens - 1)) {
      for (j in (i + 1):num_pathogens) {
        pathogen1_name <- pathogen_names[i]
        pathogen2_name <- pathogen_names[j]
        
        p1_assignment_row <- current_iter_assignments_df[current_iter_assignments_df$Pathogen_Name == pathogen1_name, ]
        p2_assignment_row <- current_iter_assignments_df[current_iter_assignments_df$Pathogen_Name == pathogen2_name, ]

        if (nrow(p1_assignment_row) > 0 && nrow(p2_assignment_row) > 0) {
          if (p1_assignment_row$Cluster_Assigned == p2_assignment_row$Cluster_Assigned) {
            coassignment_matrix[pathogen1_name, pathogen2_name] <- coassignment_matrix[pathogen1_name, pathogen2_name] + 1
            coassignment_matrix[pathogen2_name, pathogen1_name] <- coassignment_matrix[pathogen1_name, pathogen2_name] # Keep symmetric
          }
        }
      }
    }
  }
  diag(coassignment_matrix) <- num_mcmc_iterations 
  
  print("Co-assignment matrix (first 6x6):")
  print(head(coassignment_matrix[,1:min(6, ncol(coassignment_matrix))]))
  write.csv(coassignment_matrix, "Clustering/mcmc/Kmeans/main_outputs/coassignment_matrix_uniform.csv")

  # --- 9.2. Calculate Dissimilarity Matrix ---
  print("Calculating dissimilarity matrix...")
  dissimilarity_matrix <- 1 - (coassignment_matrix / N_MCMC_ITERATIONS) 
  diag(dissimilarity_matrix) <- 0 

  print("Dissimilarity matrix (first 6x6):")
  print(head(dissimilarity_matrix[,1:min(6, ncol(dissimilarity_matrix))]))
  write.csv(dissimilarity_matrix, "Clustering/mcmc/Kmeans/main_outputs/dissimilarity_matrix_uniform.csv")

  # --- 9.3. Perform Hierarchical Clustering ---
  print("Performing hierarchical clustering...")
  pathogen_dist <- as.dist(dissimilarity_matrix)
  hierarchical_clust <- hclust(pathogen_dist, method = "average")

  # --- 9.4. Loop through K values to generate outputs for each ---
  for (K_CONSENSUS in 4:8) {
    print(paste("--- Generating outputs for K =", K_CONSENSUS, "---"))

    # --- Plot Dendrogram (Revised with dendextend) ---
    print("Plotting enhanced dendrogram...")
    if (!requireNamespace("RColorBrewer", quietly = TRUE)) { install.packages("RColorBrewer", quiet = TRUE) }
    library(RColorBrewer)
    
    cluster_colors <- brewer.pal(n = max(3, K_CONSENSUS), name = "Set2") 
    if (K_CONSENSUS > length(cluster_colors)) { 
        cluster_colors <- rep(cluster_colors, length.out = K_CONSENSUS) 
    }

    dend <- as.dendrogram(hierarchical_clust)
    
    original_labels <- labels(dend)
    new_labels <- pathogen_full_name_map[original_labels]
    new_labels[is.na(new_labels)] <- original_labels[is.na(new_labels)]
    labels(dend) <- new_labels
    
    dend_colored <- color_branches(dend, k = K_CONSENSUS, col = cluster_colors[1:K_CONSENSUS])
    dend_colored <- set(dend_colored, "labels_cex", 0.7)
    dend_colored <- set(dend_colored, "branches_lwd", 3)
    
    png_filename <- paste0("Clustering/mcmc/Kmeans/figures/figure.2_k", K_CONSENSUS, "_uniform.png")
    png(png_filename, width=1200, height=800, units="px", res=100)
    par(mar = c(5, 4, 4, 10))
    plot(dend_colored, horiz = TRUE, 
         main = paste("Consensus Clustering (K=", K_CONSENSUS, ")"), 
         xlab = "Dissimilarity (1 - Proportion Co-assigned)")

    dev.off()
    print(paste("Colored consensus dendrogram saved to", png_filename))

    # --- Extract Consensus Cluster Assignments ---
    print(paste("Cutting tree to get", K_CONSENSUS, "consensus clusters..."))
    consensus_clusters <- cutree(hierarchical_clust, k = K_CONSENSUS)
    consensus_assignments_df <- data.frame(Pathogen_Name = names(consensus_clusters),
                                           Consensus_Cluster = consensus_clusters)
    print("Consensus cluster assignments:")
    print(consensus_assignments_df)
    write.csv(consensus_assignments_df, paste0("Clustering/mcmc/Kmeans/main_outputs/consensus_cluster_assignments_k", K_CONSENSUS, "_uniform.csv"))

    # --- Characterize Consensus Clusters ---
    print("Characterizing consensus clusters (means and 95% CIs from all_mcmc_samples_df)...")
    
    params_to_summarize <- c("R0_sampled", "SI_Clust_sampled", "CFR_sampled", "Presymp_Proportion_sampled",
                             "Route_resp", "Route_direct", "Route_sexual", "Route_animal", "Route_vector")

    pathogen_mcmc_samples_with_consensus_clusters <- all_mcmc_samples_df %>% 
      filter(Pathogen_Name %in% consensus_assignments_df$Pathogen_Name) %>%
      left_join(consensus_assignments_df, by = "Pathogen_Name")

    if(nrow(pathogen_mcmc_samples_with_consensus_clusters) > 0 && "Consensus_Cluster" %in% names(pathogen_mcmc_samples_with_consensus_clusters)){
      consensus_cluster_summary_list <- list()
      for (param in params_to_summarize) {
        if (param %in% names(pathogen_mcmc_samples_with_consensus_clusters)) {
          summary_for_param <- pathogen_mcmc_samples_with_consensus_clusters %>%
            group_by(Consensus_Cluster) %>%
            summarise(
              mean_val = mean(get(param), na.rm = TRUE),
              lower_ci = quantile(get(param), 0.025, na.rm = TRUE),
              upper_ci = quantile(get(param), 0.975, na.rm = TRUE),
              median_val = median(get(param), na.rm = TRUE),
              sd_val = sd(get(param), na.rm=TRUE),
              n_mcmc_samples_in_calc = sum(!is.na(get(param))),
              .groups = 'drop'
            ) %>% 
            rename_with(~paste0(param, "_", .), .cols = c(mean_val, lower_ci, upper_ci, median_val, sd_val, n_mcmc_samples_in_calc))
          consensus_cluster_summary_list[[param]] <- summary_for_param
        } else {
          warning(paste("Parameter", param, "not found in pathogen_mcmc_samples_with_consensus_clusters."))
        }
      }
      
      if(length(consensus_cluster_summary_list) > 0){
        final_consensus_summary_df <- Reduce(function(x, y) full_join(x, y, by = "Consensus_Cluster"), consensus_cluster_summary_list)
        print("Summary of Consensus Cluster Characteristics (Mean, 95% CI, Median, SD):")
        print(as.data.frame(final_consensus_summary_df))
        write_csv(final_consensus_summary_df, paste0("Clustering/mcmc/Kmeans/main_outputs/consensus_clusters_summary_k", K_CONSENSUS, "_uniform.csv"))
        print(paste0("Consensus clusters summary saved to Clustering/mcmc/Kmeans/main_outputs/consensus_clusters_summary_k", K_CONSENSUS, "_uniform.csv"))
      } else {
        print("No parameters were summarized for consensus clusters.")
      }
    } else {
      print("Could not prepare data for characterizing consensus clusters.")
    }
    
    # --- Create and Save Formatted Summary Table ---
    if(exists("final_consensus_summary_df") && nrow(final_consensus_summary_df) > 0) {
      print("Creating formatted summary table of archetype parameters...")

      # Function to format Mean (95% CI)
      format_ci <- function(mean_val, lower_ci, upper_ci) {
        # Handle potential NA values from the presymptomatic proportion calculation
        if (any(is.na(c(mean_val, lower_ci, upper_ci)))) {
          return("NA")
        }
        sprintf("%.2f (%.2f - %.2f)", mean_val, lower_ci, upper_ci)
      }

      # Use rowwise to apply the function for each cluster
      formatted_summary_table <- final_consensus_summary_df %>%
        rowwise() %>%
        transmute(
          Archetype = Consensus_Cluster,
          R0 = format_ci(R0_sampled_mean_val, R0_sampled_lower_ci, R0_sampled_upper_ci),
          `Serial Interval` = format_ci(SI_Clust_sampled_mean_val, SI_Clust_sampled_lower_ci, SI_Clust_sampled_upper_ci),
          CFR = format_ci(CFR_sampled_mean_val, CFR_sampled_lower_ci, CFR_sampled_upper_ci),
          `Presymptomatic Proportion` = format_ci(Presymp_Proportion_sampled_mean_val, Presymp_Proportion_sampled_lower_ci, Presymp_Proportion_sampled_upper_ci),
          `Respiratory Route (%)` = sprintf("%.0f", Route_resp_mean_val * 100),
          `Direct Contact Route (%)` = sprintf("%.0f", Route_direct_mean_val * 100),
          `Sexual Route (%)` = sprintf("%.0f", Route_sexual_mean_val * 100),
          `Animal Contact Route (%)` = sprintf("%.0f", Route_animal_mean_val * 100),
          `Vector-borne Route (%)` = sprintf("%.0f", Route_vector_mean_val * 100)
        ) %>%
        ungroup()

      print("Formatted Archetype Parameters:")
      print(formatted_summary_table)

      # Save the formatted table to a text file
      output_table_filename <- paste0("Clustering/mcmc/Kmeans/main_outputs/archetype_parameters_summary_table_k", K_CONSENSUS, "_uniform.txt")
      
      write.table(formatted_summary_table, file = output_table_filename, sep = "\t", row.names = FALSE, quote = FALSE)
      
      print(paste("Formatted summary table saved to", output_table_filename))
    }
    
    print(paste("--- Finished K =", K_CONSENSUS, "---"))
  }
  
  print("Ensemble/Consensus Clustering Analysis finished.")

} else {
  print("Skipping Ensemble/Consensus Clustering: combined_assignments_df or all_mcmc_samples_df not found or empty.")
}
