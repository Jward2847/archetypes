# Supplementary Figure S7 2.0 - Analysis without presymptomatic calculation

# --- 1. Load Libraries ---
library(dplyr)
library(readr)
# Add any other necessary libraries here, e.g., for specific distributions or optimization
# library(DistributionFitR) # May be useful for fitting distributions to CI


# --- 2. Configuration ---
N_MCMC_ITERATIONS <- 500 # Number of MCMC iterations 
set.seed(123) # For reproducibility


# --- 3. Load Data ---
params_df <- read_csv("Clustering/mcmc/Kmeans/data/pathogen_params.csv")
routes_df <- read_csv("Clustering/mcmc/Kmeans/data/transmission_route.csv")


# --- 4. Data Preprocessing ---

# List of pathogens for this analysis
main_pathogens <- c(
    "COVID-19_WT", "COVID-19_D", "COVID-19_O", "Ebola", "Marburg", "Mpox", 
    "H1N1_18", "H1N1_09", "SARS", "MERS"
)

# Merge and filter data using the pre-filtered params_df
data <- params_df %>%
  filter(Pathogen_Name %in% main_pathogens) %>%
  left_join(routes_df, by = "Pathogen_Name")


# --- 5. Helper Functions ---

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
      sd_val <- (upper - lower) / (2 * qnorm(0.975))
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

# Filter out rows where variance could not be calculated
data <- data %>%
  filter(!is.na(variance))

# Helper function to sample a single study with uniform probability
sample_study <- function(df_group) {
  if (nrow(df_group) == 0) return(NULL)
  sample_n(df_group, 1, replace = TRUE)
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
    
    sample_and_get_value <- function(param_name) {
      studies <- pathogen_data %>% filter(Parameter == param_name)
      if (nrow(studies) == 0) return(NA)
      
      study <- sample_study(studies)
      if (is.null(study) || is.na(study$ReportedValue) || is.na(study$variance)) return(NA)
      
      val <- rnorm(1, mean = study$ReportedValue, sd = sqrt(study$variance))
      return(max(0, val))
    }

    r0_sampled <- sample_and_get_value("R0")
    si_clust_sampled <- sample_and_get_value("SI")
    cfr_sampled <- sample_and_get_value("CFR")
    k_sampled <- sample_and_get_value("k")
    lp_sampled <- sample_and_get_value("LP")
    infp_sampled <- sample_and_get_value("InfP")
    ip_clust_sampled <- sample_and_get_value("IP")

    if(!is.na(cfr_sampled)) cfr_sampled <- pmax(0, pmin(1, cfr_sampled))
    
    pathogen_routes <- pathogen_data[1,] %>% select(starts_with("Route_"))

    pathogen_params_df <- data.frame(
      Pathogen_Name = pathogen,
      R0_sampled = r0_sampled,
      SI_Clust_sampled = si_clust_sampled,
      CFR_sampled = cfr_sampled,
      k_sampled = k_sampled,
      LatentPeriod_sampled = lp_sampled,
      InfectiousPeriod_sampled = infp_sampled,
      IncubationPeriod_sampled = ip_clust_sampled,
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
dir.create("Clustering/mcmc/Kmeans/S7_2_outputs", showWarnings = FALSE)
write_csv(all_mcmc_samples_df, "Clustering/mcmc/Kmeans/S7_2_outputs/mcmc_parameter_samples.csv")
print("MCMC parameter samples saved to Clustering/mcmc/Kmeans/S7_2_outputs/mcmc_parameter_samples.csv")


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
      k_sampled,
      LatentPeriod_sampled,
      InfectiousPeriod_sampled,
      IncubationPeriod_sampled,
      Route_resp,
      Route_direct,
      Route_sexual,
      Route_animal,
      Route_vector
    )

  # Handle NAs in the sampled parameters
  # Impute with the overall median for that parameter column
  numerical_param_cols <- c("R0_sampled", "SI_Clust_sampled", "CFR_sampled", "k_sampled", "LatentPeriod_sampled", "InfectiousPeriod_sampled", "IncubationPeriod_sampled")
  for(col in numerical_param_cols){
    if(any(is.na(clustering_data_full[[col]]))){
      na_count <- sum(is.na(clustering_data_full[[col]]))
      warning(paste("Column", col, "in all_mcmc_samples_df has", na_count, "NA values. Imputing with overall median before scaling."))
      if(is.numeric(clustering_data_full[[col]])){
        clustering_data_full[[col]][is.na(clustering_data_full[[col]])] <- median(clustering_data_full[[col]], na.rm = TRUE)
      } else {
        warning(paste("Column", col, "is not numeric or became all NAs, cannot impute with median."))
      }
    }
  }
  
  # Scale numerical features across the entire dataset
  scaled_clustering_data_full <- clustering_data_full
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

          original_data_this_iter_for_centroids <- all_mcmc_samples_df %>%
            filter(MCMC_Iteration == iter_val) %>%
            mutate(Cluster_Assigned_Iter = kmeans_iter_result$cluster)
          
          iter_unscaled_centroids <- original_data_this_iter_for_centroids %>%
            group_by(Cluster_Assigned_Iter) %>%
            summarise(
              R0_sampled = mean(R0_sampled, na.rm = TRUE),
              SI_Clust_sampled = mean(SI_Clust_sampled, na.rm = TRUE),
              CFR_sampled = mean(CFR_sampled, na.rm = TRUE),
              k_sampled = mean(k_sampled, na.rm = TRUE),
              LatentPeriod_sampled = mean(LatentPeriod_sampled, na.rm = TRUE),
              InfectiousPeriod_sampled = mean(InfectiousPeriod_sampled, na.rm = TRUE),
              IncubationPeriod_sampled = mean(IncubationPeriod_sampled, na.rm = TRUE),
              Route_resp = mean(Route_resp, na.rm=TRUE),
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
          write_csv(centroid_summary, "Clustering/mcmc/Kmeans/S7_2_outputs/cluster_centroids_summary_with_ci_k6.csv")
          print("Cluster centroids summary saved to Clustering/mcmc/Kmeans/S7_2_outputs/cluster_centroids_summary_with_ci_k6.csv")
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
          write_csv(modal_assignments, "Clustering/mcmc/Kmeans/S7_2_outputs/pathogen_modal_cluster_assignments_k6.csv")
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
              k_mean_overall = mean(k_sampled, na.rm = TRUE),
              LatentPeriod_mean_overall = mean(LatentPeriod_sampled, na.rm = TRUE),
              InfectiousPeriod_mean_overall = mean(InfectiousPeriod_sampled, na.rm = TRUE),
              IncubationPeriod_mean_overall = mean(IncubationPeriod_sampled, na.rm = TRUE),
              Route_resp_mean = mean(Route_resp, na.rm=TRUE),
              Route_direct_mean = mean(Route_direct, na.rm=TRUE),
              Route_sexual_mean = mean(Route_sexual, na.rm=TRUE),
              Route_animal_mean = mean(Route_animal, na.rm=TRUE),
              Route_vector_mean = mean(Route_vector, na.rm=TRUE),
              .groups = 'drop'
          ) %>%
          left_join(modal_assignments, by = "Pathogen_Name")

      plot_features_numerical <- c("R0_mean_overall", "SI_Clust_mean_overall", "CFR_mean_overall", "k_mean_overall", "LatentPeriod_mean_overall", "InfectiousPeriod_mean_overall", "IncubationPeriod_mean_overall")
      plot_features_routes <- c("Route_resp_mean", "Route_direct_mean", "Route_sexual_mean", "Route_animal_mean", "Route_vector_mean")

      features_for_pca_plot <- pathogen_summary_for_plot %>%
          select(all_of(plot_features_numerical), all_of(plot_features_routes))
      
      scaled_features_for_pca_plot <- features_for_pca_plot
      
      all_na_cols <- names(which(colSums(is.na(scaled_features_for_pca_plot)) == nrow(scaled_features_for_pca_plot)))
      if (length(all_na_cols) > 0) {
        warning(paste("Removing columns with all NA values before PCA:", paste(all_na_cols, collapse = ", ")))
        scaled_features_for_pca_plot <- scaled_features_for_pca_plot %>% select(-all_of(all_na_cols))
        plot_features_numerical <- plot_features_numerical[!plot_features_numerical %in% all_na_cols]
      }
      
      scaled_features_for_pca_plot[plot_features_numerical] <- scale(scaled_features_for_pca_plot[plot_features_numerical])
      
      scaled_features_for_pca_plot[!is.finite(as.matrix(scaled_features_for_pca_plot))] <- 0

      if (!requireNamespace("ggplot2", quietly = TRUE)) {
          print("Package 'ggplot2' is not installed.")
      } else if (!requireNamespace("factoextra", quietly = TRUE)){
          print("Package 'factoextra' is not installed.")
      } else if (!requireNamespace("ggrepel", quietly = TRUE)){
          print("Package 'ggrepel' is not installed. Please install it: install.packages('ggrepel')")
      } else if (nrow(scaled_features_for_pca_plot) > 1 && !any(is.na(pathogen_summary_for_plot$Modal_Cluster))) {
          
          near_zero_var_cols <- names(which(apply(scaled_features_for_pca_plot, 2, var, na.rm = TRUE) < 1e-10))
          if (length(near_zero_var_cols) > 0) {
            warning(paste("Removing columns with near-zero variance after scaling:", paste(near_zero_var_cols, collapse = ", ")))
            scaled_features_for_pca_plot <- scaled_features_for_pca_plot %>% select(-all_of(near_zero_var_cols))
          }
          
          if (ncol(scaled_features_for_pca_plot) < 2) {
            print("Skipping PCA plot because there are fewer than 2 columns with variance.")
          } else {
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
              ggsave(paste0("Clustering/mcmc/Kmeans/S7_2_outputs/pca_modal_cluster_plot_k",CHOSEN_K,".png"), plot = pca_plot_modal, width=10, height=8)
              print(paste0("PCA plot with modal cluster assignments saved to Clustering/mcmc/Kmeans/S7_2_outputs/pca_modal_cluster_plot_k",CHOSEN_K,".png"))
              
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
                            ggsave(paste0("Clustering/mcmc/Kmeans/S7_2_outputs/silhouette_plot_k",CHOSEN_K,"_iter",example_iter_val,".png"), plot = sil_plot_example_iter, width=8, height=6)
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
  library(stats) 
  library(dendextend)

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
            coassignment_matrix[pathogen2_name, pathogen1_name] <- coassignment_matrix[pathogen1_name, pathogen2_name]
          }
        }
      }
    }
  }
  diag(coassignment_matrix) <- num_mcmc_iterations 
  
  print("Co-assignment matrix (first 6x6):")
  print(head(coassignment_matrix[,1:min(6, ncol(coassignment_matrix))]))
  write.csv(coassignment_matrix, "Clustering/mcmc/Kmeans/S7_2_outputs/coassignment_matrix.csv")

  # --- 9.2. Calculate Dissimilarity Matrix ---
  print("Calculating dissimilarity matrix...")
  dissimilarity_matrix <- 1 - (coassignment_matrix / N_MCMC_ITERATIONS) 
  diag(dissimilarity_matrix) <- 0 

  print("Dissimilarity matrix (first 6x6):")
  print(head(dissimilarity_matrix[,1:min(6, ncol(dissimilarity_matrix))]))
  write.csv(dissimilarity_matrix, "Clustering/mcmc/Kmeans/S7_2_outputs/dissimilarity_matrix.csv")

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
    
    png_filename <- paste0("Clustering/mcmc/Kmeans/S7_2_outputs/figure.S7_2_k", K_CONSENSUS, ".png")
    png(png_filename, width=1200, height=800, units="px", res=100)
    par(mar = c(5, 4, 4, 10))
    plot(dend_colored, horiz = TRUE, 
         main = paste("Consensus Clustering (K=", K_CONSENSUS, ")"), 
         xlab = "Dissimilarity (1 - Proportion Co-assigned)")
  
    dev.off()
    print(paste("Colored consensus dendrogram saved to", png_filename))
  
    # --- 9.5. Extract Consensus Cluster Assignments ---
    print(paste("Cutting tree to get", K_CONSENSUS, "consensus clusters..."))
    consensus_clusters <- cutree(hierarchical_clust, k = K_CONSENSUS)
    consensus_assignments_df <- data.frame(Pathogen_Name = names(consensus_clusters),
                                           Consensus_Cluster = consensus_clusters)
    print("Consensus cluster assignments:")
    print(consensus_assignments_df)
    write.csv(consensus_assignments_df, paste0("Clustering/mcmc/Kmeans/S7_2_outputs/consensus_cluster_assignments_k", K_CONSENSUS, ".csv"))
  
    # --- 9.6. Characterize Consensus Clusters ---
    print("Characterizing consensus clusters (means and 95% CIs from all_mcmc_samples_df)...")
    
    params_to_summarize <- c("R0_sampled", "SI_Clust_sampled", "CFR_sampled",
                             "k_sampled", "LatentPeriod_sampled", "InfectiousPeriod_sampled", "IncubationPeriod_sampled",
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
        write_csv(final_consensus_summary_df, paste0("Clustering/mcmc/Kmeans/S7_2_outputs/consensus_clusters_summary_k", K_CONSENSUS, ".csv"))
        print(paste0("Consensus clusters summary saved to Clustering/mcmc/Kmeans/S7_2_outputs/consensus_clusters_summary_k", K_CONSENSUS, ".csv"))
      } else {
        print("No parameters were summarized for consensus clusters.")
      }
    } else {
      print("Could not prepare data for characterizing consensus clusters (pathogen_mcmc_samples_with_consensus_clusters is empty or missing Consensus_Cluster column).")
    }
    
    print(paste("--- Finished K =", K_CONSENSUS, "---"))
  }
  
  print("Ensemble/Consensus Clustering Analysis finished.")

} else {
  print("Skipping Ensemble/Consensus Clustering: combined_assignments_df or all_mcmc_samples_df not found or empty.")
}
