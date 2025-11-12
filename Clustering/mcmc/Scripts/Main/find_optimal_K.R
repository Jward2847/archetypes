# Script to find the optimal number of clusters (K) using the silhouette method

# --- 1. Load necessary libraries ---
# Ensure packages are installed. The main script might not have installed 'cluster'.
if (!requireNamespace("cluster", quietly = TRUE)) {
  install.packages("cluster")
}
if (!requireNamespace("ggplot2", quietly = TRUE)) {
  install.packages("ggplot2")
}
library(cluster)
library(ggplot2)

print("--- Starting Optimal K Analysis (Silhouette Method) ---")

# --- 2. Configuration & Load Data ---
# Define base directories for inputs and outputs
# This makes it easier to adapt for sensitivity analyses (e.g., S6_outputs)
main_output_dir <- "Clustering/mcmc/Kmeans/main_outputs"
figures_dir <- "Clustering/mcmc/Kmeans/figures"

# This matrix is the output from the consensus clustering step in the main analysis script.
dissimilarity_matrix_file <- file.path(main_output_dir, "dissimilarity_matrix.csv")
if (!file.exists(dissimilarity_matrix_file)) {
  stop(paste("Error: Dissimilarity matrix file not found at", dissimilarity_matrix_file, ". Please run the main mcmc_Kmeans_analysis.R script first."))
}

# Load the dissimilarity matrix.
# Using read.csv with row.names = 1 is robust for files created with write.csv.
dissim_mat <- as.matrix(read.csv(dissimilarity_matrix_file, row.names = 1))

# Convert to a 'dist' object for analysis
pathogen_dist <- as.dist(dissim_mat)

print("Successfully loaded the dissimilarity matrix.")


# --- 3. Perform Hierarchical Clustering ---
# This step is the same as in the main analysis script. We need the clustering tree to cut it at different K.
# Using 'average' linkage to be consistent with the main analysis.
hierarchical_clust <- hclust(pathogen_dist, method = "average")

print("Performed hierarchical clustering.")


# --- 4. Calculate Average Silhouette Width for a Range of K values ---
k_range <- 2:10 # Range of K to test
avg_silhouette_widths <- numeric(length(k_range))
names(avg_silhouette_widths) <- k_range

print(paste("Calculating average silhouette width for K from", min(k_range), "to", max(k_range), "..."))

for (k in k_range) {
  # Cut the tree to get k clusters
  cluster_assignments <- cutree(hierarchical_clust, k = k)
  
  # Calculate silhouette information
  # The silhouette function can take the original distance matrix directly.
  silhouette_info <- silhouette(cluster_assignments, dmatrix = dissim_mat)
  
  # Calculate and store the average silhouette width
  # We check if silhouette_info is valid before summarizing
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


# --- 5. Plot the Results ---
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
    y = "Average Silhouette Width"
  ) +
  theme_minimal(base_size = 14) +
  theme(plot.title = element_text(face = "bold"))

print(optimal_k_plot)

# --- 6. Save the Plot and Results ---
plot_filename <- file.path(figures_dir, "figure.S5.png")
ggsave(plot_filename, plot = optimal_k_plot, width = 8, height = 6)
print(paste("Optimal K plot saved to", plot_filename))

results_filename <- file.path(main_output_dir, "optimal_k_silhouette_results.csv")
write.csv(silhouette_results_df, results_filename, row.names = FALSE)
print(paste("Optimal K analysis results saved to", results_filename))

print("--- Script finished. ---") 