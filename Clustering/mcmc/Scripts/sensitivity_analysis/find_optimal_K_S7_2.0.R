# Script to find the optimal number of clusters (K) using the silhouette method for S7_2.0 Analysis

# --- 1. Load necessary libraries ---
if (!requireNamespace("cluster", quietly = TRUE)) {
  install.packages("cluster")
}
if (!requireNamespace("ggplot2", quietly = TRUE)) {
  install.packages("ggplot2")
}
library(cluster)
library(ggplot2)
library(readr)

print("--- Starting Optimal K Analysis (Silhouette Method) for S7_2.0 ---")

# --- 2. Load the Dissimilarity Matrix ---
dissimilarity_matrix_file <- "Clustering/mcmc/Kmeans/S7_2_outputs/dissimilarity_matrix.csv"
if (!file.exists(dissimilarity_matrix_file)) {
  stop(paste("Error: Dissimilarity matrix file not found at", dissimilarity_matrix_file, ". Please run the supplementary_figure_S7_2.0.R script first."))
}

dissim_df <- read_csv(dissimilarity_matrix_file, col_types = cols(.default = "d", ...1 = "c"))
pathogen_names <- dissim_df$...1
dissim_mat <- as.matrix(dissim_df[,-1])
rownames(dissim_mat) <- pathogen_names

# Convert to a 'dist' object for analysis
pathogen_dist <- as.dist(dissim_mat)

print("Successfully loaded the dissimilarity matrix.")


# --- 3. Perform Hierarchical Clustering ---
hierarchical_clust <- hclust(pathogen_dist, method = "average")

print("Performed hierarchical clustering.")


# --- 4. Calculate Average Silhouette Width for a Range of K values ---
k_range <- 2:9 # Range of K to test (for 10 pathogens)
avg_silhouette_widths <- numeric(length(k_range))
names(avg_silhouette_widths) <- k_range

print(paste("Calculating average silhouette width for K from", min(k_range), "to", max(k_range), "..."))

for (k in k_range) {
  cluster_assignments <- cutree(hierarchical_clust, k = k)
  
  silhouette_info <- silhouette(cluster_assignments, dmatrix = dissim_mat)
  
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
    title = "Optimal Number of Clusters by Silhouette Method",
    subtitle = "For S7_2.0 Analysis",
    x = "Number of Clusters (K)",
    y = "Average Silhouette Width"
  ) +
  theme_minimal(base_size = 14) +
  theme(plot.title = element_text(face = "bold"))

print(optimal_k_plot)

# --- 6. Save the Plot and Results ---
plot_filename <- "Clustering/mcmc/Kmeans/S7_2_outputs/figure_optimal_k.png"
ggsave(plot_filename, plot = optimal_k_plot, width = 8, height = 6)
print(paste("Optimal K plot saved to", plot_filename))

results_filename <- "Clustering/mcmc/Kmeans/S7_2_outputs/optimal_k_silhouette_results.csv"
write_csv(silhouette_results_df, results_filename)
print(paste("Optimal K analysis results saved to", results_filename))

print("--- Script finished. ---")
