library(ggdendro)
library(dplyr)
library(readxl)
library(ggplot2)
library(ggrepel)
library(cluster)
library(clustMixType)
library(patchwork)
library(factoextra)
library(cowplot)

rm(list = ls())

route_SI <- read_excel("Clustering/cluster_dat.xlsx", 
                       sheet = "presym", col_types = c("text", 
                                                         "numeric","numeric", "numeric", "numeric", 
                                                         "numeric", "numeric", "numeric", "numeric",
                                                         "numeric"))

route_SI <- route_SI %>%
  mutate(Pathogen = case_when(
    Pathogen == "COVID-19_WT" ~ "SARS-CoV-2 (WT)",
    Pathogen == "COVID-19_A" ~ "SARS-CoV-2 (Alpha)",
    Pathogen == "COVID-19_D" ~ "SARS-CoV-2 (Delta)",
    Pathogen == "COVID-19_O" ~ "SARS-CoV-2 (Omicron)", 
    Pathogen == "H1N1_18" ~ "A/H1N1",
    Pathogen == "H2N2" ~ "A/H2N2",
    Pathogen == "H3N2" ~ "A/H3N2",
    Pathogen == "H1N1_09" ~ "A/H1N1/09",
    Pathogen == "H5N1" ~ "A/H5N1",
    Pathogen == "Ebola" ~ "EBOV",
    Pathogen == "Marburg" ~ "MARV",
    Pathogen == "Mpox" ~ "MPV",
    Pathogen == "Lassa" ~ "LASV",
    Pathogen == "Nipah" ~ "NiV",
    Pathogen == "Zika" ~ "ZIKV",
    Pathogen == "SARS" ~ "SARS-CoV-1",
    Pathogen == "MERS" ~ "MERS-CoV",
    Pathogen == "CCHF" ~ "CCHFV",
    TRUE ~ Pathogen  # Keep other names unchanged
  ))



# Remove the non-numeric column 'pathogen'
data_numeric_5 <- route_SI[, -1]

# Scale the data (optional but recommended)
data_scaled_5 <- scale(data_numeric_5)



#Determine the Optimal Number of Clusters
# Compute and plot the WCSS (Within-cluster Sum of Squares)
fviz_nbclust(data_scaled_5, kmeans, method = "wss")

# Determine optimal clusters using the Silhouette method
fviz_nbclust(data_scaled_5, kmeans, method = "silhouette")


# Create the plots
# Create the plots with increased text size
plot_wss <- fviz_nbclust(data_scaled_5, kmeans, method = "wss") + 
  theme(
    plot.title = element_text(size = 18, face = "bold"),
    axis.title = element_text(size = 14),
    axis.text = element_text(size = 12),
    legend.title = element_text(size = 14),
    legend.text = element_text(size = 12),
    plot.background = element_rect(color = "black", size = 1)
  )

plot_silhouette <- fviz_nbclust(data_scaled_5, kmeans, method = "silhouette") + 
  theme(
    plot.title = element_text(size = 18, face = "bold"),
    axis.title = element_text(size = 14),
    axis.text = element_text(size = 12),
    legend.title = element_text(size = 14),
    legend.text = element_text(size = 12)
  )

# Combine the plots
combined_plot <- plot_grid(plot_wss, plot_silhouette, labels = c("A", "B"), label_size = 16, align = "h")

# Display the combined plot
print(combined_plot)



# Compute the Gap statistics
set.seed(123)
gap_stat <- clusGap(data_scaled_5, FUN = kmeans, nstart = 25, K.max = 10, B = 50)

# Visualize the Gap statistics
fviz_gap_stat(gap_stat)






# Apply K-means clustering
set.seed(123)  # For reproducibility
km_res_4 <- kmeans(data_scaled_5, centers = 4, nstart = 25) 
km_res_5 <- kmeans(data_scaled_5, centers = 5, nstart = 25)
km_res_6 <- kmeans(data_scaled_5, centers = 6, nstart = 25)

# View the clustering results
km_res_5$cluster

#Visualise the Clusters
#K-means figure B
clust_4 <- fviz_cluster(km_res_4, data = data_scaled_5, labelsize = 0) +  
  labs(title = "K = 4") +  
  theme(plot.title = element_text(size = 18, face = "bold"),
        legend.text = element_text(size = 12), 
        legend.title = element_text(size = 13, face = "bold"), 
        axis.text = element_text(size = 12, face = "bold"), 
        axis.title = element_text(size = 14, face = "bold"))

clust_5 <- fviz_cluster(km_res_5, data = data_scaled_5, labelsize = 0) +  
  labs(title = "K = 5") +  
  theme(plot.title = element_text(size = 18, face = "bold"),
        legend.text = element_text(size = 12), 
        legend.title = element_text(size = 13, face = "bold"), 
        axis.text = element_text(size = 12, face = "bold"), 
        axis.title = element_text(size = 14, face = "bold"))

clust_6 <- fviz_cluster(km_res_6, data = data_scaled_5, labelsize = 0) +  
  labs(title = "K = 6") +  
  theme(plot.title = element_text(size = 18, face = "bold"),
        legend.text = element_text(size = 12), 
        legend.title = element_text(size = 13, face = "bold"), 
        axis.text = element_text(size = 12, face = "bold"), 
        axis.title = element_text(size = 14, face = "bold"))

#Add Cluster Information to the Original Data
route_SI$cluster <- km_res_5$cluster

# View the updated data
head(route_SI)

# View data grouped by cluster
aggregate(data_numeric_5, by = list(cluster = route_SI$cluster), mean)

#visulasiaton 
# Create a data frame with the original pathogen names, scaled data, and clusters
data_with_clusters <- data.frame(Pathogen = route_SI$Pathogen, data_scaled_5, cluster = km_res_5$cluster)

# Visualize the clusters and label them with pathogen names 4 cluster

main_clust <- fviz_cluster(km_res_5, 
             data = data_scaled_5, 
             geom = "point", 
             labelsize = 3, 
             pointsize = 4, 
             ellipse.alpha = 0.3) +  
  geom_label_repel(aes(label = route_SI$Pathogen), 
                   size = 4,  # Increase text size
                   fontface = "bold", 
                   color = "black", 
                   fill = "white",  # Add a white background for readability
                   box.padding = 4,  # Increase padding around labels
                   point.padding = 1,  # Increase spacing from points
                   segment.color = "black",  # Make connecting lines more visible
                   segment.size = 0.5,  # Make the segment line thinner
                   force = 4,  # Increase force to spread labels
                   force_pull = 1,  # Make labels less likely to overlap
                   max.overlaps = Inf, 
                   show.legend = FALSE) +  
  labs(title = "") +
       #color = "Cluster") +  
  theme_classic(base_size = 20) +  
  theme(legend.position = "right",  
        legend.text = element_text(size = 14), 
        legend.title = element_text(size = 16, face = "bold"), 
        axis.text = element_text(size = 18), 
        axis.title = element_text(size = 20, face = "bold"),
        plot.background = element_rect(color = "black", size = 1)) +  
  scale_color_manual(values = c("red", "darkgreen", "orange", "purple", "blue")) +
  scale_fill_manual(values = c("red", "darkgreen", "orange", "purple", "blue"))


#Create figure

# Combine K-means cluster plots into one panel (a)
panel_a <- plot_grid(clust_4, clust_5, clust_6, ncol = 3, align = "hv")

# Wrap panel_a with a single black border
panel_a_boxed <- ggdraw() +
  draw_plot(panel_a) +
  annotate("rect", xmin = 0, xmax = 1, ymin = 0, ymax = 1, 
           color = "black", size = 1, fill = NA)

# Create the top panel with (a) and (b)
top_panel <- plot_grid(
  panel_a_boxed, plot_wss, 
  labels = c("(a)", "(b)"), 
  label_size = 18, 
  ncol = 2, 
  rel_widths = c(2, 1) # Make (a) wider than (b)
)

# Arrange the full figure layout
kmeans_plot <- plot_grid(
  top_panel, main_clust, 
  labels = c("", "(c)"), 
  label_size = 18, 
  ncol = 1, 
  rel_heights = c(1, 1.5) # Give more space to main_clust
)


#Save plot
ggsave("Clustering/kmeans_fig.png", kmeans_plot, width = 14, height = 10)




########
#Extract Cluster Centroids

# Extract the cluster centroids from the K-means result
cluster_centroids <- km_res_5$centers

# View the centroids
cluster_centroids


# Add the cluster assignment to the original data
route_SI$cluster <- km_res_5$cluster

# Calculate the mean of each feature (R0, CFR, SI) for each cluster
cluster_summary <- aggregate(. ~ cluster, data = route_SI[, c("R0", "CFR", "Serial", "pre_sym",
                                                              "Route_resp", "Route_direct", "Route_sexual",
                                                              "Route_animal", "Route_vector", "cluster")], mean)

# View the cluster summaries
cluster_summary


###################################################################


data <- read_excel("Clustering/cluster_dat.xlsx", 
                   sheet = "k-proto_ps", col_types = c("text", 
                                                    "numeric", "numeric","numeric","numeric", "text", "text", 
                                                    "text"))

data <- data %>%
  mutate(Pathogen = case_when(
    Pathogen == "COVID-19_WT" ~ "SARS-CoV-2 (WT)",
    Pathogen == "COVID-19_A" ~ "SARS-CoV-2 (Alpha)",
    Pathogen == "COVID-19_D" ~ "SARS-CoV-2 (Delta)",
    Pathogen == "COVID-19_O" ~ "SARS-CoV-2 (Omicron)", 
    Pathogen == "H1N1_18" ~ "A/H1N1",
    Pathogen == "H2N2" ~ "A/H2N2",
    Pathogen == "H3N2" ~ "A/H3N2",
    Pathogen == "H1N1_09" ~ "A/H1N1/09",
    Pathogen == "H5N1" ~ "A/H5N1",
    Pathogen == "Ebola" ~ "EBOV",
    Pathogen == "Marburg" ~ "MARV",
    Pathogen == "Mpox" ~ "MPV",
    Pathogen == "Lassa" ~ "LASV",
    Pathogen == "Nipah" ~ "NiV",
    Pathogen == "Zika" ~ "ZIKV",
    Pathogen == "SARS" ~ "SARS-CoV-1",
    Pathogen == "MERS" ~ "MERS-CoV",
    Pathogen == "CCHF" ~ "CCHFV",
    TRUE ~ Pathogen  # Keep other names unchanged
  ))



# Convert categorical variables into factors
data$human_human <- as.factor(data$human_human)
data$vector <- as.factor(data$vector)
data$animal <- as.factor(data$animal)

# Calculate Gower distance for mixed data types
gower_dist <- daisy(data[, c("R0", "Serial", "CFR", "human_human", "vector", "animal")], metric = "gower")

# Perform hierarchical clustering
hclust_res <- hclust(gower_dist, method = "complete")

# Plot the dendrogram
plot(hclust_res, labels = data$Pathogen, main = "Dendrogram of Pathogens (Mixed Data)", xlab = "", sub = "", cex = 0.9)


# Plot with colored clusters
plot(hclust_res, labels = data$Pathogen, main = "Dendrogram of Pathogens")
rect.hclust(hclust_res, k = 3, border = 2:5)  # k is the number of clusters, adjust as needed


##########################


# Perform hierarchical clustering
hclust_res <- hclust(gower_dist, method = "complete")

# Convert the hierarchical clustering result to a dendrogram object
dend_data <- as.dendrogram(hclust_res)

# Extract dendrogram data to use in ggplot
dendro_data <- ggdendro::dendro_data(dend_data)

# Plot the dendrogram using ggplot2
ggplot() +
  geom_segment(data = dendro_data$segments, aes(x = x, y = y, xend = xend, yend = yend)) +
  geom_text(data = dendro_data$labels, aes(x = x, y = y, label = label), hjust = 1, size = 3) +
  ggtitle("Dendrogram of Pathogens") +
  theme_minimal() +
  theme(axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank())

# Cut the dendrogram into clusters (adjust the number of clusters or height)
clusters <- cutree(hclust_res, k = 4)  # Cut into 4 clusters, for example

# Create a data frame that matches clusters with the labels
dendro_data$labels$cluster <- clusters[order.dendrogram(dend_data)]

# Plot the dendrogram with clusters colored
# Perform hierarchical clustering
hclust_res <- hclust(gower_dist, method = "complete")

# Convert the dendrogram to a dendro_data object
dend_data <- as.dendrogram(hclust_res)
dendro_data <- ggdendro::dendro_data(dend_data)

# Create clusters by cutting the dendrogram
clusters <- cutree(hclust_res, k = 5)  # Cut into 4 clusters (adjust as needed)

# Add the clusters to the dendro_data$labels dataframe
dendro_data$labels$cluster <- clusters[order.dendrogram(dend_data)]  # Match clusters to labels

# Plot using ggplot2 with cluster coloring
# Assuming 'data$Pathogen' contains the names of pathogens in the same order
# Update the label column in the dendro_data$labels data frame
dendro_data$labels$label <- data$Pathogen[order.dendrogram(dend_data)]


ggplot() +
  # Adjust line width using 'linewidth'
  geom_segment(data = dendro_data$segments, aes(x = x, y = y, xend = xend, yend = yend), 
               linewidth = 1) +  # Use 'linewidth' for line thickness
  # Add text labels with adjusted font size, angle, and color for publication
  geom_text(data = dendro_data$labels, 
            aes(x = x, y = y, label = label, color = factor(cluster)), 
            hjust = 1, angle = 90, size = 7, fontface = "bold") +  # Bold text
  # Customize color scale for cluster color (considering colorblind-friendly palette)
  scale_color_manual(values = c("blue", "red", "darkgrey", "#d4af37","green"), 
                     name = "Cluster", 
                     labels = c("Cluster 1", "Cluster 2", "Cluster 3", "Cluster 4"),
                     guide = guide_legend(override.aes = list(size = 4, shape = 16))) +  # Ensure circles in the legend
  # Add a clear and informative title
  ggtitle("") +
  # Set minimal theme with clearer font and layout
  theme_minimal(base_size = 12) +  # Increase base font size for publication
  theme(
    axis.title.x = element_blank(),  # Remove x-axis title
    axis.title.y = element_blank(),  # Remove y-axis title
    axis.text.x = element_blank(),   # Remove x-axis labels (already labeled by text)
    axis.ticks.x = element_blank(),  # Remove x-axis ticks
    panel.grid.major = element_blank(),  # Remove major grid lines for a cleaner look
    panel.grid.minor = element_blank(),  # Remove minor grid lines
    legend.position = "right",           # Position legend to the right
    legend.title = element_text(size = 20, face = "bold"),  # Bold and larger legend title
    legend.text = element_text(size = 22, face = "bold"),   # Bold and larger legend text
    axis.text.y = element_text(size = 15, face = "bold"),   # Bold y-axis labels
    plot.title = element_text(size = 14, face = "bold"),    # Title font customization
    plot.margin = margin(1, 1, 1, 1, "cm")  # Increase the plot margins to avoid cut-offs
  ) +
  # Expand limits to ensure all text is visible
  expand_limits(y = -0.5)


# Assuming hierarchical clustering has already been performed
# hclust_res contains the result of hierarchical clustering

# Cut the dendrogram into clusters (e.g., 4 clusters)
clusters <- cutree(hclust_res, k = 3)  # Adjust 'k' to your desired number of clusters

# Add the cluster assignments to the original data
data$cluster <- clusters

# Calculate the mean of each parameter (R0, CFR, etc.) by cluster
mean_params_by_cluster <- aggregate(. ~ cluster, data = data[, c("R0", "Serial", "CFR","pre_sym", "cluster")], FUN = mean)

# Display the mean parameters for each cluster
print(mean_params_by_cluster)




