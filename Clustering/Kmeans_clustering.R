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
km_res_5 <- kmeans(data_scaled_5, centers = 6, nstart = 25)
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


