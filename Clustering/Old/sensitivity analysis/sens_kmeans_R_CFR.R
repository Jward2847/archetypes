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

route_SI <- read_excel("Clustering/data/cluster_dat.xlsx", 
                       sheet = "sens_kmeans", col_types = c("text", 
                                                       "numeric","numeric", "numeric", "numeric", 
                                                       "numeric", "numeric", "numeric"))

route_SI <- route_SI %>%
  mutate(Pathogen = case_when(
    Pathogen == "COVID-19_WT" ~ "1",
    Pathogen == "COVID-19_A" ~ "2",
    Pathogen == "COVID-19_D" ~ "3",
    Pathogen == "COVID-19_O" ~ "4", 
    Pathogen == "H1N1_18" ~ "5",
    Pathogen == "H2N2" ~ "6",
    Pathogen == "H3N2" ~ "7",
    Pathogen == "H1N1_09" ~ "8",
    Pathogen == "H5N1" ~ "9",
    Pathogen == "Ebola" ~ "10",
    Pathogen == "Marburg" ~ "11",
    Pathogen == "Mpox" ~ "12",
    Pathogen == "Lassa" ~ "13",
    Pathogen == "Nipah" ~ "14",
    Pathogen == "Zika" ~ "15",
    Pathogen == "SARS" ~ "16",
    Pathogen == "MERS" ~ "17",
    Pathogen == "CCHF" ~ "18",
    Pathogen == "RVF" ~ "19",
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
  labs(x = "Number of Clusters (K)") + 
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
km_res_3 <- kmeans(data_scaled_5, centers = 3, nstart = 25) 
km_res_4 <- kmeans(data_scaled_5, centers = 4, nstart = 25)
km_res_5 <- kmeans(data_scaled_5, centers = 5, nstart = 25)

# View the clustering results
km_res_5$cluster

#Visualise the Clusters
#K-means figure B
clust_3 <- fviz_cluster(km_res_3, data = data_scaled_5, labelsize = 0) +  
  labs(title = "K = 3") +  
  theme(plot.title = element_text(size = 18, face = "bold"),
        legend.text = element_text(size = 12), 
        legend.title = element_text(size = 13, face = "bold"), 
        axis.text = element_text(size = 12, face = "bold"), 
        axis.title = element_text(size = 14, face = "bold"))

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
                           pointsize = 3, 
                           ellipse.alpha = 0.3) +  
  geom_text_repel(aes(label = route_SI$Pathogen), 
                  size = 5,  # text size
                  fontface = "bold", 
                  color = "black", 
                  vjust = -0.5,  # Moves labels above points
                  box.padding = 0.1,  # Less padding to bring labels closer
                  point.padding = 0.1,  # Closer to points
                  segment.color = NA,  # Remove segment lines
                  force = 6,  # Repulsion for spacing
                  force_pull = 1,  # Attraction to points
                  max.iter = 5000,  # Iterations
                  max.overlaps = Inf, 
                  show.legend = FALSE) +  
  labs(title = "") +
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
panel_a <- plot_grid(clust_3, clust_4, clust_5, ncol = 3, align = "hv")

# Wrap panel_a with a single black border
panel_a_boxed <- ggdraw() +
  draw_plot(panel_a) +
  annotate("rect", xmin = 0, xmax = 1, ymin = 0, ymax = 1, 
           color = "black", size = 1, fill = NA)

# Create the top panel with (a) and (b)
top_panel <- plot_grid(
  plot_wss, panel_a_boxed,
  labels = c("(a)", "(b)"), 
  label_size = 18, 
  ncol = 2, 
  rel_widths = c(1, 2) # Make (a) wider than (b)
)

# Arrange the full figure layout
sens_kmeans_plot <- plot_grid(
  top_panel, main_clust, 
  labels = c("", "(c)"), 
  label_size = 18, 
  ncol = 1, 
  rel_heights = c(1, 1.5) # Give more space to main_clust
)

sens_kmeans_plot

#Save plot
ggsave("Clustering/figs/sens_fig_kmeans.png", sens_kmeans_plot, width = 14, height = 10)