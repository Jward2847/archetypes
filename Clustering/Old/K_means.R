library(dplyr)
library(readxl)
library(ggplot2)
library(ggrepel)
library(factoextra)



rm(list=ls())
setwd("/Users/lshjw6/Documents/blueprint/epi_review/clustering")

################################################################
#R0, CFR and serial interval 

#load data
R0_CFR_SI <- read_excel("cluster_dat.xlsx", 
                          sheet = "R0_CFR_SI", col_types = c("text", 
                                                             "numeric", "numeric", "numeric"))

# Rename specific pathogens in the dataset
R0_CFR_SI <- R0_CFR_SI %>%
  mutate(Pathogen = case_when(
    Pathogen == "COVID-19_WT" ~ "COVID-19 Wild Type",
    Pathogen == "COVID-19_A" ~ "COVID-19 Alpha",
    Pathogen == "COVID-19_D" ~ "COVID-19 Delta",
    Pathogen == "COVID-19_O" ~ "COVID-19 Omicron", 
    Pathogen == "H1N1_18" ~ "H1N1", 
    Pathogen == "H1N1_09" ~ "H1N1pdm09",
    TRUE ~ Pathogen  # Keep other names unchanged
  ))

# Remove the non-numeric column 'pathogen'
data_numeric <- R0_CFR_SI[, -1]

# Scale the data (optional but recommended)
data_scaled <- scale(data_numeric)

#Determine the Optimal Number of Clusters
# Compute and plot the WCSS (Within-cluster Sum of Squares)
fviz_nbclust(data_scaled, kmeans, method = "wss")


#Apply K-means Clustering
# Set the number of clusters
k <- 4

# Apply K-means clustering
set.seed(123)  # For reproducibility
km_res <- kmeans(data_scaled, centers = k, nstart = 25)

# View the clustering results
km_res$cluster

#Visualize the Clusters
fviz_cluster(km_res, data = data_scaled)

#Add Cluster Information to the Original Data
R0_CFR_SI$cluster <- km_res$cluster

# View the updated data
head(R0_CFR_SI)

# View data grouped by cluster
aggregate(data_numeric, by = list(cluster = R0_CFR_SI$cluster), mean)

#visulasiaton 
# Create a data frame with the original pathogen names, scaled data, and clusters
data_with_clusters <- data.frame(Pathogen = R0_CFR_SI$Pathogen, data_scaled, cluster = km_res$cluster)

# Visualize the clusters and label them with pathogen names
RO_CFR_SI_plot <- fviz_cluster(km_res, data = data_scaled, geom = "point", labelsize = 3) +
  geom_text_repel(aes(label = R0_CFR_SI$Pathogen), vjust = -0.5, hjust = 0.5) +
  theme_classic()

fviz_cluster(km_res, data = data_scaled, geom = "point", labelsize = 4, 
             pointsize = 3, ellipse.alpha = 0.3) +  # Adjust point and ellipse transparency
  geom_text_repel(aes(label = R0_CFR_SI$Pathogen), 
                  size = 4, fontface = "bold", color = "black", 
                  box.padding = 0.35, point.padding = 0.5, 
                  segment.color = "grey50") +  # Customize text size and repulsion
  labs(title = "R0, CFR and serial interval", 
       x = "Dimension 1 (69.5%)", y = "Dimension 2 (21.4%)") +  # Add titles and axis labels
  theme_classic(base_size = 15) +  # Base font size for theme
  theme(legend.position = "bottom") +  # Move legend to bottom
  scale_color_manual(values = c("#E41A1C", "#377EB8", "#4DAF4A", "#984EA3"))  # Custom colors for clusters

################################################################

#R0, CFR 

#load data
R0_CFR <- read_excel("cluster_dat.xlsx", 
                        sheet = "R0_CFR", col_types = c("text", 
                                                           "numeric", "numeric"))

# Rename specific pathogens in the dataset
R0_CFR <- R0_CFR %>%
  mutate(Pathogen = case_when(
    Pathogen == "COVID-19_WT" ~ "COVID-19 Wild Type",
    Pathogen == "COVID-19_A" ~ "COVID-19 Alpha",
    Pathogen == "COVID-19_D" ~ "COVID-19 Delta",
    Pathogen == "COVID-19_O" ~ "COVID-19 Omicron", 
    Pathogen == "H1N1_18" ~ "H1N1", 
    Pathogen == "H1N1_09" ~ "H1N1pdm09",
    TRUE ~ Pathogen  # Keep other names unchanged
  ))

# Remove the non-numeric column 'pathogen'
data_numeric_2 <- R0_CFR[, -1]

# Scale the data (optional but recommended)
data_scaled_2 <- scale(data_numeric_2)

#Determine the Optimal Number of Clusters
# Compute and plot the WCSS (Within-cluster Sum of Squares)
fviz_nbclust(data_scaled_2, kmeans, method = "wss")


#Apply K-means Clustering
# Set the number of clusters
k <- 4

# Apply K-means clustering
set.seed(123)  # For reproducibility
km_res_2 <- kmeans(data_scaled_2, centers = k, nstart = 25)

# View the clustering results
km_res_2$cluster

#Visualize the Clusters
fviz_cluster(km_res_2, data = data_scaled_2)

#Add Cluster Information to the Original Data
R0_CFR$cluster <- km_res_2$cluster

# View the updated data
head(R0_CFR)

# View data grouped by cluster
aggregate(data_numeric_2, by = list(cluster = R0_CFR$cluster), mean)

#visulasiaton 
# Create a data frame with the original pathogen names, scaled data, and clusters
data_with_clusters <- data.frame(Pathogen = R0_CFR$Pathogen, data_scaled_2, cluster = km_res_2$cluster)

# Visualize the clusters and label  with pathogen names
RO_CFR_plot <- fviz_cluster(km_res, data = data_scaled_2, geom = "point", labelsize = 3) +
  geom_text_repel(aes(label = R0_CFR$Pathogen), vjust = -0.5, hjust = 0.5) +
  theme_classic()

fviz_cluster(km_res_2, data = data_scaled_2, geom = "point", labelsize = 4, 
             pointsize = 3, ellipse.alpha = 0.3) +  # Adjust point and ellipse transparency
  geom_text_repel(aes(label = R0_CFR$Pathogen), 
                  size = 4, fontface = "bold", color = "black", 
                  box.padding = 0.35, point.padding = 0.5, 
                  segment.color = "grey50") +  # Customize text size and repulsion
  labs(title = "R0 and CFR") +  # Add titles and axis labels
  theme_classic(base_size = 15) +  # Base font size for theme
  theme(legend.position = "bottom") +  # Move legend to bottom
  scale_color_manual(values = c("#E41A1C", "#377EB8", "#4DAF4A", "#984EA3"))  # Custom colors for clusters


################################################################
#R0, CFR, inc and serial interval 

#load data
inc <- read_excel("cluster_dat.xlsx", 
                        sheet = "inc", col_types = c("text", 
                                                           "numeric", "numeric", "numeric","numeric"))

# Rename specific pathogens in the dataset
inc <- inc %>%
  mutate(Pathogen = case_when(
    Pathogen == "COVID-19_WT" ~ "COVID-19 Wild Type",
    Pathogen == "COVID-19_A" ~ "COVID-19 Alpha",
    Pathogen == "COVID-19_D" ~ "COVID-19 Delta",
    Pathogen == "COVID-19_O" ~ "COVID-19 Omicron", 
    Pathogen == "H1N1_18" ~ "H1N1", 
    Pathogen == "H1N1_09" ~ "H1N1pdm09",
    TRUE ~ Pathogen  # Keep other names unchanged
  ))

# Remove the non-numeric column 'pathogen'
data_numeric_3 <- inc[, -1]

# Scale the data (optional but recommended)
data_scaled_3 <- scale(data_numeric_3)

#Determine the Optimal Number of Clusters
# Compute and plot the WCSS (Within-cluster Sum of Squares)
fviz_nbclust(data_scaled_3, kmeans, method = "wss")


#Apply K-means Clustering
# Set the number of clusters
k <- 4

# Apply K-means clustering
set.seed(123)  # For reproducibility
km_res_3 <- kmeans(data_scaled_3, centers = k, nstart = 25)

# View the clustering results
km_res_3$cluster

#Visualize the Clusters
fviz_cluster(km_res_3, data = data_scaled_3)

#Add Cluster Information to the Original Data
inc$cluster <- km_res_3$cluster

# View the updated data
head(inc)

# View data grouped by cluster
aggregate(data_numeric_3, by = list(cluster = inc$cluster), mean)

#visulasiaton 
# Create a data frame with the original pathogen names, scaled data, and clusters
data_with_clusters <- data.frame(Pathogen = inc$Pathogen, data_scaled_3, cluster = km_res_3$cluster)

# Visualize the clusters and label them with pathogen names
inc_plot <- fviz_cluster(km_res_3, data = data_scaled_3, geom = "point", labelsize = 3) +
  geom_text_repel(aes(label = inc$Pathogen), vjust = -0.5, hjust = 0.5) +
  theme_classic()

fviz_cluster(km_res_3, data = data_scaled_3, geom = "point", labelsize = 4, 
             pointsize = 3, ellipse.alpha = 0.3) +  # Adjust point and ellipse transparency
  geom_text_repel(aes(label = inc$Pathogen), 
                  size = 4, fontface = "bold", color = "black", 
                  box.padding = 0.35, point.padding = 0.5, 
                  segment.color = "grey50") +  # Customize text size and repulsion
  labs(title = "R0, incubation period, serial interval and CFR") +  # Add titles and axis labels
  theme_classic(base_size = 15) +  # Base font size for theme
  theme(legend.position = "bottom") +  # Move legend to bottom
  scale_color_manual(values = c("#E41A1C", "#377EB8", "#4DAF4A", "#984EA3"))  # Custom colors for clusters

########################################################################################

#CFR R0 and route of transmssion

setwd("/Users/lshjw6/Documents/blueprint/epi_review/clustering")

rm(list=ls())

#load data
route <- read_excel("cluster_dat.xlsx", 
                       sheet = "route", col_types = c("text", 
                                                      "numeric", "numeric", "numeric", 
                                                      "numeric", "numeric", "numeric", 
                                                      "numeric"))

route <- route %>%
  mutate(Pathogen = case_when(
    Pathogen == "COVID-19_WT" ~ "COVID-19 Wild Type",
    Pathogen == "COVID-19_A" ~ "COVID-19 Alpha",
    Pathogen == "COVID-19_D" ~ "COVID-19 Delta",
    Pathogen == "COVID-19_O" ~ "COVID-19 Omicron", 
    Pathogen == "H1N1_18" ~ "H1N1", 
    Pathogen == "H1N1_09" ~ "H1N1pdm09",
    TRUE ~ Pathogen  # Keep other names unchanged
  ))



# Remove the non-numeric column 'pathogen'
data_numeric_4 <- route[, -1]

# Scale the data (optional but recommended)
data_scaled_4 <- scale(data_numeric_4)



#Determine the Optimal Number of Clusters
# Compute and plot the WCSS (Within-cluster Sum of Squares)
fviz_nbclust(data_scaled_4, kmeans, method = "wss")


#Apply K-means Clustering
# Set the number of clusters
k <- 5

# Apply K-means clustering
set.seed(123)  # For reproducibility
km_res_4 <- kmeans(data_scaled_4, centers = k, nstart = 25)

# View the clustering results
km_res_4$cluster

#Visualize the Clusters
fviz_cluster(km_res_4, data = data_scaled_4)

#Add Cluster Information to the Original Data
route$cluster <- km_res_4$cluster

# View the updated data
head(route)

# View data grouped by cluster
aggregate(data_numeric_4, by = list(cluster = route$cluster), mean)

#visulasiaton 
# Create a data frame with the original pathogen names, scaled data, and clusters
data_with_clusters <- data.frame(Pathogen = route$Pathogen, data_scaled_4, cluster = km_res_4$cluster)

# Visualize the clusters and label them with pathogen names
inc_plot <- fviz_cluster(km_res_4, data = data_scaled_4, geom = "point", labelsize = 3) +
  geom_text_repel(aes(label = route$Pathogen), vjust = -0.5, hjust = 0.5) +
  theme_classic()

fviz_cluster(km_res_4, data = data_scaled_4, geom = "point", labelsize = 4, 
             pointsize = 3, ellipse.alpha = 0.3) +  # Adjust point and ellipse transparency
  geom_text_repel(aes(label = route$Pathogen), 
                  size = 4, fontface = "bold", color = "black", 
                  box.padding = 0.5, point.padding = 1.0, 
                  segment.color = "grey50", 
                  max.overlaps = Inf) +  # Customize text size and repulsion
  labs(title = "R0, CFR, transmission route") +  # Add titles and axis labels
  theme_classic(base_size = 15) +  # Base font size for theme
  theme(legend.position = "bottom") +  # Move legend to bottom
  scale_color_manual(values = c("#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "black"))  # Custom colors for clusters

########################################################################################
#R0, CFR, serial interval, route of transmisson 

setwd("/Users/lshjw6/Documents/blueprint/epi_review/clustering")

rm(list=ls())

#load data
route_SI <- read_excel("cluster_dat.xlsx", 
                    sheet = "route_SI", col_types = c("text", 
                                                   "numeric","numeric", "numeric", "numeric", 
                                                   "numeric", "numeric", "numeric", 
                                                   "numeric"))

route_SI <- route_SI %>%
  mutate(Pathogen = case_when(
    Pathogen == "COVID-19_WT" ~ "COVID-19 Wild Type",
    Pathogen == "COVID-19_A" ~ "COVID-19 Alpha",
    Pathogen == "COVID-19_D" ~ "COVID-19 Delta",
    Pathogen == "COVID-19_O" ~ "COVID-19 Omicron", 
    Pathogen == "H1N1_18" ~ "H1N1", 
    Pathogen == "H1N1_09" ~ "H1N1pdm09",
    TRUE ~ Pathogen  # Keep other names unchanged
  ))



# Remove the non-numeric column 'pathogen'
data_numeric_5 <- route_SI[, -1]

# Scale the data (optional but recommended)
data_scaled_5 <- scale(data_numeric_5)



#Determine the Optimal Number of Clusters
# Compute and plot the WCSS (Within-cluster Sum of Squares)
fviz_nbclust(data_scaled_5, kmeans, method = "wss")


#Apply K-means Clustering
# Set the number of clusters
k <- 4

# Apply K-means clustering
set.seed(123)  # For reproducibility
km_res_5 <- kmeans(data_scaled_5, centers = k, nstart = 25)

# View the clustering results
km_res_5$cluster

#Visualize the Clusters
fviz_cluster(km_res_5, data = data_scaled_5)

#Add Cluster Information to the Original Data
route_SI$cluster <- km_res_5$cluster

# View the updated data
head(route_SI)

# View data grouped by cluster
aggregate(data_numeric_5, by = list(cluster = route_SI$cluster), mean)

#visulasiaton 
# Create a data frame with the original pathogen names, scaled data, and clusters
data_with_clusters <- data.frame(Pathogen = route_SI$Pathogen, data_scaled_5, cluster = km_res_5$cluster)

# Visualize the clusters and label them with pathogen names
inc_plot <- fviz_cluster(km_res_5, data = data_scaled_5, geom = "point", labelsize = 3) +
  geom_text_repel(aes(label = route_SI$Pathogen), vjust = -0.5, hjust = 0.5) +
  theme_classic()

fviz_cluster(km_res_5, data = data_scaled_5, geom = "point", labelsize = 4, 
             pointsize = 3, ellipse.alpha = 0.3) +  # Adjust point and ellipse transparency
  geom_text_repel(aes(label = route_SI$Pathogen), 
                  size = 4, fontface = "bold", color = "black", 
                  box.padding = 0.5, point.padding = 1.0, 
                  segment.color = "grey50", 
                  max.overlaps = Inf) +  # Customize text size and repulsion
  labs(title = "") +  # Add titles and axis labels
  theme_classic(base_size = 15) +  # Base font size for theme
  theme(legend.position = "bottom") +  # Move legend to bottom
  scale_color_manual(values = c("#E41A1C", "#377EB8", "#4DAF4A", "#984EA3","black"))  # Custom colors for clusters

########
#Extract Cluster Centroids

# Extract the cluster centroids from the K-means result
cluster_centroids <- km_res_5$centers

# View the centroids
cluster_centroids


# Add the cluster assignment to the original data
route_SI$cluster <- km_res_5$cluster

# Calculate the mean of each feature (R0, CFR, SI) for each cluster
cluster_summary <- aggregate(. ~ cluster, data = route_SI[, c("R0", "CFR", "Serial",
                                                              "Route_resp", "Route_direct", "Route_sexual",
                                                              "Route_animal", "Route_vector", "cluster")], mean)

# View the cluster summaries
cluster_summary
