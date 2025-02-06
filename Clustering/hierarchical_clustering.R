


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


# Find optimal k (number of clusters)
nb <- NbClust(data[, c("R0", "Serial", "CFR")], diss = gower_dist, distance = NULL, method = "complete")
optimal_k <- nb$Best.nc[1]  # Retrieve optimal k
print(optimal_k) #7 optimal


# Plot with colored clusters
plot(hclust_res, labels = data$Pathogen, main = "Dendrogram of Pathogens")
rect.hclust(hclust_res, k = 6, border = 2:5)  # k is the number of clusters, adjust as needed


##########################

# Convert the hierarchical clustering result to a dendrogram object
dend_data <- as.dendrogram(hclust_res)

# Extract dendrogram data to use in ggplot
dendro_data <- ggdendro::dendro_data(dend_data)

##

# Create clusters by cutting the dendrogram
clusters <- cutree(hclust_res, k = 7)  # Cut into 4 clusters (adjust as needed)

# Add the clusters to the dendro_data$labels dataframe
dendro_data$labels$cluster <- clusters[order.dendrogram(dend_data)]  # Match clusters to labels

# Plot using ggplot2 with cluster coloring
# Assuming 'data$Pathogen' contains the names of pathogens in the same order
# Update the label column in the dendro_data$labels data frame
dendro_data$labels$label <- data$Pathogen[order.dendrogram(dend_data)]

##

########################################
#mini plots 

# K= 5 figure
dendro_5 <- ggplot() +
  geom_segment(data = dendro_data$segments, 
               aes(x = x, y = y, xend = xend, yend = yend), 
               linewidth = 1) +  
  geom_point(data = dendro_data$labels, 
             aes(x = x, y = y, color = factor(cluster)), 
             size = 4) +  
  scale_color_manual(values = c("red", "darkgreen", "orange", "purple", "blue", "darkgrey","#FF69B4"), 
                     name = "Cluster", 
                     labels = c("Cluster 1", "Cluster 2", "Cluster 3", "Cluster 4", "Cluster 5", "Cluster 6","Cluster 7"),
                     guide = guide_legend(override.aes = list(size = 4, shape = 16))) +  
  ylab("Height") +  
  ggtitle("K = 5") +  # Add plot title
  theme_minimal(base_size = 12) +  
  theme(
    axis.title.x = element_blank(),  
    axis.title.y = element_text(size = 18, face = "bold"),  
    axis.text.x = element_blank(),   
    axis.ticks.x = element_blank(),  
    panel.grid.major = element_blank(),  
    panel.grid.minor = element_blank(),  
    legend.position = "right",           
    legend.title = element_text(size = 12, face = "bold"),  
    legend.text = element_text(size = 10, face = "bold"),   
    axis.text.y = element_text(size = 8, face = "bold"),   
    plot.title = element_text(size = 16, face = "bold", hjust = 0.5)  # Centered and bold
  )  

#K=6 figure 
dendro_6 <- ggplot() +
  geom_segment(data = dendro_data$segments, 
               aes(x = x, y = y, xend = xend, yend = yend), 
               linewidth = 1) +  
  geom_point(data = dendro_data$labels, 
             aes(x = x, y = y, color = factor(cluster)), 
             size = 4) +  
  scale_color_manual(values = c("red", "darkgreen", "orange", "purple", "blue", "darkgrey","#FF69B4"), 
                     name = "Cluster", 
                     labels = c("Cluster 1", "Cluster 2", "Cluster 3", "Cluster 4", "Cluster 5", "Cluster 6","Cluster 7"),
                     guide = guide_legend(override.aes = list(size = 4, shape = 16))) +  
  ylab("") +  
  ggtitle("K = 6") +  # Add plot title
  theme_minimal(base_size = 12) +  
  theme(
    axis.title.x = element_blank(),  
    axis.title.y = element_text(size = 18, face = "bold"),  
    axis.text.x = element_blank(),   
    axis.ticks.x = element_blank(),  
    panel.grid.major = element_blank(),  
    panel.grid.minor = element_blank(),  
    legend.position = "right",           
    legend.title = element_text(size = 12, face = "bold"),  
    legend.text = element_text(size = 10, face = "bold"),   
    axis.text.y = element_text(size = 8, face = "bold"),   
    plot.title = element_text(size = 16, face = "bold", hjust = 0.5)  # Centered and bold
  )  

#K=7 figure 
dendro_7 <- ggplot() +
  geom_segment(data = dendro_data$segments, 
               aes(x = x, y = y, xend = xend, yend = yend), 
               linewidth = 1) +  
  geom_point(data = dendro_data$labels, 
             aes(x = x, y = y, color = factor(cluster)), 
             size = 4) +  
  scale_color_manual(values = c("red", "darkgreen", "orange", "purple", "blue", "darkgrey","#FF69B4"), 
                     name = "Cluster", 
                     labels = c("Cluster 1", "Cluster 2", "Cluster 3", "Cluster 4", "Cluster 5", "Cluster 6","Cluster 7"),
                     guide = guide_legend(override.aes = list(size = 4, shape = 16))) +  
  ylab("") +  
  ggtitle("K = 7") +  # Add plot title
  theme_minimal(base_size = 12) +  
  theme(
    axis.title.x = element_blank(),  
    axis.title.y = element_text(size = 18, face = "bold"),  
    axis.text.x = element_blank(),   
    axis.ticks.x = element_blank(),  
    panel.grid.major = element_blank(),  
    panel.grid.minor = element_blank(),  
    legend.position = "right",           
    legend.title = element_text(size = 12, face = "bold"),  
    legend.text = element_text(size = 10, face = "bold"),   
    axis.text.y = element_text(size = 8, face = "bold"),   
    plot.title = element_text(size = 16, face = "bold", hjust = 0.5)  # Centered and bold
  )  

# Combine dendro plots into one panel (a)
panel_a <- plot_grid(dendro_5, dendro_6, dendro_7, ncol = 3, align = "hv")
panel_a_boxed <- ggdraw() +
  draw_plot(panel_a) +
  annotate("rect", xmin = 0, xmax = 1, ymin = 0, ymax = 1, 
           color = "black", size = 1, fill = NA)
########################################

#Main plot
dendro_main <- ggplot() +
  geom_segment(data = dendro_data$segments, 
               aes(x = x, y = y, xend = xend, yend = yend), 
               linewidth = 1) +  
  # Add pathogen labels without affecting the legend
  geom_text(data = dendro_data$labels, 
            aes(x = x, y = y, label = label, color = factor(cluster)), 
            hjust = 1.2, angle = 90, size = 3.5, fontface = "bold", show.legend = FALSE) +  
  #add points for ledgend
  geom_point(data = dendro_data$labels, 
             aes(x = x, y = y, color = factor(cluster)), 
             size = 4) +  
  scale_color_manual(values = c("red", "darkgreen", "orange", "purple", "blue", "darkgrey","#FF69B4"), 
                     name = "Cluster", 
                     labels = c("Cluster 1", "Cluster 2", "Cluster 3", "Cluster 4", "Cluster 5", "Cluster 6","Cluster 7"),
                     guide = guide_legend(override.aes = list(size = 6, shape = 16))) +  
  ylab("Height") +  
  theme_minimal(base_size = 12) +  
  theme(
    axis.title.x = element_blank(),  
    axis.title.y = element_text(size = 18, face = "bold"),  
    axis.text.x = element_blank(),   
    axis.ticks.x = element_blank(),  
    panel.grid.major = element_blank(),  
    panel.grid.minor = element_blank(),  
    legend.position = "right",           
    legend.title = element_text(size = 20, face = "bold"),  
    legend.text = element_text(size = 22, face = "bold"),   
    axis.text.y = element_text(size = 15, face = "bold"),   
    plot.title = element_text(size = 14, face = "bold"),    
    plot.margin = margin(1, 1, 1, 1, "cm"),
    plot.background = element_rect(color = "black", size = 1)
  ) +
  expand_limits(y = -0.5)


# Combine plots
dendro_fig <- plot_grid(
  panel_a_boxed, dendro_main, 
  labels = c("a", "(b)"), 
  label_size = 18, 
  ncol = 1, 
  rel_heights = c(1, 1.5) # Give more space to main_clust
)
dendro_fig

ggsave("Clustering/dendro_fig.png", dendro_fig, width = 14, height = 10, dpi = 300, bg = "white")


#####################################################################################################
#Mean cluster values
# Add the cluster assignments to the original data
data$cluster <- clusters

# Calculate the mean of each parameter (R0, CFR, etc.) by cluster
mean_params_by_cluster <- aggregate(. ~ cluster, data = data[, c("R0", "Serial", "CFR","pre_sym", "cluster")], FUN = mean)

# Display the mean parameters for each cluster
print(mean_params_by_cluster)


