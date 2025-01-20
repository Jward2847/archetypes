library(dplyr)

#demo data
data("USArrests")  # Loading the data set
df <- scale(USArrests) # Scaling the data

str(df)
head(df)

setwd("/Users/lshjw6/Documents/blueprint/epi_review/clustering")

library(readxl)
test_dat <- read_excel("test_dat.xlsx", col_types = c("text", 
                                                      "numeric", "numeric", "numeric"))

# Remove the non-numeric column 'pathogen'
data_numeric <- test_dat[, -1]

# Scale the data (optional but recommended)
data_scaled <- scale(data_numeric)

fviz_nbclust(data_scaled, kmeans, method = "wss")

k <- 4

# Apply K-means clustering
set.seed(123)  # For reproducibility
km_res <- kmeans(data_scaled, centers = k, nstart = 25)

# Add cluster and pathogen names to the scaled data
data_scaled_with_labels <- data.frame(data_scaled, cluster = km_res$cluster, pathogen = test_dat$pathogen)

library(ggplot2)

fviz_cluster(km_res, data = data_scaled,
             geom = "point",
             ellipse.type = "norm", # or "convex" for convex hulls
             labelsize = 5,         # Adjust label size
             main = "K-means Clustering of Pathogens") +
  geom_text(aes(label = data_scaled_with_labels$pathogen), 
            size = 3, vjust = 1.5, hjust = 0.5)


# Add cluster information to the original data
test_dat$cluster <- km_res$cluster

# View the updated data
head(test_dat)

aggregate(data_numeric, by = list(cluster = test_dat$cluster), mean)

fviz_cluster(km_res, data = data_scaled, palette = c("#2E9FDF", "#00AFBB", "#E7B800", "#FC4E07"), ellipse.type = "euclid", # Concentration ellipse 
             star.plot = TRUE, # Add segments from centroids to items 
             repel = TRUE, # Avoid label overplotting (slow)
             ggtheme = theme_minimal() )



# View the clustering results
km_res$cluster




# View the firt 3 rows of the data
head(scaled_test_dat, n = 3)

#packages
install.packages("factoextra")
library(factoextra)

#Estimating the optimal number of clusters
#location of a bend (knee) in the plot is generally considered as an indicator of the appropriate number of clusters
#bend indicates that additional clusters beyond the fourth have little value
fviz_nbclust(df_scaled, kmeans, method = "wss") + geom_vline(xintercept = 4, linetype = 2)

#Computing k-means clustering
#Compute k-means with k = 4
set.seed(123)
km.res <- kmeans(df, 4, nstart = 25)

#Print the results
print(km.res)

#Compute the mean of each variables by clusters using the original data
aggregate(USArrests, by=list(cluster=km.res$cluster), mean)

#If you want to add the point classifications to the original data, use this
dd <- cbind(USArrests, cluster = km.res$cluster)
head(dd)

#Accessing to the results of kmeans() function

# Cluster number for each of the observations
km.res$cluster

head(km.res$cluster, 4)

# Cluster size
km.res$size

# Cluster means
km.res$centers

#Visualizing k-means clusters
fviz_cluster(km.res, data = df, palette = c("#2E9FDF", "#00AFBB", "#E7B800", "#FC4E07"), ellipse.type = "euclid", # Concentration ellipse 
             star.plot = TRUE, # Add segments from centroids to items 
             repel = TRUE, # Avoid label overplotting (slow)
             ggtheme = theme_minimal() )
             