library(ggplot2)
library(ggfittext) # For better text fitting
library(treemapify) # For treemap visualization

# Create the data frame
df <- data.frame(
  param = c("Basic reproduction number", "Dispersion parameter", 
            "Serial interval", "Incubation period", "Latent period", 
            "Infectious period", "Case fatality risk", "Infection fatality risk"),
  total = c(66, 15, 48, 68, 27, 33, 65, 22)
)


# Create treemap plot with modified text size and color scheme
ggplot(df, aes(area = total, fill = param, label = param)) +
  geom_treemap() +
  geom_treemap_text(
    grow = FALSE,  # Prevent text from expanding too much
    reflow = TRUE, 
    colour = "black",  # Ensure readability
    place = "centre",
    min.size = 10,  # Balanced text size
    fontface = "bold"  # Make text bold
  ) +
  theme_minimal() +
  scale_fill_brewer(palette = "Paired") +  # Professional color scheme
  theme(legend.position = "none")



# Create the data frame with updated labels
df <- data.frame(
  virus = c("Influenza", "CCHFV", "Ebola", "LASV", "MARV", 
            "MERS-CoV", "MPV", "NiV", "RVFV", "SARS-CoV-1", 
            "SARS-CoV-2", "ZIKV"),
  total = c(72, 7, 35, 11, 12, 21, 28, 25, 3, 21, 97, 12)
)

# Create treemap plot with modified text size and color scheme
ggplot(df, aes(area = total, fill = virus, label = virus)) +
  geom_treemap() +
  geom_treemap_text(
    grow = FALSE,  # Prevent text from expanding too much
    reflow = TRUE, 
    colour = "black",  # Ensure readability
    place = "centre",
    min.size = 10,  # Balanced text size
    fontface = "bold"  # Make text bold
  ) +
  theme_minimal() +
  scale_fill_brewer(palette = "Paired") +  # Professional color scheme
  theme(legend.position = "none")



library(RColorBrewer)
scale_fill_brewer(palette = "Set3") 






