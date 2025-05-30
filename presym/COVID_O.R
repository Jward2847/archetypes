rm(list = ls())

library(ggplot2)
library(stats)


# Parameters for Incubation Period (Log-normal)
#EpiLPS - dataset from (Data set from Backer et al. 2022)
logmean <- 1.045
sdlog <- 0.57

# Parameters for Serial Interval (Normal)
#Kremer et al. 2022
#https://pmc.ncbi.nlm.nih.gov/articles/PMC9328897/#SD1
serial_median <- 2.75
serial_sd <- 2.54
serial_mean <- serial_median # Assuming symmetry 

# Probability Density Functions
incubation_pdf <- function(x) dlnorm(x, meanlog = logmean, sdlog = sdlog)
serial_pdf <- function(x) dnorm(x, mean = serial_mean, sd = serial_sd)

# Calculate Intersection and Pre-symptomatic Transmission
intersection <- uniroot(function(x) incubation_pdf(x) - serial_pdf(x), lower = 0, upper = 3)$root
pre_symptomatic <- (pnorm(intersection, mean = serial_mean, sd = serial_sd) - 
                      plnorm(intersection, meanlog = logmean, sdlog = sdlog)) * 100

cat("Percentage of pre-symptomatic transmission:", pre_symptomatic, "%\n")

# Visualisation
time_range <- seq(-5, 15, by = 0.1)
incubation_values <- incubation_pdf(time_range)
serial_values <- serial_pdf(time_range)

# Combine data for plotting
plot_data <- data.frame(
  Time = rep(time_range, 2),
  PDF = c(incubation_values, serial_values),
  Type = rep(c("Incubation Period (Log-normal)", "Serial Interval (Normal)"), each = length(time_range))
)

# Plot
ggplot(plot_data, aes(x = Time, y = PDF, color = Type)) +
  geom_line(size = 2.0) +
  scale_color_manual(values = c("Incubation Period (Log-normal)" = "blue", 
                                "Serial Interval (Normal)" = "red")) +
  labs(
    title = "SARS-CoV-2 (Omicron)",
    x = "Time (days)",
    y = "Probability Density",
    color = ""
  ) +
  theme_minimal(base_size = 14) +
  theme(
    plot.title = element_text(face = "bold", size = 20, hjust = 0.5),
    legend.position = "top",
    axis.title = element_text(face = "bold"),
    axis.text = element_text(size = 12)
  ) +
  scale_x_continuous(expand = expansion(mult = c(0, 0.05))) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.05)))
