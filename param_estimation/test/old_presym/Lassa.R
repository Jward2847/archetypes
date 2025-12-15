rm(list = ls())

library(ggplot2)
library(stats)


# Parameters for Incubation Period (Log-Normal)
# EpiLPS, Data set from Akhmetzhanov et al. 2019
incubation_meanlog <- 2.467
incubation_sdlog <- 0.367

# Parameters for Serial Interval (Log-Normal)
serial_meanlog <- 2.0939616
serial_sdlog <- 0.8541395


# Probability Density Functions
incubation_pdf <- function(x) dlnorm(x, meanlog = incubation_meanlog, sdlog = incubation_sdlog)
serial_pdf <- function(x) dlnorm(x, meanlog = serial_meanlog, sdlog = serial_sdlog)

# Calculate Intersection and Pre-symptomatic Transmission
intersection <- uniroot(function(x) incubation_pdf(x) - serial_pdf(x), lower = 0.01, upper = 10)$root
pre_symptomatic <- (plnorm(intersection, meanlog = serial_meanlog, sdlog = serial_sdlog) - 
                      plnorm(intersection, meanlog = incubation_meanlog, sdlog = incubation_sdlog)) * 100

cat("Percentage of pre-symptomatic transmission:", pre_symptomatic, "%\n")

# Visualisation
time_range <- seq(0, 50, by = 0.1)
incubation_values <- incubation_pdf(time_range)
serial_values <- serial_pdf(time_range)

# Combine data for plotting
plot_data <- data.frame(
  Time = rep(time_range, 2),
  PDF = c(incubation_values, serial_values),
  Type = rep(c("Incubation Period (Log-Normal)", "Serial Interval (Log-Normal)"), each = length(time_range))
)

# Plot
ggplot(plot_data, aes(x = Time, y = PDF, color = Type)) +
  geom_line(size = 2.0) +
  scale_color_manual(values = c("Incubation Period (Log-Normal)" = "blue", 
                                "Serial Interval (Log-Normal)" = "red")) +
  labs(
    title = "Lassa virus",
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

