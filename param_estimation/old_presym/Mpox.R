library(ggplot2)
library(stats)

rm(list = ls())

# Parameters for Incubation Period (Log-Normal)
# EpiLPS. Data set from Miura et al. 2022
meanlog <- 2.091
sdlog <- 0.438

# Parameters for Serial Interval (Gamma)
# Ward et al. 2022
# https://www.bmj.com/content/379/bmj-2022-073153
serial_mean <- 8.0
serial_shape <- 0.8
serial_rate <- 0.1

# Probability Density Functions
incubation_pdf <- function(x) dlnorm(x, meanlog = meanlog, sdlog = sdlog)
serial_pdf <- function(x) dgamma(x, shape = serial_shape, rate = serial_rate)

# Calculate Intersection and Pre-symptomatic Transmission
intersection <- uniroot(function(x) incubation_pdf(x) - serial_pdf(x), lower = 0, upper = 50)$root
pre_symptomatic <- (pgamma(intersection, shape = serial_shape, rate = serial_rate) - 
                      plnorm(intersection, meanlog = meanlog, sdlog = sdlog)) * 100

cat("Percentage of pre-symptomatic transmission:", pre_symptomatic, "%\n")

# Visualisation
time_range <- seq(0, 40, by = 0.1)
incubation_values <- incubation_pdf(time_range)
serial_values <- serial_pdf(time_range)

# Combine data for plotting
plot_data <- data.frame(
  Time = rep(time_range, 2),
  PDF = c(incubation_values, serial_values),
  Type = rep(c("Incubation Period (Log-Normal)", "Serial Interval (Gamma)"), each = length(time_range))
)

# Plot
ggplot(plot_data, aes(x = Time, y = PDF, color = Type)) +
  geom_line(linewidth = 2.0) +
  scale_color_manual(values = c("Incubation Period (Log-Normal)" = "blue", 
                                "Serial Interval (Gamma)" = "red")) +
  labs(
    title = "Monkeypox virus",
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
