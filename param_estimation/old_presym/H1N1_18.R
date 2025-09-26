
library(ggplot2)
library(stats)
library(epiparameter)

rm(list = ls())

flu_incubation <- epiparameter_db(
  pathogen = "Influenza-A-H1N1",
  epi_name = "incubation period",
  single_epiparameter = TRUE
)

flu_incubation






# Parameters for Incubation Period 
meanlog <- 0.342
sdlog <- 0.272

# Parameters for Serial Interval (Normal)
serial_mean <- 3.0
serial_sd <- 0.9

# Probability Density Functions
incubation_pdf <- function(x) dlnorm(x, meanlog = meanlog, sdlog = sdlog)
serial_pdf <- function(x) dnorm(x, mean = serial_mean, sd = serial_sd)

# Calculate Intersection and Pre-symptomatic Transmission
intersection <- uniroot(function(x) incubation_pdf(x) - serial_pdf(x), lower = -2, upper = 1)$root
pre_symptomatic <- (pnorm(intersection, mean = serial_mean, sd = serial_sd) - 
                      plnorm(intersection, meanlog = meanlog, sdlog = sdlog)) * 100

cat("Percentage of pre-symptomatic transmission:", pre_symptomatic, "%\n")

# Visualization
time_range <- seq(-2, 6, by = 0.1)
incubation_values <- incubation_pdf(time_range)
serial_values <- serial_pdf(time_range)

# Combine data for plotting
plot_data <- data.frame(
  Time = rep(time_range, 2),
  PDF = c(incubation_values, serial_values),
  Type = rep(c("Incubation Period (Weibull)", "Serial Interval (Normal)"), each = length(time_range))
)

# Plot
ggplot(plot_data, aes(x = Time, y = PDF, color = Type)) +
  geom_line(size = 2.0) +
  scale_color_manual(values = c("Incubation Period (Weibull)" = "blue", 
                                "Serial Interval (Normal)" = "red")) +
  labs(
    title = "New Disease: Incubation Period and Serial Interval",
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


