
library(ggplot2)
library(stats)
library(epiparameter)

rm(list = ls())

#Extra incubation period gamma dist parameters
SARS_inc<- epiparameter_db(
  disease = "SARS",
  epi_name = "incubation period",
  single_epiparameter = TRUE
)

SARS_inc

# Parameters for Incubation Period (Log-normal)
# Lessler et al. 2009	https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4327893/
meanlog <- 1.386
sdlog <- 0.593

# Parameters for Serial Interval (Weibull)
serial_mean <- 8.4
serial_sd <- 3.8

# Function to calculate Weibull shape parameter numerically
weibull_shape_solver <- function(k) {
  scale <- serial_sd / sqrt(gamma(1 + 2 / k) - (gamma(1 + 1 / k))^2)
  mean_est <- scale * gamma(1 + 1 / k)
  mean_est - serial_mean
}

weibull_shape <- uniroot(weibull_shape_solver, c(0.1, 10))$root
weibull_scale <- serial_sd / sqrt(gamma(1 + 2 / weibull_shape) - (gamma(1 + 1 / weibull_shape))^2)

# Probability Density Functions
incubation_pdf <- function(x) dlnorm(x, meanlog = meanlog, sdlog = sdlog)
serial_pdf <- function(x) dweibull(x, shape = weibull_shape, scale = weibull_scale)

# Calculate Intersection and Pre-symptomatic Transmission
intersection <- uniroot(function(x) incubation_pdf(x) - serial_pdf(x), lower = 0.01, upper = 2)$root
pre_symptomatic <- (pweibull(intersection, shape = weibull_shape, scale = weibull_scale) - 
                      plnorm(intersection, meanlog = meanlog, sdlog = sdlog)) * 100

cat("Percentage of pre-symptomatic transmission:", pre_symptomatic, "%\n")

# Visualization
time_range <- seq(0, 20, by = 0.1)
incubation_values <- incubation_pdf(time_range)
serial_values <- serial_pdf(time_range)

# Combine data for plotting
plot_data <- data.frame(
  Time = rep(time_range, 2),
  PDF = c(incubation_values, serial_values),
  Type = rep(c("Incubation Period (Log-normal)", "Serial Interval (Weibull)"), each = length(time_range))
)

# Plot
ggplot(plot_data, aes(x = Time, y = PDF, color = Type)) +
  geom_line(size = 2.0) +
  scale_color_manual(values = c("Incubation Period (Log-normal)" = "blue", 
                                "Serial Interval (Weibull)" = "red")) +
  labs(
    title = "SARS-CoV-1",
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