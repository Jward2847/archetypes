
library(ggplot2)
library(stats)

rm(list = ls())
# Parameters for Incubation Period (Weibull)
incubation_mean <- 3.30
incubation_sd <- 1.50
incubation_95th <- 6.0

# Solve for Weibull shape and scale parameters
weibull_solver <- function(shape) {
  scale <- incubation_sd / sqrt(gamma(1 + 2 / shape) - (gamma(1 + 1 / shape))^2)
  mean_est <- scale * gamma(1 + 1 / shape)
  mean_est - incubation_mean
}

incubation_shape <- uniroot(weibull_solver, c(0.1, 10))$root
incubation_scale <- incubation_sd / sqrt(gamma(1 + 2 / incubation_shape) - (gamma(1 + 1 / incubation_shape))^2)

# Validate 95th percentile
incubation_95th_check <- qweibull(0.95, shape = incubation_shape, scale = incubation_scale)
cat("Fitted Weibull Parameters:\n")
cat("Shape:", incubation_shape, "\n")
cat("Scale:", incubation_scale, "\n")
cat("95th Percentile Check:", incubation_95th_check, "\n")

# Parameters for Serial Interval (Gamma)
serial_mean <- 8.0
serial_median <- 6.8

# Solve for Gamma shape and rate parameters
gamma_solver <- function(shape) {
  rate <- shape / serial_mean
  median_est <- qgamma(0.5, shape = shape, rate = rate)
  median_est - serial_median
}

serial_shape <- uniroot(gamma_solver, c(0.1, 10))$root
serial_rate <- serial_shape / serial_mean

cat("Fitted Gamma Parameters:\n")
cat("Shape:", serial_shape, "\n")
cat("Rate:", serial_rate, "\n")

# Probability Density Functions
incubation_pdf <- function(x) dweibull(x, shape = incubation_shape, scale = incubation_scale)
serial_pdf <- function(x) dgamma(x, shape = serial_shape, rate = serial_rate)

# Calculate Intersection and Pre-symptomatic Transmission
intersection <- uniroot(function(x) incubation_pdf(x) - serial_pdf(x), lower = 0, upper = 5)$root
pre_symptomatic <- (pgamma(intersection, shape = serial_shape, rate = serial_rate) - 
                      pweibull(intersection, shape = incubation_shape, scale = incubation_scale)) * 100

cat("Percentage of pre-symptomatic transmission:", pre_symptomatic, "%\n")

# Visualisation
time_range <- seq(0, 30, by = 0.1)
incubation_values <- incubation_pdf(time_range)
serial_values <- serial_pdf(time_range)

# Combine data for plotting
plot_data <- data.frame(
  Time = rep(time_range, 2),
  PDF = c(incubation_values, serial_values),
  Type = rep(c("Incubation Period (Weibull)", "Serial Interval (Gamma)"), each = length(time_range))
)

# Plot
ggplot(plot_data, aes(x = Time, y = PDF, color = Type)) +
  geom_line(size = 2.0) +
  scale_color_manual(values = c("Incubation Period (Weibull)" = "blue", 
                                "Serial Interval (Gamma)" = "red")) +
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
