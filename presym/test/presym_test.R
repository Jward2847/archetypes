# Load necessary libraries
library(ggplot2)
library(dplyr)
library(stats)   


# Define parameters for the distributions
# Incubation period (Lognormal): mean = 4.4, median = 4.1, 25th = 3.2, 75th = 5.3
# Approximate the lognormal parameters
incubation_meanlog <- log(4.4)  # Approximate mean of the log-transformed data
incubation_sdlog <- (log(5.3) - log(3.2)) / (2 * qnorm(0.75)) # Approximate SD of log-transformed data

# Serial interval (Normal): mean = 4.6, SD = 4.4
serial_mean <- 4.6
serial_sd <- 4.4

# Define the probability density functions (PDFs)
incubation_pdf <- function(x) dlnorm(x, meanlog = incubation_meanlog, sdlog = incubation_sdlog)
serial_pdf <- function(x) dnorm(x, mean = serial_mean, sd = serial_sd)

# Calculate the AUC for pre-symptomatic transmission (where time < 0 for serial interval)
pre_symptomatic_auc <- integrate(serial_pdf, lower = -Inf, upper = 0)$value
total_serial_auc <- integrate(serial_pdf, lower = -Inf, upper = Inf)$value

# Proportion of pre-symptomatic transmission
proportion_pre_symptomatic <- pre_symptomatic_auc / total_serial_auc
cat("Proportion of pre-symptomatic transmission:", proportion_pre_symptomatic, "\n")

# Visualization
time_range <- seq(-10, 20, by = 0.1) # Define time range for plotting
incubation_values <- incubation_pdf(time_range)
serial_values <- serial_pdf(time_range)

# Combine data for ggplot
plot_data <- data.frame(
  Time = rep(time_range, 2),
  PDF = c(incubation_values, serial_values),
  Type = rep(c("Incubation Period (Lognormal)", "Serial Interval (Normal)"), each = length(time_range))
)

# Plot
ggplot(plot_data, aes(x = Time, y = PDF, color = Type)) +
  geom_line(size = 1) +
  geom_area(data = subset(plot_data, Type == "Serial Interval (Normal)" & Time < 0),
            aes(y = PDF), fill = "blue", alpha = 0.3) +
  labs(
    title = "Comparison of Serial Interval and Incubation Period",
    x = "Time (days)",
    y = "Probability Density",
    color = "Distribution"
  ) +
  theme_minimal()
