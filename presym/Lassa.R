rm(list = ls())

library(ggplot2)
library(stats)

# Parameters for Incubation Period (Log-normal)
# EpiLPS, Data set from Akhmetzhanov et al. 2019
meanlog <- 2.467
sdlog <- 0.367

# Parameters for Serial Interval (Gamma)
# Zhao et al. 2020	https://pmc.ncbi.nlm.nih.gov/articles/PMC7019145/
serial_mean <- 7.8
serial_sd <- 10.7
serial_variance <- serial_sd^2
serial_shape <- (serial_mean^2) / serial_variance
serial_rate <- serial_mean / serial_variance

# Probability Density Functions
incubation_pdf <- function(x) dlnorm(x, meanlog = meanlog, sdlog = sdlog)
serial_pdf <- function(x) dgamma(x, shape = serial_shape, rate = serial_rate)

# Calculate Intersection and Pre-symptomatic Transmission
intersection <- uniroot(function(x) incubation_pdf(x) - serial_pdf(x), lower = 0, upper = 10)$root
pre_symptomatic <- (pgamma(intersection, shape = serial_shape, rate = serial_rate) - 
                      plnorm(intersection, meanlog = meanlog, sdlog = sdlog)) * 100

cat("Percentage of pre-symptomatic transmission:", pre_symptomatic, "%\n")

# Visualisation
time_range <- seq(0, 50, by = 0.1)
incubation_values <- incubation_pdf(time_range)
serial_values <- serial_pdf(time_range)

# Combine data for plotting
plot_data <- data.frame(
  Time = rep(time_range, 2),
  PDF = c(incubation_values, serial_values),
  Type = rep(c("Incubation Period (Log-normal)", "Serial Interval (Gamma)"), each = length(time_range))
)

# Plot
ggplot(plot_data, aes(x = Time, y = PDF, color = Type)) +
  geom_line(size = 2.0) +
  scale_color_manual(values = c("Incubation Period (Log-normal)" = "blue", 
                                "Serial Interval (Gamma)" = "red")) +
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


#######################################################################

#Lassa Lo locarno 

# Parameters for Incubation Period (Log-Normal)
incubation_meanlog <- 2.467
incubation_sdlog <- 0.367

# Parameters for Serial Interval (Gamma)
serial_shape <- 2.144
serial_rate <- 0.358

# Probability Density Functions
incubation_pdf <- function(x) dlnorm(x, meanlog = incubation_meanlog, sdlog = incubation_sdlog)
serial_pdf <- function(x) dgamma(x, shape = serial_shape, rate = serial_rate)

# Calculate Intersection and Pre-symptomatic Transmission
intersection <- uniroot(function(x) incubation_pdf(x) - serial_pdf(x), lower = 0.01, upper = 10)$root
pre_symptomatic <- (pgamma(intersection, shape = serial_shape, rate = serial_rate) - 
                      plnorm(intersection, meanlog = incubation_meanlog, sdlog = incubation_sdlog)) * 100

cat("Percentage of pre-symptomatic transmission:", pre_symptomatic, "%\n")

# Visualization
time_range <- seq(0, 35, by = 0.1)
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
  geom_line(size = 2.0) +
  scale_color_manual(values = c("Incubation Period (Log-Normal)" = "blue", 
                                "Serial Interval (Gamma)" = "red")) +
  labs(
    title = "Lo Iacono data set",
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

#################################################################################
#Outbreak data set

library(ggplot2)

# Parameters for Incubation Period (Log-Normal)
incubation_meanlog <- 2.467
incubation_sdlog <- 0.367

# Parameters for Serial Interval (Log-Normal)
serial_meanlog <- 2.667
serial_sdlog <- 0.360

# Probability Density Functions
incubation_pdf <- function(x) dlnorm(x, meanlog = incubation_meanlog, sdlog = incubation_sdlog)
serial_pdf <- function(x) dlnorm(x, meanlog = serial_meanlog, sdlog = serial_sdlog)

# Calculate Intersection and Pre-symptomatic Transmission
intersection <- uniroot(function(x) incubation_pdf(x) - serial_pdf(x), lower = 10, upper = 15)$root
pre_symptomatic <- (plnorm(intersection, meanlog = serial_meanlog, sdlog = serial_sdlog) - 
                      plnorm(intersection, meanlog = incubation_meanlog, sdlog = incubation_sdlog)) * 100

cat("Percentage of pre-symptomatic transmission:", pre_symptomatic, "%\n")

# Visualization
time_range <- seq(0, 40, by = 0.1)
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
    title = "Collected data",
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

#############################################################################

#Combined 

# Parameters for Incubation Period (Log-Normal)
incubation_meanlog <- 2.467
incubation_sdlog <- 0.367

# Parameters for Serial Interval (Gamma)
serial_shape <- 2.062
serial_rate <- 0.207

# Probability Density Functions
incubation_pdf <- function(x) dlnorm(x, meanlog = incubation_meanlog, sdlog = incubation_sdlog)
serial_pdf <- function(x) dgamma(x, shape = serial_shape, rate = serial_rate)

# Calculate Intersection and Pre-symptomatic Transmission
intersection <- uniroot(function(x) incubation_pdf(x) - serial_pdf(x), lower = 0.01, upper = 10)$root
pre_symptomatic <- (pgamma(intersection, shape = serial_shape, rate = serial_rate) - 
                      plnorm(intersection, meanlog = incubation_meanlog, sdlog = incubation_sdlog)) * 100

cat("Percentage of pre-symptomatic transmission:", pre_symptomatic, "%\n")

# Visualization
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
  geom_line(size = 2.0) +
  scale_color_manual(values = c("Incubation Period (Log-Normal)" = "blue", 
                                "Serial Interval (Gamma)" = "red")) +
  labs(
    title = "Combined data set",
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


