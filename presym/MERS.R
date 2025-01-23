
library(ggplot2)
library(stats)
library(epiparameter)

rm(list = ls())

#Extra SI gamma dist paramters
MERS_SI<- epiparameter_db(
  disease = "MERS",
  epi_name = "serial interval",
  single_epiparameter = TRUE
)



# Parameters for Incubation Period (Log-normal)
# EpiLPS, Data set from Cauchemez et al. 2014
meanlog <- 1.567 
sdlog <- 0.486

# Parameters for Serial Interval (Gamma)
# Cowling et al. 2015	
# https://www.eurosurveillance.org/content/10.2807/1560-7917.ES2015.20.25.21163?crawler=true
serial_shape <- 20.250
serial_scale <- 0.622

# Probability Density Functions
incubation_pdf <- function(x) dlnorm(x, meanlog = meanlog, sdlog = sdlog)
serial_pdf <- function(x) dgamma(x, shape = serial_shape, scale = serial_scale)

# Calculate Intersection and Pre-symptomatic Transmission
intersection <- uniroot(function(x) incubation_pdf(x) - serial_pdf(x), lower = 5, upper = 10)$root
pre_symptomatic <- (pgamma(intersection, shape = serial_shape, scale = serial_scale) - 
                      plnorm(intersection, meanlog = meanlog, sdlog = sdlog)) * 100

cat("Percentage of pre-symptomatic transmission:", pre_symptomatic, "%\n")

# Visualisation
time_range <- seq(0, 25, by = 0.1)
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
    title = "MERS-CoV",
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
