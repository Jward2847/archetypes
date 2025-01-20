library(ggplot2)
library(stats)
library(MASS) # For fitting distributions


# Define parameters for the Weibull distribution (Incubation Period)
# EpiLPS using data set from Pavlin 2014
weibull_percentiles <- c(4.50, 9.0)
incubation_mean <- 6.9
weibull_shape <- 2.302
weibull_scale <- 7.804

# Gamma distribution (Serial Interval)
# Qian et al. 2023
# https://pubmed.ncbi.nlm.nih.gov/37964296/#:~:text=Of%20six%20vaccination%20strategies%20explored,CI%200.90%2D0.91)%2C%200.89
serial_mean <- 9.20
serial_sd <- 4.40
gamma_shape <- (serial_mean / serial_sd)^2
gamma_rate <- serial_mean / (serial_sd^2)

# Define Probability Density Functions
incubation_pdf <- function(x) dweibull(x, shape = weibull_shape, scale = weibull_scale)
serial_pdf <- function(x) dgamma(x, shape = gamma_shape, rate = gamma_rate)

# Estimate intersection and pre-symptomatic transmission
intersection <- uniroot(function(x) incubation_pdf(x) - serial_pdf(x), lower = 0, upper = 15)$root
pre_symptomatic <- (pgamma(intersection, shape = gamma_shape, rate = gamma_rate) - 
                      pweibull(intersection, shape = weibull_shape, scale = weibull_scale)) * 100

cat("Percentage of pre-symptomatic transmission:", pre_symptomatic, "%\n")

# Visualization
time_range <- seq(0, 20, by = 0.1)
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
  geom_line(linewidth = 2.0) +
  scale_color_manual(values = c("Incubation Period (Weibull)" = "blue", "Serial Interval (Gamma)" = "red")) +
  labs(
    title = "Marburg virus",
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

