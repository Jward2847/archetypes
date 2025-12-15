library(epiparameter)
library(ggplot2)
library(stats)

rm(list = ls())

ebola_incubation <- epiparameter_db(
  disease = "ebola",
  epi_name = "incubation period",
  single_epiparameter = TRUE
)

ebola_SI <- epiparameter_db(
  disease = "ebola",
  epi_name = "serial interval",
  single_epiparameter = TRUE
)

ebola_incubation 

ebola_SI




# Parameters for Incubation Period (Gamma)
# WHO Ebola Response Team	https://www.nejm.org/doi/full/10.1056/NEJMc1414992

incubation_shape <- 1.578
incubation_scale <- 6.528

# Parameters for Serial Interval (Gamma)
# WHO Ebola Response Team	https://www.nejm.org/doi/full/10.1056/NEJMc1414992
serial_shape <- 2.188
serial_scale <- 6.490

# Probability Density Functions
incubation_pdf <- function(x) dgamma(x, shape = incubation_shape, scale = incubation_scale)
serial_pdf <- function(x) dgamma(x, shape = serial_shape, scale = serial_scale)

# Calculate Intersection and Pre-symptomatic Transmission
intersection <- uniroot(function(x) incubation_pdf(x) - serial_pdf(x), lower = 5, upper = 10)$root
pre_symptomatic <- (pgamma(intersection, shape = serial_shape, scale = serial_scale) - 
                      pgamma(intersection, shape = incubation_shape, scale = incubation_scale)) * 100

cat("Percentage of pre-symptomatic transmission:", pre_symptomatic, "%\n")

# Visualisation
time_range <- seq(0, 50, by = 0.1)
incubation_values <- incubation_pdf(time_range)
serial_values <- serial_pdf(time_range)

# Combine data for plotting
plot_data <- data.frame(
  Time = rep(time_range, 2),
  PDF = c(incubation_values, serial_values),
  Type = rep(c("Incubation Period (Gamma)", "Serial Interval (Gamma)"), each = length(time_range))
)

# Plot
ggplot(plot_data, aes(x = Time, y = PDF, color = Type)) +
  geom_line(size = 2.0) +
  scale_color_manual(values = c("Incubation Period (Gamma)" = "blue", 
                                "Serial Interval (Gamma)" = "red")) +
  labs(
    title = "Ebola virus",
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
