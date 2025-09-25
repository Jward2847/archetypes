library(ggplot2)
library(stats)


rm(list = ls())



# Parameters for Incubation Period (Gamma)
# Modelled with EpiLPS. Data set from Nikolay et al. 2019
incubation_shape <- 16.123
incubation_rate <- 1.714

# Parameters for Serial Interval (Gamma)
# Nikolay et al. 2019	
# https://www.nejm.org/doi/10.1056/NEJMoa1805376?url_ver=Z39.88-2003&rfr_id=ori:rid:crossref.org&rfr_dat=cr_pub%20%200pubmed#f1
serial_mean <- 12.7
serial_sd <- 3.0
serial_variance <- serial_sd^2
serial_shape <- (serial_mean^2) / serial_variance
serial_rate <- serial_mean / serial_variance

# Probability Density Functions
incubation_pdf <- function(x) dgamma(x, shape = incubation_shape, rate = incubation_rate)
serial_pdf <- function(x) dgamma(x, shape = serial_shape, rate = serial_rate)

# Calculate Intersection and Pre-symptomatic Transmission
intersection <- uniroot(function(x) incubation_pdf(x) - serial_pdf(x), lower = 10, upper = 15)$root
pre_symptomatic <- (pgamma(intersection, shape = serial_shape, rate = serial_rate) - 
                      pgamma(intersection, shape = incubation_shape, rate = incubation_rate)) * 100

cat("Percentage of pre-symptomatic transmission:", pre_symptomatic, "%\n")

# Visualisation
time_range <- seq(0, 25, by = 0.1)
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
    title = "Nipah virus",
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

