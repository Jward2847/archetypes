library(ggplot2)

# Parameters for Incubation Period (Gamma)
incubation_mean <- 5.7
incubation_sd <- 2.1
incubation_variance <- incubation_sd^2
incubation_shape <- (incubation_mean^2) / incubation_variance
incubation_rate <- incubation_mean / incubation_variance

# Parameters for Serial Interval (Gamma)
serial_shape <- 3.9227199
serial_rate <- 0.3298345

# Print fitted parameters
cat("Incubation Period Parameters:\n")
cat("Shape:", incubation_shape, "\n")
cat("Rate:", incubation_rate, "\n\n")

# Probability Density Functions
incubation_pdf <- function(x) dgamma(x, shape = incubation_shape, rate = incubation_rate)
serial_pdf <- function(x) dgamma(x, shape = serial_shape, rate = serial_rate)

# Calculate Intersection and Pre-symptomatic Transmission
intersection <- uniroot(function(x) incubation_pdf(x) - serial_pdf(x), lower = 0.01, upper = 5)$root
pre_symptomatic <- (pgamma(intersection, shape = serial_shape, rate = serial_rate) - 
                      pgamma(intersection, shape = incubation_shape, rate = incubation_rate)) * 100

cat("Percentage of pre-symptomatic transmission:", pre_symptomatic, "%\n")

# Visualization
time_range <- seq(0, 35, by = 0.1)
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
    title = "CCHVF",
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
