library(ggplot2)
library(stats)   

# Parameters for the incubation period (Gamma)
# Galmiche et al. 2023 
# https://www.thelancet.com/journals/lanmic/article/PIIS2666-5247(23)00005-8/fulltext
incubation_mean <- 4.96
incubation_sd <- 2.32
incubation_shape <- (incubation_mean / incubation_sd)^2
incubation_rate <- incubation_mean / (incubation_sd^2)

# Parameters for the serial interval (Gamma)
# Heiden and Buchholz 2022
# https://www.ncbi.nlm.nih.gov/pmc/articles/PMC9530378/#:~:text=The%20mean%20serial%20interval%20dropped,Omicron%20circulation%20(Table%201).
serial_mean <- 4.5
serial_sd <- (4.54 - 4.46) / (2 * qnorm(0.975))  # Calculate SD from 95% CI
serial_shape <- (serial_mean / serial_sd)^2
serial_rate <- serial_mean / (serial_sd^2)

# Define Probability Density Functions
incubation_pdf <- function(x) dgamma(x, shape = incubation_shape, rate = incubation_rate)
serial_pdf <- function(x) dgamma(x, shape = serial_shape, rate = serial_rate)

# Calculate intersection point and pre-symptomatic transmission
intersection <- uniroot(function(x) incubation_pdf(x) - serial_pdf(x), lower = 0, upper = 5)$root
pre_symptomatic <- (pgamma(intersection, shape = serial_shape, rate = serial_rate) - 
                      pgamma(intersection, shape = incubation_shape, rate = incubation_rate)) * 100

cat("Percentage of pre-symptomatic transmission:", pre_symptomatic, "%\n")

# Visualization
time_range <- seq(0, 20, by = 0.1)
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
  scale_color_manual(values = c("Incubation Period (Gamma)" = "blue", "Serial Interval (Gamma)" = "red")) +
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


library(ggplot2)
library(stats)

# Parameters for Incubation Period (Gamma)
incubation_mean <- 4.96
incubation_sd <- 2.32
incubation_variance <- incubation_sd^2
incubation_shape <- (incubation_mean^2) / incubation_variance
incubation_rate <- incubation_mean / incubation_variance

# Parameters for Serial Interval (Gamma)
serial_mean <- 4.5
serial_shape <- 2.0
serial_rate <- 0.44

# Probability Density Functions
incubation_pdf <- function(x) dgamma(x, shape = incubation_shape, rate = incubation_rate)
serial_pdf <- function(x) dgamma(x, shape = serial_shape, rate = serial_rate)

# Calculate Intersection and Pre-symptomatic Transmission
v <- uniroot(\(x)incubation_pdf(x) - serial_pdf(x), lower = 0, upper = 2)$root
pre_symptomatic <-( pnorm(v, serial_mean, serial_sd) - plnorm(v, incubation_meanlog, incubation_sdlog)) * 100

# Percentage of pre-symptomatic transmission
cat("Percentage of pre-symptomatic transmission:",pre_symptomatic, "%")


# Visualization
time_range <- seq(0, 15, by = 0.1)
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
  scale_color_manual(values = c("Incubation Period (Gamma)" = "blue", "Serial Interval (Gamma)" = "red")) +
  labs(
    title = "SARS-CoV-2 (Alpha)",
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
