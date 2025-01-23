library(epiparameter)
library(ggplot2)
library(stats)

rm(list = ls())

#Extract Zika incubation period paramters
zika_incubation <- epiparameter_db(
  disease = "zika",
  epi_name = "incubation period",
  single_epiparameter = TRUE
)

zika_incubation



# Parameters for Incubation Period (Log-Normal)
# Lessler et al. 2016	
# https://pmc.ncbi.nlm.nih.gov/articles/PMC5096355/#:~:text=The%20virus%20was%20detectable%20in,incidence%20of%20Zika%20virus%20infection.
meanlog <- 1.775
sdlog <- 0.405

# Serial Interval Parameters (Fitting Gamma Distribution)
serial_median <- 12
serial_iqr <- c(10, 14.5)
q25 <- serial_iqr[1]
q75 <- serial_iqr[2]

# gamma solver
gamma_solver <- function(k) {
  scale <- (q75 - q25) / (qgamma(0.75, shape = k) - qgamma(0.25, shape = k))
  rate <- 1 / scale
  median_est <- qgamma(0.5, shape = k, rate = rate)
  median_est - serial_median
}

# Test gamma_solver at the bounds
cat("gamma_solver(0.1):", gamma_solver(0.1), "\n")
cat("gamma_solver(10):", gamma_solver(10), "\n")

# Adjust bounds if necessary
serial_shape <- uniroot(gamma_solver, c(0.01, 20))$root
serial_rate <- 1 / ((q75 - q25) / (qgamma(0.75, shape = serial_shape) - qgamma(0.25, shape = serial_shape)))

cat("Fitted Gamma Parameters:\n")
cat("Shape:", serial_shape, "\n")
cat("Rate:", serial_rate, "\n")

# Serial Interval Parameters (Fitted Gamma Distribution)
serial_shape <- 13.24165
serial_rate <- 1.075821

# Probability Density Functions
incubation_pdf <- function(x) dlnorm(x, meanlog = meanlog, sdlog = sdlog)
serial_pdf <- function(x) dgamma(x, shape = serial_shape, rate = serial_rate)

# Calculate Intersection and Pre-symptomatic Transmission
intersection <- uniroot(function(x) incubation_pdf(x) - serial_pdf(x), lower = 5, upper = 10)$root
pre_symptomatic <- (pgamma(intersection, shape = serial_shape, rate = serial_rate) - 
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
  Type = rep(c("Incubation Period (Log-Normal)", "Serial Interval (Gamma-Fitted)"), each = length(time_range))
)

# Plot
ggplot(plot_data, aes(x = Time, y = PDF, color = Type)) +
  geom_line(size = 2.0) +
  scale_color_manual(values = c("Incubation Period (Log-Normal)" = "blue", 
                                "Serial Interval (Gamma-Fitted)" = "red")) +
  labs(
    title = "Zika virus",
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