library(ggplot2)
library(stats)   

# Define parameters for the distributions
# Incubation period (Lognormal): mean = 4.4, median = 4.1, 25th = 3.2, 75th = 5.3 
# source EpiLPS using data set from Backer et al. 2020
# Approximate the lognormal parameters
incubation_meanlog <- 1.417
incubation_sdlog <- 0.381

# Serial interval (Normal): mean = 4.6, SD = 4.4
# source Yang et al. 2020 
# https://www.ncbi.nlm.nih.gov/pmc/articles/PMC7324649/#:~:text=Our%20estimated%20median%20incubation%20period,has%20a%20median%20of%204.6

#Serial interval parameters
serial_mean <- 4.6
serial_sd <- 4.4

# Define the probability density functions (PDFs)
incubation_pdf <- function(x) dlnorm(x, meanlog = incubation_meanlog, sdlog = incubation_sdlog) 
serial_pdf <- function(x) dnorm(x, mean = serial_mean, sd = serial_sd)

# Estimate area under curve between incubation period and serial interval 
v <- uniroot(\(x)incubation_pdf(x) - serial_pdf(x), lower = 0, upper = 5)$root
pre_symptomatic <-( pnorm(v, serial_mean, serial_sd) - plnorm(v, incubation_meanlog, incubation_sdlog)) * 100

# Percentage of pre-symptomatic transmission
cat("Percentage of pre-symptomatic transmission:",pre_symptomatic, "%")


# Figure
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
  geom_line(size = 2.0) +  
  scale_color_manual(values = c("Incubation Period (Lognormal)" = "blue", 
                                "Serial Interval (Normal)" = "red")) +
  labs(
    title = "SARS-CoV-2 (Wild type)",
    subtitle = "",
    x = "Time (days)",
    y = "Probability Density",
    color = ""
  ) +
  theme_minimal(base_size = 14) +  # Increase base font size
  theme(
    plot.title = element_text(face = "bold", size = 30, hjust = 0.5),  # Title styling
    plot.subtitle = element_text(size = 12, hjust = 0.5),  # Subtitle styling
    legend.position = "top",  # Move legend to the top
    legend.title = element_text(face = "bold"),
    legend.text = element_text(size = 20),
    axis.title = element_text(face = "bold"),
    axis.text = element_text(size = 18)
  ) +
  scale_x_continuous(expand = expansion(mult = c(0, 0.05))) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.05)))