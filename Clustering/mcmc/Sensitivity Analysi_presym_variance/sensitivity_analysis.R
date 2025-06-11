# Sensitivity Analysis for Presymptomatic Transmission of Ebola
#
# This script investigates how the variance of the Serial Interval (SI) and
# Incubation Period (IP) distributions affects the calculated proportion of
# presymptomatic transmission, P(SI < IP).
#
# We use the original parameters for Ebola from the main analysis and compare
# the result with a second set of "low-variance" parameters that share the
# same mean but are less dispersed.

# --- 1. Load Libraries ---
library(ggplot2)
library(dplyr)
set.seed(456) # For reproducibility of this specific analysis

# --- 2. Define Parameters ---

# Original High-Variance Parameters (from your data source)
ip_shape_orig <- 1.578
ip_scale_orig <- 6.528
si_shape_orig <- 2.188
si_scale_orig <- 6.490

# Calculate original means (for verification)
mean_ip_orig <- ip_shape_orig * ip_scale_orig
mean_si_orig <- si_shape_orig * si_scale_orig

# Create new Low-Variance Parameters
# We reduce variance while keeping the mean constant.
# For a Gamma distribution: Mean = shape * scale, Variance = shape * scale^2.
# To reduce variance by a factor of F^2, we divide scale by F and multiply shape by F.
# Let's use F=10, for a final, very aggressive reduction to demonstrate the point.
variance_reduction_factor <- 10

ip_shape_low_var <- ip_shape_orig * variance_reduction_factor
ip_scale_low_var <- ip_scale_orig / variance_reduction_factor
si_shape_low_var <- si_shape_orig * variance_reduction_factor
si_scale_low_var <- si_scale_orig / variance_reduction_factor

# --- 3. Simulation Function ---

estimate_presymp_proportion <- function(si_shape, si_scale, ip_shape, ip_scale, n_samples = 50000) {
  # Generate samples from the specified Gamma distributions
  si_samples <- rgamma(n_samples, shape = si_shape, scale = si_scale)
  ip_samples <- rgamma(n_samples, shape = ip_shape, scale = ip_scale)
  
  # Calculate the proportion of times SI < IP
  presymp_proportion <- mean(si_samples < ip_samples, na.rm = TRUE)
  return(presymp_proportion)
}

# --- 4. Run Simulations ---

prop_high_variance <- estimate_presymp_proportion(si_shape_orig, si_scale_orig, ip_shape_orig, ip_scale_orig)
prop_low_variance <- estimate_presymp_proportion(si_shape_low_var, si_scale_low_var, ip_shape_low_var, ip_scale_low_var)

# --- 5. Print Results ---

cat("--- Sensitivity Analysis Results ---\n\n")
cat(sprintf("Original (High Variance) Parameters:\n"))
cat(sprintf("  - Mean IP: %.2f days, Mean SI: %.2f days\n", mean_ip_orig, mean_si_orig))
cat(sprintf("  - Calculated Presymptomatic Proportion: %.1f%%\n\n", prop_high_variance * 100))

cat(sprintf("Synthesized (Low Variance) Parameters:\n"))
cat(sprintf("  - Mean IP: %.2f days, Mean SI: %.2f days (Means are preserved)\n", 
            (ip_shape_low_var * ip_scale_low_var), (si_shape_low_var * si_scale_low_var)))
cat(sprintf("  - Calculated Presymptomatic Proportion: %.1f%%\n", prop_low_variance * 100))


# --- 6. Visualization ---

# Create a data frame for ggplot
time_range <- seq(0, 50, length.out = 500)
plot_data <- tibble(
  Time = time_range,
  # High-variance densities
  IP_High_Var = dgamma(Time, shape = ip_shape_orig, scale = ip_scale_orig),
  SI_High_Var = dgamma(Time, shape = si_shape_orig, scale = si_scale_orig),
  # Low-variance densities
  IP_Low_Var = dgamma(Time, shape = ip_shape_low_var, scale = ip_scale_low_var),
  SI_Low_Var = dgamma(Time, shape = si_shape_low_var, scale = si_scale_low_var)
) %>%
  tidyr::pivot_longer(
    cols = -Time,
    names_to = c("Distribution", "Variance"),
    names_pattern = "(.*)_(.*)_Var",
    values_to = "Density"
  )

# Generate the plot
sensitivity_plot <- ggplot(plot_data, aes(x = Time, y = Density, color = Distribution, linetype = Variance)) +
  geom_line(linewidth = 1.1) +
  scale_color_manual(values = c("IP" = "#0072B2", "SI" = "#D55E00")) +
  scale_linetype_manual(values = c("High" = "solid", "Low" = "dashed")) +
  labs(
    title = "Sensitivity of Presymptomatic Transmission to Distribution Variance",
    subtitle = "Comparing original high-variance parameters with synthesized low-variance parameters for Ebola",
    x = "Time (days)",
    y = "Probability Density",
    color = "Distribution",
    linetype = "Variance"
  ) +
  annotate(
    "text", x = 35, y = 0.08, hjust = 0, size = 4,
    label = paste0("High-Variance P(SI < IP) = ", round(prop_high_variance * 100, 1), "%")
  ) +
    annotate(
    "text", x = 35, y = 0.07, hjust = 0, size = 4,
    label = paste0("Low-Variance P(SI < IP) = ", round(prop_low_variance * 100, 1), "%")
  ) +
  theme_minimal(base_size = 14) +
  theme(legend.position = "bottom")

print(sensitivity_plot)

# Save the plot
ggsave("Clustering/mcmc/ebola_variance_sensitivity_plot.png", plot = sensitivity_plot, width = 12, height = 8)

cat("\nPlot saved to 'Clustering/mcmc/ebola_variance_sensitivity_plot.png'\n")

print("Script finished.") 