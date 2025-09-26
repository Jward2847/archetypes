library(ggplot2)
library(dplyr)
if (!require(fitdistrplus)) install.packages("fitdistrplus", dependencies=TRUE)
library(ggplot2)
library(fitdistrplus)

# Define the centile values
gamma_shape <- 0.8
gamma_rate <- 0.1

# Define the Weibull parameters from the second table
weibull_shape <- 1.5
weibull_scale <- 8.4

# Generate the Gamma distribution
set.seed(123)  # for reproducibility
n <- 1000  # number of samples
gamma_data <- rgamma(n, shape = gamma_shape, rate = gamma_rate)

# Generate the Weibull distribution
weibull_data <- rweibull(n, shape = weibull_shape, scale = weibull_scale)

# Create data frames for plotting
gamma_df <- data.frame(value = gamma_data, distribution = "Gamma")
weibull_df <- data.frame(value = weibull_data, distribution = "Weibull")

# Combine the data frames
combined_df <- rbind(gamma_df, weibull_df)

# Plot both distributions on the same figure
ggplot(combined_df, aes(x = value, fill = distribution)) +
  geom_density(alpha = 0.5) +
  labs(title = "Gamma and Weibull Distributions",
       x = "Value",
       y = "Density",
       fill = "Distribution") +
  theme_minimal()

set.seed(123)  # for reproducibility
n <- 100000  # number of samples for Monte Carlo simulation
gamma_samples <- rgamma(n, shape = gamma_shape, rate = gamma_rate)
weibull_samples <- rweibull(n, shape = weibull_shape, scale = weibull_scale)

# Calculate the probability that the Gamma sample is less than the Weibull sample
probability <- mean(gamma_samples < weibull_samples)

# Convert the probability to a percentage
probability_percent <- probability * 100
probability_percent

###

set.seed(123)  # for reproducibility
n <- 100000  # number of samples for Monte Carlo simulation
gamma_samples <- rgamma(n, shape = gamma_shape, rate = gamma_rate)
weibull_samples <- rweibull(n, shape = weibull_shape, scale = weibull_scale)

# Calculate the probability that the Gamma sample is less than the Weibull sample
probability <- mean(gamma_samples < weibull_samples)

# Convert the probability to a percentage
probability_percent <- probability * 100

# Print the probability in percentage
print(paste("The probability that the serial interval is shorter than the incubation period is", round(probability_percent, 2), "%"))


###################################
#Ebola

# Given information for the serial interval
mean_si <- 15.40
ci_lower_si <- 13.20
ci_upper_si <- 17.50

# Given information for the incubation period
mean_ip <- 8.50
ci_lower_ip <- 7.70
ci_upper_ip <- 9.20

# Calculate standard error for both
se_si <- (ci_upper_si - mean_si) / 1.96
se_ip <- (ci_upper_ip - mean_ip) / 1.96

# Estimate the variance for both
variance_si <- se_si^2
variance_ip <- se_ip^2

# Calculate shape (k) and scale (θ) parameters for both
theta_si <- variance_si / mean_si
k_si <- mean_si^2 / variance_si

theta_ip <- variance_ip / mean_ip
k_ip <- mean_ip^2 / variance_ip

# Print the estimated parameters
cat("Serial Interval - Shape (k):", k_si, "Scale (θ):", theta_si, "\n")
cat("Incubation Period - Shape (k):", k_ip, "Scale (θ):", theta_ip, "\n")

# Monte Carlo simulation
set.seed(123)
n <- 100000  # number of samples for simulation

# Generate samples for both distributions
samples_si <- rgamma(n, shape = k_si, scale = theta_si)
samples_ip <- rgamma(n, shape = k_ip, scale = theta_ip)

# Calculate the probability that a serial interval sample is less than an incubation period sample
comparison <- samples_si < samples_ip
prob <- mean(comparison)

# Calculate 95% CI for the probability
alpha <- 0.05
z <- qnorm(1 - alpha / 2)
se_prob <- sqrt(prob * (1 - prob) / n)
ci_lower <- prob - z * se_prob
ci_upper <- prob + z * se_prob

cat("Probability that the serial interval is shorter than the incubation period:", prob, "\n")
cat("95% Confidence Interval:", ci_lower, "to", ci_upper, "\n")

# Plot the histograms of the samples with fitted Gamma distributions
par(mfrow = c(2, 1))

hist(samples_si, breaks = 30, probability = TRUE, col = "blue", xlim = c(0, max(samples_si)),
     main = "Histogram of Serial Interval Samples with Fitted Gamma Distribution",
     xlab = "Serial Interval (days)", ylab = "Density")
curve(dgamma(x, shape = k_si, scale = theta_si), col = "red", lwd = 2, add = TRUE)
legend("topright", legend = c("Serial Interval Samples", "Fitted Gamma Distribution"),
col = c("blue", "red"), lwd = 2)

hist(samples_ip, breaks = 30, probability = TRUE, col = "green", xlim = c(0, max(samples_ip)),
     main = "Histogram of Incubation Period Samples with Fitted Gamma Distribution",
     xlab = "Incubation Period (days)", ylab = "Density")
curve(dgamma(x, shape = k_ip, scale = theta_ip), col = "red", lwd = 2, add = TRUE)
# Given information for the serial interval
mean_si <- 15.40
ci_lower_si <- 13.20
ci_upper_si <- 17.50

# Given information for the incubation period
mean_ip <- 8.50
ci_lower_ip <- 7.70
ci_upper_ip <- 9.20

# Calculate standard error for both
se_si <- (ci_upper_si - mean_si) / 1.96
se_ip <- (ci_upper_ip - mean_ip) / 1.96

# Estimate the variance for both
variance_si <- se_si^2
variance_ip <- se_ip^2

# Calculate shape (k) and scale (θ) parameters for both
theta_si <- variance_si / mean_si
k_si <- mean_si^2 / variance_si

theta_ip <- variance_ip / mean_ip
k_ip <- mean_ip^2 / variance_ip

# Print the estimated parameters
cat("Serial Interval - Shape (k):", k_si, "Scale (θ):", theta_si, "\n")
cat("Incubation Period - Shape (k):", k_ip, "Scale (θ):", theta_ip, "\n")

# Monte Carlo simulation
set.seed(123)
n <- 100000  # number of samples for simulation

# Generate samples for both distributions
samples_si <- rgamma(n, shape = k_si, scale = theta_si)
samples_ip <- rgamma(n, shape = k_ip, scale = theta_ip)

# Calculate the probability that a serial interval sample is less than an incubation period sample
comparison <- samples_si < samples_ip
prob <- mean(comparison)

# Calculate 95% CI for the probability
alpha <- 0.05
z <- qnorm(1 - alpha / 2)
se_prob <- sqrt(prob * (1 - prob) / n)
ci_lower <- prob - z * se_prob
ci_upper <- prob + z * se_prob

cat("Probability that the serial interval is shorter than the incubation period:", prob, "\n")
cat("95% Confidence Interval:", ci_lower, "to", ci_upper, "\n")

# Plot the histograms of the samples with fitted Gamma distributions
par(mfrow = c(2, 1))

hist(samples_si, breaks = 30, probability = TRUE, col = "blue", xlim = c(0, max(samples_si)),
     main = "Histogram of Serial Interval Samples with Fitted Gamma Distribution",
     xlab = "Serial Interval (days)", ylab = "Density")
curve(dgamma(x, shape = k_si, scale = theta_si), col = "red", lwd = 2, add = TRUE)
legend("topright", legend = c("Serial Interval Samples", "Fitted Gamma Distribution"),
col = c("blue", "red"), lwd = 2)

hist(samples_ip, breaks = 30, probability = TRUE, col = "green", xlim = c(0, max(samples_ip)),
     main = "Histogram of Incubation Period Samples with Fitted Gamma Distribution",
     xlab = "Incubation Period (days)", ylab = "Density")
curve(dgamma(x, shape = k_ip, scale = theta_ip), col = "red", lwd = 2, add = TRUE)
legend("topright", legend = c("Incubation Period Samples", "Fitted Gamma Distribution"),
       col = c("green", "red"), lwd = 2)

###################################
#MArburg 

# Define the Gamma parameters for the serial interval
gamma_mean <- 9.2
gamma_sd <- 4.4

# Convert mean and standard deviation to shape and rate for the Gamma distribution
gamma_shape <- (gamma_mean / gamma_sd)^2
gamma_rate <- gamma_mean / (gamma_sd^2)

# Define the Weibull parameters for the incubation period
weibull_mean <- 6.9
weibull_sd <- 3.2

# Approximate the shape and scale for Weibull distribution using method of moments
weibull_shape <- (weibull_sd / weibull_mean)^(-1.086)
weibull_scale <- weibull_mean / gamma(1 + 1 / weibull_shape)

# Generate samples from the Gamma and Weibull distributions
set.seed(123)  # for reproducibility
n <- 100000  # number of samples for Monte Carlo simulation
gamma_samples <- rgamma(n, shape = gamma_shape, rate = gamma_rate)
weibull_samples <- rweibull(n, shape = weibull_shape, scale = weibull_scale)

# Calculate the probability that the Gamma sample is less than the Weibull sample
probability <- mean(gamma_samples < weibull_samples)

# Convert the probability to a percentage
probability_percent <- probability * 100

# Function to calculate the probability for a bootstrap sample
calc_prob <- function(data, indices) {
  gamma_samples_boot <- gamma_samples[indices]
  weibull_samples_boot <- weibull_samples[indices]
  return(mean(gamma_samples_boot < weibull_samples_boot))
}

# Perform bootstrap resampling
set.seed(123)
n_bootstrap <- 1000
bootstrap_results <- boot(data = 1:n, statistic = calc_prob, R = n_bootstrap)

# Calculate the 95% CI
ci <- boot.ci(bootstrap_results, type = "perc")

# Convert the CI to percentages
lower_ci <- ci$percent[4] * 100
upper_ci <- ci$percent[5] * 100

# Print the probability and 95% CI
print(paste("The probability that the serial interval is shorter than the incubation period is", round(probability_percent, 2), "%"))
print(paste("The 95% CI is", round(lower_ci, 2), "% to", round(upper_ci, 2), "%"))

############

#Lassa fever 
# Install and load the necessary packages
if (!require(ggplot2)) install.packages("ggplot2", dependencies=TRUE)
if (!require(boot)) install.packages("boot", dependencies=TRUE)
library(ggplot2)
library(boot)

# Define the Gamma parameters for the serial interval
gamma_mean <- 7.80
gamma_sd <- 10.7

# Convert mean and standard deviation to shape and rate for the Gamma distribution
gamma_shape <- (gamma_mean / gamma_sd)^2
gamma_rate <- gamma_mean / (gamma_sd^2)

# Define the Log-normal parameters for the incubation period
lognorm_meanlog <- log(12.80^2 / sqrt(4.80^2 + 12.80^2))
lognorm_sdlog <- sqrt(log(1 + (4.80^2 / 12.80^2)))

# Function to calculate the probability for a bootstrap sample
calc_prob <- function(data, indices) {
  gamma_samples <- rgamma(length(indices), shape = gamma_shape, rate = gamma_rate)
  lognorm_samples <- rlnorm(length(indices), meanlog = lognorm_meanlog, sdlog = lognorm_sdlog)
  return(mean(gamma_samples < lognorm_samples))
}

# Number of bootstrap resamples
n_bootstrap <- 1000

# Perform the bootstrap
set.seed(123)
bootstrap_results <- boot(data = 1:n_bootstrap, statistic = calc_prob, R = n_bootstrap)

# Calculate the 95% CI
ci <- boot.ci(bootstrap_results, type = "perc")

# Calculate the mean probability
probability_percent <- mean(bootstrap_results$t) * 100
lower_ci <- ci$percent[4] * 100
upper_ci <- ci$percent[5] * 100

# Print the probability and 95% CI
print(paste("The probability that the serial interval is shorter than the incubation period is", round(probability_percent, 2), "%"))
print(paste("The 95% CI is", round(lower_ci, 2), "% to", round(upper_ci, 2), "%"))


#########
#MERS
# Install and load the necessary packages
if (!require(ggplot2)) install.packages("ggplot2", dependencies=TRUE)
if (!require(boot)) install.packages("boot", dependencies=TRUE)
library(ggplot2)
library(boot)

# Calculate parameters for the log-normal distribution of the serial interval
serial_meanlog <- log(7.60)
serial_sdlog <- (log(23.10) - log(2.50)) / (2 * 1.96)

# Calculate parameters for the log-normal distribution of the incubation period
incubation_median <- 5.2
incubation_sdlog <- (log(14.7) - log(1.9)) / (2 * 1.96)
incubation_meanlog <- log(incubation_median)

# Function to calculate the probability for a bootstrap sample
calc_prob <- function(data, indices) {
  serial_samples <- rlnorm(length(indices), meanlog = serial_meanlog, sdlog = serial_sdlog)
  incubation_samples <- rlnorm(length(indices), meanlog = incubation_meanlog, sdlog = incubation_sdlog)
  return(mean(serial_samples < incubation_samples))
}

# Number of bootstrap resamples
n_bootstrap <- 1000

# Perform the bootstrap
set.seed(123)
bootstrap_results <- boot(data = 1:n_bootstrap, statistic = calc_prob, R = n_bootstrap)

# Calculate the 95% CI
ci <- boot.ci(bootstrap_results, type = "perc")

# Calculate the mean probability
probability_percent <- mean(bootstrap_results$t) * 100
lower_ci <- ci$percent[4] * 100
upper_ci <- ci$percent[5] * 100

# Print the probability and 95% CI
print(paste("The probability that the serial interval is shorter than the incubation period is", round(probability_percent, 2), "%"))
print(paste("The 95% CI is", round(lower_ci, 2), "% to", round(upper_ci, 2), "%"))

# Generate samples for plotting
set.seed(123)
n <- 10000  # number of samples for plotting
serial_samples <- rlnorm(n, meanlog = serial_meanlog, sdlog = serial_sdlog)
incubation_samples <- rlnorm(n, meanlog = incubation_meanlog, sdlog = incubation_sdlog)

# Create data frames for plotting
serial_df <- data.frame(value = serial_samples, distribution = "Serial Interval (Log-normal)")
incubation_df <- data.frame(value = incubation_samples, distribution = "Incubation Period (Log-normal)")

# Combine the data frames
combined_df <- rbind(serial_df, incubation_df)

# Plot the distributions
p <- ggplot(combined_df, aes(x = value, fill = distribution)) +
  geom_density(alpha = 0.5) +
  labs(title = "Comparison of Serial Interval and Incubation Period Distributions",
       x = "Days",
       y = "Density",
       fill = "Distribution") +
  theme_minimal()

# Add text to the plot for the estimated probability
p + annotate("text", x = 40, y = 0.03, 
             label = paste("P(Serial Interval < Incubation Period) =", round(probability_percent, 2), "%"), 
             size = 5, color = "black")

############
#SARS

# Install and load the necessary packages
if (!require(ggplot2)) install.packages("ggplot2", dependencies=TRUE)
if (!require(boot)) install.packages("boot", dependencies=TRUE)
library(ggplot2)
library(boot)

# Calculate parameters for the Weibull distribution of the serial interval
weibull_shape <- (3.8 / 8.4)^(-1.086)
weibull_scale <- 8.4 / gamma(1 + 1 / weibull_shape)

# Calculate parameters for the log-normal distribution of the incubation period
incubation_median <- 4.00
incubation_sdlog <- (log(4.40) - log(3.60)) / (2 * 1.96)
incubation_meanlog <- log(incubation_median)

# Function to calculate the probability for a bootstrap sample
calc_prob <- function(data, indices) {
  weibull_samples <- rweibull(length(indices), shape = weibull_shape, scale = weibull_scale)
  lognorm_samples <- rlnorm(length(indices), meanlog = incubation_meanlog, sdlog = incubation_sdlog)
  return(mean(weibull_samples < lognorm_samples))
}

# Number of bootstrap resamples
n_bootstrap <- 1000

# Perform the bootstrap
set.seed(123)
bootstrap_results <- boot(data = 1:n_bootstrap, statistic = calc_prob, R = n_bootstrap)

# Calculate the 95% CI
ci <- boot.ci(bootstrap_results, type = "perc")

# Calculate the mean probability
probability_percent <- mean(bootstrap_results$t) * 100
lower_ci <- ci$percent[4] * 100
upper_ci <- ci$percent[5] * 100

# Print the probability and 95% CI
print(paste("The probability that the serial interval is shorter than the incubation period is", round(probability_percent, 2), "%"))
print(paste("The 95% CI is", round(lower_ci, 2), "% to", round(upper_ci, 2), "%"))

# Generate samples for plotting
set.seed(123)
n <- 10000  # number of samples for plotting
weibull_samples <- rweibull(n, shape = weibull_shape, scale = weibull_scale)
lognorm_samples <- rlnorm(n, meanlog = incubation_meanlog, sdlog = incubation_sdlog)

# Create data frames for plotting
weibull_df <- data.frame(value = weibull_samples, distribution = "Serial Interval (Weibull)")
lognorm_df <- data.frame(value = lognorm_samples, distribution = "Incubation Period (Log-normal)")

# Combine the data frames
combined_df <- rbind(weibull_df, lognorm_df)

# Plot the distributions
p <- ggplot(combined_df, aes(x = value, fill = distribution)) +
  geom_density(alpha = 0.5) +
  labs(title = "Comparison of Serial Interval and Incubation Period Distributions",
       x = "Days",
       y = "Density",
       fill = "Distribution") +
  theme_minimal()

# Add text to the plot for the estimated probability
p + annotate("text", x = 20, y = 0.05, 
             label = paste("P(Serial Interval < Incubation Period) =", round(probability_percent, 2), "%"), 
             size = 5, color = "black")

####
#Zika 

# Install and load the necessary packages
if (!require(ggplot2)) install.packages("ggplot2", dependencies=TRUE)
if (!require(boot)) install.packages("boot", dependencies=TRUE)
library(ggplot2)
library(boot)

# Calculate parameters for the log-normal distribution of the serial interval
serial_median <- 12
serial_iqr <- c(10, 14.5)
serial_sdlog <- (log(serial_iqr[2]) - log(serial_iqr[1])) / (2 * 0.674)
serial_meanlog <- log(serial_median)

# Calculate parameters for the log-normal distribution of the incubation period
incubation_median <- 5.9
incubation_sdlog <- (log(7.6) - log(4.4)) / (2 * 1.96)
incubation_meanlog <- log(incubation_median)

# Function to calculate the probability for a bootstrap sample
calc_prob <- function(data, indices) {
  serial_samples <- rlnorm(length(indices), meanlog = serial_meanlog, sdlog = serial_sdlog)
  lognorm_samples <- rlnorm(length(indices), meanlog = incubation_meanlog, sdlog = incubation_sdlog)
  return(mean(serial_samples < lognorm_samples))
}

# Number of bootstrap resamples
n_bootstrap <- 1000

# Perform the bootstrap
set.seed(123)
bootstrap_results <- boot(data = 1:n_bootstrap, statistic = calc_prob, R = n_bootstrap)

# Calculate the 95% CI
ci <- boot.ci(bootstrap_results, type = "perc")

# Calculate the mean probability
probability_percent <- mean(bootstrap_results$t) * 100
lower_ci <- ci$percent[4] * 100
upper_ci <- ci$percent[5] * 100

# Print the probability and 95% CI
print(paste("The probability that the serial interval is shorter than the incubation period is", round(probability_percent, 2), "%"))
print(paste("The 95% CI is", round(lower_ci, 2), "% to", round(upper_ci, 2), "%"))

# Generate samples for plotting
set.seed(123)
n <- 10000  # number of samples for plotting
serial_samples <- rlnorm(n, meanlog = serial_meanlog, sdlog = serial_sdlog)
lognorm_samples <- rlnorm(n, meanlog = incubation_meanlog, sdlog = incubation_sdlog)

# Create data frames for plotting
serial_df <- data.frame(value = serial_samples, distribution = "Serial Interval (Log-normal)")
lognorm_df <- data.frame(value = lognorm_samples, distribution = "Incubation Period (Log-normal)")

# Combine the data frames
combined_df <- rbind(serial_df, lognorm_df)

# Plot the distributions
p <- ggplot(combined_df, aes(x = value, fill = distribution)) +
  geom_density(alpha = 0.5) +
  labs(title = "Comparison of Serial Interval and Incubation Period Distributions",
       x = "Days",
       y = "Density",
       fill = "Distribution") +
  theme_minimal()

# Add text to the plot for the estimated probability
p + annotate("text", x = 30, y = 0.05, 
             label = paste("P(Serial Interval < Incubation Period) =", round(probability_percent, 2), "%"), 
             size = 5, color = "black")


#H1N1
# Install and load the necessary packages
if (!require(ggplot2)) install.packages("ggplot2", dependencies=TRUE)
if (!require(boot)) install.packages("boot", dependencies=TRUE)
library(ggplot2)
library(boot)

# Calculate parameters for the Gamma distribution of the serial interval
serial_mean <- 2.4
serial_sd <- 1  # We don't have the SD, assume SD = 1 for estimation
serial_shape <- (serial_mean / serial_sd)^2
serial_rate <- serial_mean / (serial_sd^2)

# Calculate parameters for the Gamma distribution of the incubation period
incubation_mean <- 1.65
incubation_ci <- c(1.41, 1.89)
incubation_sd <- (incubation_ci[2] - incubation_ci[1]) / (2 * 1.96)
incubation_shape <- (incubation_mean / incubation_sd)^2
incubation_rate <- incubation_mean / (incubation_sd^2)

# Function to calculate the probability for a bootstrap sample
calc_prob <- function(data, indices) {
  serial_samples <- rgamma(length(indices), shape = serial_shape, rate = serial_rate)
  incubation_samples <- rgamma(length(indices), shape = incubation_shape, rate = incubation_rate)
  return(mean(serial_samples < incubation_samples))
}

# Number of bootstrap resamples
n_bootstrap <- 1000

# Perform the bootstrap
set.seed(123)
bootstrap_results <- boot(data = 1:n_bootstrap, statistic = calc_prob, R = n_bootstrap)

# Calculate the 95% CI
ci <- boot.ci(bootstrap_results, type = "perc")

# Calculate the mean probability
probability_percent <- mean(bootstrap_results$t) * 100
lower_ci <- ci$percent[4] * 100
upper_ci <- ci$percent[5] * 100

# Print the probability and 95% CI
print(paste("The probability that the serial interval is shorter than the incubation period is", round(probability_percent, 2), "%"))
print(paste("The 95% CI is", round(lower_ci, 2), "% to", round(upper_ci, 2), "%"))

# Generate samples for plotting
set.seed(123)
n <- 10000  # number of samples for plotting
serial_samples <- rgamma(n, shape = serial_shape, rate = serial_rate)
incubation_samples <- rgamma(n, shape = incubation_shape, rate = incubation_rate)

# Create data frames for plotting
serial_df <- data.frame(value = serial_samples, distribution = "Serial Interval (Gamma)")
incubation_df <- data.frame(value = incubation_samples, distribution = "Incubation Period (Gamma)")

# Combine the data frames
combined_df <- rbind(serial_df, incubation_df)

# Plot the distributions
p <- ggplot(combined_df, aes(x = value, fill = distribution)) +
  geom_density(alpha = 0.5) +
  labs(title = "Comparison of Serial Interval and Incubation Period Distributions",
       x = "Days",
       y = "Density",
       fill = "Distribution") +
  theme_minimal()

# Add text to the plot for the estimated probability
p + annotate("text", x = 4, y = 0.6, 
             label = paste("P(Serial Interval < Incubation Period) =", round(probability_percent, 2), "%"), 
             size = 5, color = "black")

#################
#H3N2
# Install and load the necessary packages
if (!require(ggplot2)) install.packages("ggplot2", dependencies=TRUE)
if (!require(boot)) install.packages("boot", dependencies=TRUE)
library(ggplot2)
library(boot)

# Calculate parameters for the Weibull distribution of the serial interval
serial_mean <- 3.60
serial_sd <- 1.60
serial_shape <- (serial_sd / serial_mean)^(-1.086)
serial_scale <- serial_mean / gamma(1 + 1 / serial_shape)

# Calculate parameters for the Weibull distribution of the incubation period
incubation_mean <- 1.48
incubation_sd <- 0.47
incubation_shape <- (incubation_sd / incubation_mean)^(-1.086)
incubation_scale <- incubation_mean / gamma(1 + 1 / incubation_shape)

# Function to calculate the probability for a bootstrap sample
calc_prob <- function(data, indices) {
  serial_samples <- rweibull(length(indices), shape = serial_shape, scale = serial_scale)
  incubation_samples <- rweibull(length(indices), shape = incubation_shape, scale = incubation_scale)
  return(mean(serial_samples < incubation_samples))
}

# Number of bootstrap resamples
n_bootstrap <- 1000

# Perform the bootstrap
set.seed(123)
bootstrap_results <- boot(data = 1:n_bootstrap, statistic = calc_prob, R = n_bootstrap)

# Calculate the 95% CI
ci <- boot.ci(bootstrap_results, type = "perc")

# Calculate the mean probability
probability_percent <- mean(bootstrap_results$t) * 100
lower_ci <- ci$percent[4] * 100
upper_ci <- ci$percent[5] * 100

# Print the probability and 95% CI
print(paste("The probability that the serial interval is shorter than the incubation period is", round(probability_percent, 2), "%"))
print(paste("The 95% CI is", round(lower_ci, 2), "% to", round(upper_ci, 2), "%"))

# Generate samples for plotting
set.seed(123)
n <- 10000  # number of samples for plotting
serial_samples <- rweibull(n, shape = serial_shape, scale = serial_scale)
incubation_samples <- rweibull(n, shape = incubation_shape, scale = incubation_scale)

# Create data frames for plotting
serial_df <- data.frame(value = serial_samples, distribution = "Serial Interval (Weibull)")
incubation_df <- data.frame(value = incubation_samples, distribution = "Incubation Period (Weibull)")

# Combine the data frames
combined_df <- rbind(serial_df, incubation_df)

# Plot the distributions
p <- ggplot(combined_df, aes(x = value, fill = distribution)) +
  geom_density(alpha = 0.5) +
  labs(title = "Comparison of Serial Interval and Incubation Period Distributions",
       x = "Days",
       y = "Density",
       fill = "Distribution") +
  theme_minimal()

# Add text to the plot for the estimated probability
p + annotate("text", x = 5, y = 0.5, 
             label = paste("P(Serial Interval < Incubation Period) =", round(probability_percent, 2), "%"), 
             size = 5, color = "black")

####################
#COVID WT

# Install and load the necessary packages
if (!require(ggplot2)) install.packages("ggplot2", dependencies=TRUE)
if (!require(boot)) install.packages("boot", dependencies=TRUE)
library(ggplot2)
library(boot)

# Define the normal distribution parameters for the serial interval
serial_mean <- 4.6
serial_sd <- 4.4

# Define the log-normal distribution parameters for the incubation period
incubation_meanlog <- log(4.4^2 / sqrt(1.8^2 + 4.4^2))
incubation_sdlog <- sqrt(log(1 + (1.8^2 / 4.4^2)))

# Function to calculate the probability for a bootstrap sample
calc_prob <- function(data, indices) {
  serial_samples <- rnorm(length(indices), mean = serial_mean, sd = serial_sd)
  lognorm_samples <- rlnorm(length(indices), meanlog = incubation_meanlog, sdlog = incubation_sdlog)
  return(mean(serial_samples < lognorm_samples))
}

# Number of bootstrap resamples
n_bootstrap <- 1000

# Perform the bootstrap
set.seed(123)
bootstrap_results <- boot(data = 1:n_bootstrap, statistic = calc_prob, R = n_bootstrap)

# Calculate the 95% CI
ci <- boot.ci(bootstrap_results, type = "perc")

# Calculate the mean probability
probability_percent <- mean(bootstrap_results$t) * 100
lower_ci <- ci$percent[4] * 100
upper_ci <- ci$percent[5] * 100

# Print the probability and 95% CI
print(paste("The probability that the serial interval is shorter than the incubation period is", round(probability_percent, 2), "%"))
print(paste("The 95% CI is", round(lower_ci, 2), "% to", round(upper_ci, 2), "%"))

# Generate samples for plotting
set.seed(123)
n <- 10000  # number of samples for plotting
serial_samples <- rnorm(n, mean = serial_mean, sd = serial_sd)
lognorm_samples <- rlnorm(n, meanlog = incubation_meanlog, sdlog = incubation_sdlog)

# Create data frames for plotting
serial_df <- data.frame(value = serial_samples, distribution = "Serial Interval (Normal)")
lognorm_df <- data.frame(value = lognorm_samples, distribution = "Incubation Period (Log-normal)")

# Combine the data frames
combined_df <- rbind(serial_df, lognorm_df)

# Plot the distributions
p <- ggplot(combined_df, aes(x = value, fill = distribution)) +
  geom_density(alpha = 0.5) +
  labs(title = "Comparison of Serial Interval and Incubation Period Distributions",
       x = "Days",
       y = "Density",
       fill = "Distribution") +
  theme_minimal()

# Add text to the plot for the estimated probability
p + annotate("text", x = 10, y = 0.1, 
             label = paste("P(Serial Interval < Incubation Period) =", round(probability_percent, 2), "%"), 
             size = 5, color = "black")

########################
#COVID Alpha 

# Install and load the necessary packages
if (!require(ggplot2)) install.packages("ggplot2", dependencies=TRUE)
if (!require(boot)) install.packages("boot", dependencies=TRUE)
library(ggplot2)
library(boot)

# Calculate parameters for the Gamma distribution of the serial interval
serial_mean <- 4.5
serial_ci <- c(4.46, 4.54)
serial_sd <- (serial_ci[2] - serial_ci[1]) / (2 * 1.96)  # Approximate SD from CI
serial_shape <- (serial_mean / serial_sd)^2
serial_rate <- serial_mean / (serial_sd^2)

# Calculate parameters for the Gamma distribution of the incubation period
incubation_mean <- 4.96
incubation_sd <- 2.32
incubation_shape <- (incubation_mean / incubation_sd)^2
incubation_rate <- incubation_mean / (incubation_sd^2)

# Function to calculate the probability for a bootstrap sample
calc_prob <- function(data, indices) {
  serial_samples <- rgamma(length(indices), shape = serial_shape, rate = serial_rate)
  incubation_samples <- rgamma(length(indices), shape = incubation_shape, rate = incubation_rate)
  return(mean(serial_samples < incubation_samples))
}

# Number of bootstrap resamples
n_bootstrap <- 1000

# Perform the bootstrap
set.seed(123)
bootstrap_results <- boot(data = 1:n_bootstrap, statistic = calc_prob, R = n_bootstrap)

# Calculate the 95% CI
ci <- boot.ci(bootstrap_results, type = "perc")

# Calculate the mean probability
probability_percent <- mean(bootstrap_results$t) * 100
lower_ci <- ci$percent[4] * 100
upper_ci <- ci$percent[5] * 100

# Print the probability and 95% CI
print(paste("The probability that the serial interval is shorter than the incubation period is", round(probability_percent, 2), "%"))
print(paste("The 95% CI is", round(lower_ci, 2), "% to", round(upper_ci, 2), "%"))

# Generate samples for plotting
set.seed(123)
n <- 10000  # number of samples for plotting
serial_samples <- rgamma(n, shape = serial_shape, rate = serial_rate)
incubation_samples <- rgamma(n, shape = incubation_shape, rate = incubation_rate)

# Create data frames for plotting
serial_df <- data.frame(value = serial_samples, distribution = "Serial Interval (Gamma)")
incubation_df <- data.frame(value = incubation_samples, distribution = "Incubation Period (Gamma)")

# Combine the data frames
combined_df <- rbind(serial_df, incubation_df)

# Plot the distributions
p <- ggplot(combined_df, aes(x = value, fill = distribution)) +
  geom_density(alpha = 0.5) +
  labs(title = "Comparison of Serial Interval and Incubation Period Distributions",
       x = "Days",
       y = "Density",
       fill = "Distribution") +
  theme_minimal()

# Add text to the plot for the estimated probability
p + annotate("text", x = 10, y = 0.1, 
             label = paste("P(Serial Interval < Incubation Period) =", round(probability_percent, 2), "%"), 
             size = 5, color = "black")

##########
#COVID-19 detla 

# Install and load the necessary packages
if (!require(ggplot2)) install.packages("ggplot2", dependencies=TRUE)
if (!require(boot)) install.packages("boot", dependencies=TRUE)
library(ggplot2)
library(boot)

# Calculate parameters for the Gamma distribution of the serial interval
serial_mean <- 4.19
serial_ci <- c(4.16, 4.22)
serial_sd <- (serial_ci[2] - serial_ci[1]) / (2 * 1.96)  # Approximate SD from CI
serial_shape <- (serial_mean / serial_sd)^2
serial_rate <- serial_mean / (serial_sd^2)

# Calculate parameters for the Gamma distribution of the incubation period
incubation_mean <- 4.30
incubation_sd <- 2.40
incubation_shape <- (incubation_mean / incubation_sd)^2
incubation_rate <- incubation_mean / (incubation_sd^2)

# Function to calculate the probability for a bootstrap sample
calc_prob <- function(data, indices) {
  serial_samples <- rgamma(length(indices), shape = serial_shape, rate = serial_rate)
  incubation_samples <- rgamma(length(indices), shape = incubation_shape, rate = incubation_rate)
  return(mean(serial_samples < incubation_samples))
}

# Number of bootstrap resamples
n_bootstrap <- 1000

# Perform the bootstrap
set.seed(123)
bootstrap_results <- boot(data = 1:n_bootstrap, statistic = calc_prob, R = n_bootstrap)

# Calculate the 95% CI
ci <- boot.ci(bootstrap_results, type = "perc")

# Calculate the mean probability
probability_percent <- mean(bootstrap_results$t) * 100
lower_ci <- ci$percent[4] * 100
upper_ci <- ci$percent[5] * 100

# Print the probability and 95% CI
print(paste("The probability that the serial interval is shorter than the incubation period is", round(probability_percent, 2), "%"))
print(paste("The 95% CI is", round(lower_ci, 2), "% to", round(upper_ci, 2), "%"))

# Generate samples for plotting
set.seed(123)
n <- 10000  # number of samples for plotting
serial_samples <- rgamma(n, shape = serial_shape, rate = serial_rate)
incubation_samples <- rgamma(n, shape = incubation_shape, rate = incubation_rate)

# Create data frames for plotting
serial_df <- data.frame(value = serial_samples, distribution = "Serial Interval (Gamma)")
incubation_df <- data.frame(value = incubation_samples, distribution = "Incubation Period (Gamma)")

# Combine the data frames
combined_df <- rbind(serial_df, incubation_df)

# Plot the distributions
p <- ggplot(combined_df, aes(x = value, fill = distribution)) +
  geom_density(alpha = 0.5) +
  labs(title = "Comparison of Serial Interval and Incubation Period Distributions",
       x = "Days",
       y = "Density",
       fill = "Distribution") +
  theme_minimal()

# Add text to the plot for the estimated probability
p + annotate("text", x = 10, y = 0.1, 
             label = paste("P(Serial Interval < Incubation Period) =", round(probability_percent, 2), "%"), 
             size = 5, color = "black")

#########
#COVID-19 Omicron 

# Install and load the necessary packages
if (!require(ggplot2)) install.packages("ggplot2", dependencies=TRUE)
if (!require(boot)) install.packages("boot", dependencies=TRUE)
library(ggplot2)
library(boot)

# Calculate parameters for the Gamma distribution of the serial interval
serial_mean <- 2.58
serial_ci <- c(0.63, 6.76)
serial_sd <- (serial_ci[2] - serial_ci[1]) / (2 * 1.96)  # Approximate SD from CI
serial_shape <- (serial_mean / serial_sd)^2
serial_rate <- serial_mean / (serial_sd^2)

# Calculate parameters for the log-normal distribution of the incubation period
incubation_mean <- 3.30
incubation_sd <- 2.10
incubation_meanlog <- log(incubation_mean^2 / sqrt(incubation_sd^2 + incubation_mean^2))
incubation_sdlog <- sqrt(log(1 + (incubation_sd^2 / incubation_mean^2)))

# Function to calculate the probability for a bootstrap sample
calc_prob <- function(data, indices) {
  serial_samples <- rgamma(length(indices), shape = serial_shape, rate = serial_rate)
  lognorm_samples <- rlnorm(length(indices), meanlog = incubation_meanlog, sdlog = incubation_sdlog)
  return(mean(serial_samples < lognorm_samples))
}

# Number of bootstrap resamples
n_bootstrap <- 1000

# Perform the bootstrap
set.seed(123)
bootstrap_results <- boot(data = 1:n_bootstrap, statistic = calc_prob, R = n_bootstrap)

# Calculate the 95% CI
ci <- boot.ci(bootstrap_results, type = "perc")

# Calculate the mean probability
probability_percent <- mean(bootstrap_results$t) * 100
lower_ci <- ci$percent[4] * 100
upper_ci <- ci$percent[5] * 100

# Print the probability and 95% CI
print(paste("The probability that the serial interval is shorter than the incubation period is", round(probability_percent, 2), "%"))
print(paste("The 95% CI is", round(lower_ci, 2), "% to", round(upper_ci, 2), "%"))

# Generate samples for plotting
set.seed(123)
n <- 10000  # number of samples for plotting
serial_samples <- rgamma(n, shape = serial_shape, rate = serial_rate)
lognorm_samples <- rlnorm(n, meanlog = incubation_meanlog, sdlog = incubation_sdlog)

# Create data frames for plotting
serial_df <- data.frame(value = serial_samples, distribution = "Serial Interval (Gamma)")
lognorm_df <- data.frame(value = lognorm_samples, distribution = "Incubation Period (Log-normal)")

# Combine the data frames
combined_df <- rbind(serial_df, lognorm_df)

# Plot the distributions
p <- ggplot(combined_df, aes(x = value, fill = distribution)) +
  geom_density(alpha = 0.5) +
  labs(title = "Comparison of Serial Interval and Incubation Period Distributions",
       x = "Days",
       y = "Density",
       fill = "Distribution") +
  theme_minimal()

# Add text to the plot for the estimated probability
p + annotate("text", x = 10, y = 0.1, 
             label = paste("P(Serial Interval < Incubation Period) =", round(probability_percent, 2), "%"), 
             size = 5, color = "black")

###########
#H5N1


# Install and load the necessary packages
if (!require(ggplot2)) install.packages("ggplot2", dependencies=TRUE)
if (!require(boot)) install.packages("boot", dependencies=TRUE)
library(ggplot2)
library(boot)

# Define the Log-normal parameters for the serial interval
serial_mean <- 5.648005
serial_median <- 4.028020

# Calculate meanlog and sdlog for the Log-normal distribution
serial_meanlog <- log(serial_median)
serial_sdlog <- sqrt(2 * (log(serial_mean) - serial_meanlog))

# Define the Weibull parameters for the incubation period
weibull_mean <- 3.30
weibull_ci <- c(2.70, 3.90)
weibull_sd <- 1.5

# Approximate the shape and scale for Weibull distribution using method of moments
weibull_shape <- (weibull_sd / weibull_mean)^(-1.086)
weibull_scale <- weibull_mean / gamma(1 + 1 / weibull_shape)

# Generate samples from the Log-normal and Weibull distributions
set.seed(123)  # for reproducibility
n <- 100000  # number of samples for Monte Carlo simulation
lognorm_samples <- rlnorm(n, meanlog = serial_meanlog, sdlog = serial_sdlog)
weibull_samples <- rweibull(n, shape = weibull_shape, scale = weibull_scale)

# Calculate the probability that the Log-normal sample is less than the Weibull sample
probability <- mean(lognorm_samples < weibull_samples)

# Convert the probability to a percentage
probability_percent <- probability * 100

# Function to calculate the probability for a bootstrap sample
calc_prob <- function(data, indices) {
  lognorm_samples_boot <- lognorm_samples[indices]
  weibull_samples_boot <- weibull_samples[indices]
  return(mean(lognorm_samples_boot < weibull_samples_boot))
}

# Perform bootstrap resampling
set.seed(123)
n_bootstrap <- 1000
bootstrap_results <- boot(data = 1:n, statistic = calc_prob, R = n_bootstrap)

# Calculate the 95% CI
ci <- boot.ci(bootstrap_results, type = "perc")

# Convert the CI to percentages
lower_ci <- ci$percent[4] * 100
upper_ci <- ci$percent[5] * 100

# Print the probability and 95% CI
print(paste("The probability that the serial interval is shorter than the incubation period is", round(probability_percent, 2), "%"))
print(paste("The 95% CI is", round(lower_ci, 2), "% to", round(upper_ci, 2), "%"))

# Generate samples for plotting
set.seed(123)
n <- 10000  # number of samples for plotting
lognorm_samples <- rlnorm(n, meanlog = serial_meanlog, sdlog = serial_sdlog)
weibull_samples <- rweibull(n, shape = weibull_shape, scale = weibull_scale)

# Create data frames for plotting
serial_df <- data.frame(value = lognorm_samples, distribution = "Serial Interval (Log-normal)")
incubation_df <- data.frame(value = weibull_samples, distribution = "Incubation Period (Weibull)")

# Combine the data frames
combined_df <- rbind(serial_df, incubation_df)

# Plot the distributions
p <- ggplot(combined_df, aes(x = value, fill = distribution)) +
  geom_density(alpha = 0.5) +
  labs(title = "Comparison of Serial Interval and Incubation Period Distributions",
       x = "Days",
       y = "Density",
       fill = "Distribution") +
  theme_minimal()

# Add text to the plot for the estimated probability and confidence interval
p + annotate("text", x = 15, y = 0.05, 
             label = paste("P(Serial Interval < Incubation Period) =", round(probability_percent, 2), "%\n95% CI:", round(lower_ci, 2), "% to", round(upper_ci, 2), "%"), 
             size = 5, color = "black")

########################################
#Nipah 

# Install and load the necessary packages
if (!require(ggplot2)) install.packages("ggplot2", dependencies=TRUE)
if (!require(boot)) install.packages("boot", dependencies=TRUE)
library(ggplot2)
library(boot)

# Define the Gamma parameters for the serial interval
serial_mean <- 12.7
serial_sd <- 3

# Convert mean and standard deviation to shape and rate for the Gamma distribution
serial_shape <- (serial_mean / serial_sd)^2
serial_rate <- serial_mean / (serial_sd^2)

# Define the Gamma parameters for the incubation period
incubation_mean <- 9.7
incubation_sd <- 2.2

# Convert mean and standard deviation to shape and rate for the Gamma distribution
incubation_shape <- (incubation_mean / incubation_sd)^2
incubation_rate <- incubation_mean / (incubation_sd^2)

# Generate samples from the Gamma distributions
set.seed(123)  # for reproducibility
n <- 100000  # number of samples for Monte Carlo simulation
serial_samples <- rgamma(n, shape = serial_shape, rate = serial_rate)
incubation_samples <- rgamma(n, shape = incubation_shape, rate = incubation_rate)

# Calculate the probability that the serial interval is shorter than the incubation period
probability <- mean(serial_samples < incubation_samples)

# Convert the probability to a percentage
probability_percent <- probability * 100

# Function to calculate the probability for a bootstrap sample
calc_prob <- function(data, indices) {
  serial_samples_boot <- serial_samples[indices]
  incubation_samples_boot <- incubation_samples[indices]
  return(mean(serial_samples_boot < incubation_samples_boot))
}

# Perform bootstrap resampling
set.seed(123)
n_bootstrap <- 1000
bootstrap_results <- boot(data = 1:n, statistic = calc_prob, R = n_bootstrap)

# Calculate the 95% CI
ci <- boot.ci(bootstrap_results, type = "perc")

# Convert the CI to percentages
lower_ci <- ci$percent[4] * 100
upper_ci <- ci$percent[5] * 100

# Print the probability and 95% CI
print(paste("The probability that the serial interval is shorter than the incubation period is", round(probability_percent, 2), "%"))
print(paste("The 95% CI is", round(lower_ci, 2), "% to", round(upper_ci, 2), "%"))

# Generate samples for plotting
set.seed(123)
n <- 10000  # number of samples for plotting
serial_samples <- rgamma(n, shape = serial_shape, rate = serial_rate)
incubation_samples <- rgamma(n, shape = incubation_shape, rate = incubation_rate)

# Create data frames for plotting
serial_df <- data.frame(value = serial_samples, distribution = "Serial Interval (Gamma)")
incubation_df <- data.frame(value = incubation_samples, distribution = "Incubation Period (Gamma)")

# Combine the data frames
combined_df <- rbind(serial_df, incubation_df)

# Plot the distributions
p <- ggplot(combined_df, aes(x = value, fill = distribution)) +
  geom_density(alpha = 0.5) +
  labs(title = "Comparison of Serial Interval and Incubation Period Distributions",
       x = "Days",
       y = "Density",
       fill = "Distribution") +
  theme_minimal()

# Add text to the plot for the estimated probability and confidence interval
p + annotate("text", x = 25, y = 0.03, 
             label = paste("P(Serial Interval < Incubation Period) =", round(probability_percent, 2), "%\n95% CI:", round(lower_ci, 2), "% to", round(upper_ci, 2), "%"), 
             size = 5, color = "black")

################################################################
#H2N2

#fit inc to a dist
# Given data
mean <- 2.0
median <- 2.0
q25 <- 2.0
q75 <- 2.5

# Load necessary library
library(MASS)

# Generate a grid of values for plotting
x_values <- seq(1.5, 3, length.out = 1000)

# Fit Normal distribution based on quantiles
std_normal <- (q75 - q25) / (qnorm(0.75) - qnorm(0.25))
normal_fit <- dnorm(x_values, mean = mean, sd = std_normal)

# Fit Lognormal distribution
lognorm_fit <- fitdistr(c(q25, median, q75), "lognormal")
lognorm_density <- dlnorm(x_values, meanlog = lognorm_fit$estimate["meanlog"], sdlog = lognorm_fit$estimate["sdlog"])

# Fit Gamma distribution
gamma_fit <- fitdistr(c(q25, median, q75), "gamma")
gamma_density <- dgamma(x_values, shape = gamma_fit$estimate["shape"], rate = gamma_fit$estimate["rate"])

# Plot the fitted distributions
plot(x_values, normal_fit, type = "l", col = "blue", lwd = 2, ylab = "Density", xlab = "Value", main = "Fitted Distributions")
lines(x_values, lognorm_density, col = "green", lwd = 2)


# Install and load the necessary packages
library(ggplot2)
library(boot)

# Define the Gamma parameters for the serial interval
serial_mean <- 3.66
serial_median <- 3.40
serial_q1 <- 3.0
serial_q3 <- 4.0

# Approximate the shape and rate for the Gamma distribution using method of moments
serial_sd <- (serial_q3 - serial_q1) / 1.35  # Approximate SD from IQR (Interquartile Range)
serial_shape <- (serial_mean / serial_sd)^2
serial_rate <- serial_mean / (serial_sd^2)

# Define the Log-normal parameters for the incubation period
incubation_mean <- 1.34

# Calculate the meanlog and sdlog for the Log-normal distribution using the mean
incubation_meanlog <- log(incubation_mean)  # Assuming log transformation
incubation_sdlog <- 0.25  # A common assumption for log-normal distributions (adjustable if needed)

# Generate samples from the Gamma and Log-normal distributions
set.seed(123)  # for reproducibility
n <- 100000  # number of samples for Monte Carlo simulation
gamma_samples <- rgamma(n, shape = serial_shape, rate = serial_rate)
lognorm_samples <- rlnorm(n, meanlog = incubation_meanlog, sdlog = incubation_sdlog)

# Calculate the probability that the Gamma sample is less than the Log-normal sample
probability <- mean(gamma_samples < lognorm_samples)

# Convert the probability to a percentage
probability_percent <- probability * 100

# Function to calculate the probability for a bootstrap sample
calc_prob <- function(data, indices) {
  gamma_samples_boot <- gamma_samples[indices]
  lognorm_samples_boot <- lognorm_samples[indices]
  return(mean(gamma_samples_boot < lognorm_samples_boot))
}

# Perform bootstrap resampling
set.seed(123)
n_bootstrap <- 1000
bootstrap_results <- boot(data = 1:n, statistic = calc_prob, R = n_bootstrap)

# Calculate the 95% CI
ci <- boot.ci(bootstrap_results, type = "perc")

# Convert the CI to percentages
lower_ci <- ci$percent[4] * 100
upper_ci <- ci$percent[5] * 100

# Print the probability and 95% CI
print(paste("The probability that the serial interval is shorter than the incubation period is", round(probability_percent, 2), "%"))
print(paste("The 95% CI is", round(lower_ci, 2), "% to", round(upper_ci, 2), "%"))

# Generate samples for plotting
set.seed(123)
n <- 10000  # number of samples for plotting
gamma_samples <- rgamma(n, shape = serial_shape, rate = serial_rate)
lognorm_samples <- rlnorm(n, meanlog = incubation_meanlog, sdlog = incubation_sdlog)

# Create data frames for plotting
serial_df <- data.frame(value = gamma_samples, distribution = "Serial Interval (Gamma)")
incubation_df <- data.frame(value = lognorm_samples, distribution = "Incubation Period (Log-normal)")

# Combine the data frames
combined_df <- rbind(serial_df, incubation_df)

# Plot the distributions
p <- ggplot(combined_df, aes(x = value, fill = distribution)) +
  geom_density(alpha = 0.5) +
  labs(title = "Comparison of Serial Interval and Incubation Period Distributions",
       x = "Days",
       y = "Density",
       fill = "Distribution") +
  theme_minimal()

# Add text to the plot for the estimated probability and confidence interval
p + annotate("text", x = 8, y = 0.15, 
             label = paste("P(Serial Interval < Incubation Period) =", round(probability_percent, 2), "%\n95% CI:", round(lower_ci, 2), "% to", round(upper_ci, 2), "%"), 
             size = 5, color = "black")


#######################################
#H2N2

# Install and load the necessary packages
library(ggplot2)
library(boot)

# Define the Normal distribution parameters for the serial interval
serial_mean <- 3.45
serial_median <- 3.50
serial_q1 <- 3.0
serial_q3 <- 4.0

# Approximate the standard deviation from the interquartile range (IQR)
serial_sd <- (serial_q3 - serial_q1) / 1.35  # Approximate SD from IQR for Normal distribution

# Define the Log-normal parameters for the incubation period
incubation_mean <- 2.0
incubation_median <- 2.0
incubation_q1 <- 2.0
incubation_q3 <- 2.50

# Calculate the meanlog and sdlog for the Log-normal distribution
incubation_meanlog <- log(incubation_median)  # Assuming log transformation of the median
incubation_sdlog <- (log(incubation_q3) - log(incubation_q1)) / (2 * 0.674)  # Approximate sdlog using IQR for Log-normal

# Generate samples from the Normal and Log-normal distributions
set.seed(123)  # for reproducibility
n <- 100000  # number of samples for Monte Carlo simulation
normal_samples <- rnorm(n, mean = serial_mean, sd = serial_sd)
lognorm_samples <- rlnorm(n, meanlog = incubation_meanlog, sdlog = incubation_sdlog)

# Calculate the probability that the Normal sample is less than the Log-normal sample
probability <- mean(normal_samples < lognorm_samples)

# Convert the probability to a percentage
probability_percent <- probability * 100

# Function to calculate the probability for a bootstrap sample
calc_prob <- function(data, indices) {
  normal_samples_boot <- normal_samples[indices]
  lognorm_samples_boot <- lognorm_samples[indices]
  return(mean(normal_samples_boot < lognorm_samples_boot))
}

# Perform bootstrap resampling
set.seed(123)
n_bootstrap <- 1000
bootstrap_results <- boot(data = 1:n, statistic = calc_prob, R = n_bootstrap)

# Calculate the 95% CI
ci <- boot.ci(bootstrap_results, type = "perc")

# Convert the CI to percentages
lower_ci <- ci$percent[4] * 100
upper_ci <- ci$percent[5] * 100

# Print the probability and 95% CI
print(paste("The probability that the serial interval is shorter than the incubation period is", round(probability_percent, 2), "%"))
print(paste("The 95% CI is", round(lower_ci, 2), "% to", round(upper_ci, 2), "%"))

# Generate samples for plotting
set.seed(123)
n <- 10000  # number of samples for plotting
normal_samples <- rnorm(n, mean = serial_mean, sd = serial_sd)
lognorm_samples <- rlnorm(n, meanlog = incubation_meanlog, sdlog = incubation_sdlog)

# Create data frames for plotting
serial_df <- data.frame(value = normal_samples, distribution = "Serial Interval (Normal)")
incubation_df <- data.frame(value = lognorm_samples, distribution = "Incubation Period (Log-normal)")

# Combine the data frames
combined_df <- rbind(serial_df, incubation_df)

# Plot the distributions
p <- ggplot(combined_df, aes(x = value, fill = distribution)) +
  geom_density(alpha = 0.5) +
  labs(title = "Comparison of Serial Interval and Incubation Period Distributions",
       x = "Days",
       y = "Density",
       fill = "Distribution") +
  theme_minimal()

# Add text to the plot for the estimated probability and confidence interval
p + annotate("text", x = 5, y = 0.5, 
             label = paste("P(Serial Interval < Incubation Period) =", round(probability_percent, 2), "%\n95% CI:", round(lower_ci, 2), "% to", round(upper_ci, 2), "%"), 
             size = 5, color = "black")



#######################################

#Figure
library(readxl)
setwd("/Users/lshjw6/Documents/blueprint/epi_review/params")
asym <- read_excel("bigboi.xlsx", sheet = "asym", 
                     col_types = c("text", "text", "text", 
                                   "text", "text", "text", "numeric", 
                                   "numeric", "numeric", "text", "text", 
                                   "text", "text"))
View(asym)

library(ggplot2)

# Create plot
ggplot(asym, aes(x = pathogen)) +
  geom_point(aes(y = prop_presym), color = "red", size = 3) +  # Dot for prop_presym
  geom_errorbar(aes(ymin = lower_CI, ymax = upper_CI), width = 0.2, color = "black") +  # Error bars for CI
  labs(x = "Pathogen", y = "Proportion Presymptomatic (%)") +
  theme_classic() +
  coord_flip() 

ggplot(asym, aes(x = pathogen)) +
  geom_point(aes(y = prop_presym), color = "red", size = 3) +  # Dot for prop_presym
  geom_errorbar(aes(ymin = lower_CI, ymax = upper_CI, linetype = "95% CI"), width = 0.2, color = "black") +  # Error bars for CI
  labs(x = "Pathogen", y = "Proportion Presymptomatic (%)", linetype = "Legend") +
  theme_classic() +  
  coord_flip() +  
  theme(
    axis.title.x = element_text(size = 14),  # Increase X axis label size
    axis.title.y = element_text(size = 14),  # Increase Y axis label size
    axis.text = element_text(size = 12),  # Increase axis text size
    legend.title = element_text(size = 14),  # Increase legend title size
    legend.text = element_text(size = 12)  # Increase legend text size
  )

#table
library(dplyr)
install.packages("gt")
library(gt)

asym_cleaned <- asym %>%
  mutate_all(~ifelse(is.na(.), "-", .))

selected_columns <- asym_cleaned %>%
  select(pathogen, asym_source, inc_source, SI_source, inc_dist, SI_dist, Notes) %>%
  rename(
    Pathogen = pathogen,
    `Asymptomatic data source` = asym_source,
    `Incubation period source` = inc_source,
    `Serial interval source` = SI_source,
    `Incubation distribution` = inc_dist,
    `Serial interval distribution` = SI_dist,
    Notes = Notes
  )

table_figure <- selected_columns %>%
  gt() %>%
  tab_header(
    title = "pre-symptomatic transmission proportion") %>%
  fmt_missing(columns = everything(), missing_text = "N/A") %>%
  cols_width(
    everything() ~ px(200)  # Adjust column widths if necessary
  ) %>%
  tab_options(
    table.font.size = 12,  # Adjust the font size
    heading.title.font.size = 14,
    heading.subtitle.font.size = 12
  )

# Display the table
print(table_figure)
gtsave(table_figure, "presym_table.html")




