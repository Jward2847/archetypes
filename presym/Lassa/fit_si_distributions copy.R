# if primarycensored is not installed, install the development version from the
# repository using this command:
# devtools::install_github("epinowcast/primarycensored")
library(data.table)
library(primarycensored) # devtools::install_github("epinowcast/primarycensored")
library(ggplot2)
library(loo)
library(tidybayes)
library(cmdstanr)

# Loading in data
data_files <- list(
  "lassa" = "presym/Lassa/lassa_ob.csv")

#Load both sets in
data_files <- list(
  "lassa_1" = "presym/Lassa/lassa_ob.csv",
  "lassa_2" = "presym/Lassa/lassa_ll.csv"
)



data_tables <- lapply(data_files, fread)

# Combining data.tables for panels from both lassa datasets
# Ensure merging works correctly by explicitly joining on "days"
dt_onsets <- merge(
  data_tables$lassa_1[, .(days, onsets_lassa_1 = count)],
  data_tables$lassa_2[, .(days, onsets_lassa_2 = count)],
  by = "days",      # Merge on the "days" column
  all = TRUE        # Keep all rows even if one table lacks certain days
)[
  , onsets := fifelse(is.na(onsets_lassa_1), 0, onsets_lassa_1) + 
    fifelse(is.na(onsets_lassa_2), 0, onsets_lassa_2)
][
  , .(days, onsets) # Select final columns
]



# one single data set dataset creation
dt_onsets <- data_tables$lassa[, .(days, onsets = count)]

## upper limits
dt_onsets[, days_upper := days + 1]
## primary window
dt_onsets[, pwindow := 1]
## observation time
dt_onsets[, relative_obs_time := Inf]

# Calculate the intervals for fitting the serial interval distributions where
# Interpret the “onsets” as frequency counts. I.e., the raw data provides counts
# of the number of cases that showed symptoms on specific days.
# We convert this into a format representing each interval explicitly
dt_onsets_intervals <- data.table(
  time =  0:length(dt_onsets$days),
  onsets = rep(dt_onsets$days, dt_onsets$onsets))

# This compiles the primarycensored Stan model
pcd_model <- pcd_cmdstan_model()

# Fit the model using fitdistcens with censored gamma distribution
# See documentation for primarycensored to find which dist_ids refer to which parametric distributions
# Here: https://cran.r-project.org/web/packages/primarycensored/primarycensored.pdf

# dist_id = 2 -> Gamma distribution
gamma_data <- pcd_as_stan_data(
  dt_onsets,
  delay = "days",
  delay_upper = "days_upper",
  n = "onsets",
  dist_id = 2, # This is where the parametric form of the distribution being fit is specified. I.e., 2 = Gamma
  primary_id = 1,
  param_bounds = list(lower = c(0, 0), upper = c(Inf, Inf)),
  primary_param_bounds = list(lower = numeric(0), upper = numeric(0)),
  priors = list(location = c(1, 1), scale = c(1, 1)),
  primary_priors = list(location = numeric(0), scale = numeric(0)),
  compute_log_lik = TRUE
)

# dist_id = 1 -> Lognormal distribution
ln_data <- pcd_as_stan_data(
  dt_onsets,
  delay = "days",
  delay_upper = "days_upper",
  n = "onsets",
  dist_id = 1, # This is where the parametric form of the distribution being fit is specified. I.e., 1 = Lognormal
  primary_id = 1,
  param_bounds = list(lower = c(-Inf, 0), upper = c(Inf, Inf)),
  primary_param_bounds = list(lower = numeric(0), upper = numeric(0)),
  priors = list(location = c(1, 1), scale = c(1, 1)),
  primary_priors = list(location = numeric(0), scale = numeric(0)),
  compute_log_lik = TRUE
)

# Fit models to data
fit_lognormal <- pcd_model$sample(ln_data, chains = 4, parallel_chains = 4)
fit_gamma <- pcd_model$sample(gamma_data, chains = 4, parallel_chains = 4)

# Extract posterior samples
dt_lognormal <- as.data.table(spread_draws(fit_lognormal, params[id]))
dt_lognormal <- dcast(
  dt_lognormal, .chain + .iteration + .draw ~ id, 
  value.var = "params")

dt_lognormal[, value := rlnorm(n = .N, meanlog = `1`, sdlog = `2`)]

dt_gamma <- as.data.table(spread_draws(fit_gamma, params[id]))
dt_gamma <- dcast(dt_gamma, .chain + .iteration + .draw ~ id, value.var = "params")
dt_gamma[, value := rgamma(n = .N, shape = `1`, rate = `2`)]

# Combine distribution type
dt_si_posteriors <- rbind(
  dt_lognormal[, type := "Lognormal"],
  dt_gamma[, type := "Gamma"])

# Summarise distributions
dt_summary_median <- dt_si_posteriors[, .(
  me = median(value),
  lo = median(value) - IQR(value),
  hi = median(value) + IQR(value)),
  by = "type"][, average := "Median"]

# Summarise distributions
dt_summary_mean <- dt_si_posteriors[, .(
  me = mean(value), 
  lo = mean(value) - IQR(value),
  hi = mean(value) + IQR(value)), 
  by = "type"][, average := "Mean"]

# Creating PMF of underlying data for plotting behind fitted distributions
dt_onsets_plot <- rbind(
  data.table(time = dt_onsets_intervals$onsets, type = "Lognormal"),
  data.table(time = dt_onsets_intervals$onsets, type = "Gamma")
)

# Create the plot with observed data
p_si <- ggplot() +
  geom_histogram(
    data = dt_onsets_plot,
    aes(x = time, y = after_stat(density)),
    binwidth = 1,
    fill = "gray",
    alpha = 1,
    boundary = 0,
    closed = "left"
  ) +
  geom_density(
    data = dt_si_posteriors,
    aes(x = value, fill = type),
    alpha = 0.8
  ) +
  geom_vline(
    data = dt_summary_median,
    aes(xintercept = me),
    linetype = "dashed"
  ) +
  facet_grid(~type) +
  theme_linedraw() +
  theme(legend.position = "none") +
  labs(x = "Time (days since exposure)", y = "Density") +
  lims(x = c(0, 30))

ggsave("plots/serial_interval.png", p_si, width = 8, height = 4)

knitr::kable(
  rbind(dt_summary_median, dt_summary_mean)[order(type, average)])

#--- LOO analysis 
knitr::kable(loo_compare(fit_lognormal$loo(), fit_gamma$loo()))

# After the LOO analysis and before saving the CSV
# Summarize the parameters for both distributions
param_summary <- rbind(
  dt_lognormal[, .(
    param = c("meanlog", "sdlog"),
    median = c(median(`1`), median(`2`)),
    mean = c(mean(`1`), mean(`2`)),
    q2.5 = c(quantile(`1`, 0.025), quantile(`2`, 0.025)),
    q97.5 = c(quantile(`1`, 0.975), quantile(`2`, 0.975))
  )][, dist := "Lognormal"],
  
  dt_gamma[, .(
    param = c("shape", "rate"),
    median = c(median(`1`), median(`2`)),
    mean = c(mean(`1`), mean(`2`)),
    q2.5 = c(quantile(`1`, 0.025), quantile(`2`, 0.025)),
    q97.5 = c(quantile(`1`, 0.975), quantile(`2`, 0.975))
  )][, dist := "Gamma"]
)

# Print the parameter summary table
knitr::kable(param_summary, digits = 3)

# Save the parameter summary as CSV
fwrite(param_summary, "data/si_param_summary.csv")

# Save the serial interval posteriors as CSV
fwrite(dt_si_posteriors, "data/si_posteriors.csv")




