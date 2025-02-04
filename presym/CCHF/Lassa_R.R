---
  title: "How to estimate R and k from cluster size data using Markov chain Monte Carlo?"
format: 
  html:
  code-link: true
editor: source
editor_options: 
  chunk_output_type: console
date: last-modified
toc: true
toc_float: true
---
  
  ## Ingredients
  
  - Use Bayesian estimation methods to estimate the reproduction number ($R$) and extent of superspreading, represented by the dispersion of a negative binomial distribution for individual-level seconday cases ($k$), from data on MERS-CoV- outbreak clusters.

- We will use [Markov chain Monte Carlo (MCMC)](https://en.wikipedia.org/wiki/Markov_chain_Monte_Carlo), specifically a simple Metropolis-Hastings algorithm to estimate the parameters.


## Steps in code

```{r}
#| warning: false

install.packages("epichains")
install.packages("MCMCpack")
install.packages("coda")

# check whether {pak} is installed
if(!require("pak")) install.packages("pak")
pak::pak("epiverse-trace/epiparameter")

# Load required packages
library(epichains)
library(MCMCpack)
library(epiparameter)
library(coda)
library(ggplot2)
library(dplyr)
```

```{r}
# Define data
h5_clusters = c(rep(1,13),2)

scenario_2 <- c(rep(1, 56), 3, 2)

# Get prior for R based on Aditama et al, PLOS ONE, 2012
get_prior <- extract_param(
  type = "percentiles",
  values = c(0.009, 0.315),
  distribution = "gamma",
  percentiles = c(0.025,0.975)
)


h5_prior_r <- function(x){dgamma(x,shape = get_prior[["shape"]], scale = get_prior[["scale"]])}

# Show summary table of frequencies
freq_df <- as.data.frame(table(h5_clusters)); names(freq_df) <- c("Cluster size", "Frequency")

# Create a table for the HTML document
knitr::kable(freq_df, caption = "Frequencies of H5 Clusters")

# Define likelihood function
lik_function <- function(param) {
  if (any(param <= 0)) return(-Inf) # Ensure positive parameters
  
  # Extract values of R and k
  r_val <- as.numeric(param[1])
  k_val <- as.numeric(param[2])
  
  # Define likelihood
  log_likelihood <- likelihood(
    chains = h5_clusters,
    statistic = "size",
    offspring_dist = rnbinom,
    size = k_val,
    mu = r_val
  )
  
  # Assume non-informative priors for R and k
  log_prior <- h5_prior_r(r_val)
  
  # Return log-posterior (log-likelihood + log-prior)
  return(log_likelihood + log_prior)
}

# Define number of MCMC iterations
n_iter <- 1e4

# Define 'burn in' period for fitting, to be discarded
n_burn <- 1e3

# Initial guess for c(R,k):
init_param <- c(R=0.1, k=0.5)

# Run MCMC to estimate parameters
result_mcmcpack <- MCMCmetrop1R(lik_function, 
                                theta.init = init_param, 
                                burnin = n_burn, 
                                mcmc = n_iter, 
                                thin = 1)

# Calculate effective sample size (i.e. measure of MCMC mixing)
ess_mcmcpack <- effectiveSize(result_mcmcpack)

# Plot posterior estimates
plot(result_mcmcpack)

# Define helper function to calculate median and 95% credible interval from data.frame of MCMC samples
get_param <- function(x){
  apply(x,2,function(y){val = signif(quantile(y,c(0.5,0.025,0.975)),3);
  val_text <- paste0(val[1]," (95%: CrI: ",val[2],"-",val[3],")")})
}

# Get posterior median and 95% CrI
posterior_estimates <- get_param(result_mcmcpack)

# Compile table
results_table <- data.frame(
  Package = "MCMCpack",
  Posterior_R = posterior_estimates[1],
  Posterior_k = posterior_estimates[2],
  ESS_R = ess_mcmcpack[1],
  ESS_k = ess_mcmcpack[2]
)

# Output the table with kable
knitr::kable(results_table, caption = "MCMC Comparison Table", align = 'c')

# Compare prior and posterior
posterior_samples_R <- result_mcmcpack[, 1] 

# Set up the plotting range based on the data
x_range <- seq(min(posterior_samples_R), max(posterior_samples_R), length.out = 1000)

plot(density(posterior_samples_R), col = "blue", lwd = 2, 
     main = "Posterior and Prior for R", xlab = "R", ylab = "Density")

# Add prior distribution curve
lines(x_range, h5_prior_r(x_range), col = "red", lwd = 2, lty = 2)

# Add a legend
legend("topright", legend = c("Posterior", "Prior"), col = c("blue", "red"), lwd = 2, lty = c(1, 2))



###########################################################
#GGPlot

# Load ggplot2 library
library(ggplot2)

# Create a data frame for ggplot
posterior_density <- density(posterior_samples_R)
prior_density <- data.frame(
  x = x_range,
  y = h5_prior_r(x_range)
)

posterior_data <- data.frame(
  x = posterior_density$x,
  y = posterior_density$y
)

# Combine both densities into one data frame
density_df <- rbind(
  data.frame(x = posterior_density$x, y = posterior_density$y, Type = "Posterior"),
  data.frame(x = prior_density$x, y = prior_density$y, Type = "Prior")
)

# Plot using ggplot
ggplot(density_df, aes(x = x, y = y, color = Type, linetype = Type)) +
  geom_line(size = 1) +
  labs(
    title = "Prior and Posterior Distributions for R",
    x = "R",
    y = "Density",
    color = "Distribution",
    linetype = "Distribution"
  ) + scale_x_continuous(
    limits = c(0, NA),
    expand = c(0, 0)
  ) +
  scale_y_continuous(
    limits = c(0, NA),
    expand = c(0, 0)
  ) + 
  scale_color_manual(values = c("Posterior" = "blue", "Prior" = "red")) +
  theme_classic() +
  theme(
    text = element_text(size = 16, face = "bold"),
    axis.text.y = element_text(size = 16, face = "bold"), # Y-axis text size
    axis.text.x = element_text(size = 16, face = "bold"), # Bold X-axis text
    legend.text = element_text(size = 16, face = "bold"),
    legend.title = element_text(size = 18, face = "bold"),
    legend.position = "top"
  )




# summary statistics
mean_posterior <- mean(posterior_samples_R)
median_posterior <- median(posterior_samples_R)
credible_interval <- quantile(posterior_samples_R, c(0.025, 0.975))

# Print the results
cat("Posterior Summary for R:\n")
cat("Mean:", signif(mean_posterior, 3), "\n")
cat("Median:", signif(median_posterior, 3), "\n")
cat("95% Credible Interval:", 
    paste0("(", signif(credible_interval[1], 3), ", ", signif(credible_interval[2], 3), ")"), "\n")

################################
# Define scenarios
scenario_1 <- c(rep(1, 56), c(2,2)) #53 single spillover, 2xcluster of 2(cali and source, Miss and source)
scenario_2 <- c(rep(1, 56), 3, 2) #53 single spillover , 1x cluster of 3 (Miss, source and household contact), 1x cluster of 2 (Cali and source)

# Function to perform MCMC and summarize results
run_mcmc <- function(h5_clusters, n_iter = 1e4, n_burn = 1e3) {
  # Likelihood function (as before)
  lik_function <- function(param) {
    if (any(param <= 0)) return(-Inf)
    r_val <- as.numeric(param[1])
    k_val <- as.numeric(param[2])
    log_likelihood <- likelihood(
      chains = h5_clusters,
      statistic = "size",
      offspring_dist = rnbinom,
      size = k_val,
      mu = r_val
    )
    log_prior <- h5_prior_r(r_val)
    return(log_likelihood + log_prior)
  }
  
  # Run MCMC
  result <- MCMCmetrop1R(
    lik_function,
    theta.init = c(R = 0.1, k = 1.0),
    burnin = n_burn,
    mcmc = n_iter,
    thin = 1
  )
  
  # Summarize posterior
  posterior_samples_R <- result[, 1]
  summary_stats <- list(
    mean = mean(posterior_samples_R),
    median = median(posterior_samples_R),
    ci95 = quantile(posterior_samples_R, c(0.025, 0.975))
  )
  
  return(list(result = result, posterior_samples_R = posterior_samples_R, summary_stats = summary_stats))
}

# Run MCMC for both scenarios
mcmc_scenario_1 <- run_mcmc(scenario_1)
mcmc_scenario_2 <- run_mcmc(scenario_2)

# Create prior density data
prior_density <- data.frame(
  x = seq(0, 5, length.out = 1000),  # Adjust range as needed
  y = h5_prior_r(seq(0, 5, length.out = 1000)),
  Scenario = "Prior"
)

# Create prior density data
prior_density <- data.frame(
  x = seq(0, 5, length.out = 1000),  # Adjust range as needed
  y = h5_prior_r(seq(0, 5, length.out = 1000)),
  Scenario = "Prior"
)

# Combine posterior densities with prior density
posterior_density_1 <- density(mcmc_scenario_1$posterior_samples_R)
posterior_density_2 <- density(mcmc_scenario_2$posterior_samples_R)

density_df <- rbind(
  data.frame(x = posterior_density_1$x, y = posterior_density_1$y, Scenario = "Scenario 1"),
  data.frame(x = posterior_density_2$x, y = posterior_density_2$y, Scenario = "Scenario 2"),
  prior_density
)

# Plot prior and posterior distributions
R0_US <- ggplot(density_df, aes(x = x, y = y, color = Scenario, linetype = Scenario)) +
  geom_line(size = 1) +
  labs(
    title = "",
    x = "R",
    y = "Density",
    color = "Distribution",
    linetype = "Distribution"
  ) +
  scale_color_manual(values = c("Scenario 1" = "blue", "Scenario 2" = "darkgreen", "Prior" = "red")) +
  scale_linetype_manual(values = c("Scenario 1" = "solid", "Scenario 2" = "solid", "Prior" = "dashed")) +
  scale_x_continuous(
    limits = c(0, 0.15),
    expand = c(0, 0)
  )+
  scale_y_continuous(
    limits = c(0, NA),
    expand = c(0, 0)
  ) +
  theme_classic() +
theme(
  text = element_text(size = 16, face = "bold"),
  axis.text.y = element_text(size = 16, face = "bold"), # Y-axis text size
  axis.text.x = element_text(size = 16, face = "bold"), # Bold X-axis text
  legend.text = element_text(size = 16, face = "bold"),
  legend.title = element_text(size = 18, face = "bold"),
  legend.position = "top"
)




# Summarize results for both scenarios
summary_df <- data.frame(
  Scenario = c("Scenario 1", "Scenario 2"),
  Mean_R = c(mcmc_scenario_1$summary_stats$mean, mcmc_scenario_2$summary_stats$mean),
  Median_R = c(mcmc_scenario_1$summary_stats$median, mcmc_scenario_2$summary_stats$median),
  CI95 = c(
    paste0("(", signif(mcmc_scenario_1$summary_stats$ci95[1], 3), ", ",
           signif(mcmc_scenario_1$summary_stats$ci95[2], 3), ")"),
    paste0("(", signif(mcmc_scenario_2$summary_stats$ci95[1], 3), ", ",
           signif(mcmc_scenario_2$summary_stats$ci95[2], 3), ")")
  )
)

# Display summary table
knitr::kable(summary_df, caption = "Summary of Posterior Estimates for Both Scenarios")


posterior_samples_k_scenario_1 <- mcmc_scenario_1$result[, 2]
posterior_samples_k_scenario_2 <- mcmc_scenario_2$result[, 2]

summary_k_scenario_1 <- list(
  mean = mean(posterior_samples_k_scenario_1),
  median = median(posterior_samples_k_scenario_1),
  ci95 = quantile(posterior_samples_k_scenario_1, c(0.025, 0.975))
)

summary_k_scenario_2 <- list(
  mean = mean(posterior_samples_k_scenario_2),
  median = median(posterior_samples_k_scenario_2),
  ci95 = quantile(posterior_samples_k_scenario_2, c(0.025, 0.975))
)

results_table <- data.frame(
  Scenario = c("Scenario 1", "Scenario 2"),
  Mean_R = c(mcmc_scenario_1$summary_stats$mean, mcmc_scenario_2$summary_stats$mean),
  Median_R = c(mcmc_scenario_1$summary_stats$median, mcmc_scenario_2$summary_stats$median),
  CI95_R = c(
    paste0("(", signif(mcmc_scenario_1$summary_stats$ci95[1], 3), ", ", 
           signif(mcmc_scenario_1$summary_stats$ci95[2], 3), ")"),
    paste0("(", signif(mcmc_scenario_2$summary_stats$ci95[1], 3), ", ", 
           signif(mcmc_scenario_2$summary_stats$ci95[2], 3), ")")
  ),
  Mean_k = c(summary_k_scenario_1$mean, summary_k_scenario_2$mean),
  Median_k = c(summary_k_scenario_1$median, summary_k_scenario_2$median),
  CI95_k = c(
    paste0("(", signif(summary_k_scenario_1$ci95[1], 3), ", ", 
           signif(summary_k_scenario_1$ci95[2], 3), ")"),
    paste0("(", signif(summary_k_scenario_2$ci95[1], 3), ", ", 
           signif(summary_k_scenario_2$ci95[2], 3), ")")
  )
)

knitr::kable(results_table, caption = "Posterior Estimates for R and k Across Scenarios")

barplot(table(scenario_1), main = "Cluster Sizes for Scenario 1")

