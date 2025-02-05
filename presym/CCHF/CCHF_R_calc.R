
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


# Define CCHF cluster 
# Cases in the EU from 2013-2024
#https://www.ecdc.europa.eu/en/crimean-congo-haemorrhagic-fever/surveillance/cases-eu-since-2013
CCHF_clusters = c(rep(1, 59), 2)

# Show summary table of frequencies
freq_df <- as.data.frame(table(CCHF_clusters)); names(freq_df) <- c("Cluster size", "Frequency")

# Create a table for the HTML document
knitr::kable(freq_df, caption = "Frequencies of MERS Clusters")

# Define likelihood function
lik_function <- function(param) {
  if (any(param <= 0)) return(-Inf) # Ensure positive parameters
  
  # Extract values of R and k
  r_val <- as.numeric(param[1])
  k_val <- as.numeric(param[2])
  
  # Define likelihood
  log_likelihood <- likelihood(
    chains = CCHF_clusters,
    statistic = "size",
    offspring_dist = rnbinom,
    size = k_val,
    mu = r_val
  )
  
  # Assume non-informative priors for R and k
  log_prior <- 0 # But could add informative priors here if required
  
  # Return log-posterior (log-likelihood + log-prior)
  return(log_likelihood + log_prior)
}

# Define number of MCMC iterations
n_iter <- 1e4

# Define 'burn in' period for fitting, to be discarded
n_burn <- 1e3

# Initial guess for c(R,k):
init_param <- c(R=0.2, k=2)

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

library(ggplot2)
library(dplyr)

# Convert MCMC samples to a dataframe
posterior_df <- as.data.frame(result_mcmcpack)
colnames(posterior_df) <- c("R0", "k")  # Rename for clarity

# Calculate median and 95% CrI
R0_median <- median(posterior_df$R0)
R0_CrI <- quantile(posterior_df$R0, c(0.025, 0.975))

# Create density plot with ggplot2
R_1 <- ggplot(posterior_df, aes(x = R0)) +
  geom_density(fill = "lightblue", alpha = 0.6, color = "black") +
  geom_vline(xintercept = R0_median, color = "red", linetype = "dashed", linewidth = 1.2) +
  labs(
    title = "R0 posterior Distribution for CCHFV (EU cases 2013-2024)",
    x = expression(R[0]),
    y = "Density"
  ) +
  theme_classic() +
  scale_x_continuous(
    limits = c(0, 0.15),
    expand = c(0, 0)
  ) +
  scale_y_continuous(
    limits = c(0, NA),
    expand = c(0, 0)
  ) 

ggsave("presym/CCHF/CCHF_R.png", R_1, width = 8, height = 4)