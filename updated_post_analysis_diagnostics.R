library(coda)
library(bayesplot)
library(ggplot2)

# Load the saved results
fit <- readRDS("jags_comprehensive_model_results.rds")

# Convert to mcmc object for diagnostics
mcmc_samples <- as.mcmc(fit)

# Extract posterior samples and observed data
y_rep <- fit$BUGSoutput$sims.list$y_rep
y <- fit$BUGSoutput$data$y

# Convergence diagnostic: Gelman-Rubin statistic
print(gelman.diag(mcmc_samples))

# Posterior predictive check: density overlay
ppc_dens_overlay(y, y_rep[1:100,]) + ggtitle("Posterior Predictive Density Check")
ggsave("ppc_density.pdf", width = 10, height = 6)

# Calculate Bayesian p-value
discrepancy_function <- function(y, theta) {
  residuals <- y - theta
  sum(residuals^2)
}

T_obs <- apply(y_rep, 1, function(y_rep_i) discrepancy_function(y, y_rep_i))
T_rep <- apply(y_rep, 1, function(y_rep_i) discrepancy_function(y_rep_i, y_rep_i))

bayesian_p_value <- mean(T_rep > T_obs)

print(paste("Bayesian p-value:", bayesian_p_value))

# Plot the discrepancy distribution
pdf("discrepancy_distribution.pdf")
hist(T_rep, main="Posterior Predictive Check", xlab="Discrepancy")
abline(v=mean(T_obs), col="red", lwd=2)
dev.off()

# Summary of posterior distributions
summary_stats <- summary(mcmc_samples)
write.csv(summary_stats$statistics, "summary_statistics.csv")

print("Simplified post-analysis diagnostics completed. Check the generated PDF and CSV files for results.")