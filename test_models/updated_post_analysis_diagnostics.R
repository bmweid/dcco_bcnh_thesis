# Updated post-analysis diagnostics

library(coda)
library(bayesplot)
library(ggplot2)
library(loo)  # for WAIC calculation

# Load the saved results
fit <- readRDS("jags_comprehensive_model_results.rds")

# Convert to mcmc object for diagnostics
mcmc_samples <- as.mcmc(fit)

# Print DIC
print(paste("DIC:", fit$BUGSoutput$DIC))

# Calculate WAIC
log_lik <- fit$BUGSoutput$sims.list$loglik
waic <- waic(log_lik)
print(waic)

# Plot WAIC weights
plot(waic, main="WAIC plot")

# Posterior predictive checks
y_rep <- fit$BUGSoutput$sims.list$y_rep
y <- fit$BUGSoutput$data$y

ppc_dens_overlay(y, y_rep[1:100,]) + ggtitle("Posterior Predictive Density Check")
ggsave("ppc_density.pdf", width = 10, height = 6)

ppc_stat(y, y_rep, stat = "mean") + ggtitle("Posterior Predictive Check: Mean")
ggsave("ppc_mean.pdf", width = 10, height = 6)

ppc_stat(y, y_rep, stat = "sd") + ggtitle("Posterior Predictive Check: SD")
ggsave("ppc_sd.pdf", width = 10, height = 6)

# Summary of posterior distributions
summary_stats <- summary(mcmc_samples)
write.csv(summary_stats$statistics, "summary_statistics.csv")
write.csv(summary_stats$quantiles, "summary_quantiles.csv")

print("Post-analysis diagnostics completed. Check the generated PDF and CSV files for results.")