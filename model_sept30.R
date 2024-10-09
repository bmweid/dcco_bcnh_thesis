# add t distribution for response variable to handle outliers better than normal distribution
# quadratic term for nest density to cpature nonlinear effects
# AR(1) term to account for temporal autocorrelation
# implement variable selection using spike-and-slab priors with a horseshoe prior for regularization
#code for k-fold cross validation to assess performance 

# Load required libraries
library(tidyverse)
library(R2jags)
library(bayesplot)

# Load the data
model_data <- read_csv("combined_dataset.csv")

# Prepare the data
prepared_data <- model_data %>%
  select(year, bcnh_growthindex, total_nest_density, bcnh_nest_density, 
         dcco_growthindex, bcnh_nest_success, dcco_usurpation, 
         winter_nestremoval, deterrence_activenestremoval, deterrence_night,
         bcnh_road_proximity, raccoon_predation) %>%
  mutate(across(everything(), ~ifelse(is.na(.), mean(., na.rm = TRUE), .))) %>%
  mutate(across(everything(), as.numeric)) %>%
  mutate(across(-year, scale)) %>%
  arrange(year)  # Ensure data is in chronological order

# Create lagged bcnh_growthindex
prepared_data <- prepared_data %>%
  mutate(bcnh_growthindex_lag = lag(bcnh_growthindex))

# Define the number of observations (N) and unique years (N_years)
N <- nrow(prepared_data)
year <- as.numeric(factor(prepared_data$year))
N_years <- length(unique(year))

# Prepare data for JAGS
jags_data <- list(
  N = N,
  N_years = N_years,
  bcnh_growth_index = as.vector(prepared_data$bcnh_growthindex),
  bcnh_growth_index_lag = as.vector(prepared_data$bcnh_growthindex_lag),
  nest_density = as.vector(prepared_data$total_nest_density),
  dcco_growth_index = as.vector(prepared_data$dcco_growthindex),
  year = year,
  c = as.matrix(prepared_data[, c("total_nest_density", "bcnh_nest_density", 
                                  "dcco_growthindex", "bcnh_nest_success", 
                                  "winter_nestremoval", "deterrence_activenestremoval", 
                                  "deterrence_night", "dcco_usurpation", 
                                  "bcnh_road_proximity", "raccoon_predation")])
)

# JAGS model specification (updated)
model_string <- "
model {
  for (i in 1:N) {
    # BCNH growth index as response variable (observed data)
    bcnh_growth_index[i] ~ dt(mu[i], tau, nu)
    
    # Mean of the response variable
    mu[i] <- inprod(X[i,], beta[]) + 
             beta[11] * pow(nest_density[i], 2) +  # Quadratic term
             beta[12] * nest_density[i] * dcco_growth_index[i] + 
             alpha_year[year[i]] +
             phi * bcnh_growth_index_lag[i]  # AR(1) term using lagged variable
    
    # Design matrix X for variable selection
    for (j in 1:10) {
      X[i,j] <- equals(z[j], 1) * c[i,j]
    }
  }
  
  # Priors for fixed effects with variable selection
  for (j in 1:12) {
    beta[j] ~ dnorm(0, tau_beta * lambda[j])
    lambda[j] ~ dt(0, 1, 1)T(0,)  # Half-t prior for horseshoe
  }
  for (j in 1:10) {
    z[j] ~ dbern(0.5)  # Binary indicator for variable selection
  }
  
  # Global shrinkage parameter
  tau_beta ~ dt(0, 1, 1)T(0,)
  
  # Prior for residual precision and degrees of freedom
  tau <- pow(sigma, -2)
  sigma ~ dunif(0, 10)
  nu ~ dgamma(2, 0.1)  # Degrees of freedom for t-distribution
  
  # Priors for random year effect
  for (t in 1:N_years) {
    alpha_year[t] ~ dnorm(0, tau_year)
  }
  tau_year <- pow(sigma_year, -2)
  sigma_year ~ dunif(0, 10)
  
  # Prior for AR(1) coefficient
  phi ~ dunif(-1, 1)
}
"

# Initial values function (unchanged)
inits <- function() {
  list(
    beta = rnorm(12, 0, 0.1),
    sigma = runif(1, 0, 1),
    sigma_year = runif(1, 0, 1),
    alpha_year = rnorm(N_years, 0, 0.1),
    nu = 5,
    phi = 0,
    z = rbinom(10, 1, 0.5),
    lambda = rexp(12, 1)
  )
}

# Parameters to monitor (unchanged)
parameters <- c("beta", "sigma", "sigma_year", "alpha_year", "nu", "phi", "z", "lambda")

# Fit the model
fit <- jags(
  data = jags_data,
  inits = inits,
  parameters.to.save = parameters,
  model.file = textConnection(model_string),
  n.chains = 3,
  n.iter = 50000,
  n.burnin = 2000,
  n.thin = 10
)

# Convert to MCMC for diagnostics
fit_mcmc <- as.mcmc(fit)

# Convert JAGS output to a format suitable for bayesplot
posterior_samples <- as.array(fit$BUGSoutput$sims.array)

# Create trace plots
pdf("trace_plots.pdf", width = 12, height = 8)
mcmc_trace(posterior_samples, 
           pars = c("beta[1]", "beta[2]", "beta[3]", "beta[4]", "beta[5]",
                    "beta[6]", "beta[7]", "beta[8]", "beta[9]", "beta[10]",
                    "beta[11]", "beta[12]", "beta[13]", "sigma", "sigma_year", "nu", "phi"),
           facet_args = list(ncol = 4, strip.position = "top"))
dev.off()

# Create density plots
pdf("density_plots.pdf", width = 12, height = 8)
mcmc_dens(posterior_samples,
          pars = c("beta[1]", "beta[2]", "beta[3]", "beta[4]", "beta[5]",
                   "beta[6]", "beta[7]", "beta[8]", "beta[9]", "beta[10]",
                   "beta[11]", "beta[12]", "beta[13]", "sigma", "sigma_year", "nu", "phi"),
          facet_args = list(ncol = 4))
dev.off()

# Create a combined plot with both traces and densities
pdf("combined_diagnostic_plots.pdf", width = 15, height = 12)
mcmc_combo(posterior_samples,
           pars = c("beta[1]", "beta[2]", "beta[3]", "beta[4]", "beta[5]",
                    "beta[6]", "beta[7]", "beta[8]", "beta[9]", "beta[10]",
                    "beta[11]", "beta[12]", "beta[13]", "sigma", "sigma_year", "nu", "phi"),
           combo = c("dens", "trace"),
           facet_args = list(ncol = 4))
dev.off()

# Summary of posterior distributions
summary(fit_mcmc)

# Save the summary results as a CSV file
posterior_samples <- as.matrix(fit_mcmc)
results <- data.frame(
  parameter = colnames(posterior_samples),
  mean = colMeans(posterior_samples),
  lower_ci = apply(posterior_samples, 2, quantile, probs = 0.025),
  upper_ci = apply(posterior_samples, 2, quantile, probs = 0.975)
)

# Print and save results table
print(results)
write.csv(results, "jags_comprehensive_model_results_table.csv", row.names = FALSE)

# Calculate residuals
y_obs <- jags_data$bcnh_growth_index
mu <- fit$BUGSoutput$sims.list$mu
mu_mean <- apply(mu, 2, mean)
residuals <- y_obs - mu_mean

# Check residuals
print(residuals)

# Visualize residuals
pdf("residual_plots.pdf")
hist(residuals, main="Histogram of Residuals", xlab="Residuals")
plot(y_obs, residuals, main="Residuals vs Observed", xlab="Observed Values", ylab="Residuals")
abline(h=0, col="red")
dev.off()

# Bayesian P-value calculation
deviance_samples <- fit$BUGSoutput$sims.list$deviance
y_obs <- jags_data$bcnh_growth_index
y_rep <- fit$BUGSoutput$sims.list$y_rep
sigma <- fit$BUGSoutput$sims.list$sigma

calc_deviance <- function(y, y_pred, sigma) {
  -2 * sum(dnorm(y, mean = y_pred, sd = sigma, log = TRUE))
}

deviance_obs <- apply(y_rep, 1, function(y_pred) calc_deviance(y_obs, y_pred, sigma))
deviance_rep <- apply(y_rep, 1, function(y_pred) calc_deviance(y_pred, y_pred, sigma))

bayes_p_value <- mean(deviance_rep > deviance_obs)
print(paste("Bayesian P-value: ", bayes_p_value))

# Compute DIC
dic <- dic.samples(fit$model, n.iter = 1000)
print(dic)

# Perform 5-fold cross-validation
k <- 5
folds <- cut(seq(1, N), breaks = k, labels = FALSE)
cv_results <- vector("list", k)

for (i in 1:k) {
  test_indices <- which(folds == i)
  train_indices <- which(folds != i)
  
  # Prepare training data
  train_data <- jags_data
  train_data$N <- length(train_indices)
  train_data$bcnh_growth_index <- jags_data$bcnh_growth_index[train_indices]
  train_data$c <- jags_data$c[train_indices,]
  train_data$year <- jags_data$year[train_indices]
  
  # Fit model on training data
  cv_fit <- jags(data = train_data, inits = inits, parameters.to.save = parameters,
                 model.file = textConnection(model_string),
                 n.chains = 3, n.iter = 10000, n.burnin = 1000, n.thin = 10)
  
  # Predict on test data
  test_data <- jags_data
  test_data$N <- length(test_indices)
  test_data$bcnh_growth_index <- jags_data$bcnh_growth_index[test_indices]
  test_data$c <- jags_data$c[test_indices,]
  test_data$year <- jags_data$year[test_indices]
  
  pred <- jags.predict(cv_fit, test_data)
  
  # Calculate RMSE
  rmse <- sqrt(mean((jags_data$bcnh_growth_index[test_indices] - colMeans(pred))^2))
  cv_results[[i]] <- rmse
}

mean_cv_rmse <- mean(unlist(cv_results))
print(paste("Mean cross-validated RMSE:", mean_cv_rmse))