# Load required libraries
library(tidyverse)
library(R2jags)

# Load the data
model_data <- read_csv("C:/Users/baill/OneDrive/r-projects/cormorant-colony-analysis/current working files/model datasets/combined_dataset.csv")

# Prepare the data
prepared_data <- model_data %>%
  select(year, growth_index_BCNH, total_nest_density, bcnh_nest_density, 
         growth_index_DCCO, bcnh_nest_success, dcco_usurpation, 
         winter_nestremoval, deterrence_activenestremoval, deterrence_night,
         bcnh_road_proximity)

# Handle missing values
prepared_data <- prepared_data %>%
  mutate(across(everything(), ~ifelse(is.na(.), mean(., na.rm = TRUE), .)))

# Ensure all variables are numeric
prepared_data <- prepared_data %>%
  mutate(across(everything(), as.numeric))

# Standardize all variables except year
prepared_data <- prepared_data %>%
  mutate(across(-year, scale))

# Calculate number of observations (N) and unique years (N_years)
N <- nrow(prepared_data)
year <- as.numeric(factor(prepared_data$year))
N_years <- length(unique(year))

# Prepare data for JAGS
jags_data <- list(
  N = nrow(prepared_data),
  N_years = length(unique(prepared_data$year)),
  bcnh_growth_index = as.vector(prepared_data$growth_index_BCNH),
  nest_density = as.vector(prepared_data$total_nest_density),
  bcnh_nest_density = as.vector(prepared_data$bcnh_nest_density),
  dcco_growth_index = as.vector(prepared_data$growth_index_DCCO),
  nest_success = as.vector(prepared_data$bcnh_nest_success),
  winter_nest_removal = as.vector(prepared_data$winter_nestremoval),
  deterrence_active_nest = as.vector(prepared_data$deterrence_activenestremoval),
  deterrence_at_night = as.vector(prepared_data$deterrence_night),
  usurpation = as.vector(prepared_data$dcco_usurpation),
  road_proximity = as.vector(prepared_data$bcnh_road_proximity),
  year = as.vector(as.numeric(factor(prepared_data$year)))
)

# Verify that all variables are vectors
for (var in names(jags_data)) {
  if (!is.null(dim(jags_data[[var]]))) {
    print(paste(var, "is not a vector. Converting to vector."))
    jags_data[[var]] <- as.vector(jags_data[[var]])
  }
}

# Print dimensions of key variables for verification
print("Dimensions of key variables:")
for (var in names(jags_data)) {
  if (is.vector(jags_data[[var]])) {
    print(paste(var, ":", length(jags_data[[var]])))
  } else {
    print(paste(var, ":", paste(dim(jags_data[[var]]), collapse="x")))
  }
}
# Print structure of jags_data for debugging
print("Structure of jags_data:")
print(str(jags_data))

# Print first few rows of each variable
print("First few rows of each variable:")
for (var in names(jags_data)) {
  if (is.vector(jags_data[[var]])) {
    print(paste(var, ":", paste(head(jags_data[[var]]), collapse = ", ")))
  }
}

# Write the updated model to a file
cat("
model {
  for (i in 1:N) {
    # BCNH growth index as response variable
    bcnh_growth_index[i] ~ dnorm(mu[i], tau)
    mu[i] <- beta[1] + beta[2] * nest_density[i] + beta[3] * bcnh_nest_density[i] + 
             beta[4] * dcco_growth_index[i] + beta[5] * nest_success[i] + 
             beta[6] * winter_nest_removal[i] + beta[7] * deterrence_active_nest[i] + 
             beta[8] * deterrence_at_night[i] + beta[9] * usurpation[i] + 
             beta[10] * road_proximity[i] + alpha_year[year[i]]
  }
  
  # Priors for fixed effects
  for (j in 1:10) {
    beta[j] ~ dnorm(0, 0.01)
  }
  
  # Prior for residual precision
  tau <- pow(sigma, -2)
  sigma ~ dunif(0, 10)
  
  # Priors for random year effect
  for (t in 1:N_years) {
    alpha_year[t] ~ dnorm(0, tau_year)
  }
  tau_year <- pow(sigma_year, -2)
  sigma_year ~ dunif(0, 10)
}", file = "comprehensive_model.txt")

# Initial values function
inits <- function() {
  list(
    beta = rnorm(10, 0, 0.1),
    sigma = runif(1, 0, 1),
    sigma_year = runif(1, 0, 1),
    alpha_year = rnorm(N_years, 0, 0.1)
  )
}

# Parameters to monitor
parameters <- c("beta", "sigma", "sigma_year", "alpha_year")

# Fit the model
fit <- jags(
  data = jags_data,
  inits = inits,
  parameters.to.save = parameters,
  model.file = "comprehensive_model.txt",
  n.chains = 3,
  n.iter = 50000,
  n.burnin = 2000,
  n.thin = 10
)


# Print results
print(fit)

# Plot diagnostics
plot(fit)

# Get summary of the posterior distributions
fit_mcmc <- as.mcmc(fit)
summary(fit_mcmc)

# Save results
saveRDS(fit, "jags_comprehensive_model_results_testAug25.rds")

# Create diagnostic plots
traceplot(fit)
densplot(fit)

# Extract posterior samples

# Calculate mean and 95% credible intervals
results <- data.frame(
  parameter = colnames(posterior_samples),
  mean = colMeans(posterior_samples),
  lower_ci = apply(posterior_samples, 2, quantile, probs = 0.025),
  upper_ci = apply(posterior_samples, 2, quantile, probs = 0.975)
)

# Print results table
print(results)

# Save results table
write.csv(results, "jags_comprehensive_model_results_table.csv", row.names = FALSE)