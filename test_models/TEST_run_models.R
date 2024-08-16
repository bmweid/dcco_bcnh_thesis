# File: TEST_run_models.R
# Purpose: Debug and run the comprehensive test model using R2jags (Peninsula A, 2010-2013)

library(tidyverse)
library(R2jags)

# Resolve the filter conflict
conflicts_prefer(dplyr::filter)

# Load the data
model_data <- read_csv("C:/Users/baill/OneDrive/r-projects/cormorant-colony-analysis/current working files/model datasets/combined_model_data.csv")

# Filter data for Peninsula A and time range 2010-2013
filtered_data <- model_data %>%
  dplyr::filter(peninsula == "A", year >= 2010, year <= 2013)

# Rename columns to match JAGS model
filtered_data <- filtered_data %>%
  rename(
    nest_density = density,
    road_proximity = proximity,
    bcnh_nest_success = success,
    dcco_usurpation = usurpation,
    dcco_management = management
  )

N <- nrow(filtered_data)
year <- as.numeric(factor(filtered_data$year))
N_years <- length(unique(year))

print(paste("Number of observations:", N))
print(paste("Number of years:", N_years))

# Prepare data for JAGS
jags_data <- list(
  N = N,
  y = filtered_data$growth_index,
  nest_density = filtered_data$nest_density,
  road_proximity = filtered_data$road_proximity,
  bcnh_nest_success = filtered_data$bcnh_nest_success,
  dcco_usurpation = filtered_data$dcco_usurpation,
  dcco_management = filtered_data$dcco_management,
  N_years = N_years,
  year = year
)

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

# Rest of the script (inits, parameters, model running) remains the same...

# Initial values function
inits <- function() {
  list(
    beta = rnorm(6, 0, 1),
    sigma = runif(1, 0, 10),
    sigma_year = runif(1, 0, 10),
    alpha_year = rnorm(N_years, 0, 1)
  )
}

# Parameters to monitor
parameters <- c("beta", "sigma", "sigma_year", "alpha_year")

# Run the model with R2jags
fit <- jags(
  data = jags_data,
  inits = inits,
  parameters.to.save = parameters,
  model.file = "C:/Users/baill/OneDrive/r-projects/cormorant-colony-analysis/current working files/winbugs model scripts/test_models/comprehensive_model.txt",
  n.chains = 3,
  n.iter = 50000,
  n.burnin = 5000,
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
saveRDS(fit, "jags_comprehensive_model_results.rds")

# Create diagnostic plots
traceplot(fit)
densplot(fit)

# Extract posterior samples
posterior_samples <- as.matrix(fit_mcmc)

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