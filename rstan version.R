library (tidyverse)
library (dplyr)
library (tidyr)
library (rstan)

#prepare data list for Stan
#create data frame named df with the required variables
df <- data.frame(
  colony_growth_index_dcco = rnorm(100),
  site = sample(1:10, 100, replace = TRUE),
  survey = sample(1:4, 100, replace = TRUE),
  prior_growth_rate = rnorm(100),
  other_disturbance = rnorm(100),
  DCCO_abundance = rnorm(100),
  dcco_population_growth = rnorm(100),
  dcco_bcnh_nest_density = rnorm(100),
  dcco_bcnh_clustering = rnorm(100),
  bcnh_nest_success = rnorm(100),
  dcco_nest_usurpation = rnorm(100),
  dcco_management = rnorm(100),
  linear_features = rnorm(100)
)

# Create the list to be passed to Stan
data_list <- list(
  N = nrow(df),
  J = length(unique(df$site)),
  L = length(unique(df$survey)),
  site = df$site,
  survey = df$survey,
  prior_growth_rate = df$prior_growth_rate,
  other_disturbance = df$other_disturbance,
  DCCO_abundance = df$DCCO_abundance,
  dcco_population_growth = df$dcco_population_growth,
  dcco_bcnh_nest_density = df$dcco_bcnh_nest_density,
  dcco_bcnh_clustering = df$dcco_bcnh_clustering,
  bcnh_nest_success = df$bcnh_nest_success,
  dcco_nest_usurpation = df$dcco_nest_usurpation,
  dcco_management = df$dcco_management,
  linear_features = df$linear_features,
  colony_growth_index = df$colony_growth_index
)

# Verify the data list structure
str(data_list)

# Define the Stan model
stan_model_code <- "
data {
  int<lower=0> N;  // number of observations
  int<lower=0> J;  // number of sites
  int<lower=0> L;  // number of survey periods
  int<lower=1,upper=J> site[N];  // site indicator
  int<lower=1,upper=L> survey[N];  // survey period indicator
  vector[N] dcco_population_growth;
  vector[N] dcco_bcnh_nest_density;
  vector[N] dcco_bcnh_clustering;
  vector[N] bcnh_nest_success;
  vector[N] dcco_nest_usurpation;
  vector[N] dcco_management;
  vector[N] linear_features;
  vector[N] colony_growth_index;  // dependent variable
}

parameters {
  real beta_0;
  real beta_1;
  real beta_2;
  real beta_3;
  real beta_4;
  real beta_5;
  real beta_6;
  real beta_7;
  vector[J] a_site;
  vector[L] a_survey;
  real<lower=0> sigma;
  real<lower=0> sigma_site;
  real<lower=0> sigma_survey;
}

model {
  vector[N] mu;

  for (i in 1:N) {
    mu[i] = beta_0 + beta_1 * dcco_population_growth[i] + beta_2 * dcco_bcnh_nest_density[i] +
            beta_3 * dcco_bcnh_clustering[i] + beta_4 * bcnh_nest_success[i] +
            beta_5 * dcco_nest_usurpation[i] + beta_6 * dcco_management[i] + 
            beta_7 + linear_features[i] +
            a_site[site[i]] + a_survey[survey[i]];
  }

  // Priors
  beta_0 ~ normal(0, 10);
  beta_1 ~ normal(0, 10);
  beta_2 ~ normal(0, 10);
  beta_3 ~ normal(0, 10);
  beta_4 ~ normal(0, 10);
  beta_5 ~ normal(0, 10);
  beta_6 ~ normal(0, 10);
  beta_7 ~ normal(0, 10);
  a_site ~ normal(0, sigma_site);
  a_survey ~ normal(0, sigma_survey);
  sigma ~ normal(0, 10);
  sigma_site ~ normal(0, 10);
  sigma_lake ~ normal(0, 10);
  sigma_survey ~ normal(0, 10);

  // Likelihood
  colony_growth_index ~ normal(mu, sigma);
}
"