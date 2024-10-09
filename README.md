# dcco_bcnh_thesis
Thesis files for analysis of colony at TTP


[following Wyman et al]
Used linear mixed-effects models to relate cormorant abundance and management to bcnh covariates
Model specification:
Fit model for each peninsula (3 models total?)
All models to include nest abundance, growth index, nest density, dcco-bcnh clustering, bcnh nest success, and management
Scale all numerical covariates to have a standard deviation of 1 (to be completed)
Significant missing data depending on covariate; fit models using Bayesian approach to compensate
Percentages of missing data ranges: [insert]
Constructed distributions for missing data using means and standard deviations of observed data 
[insert how distributions were constructed for which variable; empirical, exponential, Bernouilli, uniform?]
Choose priors:
Used diffuse normal priors for effects of linear predictors and diffuse uniform priors for variance parameters
[Sampling and Inference, convergence diagnostics]
R; used STAN, Hamilton Monte Carlo (MCMC method) for efficiency, ran 3 simultaneous chains for [30,000-50,000 iterations,] as determined by number of iterations required to achieve Monte Carlo standard <2% of posterior distribution standard deviation and R-hat statistic [measure used in Bayesian statistics to assess convergence of MCMC simulations]
Discard a burn-in period of [2,000?] and thinned chains by [10?] to reduce computational memory limitations
Evaluate model fit:
Using Bayesian P-values (King et al., 2010)
Discrepancy function = sum of squared residuals 
P-value of 0.5 indicates good model fit (that the fit model to the observed data is no better or worse than the fit of the model to data simulated from the model itself
Bayesian P-values range from 0 to 1 (King et al. 2010)

