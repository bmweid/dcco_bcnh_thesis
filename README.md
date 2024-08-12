# dcco_bcnh_thesis
Thesis files for analysis of colony at TTP

Contains prepped data for analysis in R using linear mixed models.


Data/database
Has been very challenging to clean and set up properly for R and ArcGIS
Started with integrating new counts for 2015-2023 using excel and VLOOKUP formulas 
Very challenging as sets varied significantly across years- eg. some years they did not use ‘new’ to indicate new tags, and sometimes they did. Originally I had been adding new tags using this indication and had to go back over everything with VLOOKUP formulas when I found massive variances in nest counts.
Methods of coordinate collection varied year to year- different geographic coordinate systems. Had to convert them all to X & Y in UTM 17N to get them to show in ArcGIS, which took time. 
Coordinates are still missing for many tree tags but I searched all provided datasets multiple times using VLOOKUP- as complete as likely to get.
QUESTION: My impression is that I’ll need to omit trees with missing coordinates for use in the model as the model will look at each peninsula separately for trends. Wanted to confirm this makes sense.
In github: https://github.com/bmweid/dcco_bcnh_thesis/blob/main/dcco_bcnh_database_%201992_2023_july30th.xlsx
Database contains all treeids (tags vary if replace, treeid used as common identifier) and counts yearly, including NA those with coordinates

Working with the database in ArcGIS and R:
Data had to be ‘pivoted longer’ so that each unifying tree tag “treeid” and year had an associated dcco & bcnh count, significantly increasing the number of rows in the data
Code to do that here:
To prep data for Arcgis (can have no NAs, must have no characters in numeric columns, etc), Data was originally in two separate layers, one for dcco and one for bcnh, for easier database management since pivoting increases rows significantly, as well as I was unsure if I would need them separate for analysis in ArcGIS
o https://github.com/bmweid/dcco_bcnh_thesis/blob/main/prep_forarcgis_databases.rmd
After uploading to Arcgis Pro and running some analysis, code to merge separate layers into one database:
https://github.com/bmweid/dcco_bcnh_thesis/blob/86791faa05fb1c74c50e4ce09f5320bef24956a4/steps_clean_postarcgis_data.R
Resulting data after uploading and merging in Arcgis, In github: https://github.com/bmweid/dcco_bcnh_thesis/blob/main/arcgismerged_cleanedbcnhdcco.csv
Had done significant amount of work in QGIS but was unable to determine how to run Moran’s I in QGIS; recently switched to ArcGIS, resulting in some delays


Response variable- partially complete
Colony growth index (following Wyman et al 2018, Guillamet et al 2004)
DCCO population growth by subcolony

to be completed- have code and had run it, but not divided by sub colony. Had to go back to add missing peninsula data in ArcGIS.
Code is here, will rerun with current database this week:
https://github.com/bmweid/dcco_bcnh_thesis/blob/main/colony_growth_index.R
Data predictor variables
Nest Density - not complete
Looking to calculate mean nest density in ArcGIS Pro, per peninsula per year
This number will be used in the model, with all covariates scaled to standard deviation of 1
Next steps: having been looking to determine best tool for this, and how to determine a “tipping point” per tree- a number of dcco nests by which bcnh will no longer nest in the tree
Had challenges formatting data correctly to determine this in ArcGIS
Next steps:
Complete nest density analysis
Create heat maps for inclusion in thesis from points
QUESTION/LOOKING FOR INPUT:
Any advice on how to approach this?

Linear Features- partially complete
In github: https://github.com/bmweid/dcco_bcnh_thesis/blob/main/bcnh_roadproximity.csv
Dataset contains count of bcnh nests in proximity to roads yearly
Determined creating linear features for the roads, setting an average road width, and a buffer of 40m, selecting all trees that contain bcnh in that buffer yearly
Next steps: 
Looking to verify that a 40 m buffer is adequate to determine a distance that may impact bcnh nest success- I choose this rather arbitrarily to test out the tool
Need to verify that a simple count is suitable for use in the model 
LOOKING FOR INPUT:
Looking for papers that may determine an ‘impact distance’ for roads for nest success, let me know if you have any ideas

Understory vegetation- scrapped for time
Unfortunately I don’t think I’ll be able to get to this considering current time demands

BCNH Nest Success- complete
https://github.com/bmweid/dcco_bcnh_thesis/blob/main/bcnh_nest_success.csv
Using Gail’s provided data, nest success from 2009-2021 as percent per year of success nests
Will capture impact of raccoon predation

DCCO Usurpation-  complete
Using Gail’s provided data for # of BCNH nests usurped annually 2008-2021
https://github.com/bmweid/dcco_bcnh_thesis/blob/main/dcco_usurpation.csv

DCCO Management- partially complete
Using McDonald et al 2018 deterrence activities time summary, 2008 to 2016. Requires update from TRCA from 2016-present
Decided on a binary (0/1) indicator for management activity (aka no management from 1992-2008). Management activities divided into categorical variables and binary yes/no; winter nest removal, deterrence with active nest removal, nighttime deterrence, ground nest enhancements.
https://github.com/bmweid/dcco_bcnh_thesis/blob/main/dcco_management.csv
Next Steps:
Need to reach out to TRCA to see if any more can be added to this from 2016-present
QUESTION/INPUT REQUESTED:
Should this be simplified to simple yes/no for management yearly or should I try to ‘weight’ management activities? Unfortunately in discussion with the TRCA it didn’t seem like there would be more complex data available (eg were bcnh in proximity to management activities, how many nests well removed, intensity of management activities)


Model Notes
Draft code included here, based on Wyman et al’s paper- needs to be tested and cleaned, certainly happy to get suggestions and feedback
Using STAN instead of Winbugs at the moment- it has more active development, HMC methods of MCMC that are more efficient, more flexibility and faster performance because of the MCMC method (may revert to Winbugs after testing, since STAN has a steeper learning curve)
https://github.com/bmweid/dcco_bcnh_thesis/blob/main/rstan%20version.R
This hasn’t been run or tested yet, planning to test with sample data this week
Further notes on model in Thesis-Outline July 20th.docx
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
