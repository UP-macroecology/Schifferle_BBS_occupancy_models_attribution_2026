# Fit dynamic occupancy models to the North American breeding bird survey data

We (try to) fit dynamic occupancy models (DOMs) to predict probabilities of initial occupancy, local extinction and colonisation for North American breeding bird species based on climate and land use variables.
We focus on the conterminous US and the time period 1991 to 2015. Models are fit for single species using a Bayesian approach.

In a later step, we aim to determine the relative importance of climate and land use change on changes in occupancy, extinction and colonisation across species. Therefore we use the same set of environmental covariates for all species.

## General workflow:

### Data preparation:

#### Species data:

We use North American breeding bird survey (BBS) data that had been further processed [here](https://github.com/UP-macroecology/Schifferle_BBS_explorations_2023/blob/main/BBS_data_prep.R). We reformatted them for occupancy modelling (presences and absences for each route, year and species) (*1_0_reformat_BBS_data.R*).

We further filtered the BBS routes (*1_1_route_selection.R*): we discarded routes that were rarely surveyed within our focal time period (routes required to be sampled in the first and last year and maximum 5 years missing in between). We spatially thinned the routes to avoid that the same environmental data cell contains multiple routes and we retained only 30 routes per Bird Conservation Region to get a more equal coverage of the conterminous US. We ended up with 476 routes.

We also filtered the species for which we fit DOMs (*1_2_species_selection.R*): we discarded species for which the BBS methodology is not suited (nocturnal and water-related species), as well as particularly widespread or rare species and ended up with 174 species.

#### Environmental data:

As environmental data we use climate and land use data.
For climate, we use the ISIMIP Chelsa-W5E5v1.0 dataset and calculated 19 bioclimatic variables as well as mean temperature and precipitation sums
of each season (spring, summer, autumn, winter) (*1_0_prepare_bioclim_data.R*).
For land use, we use the land use states of the LUH2 dataset (*1_0_prepare_landuse_data.R*).
We use all variables at a resolution of 0.5°.

#### Variable selection for DOMs:

As covariates in the DOMs we chose a subset of the climatic and land use variables. We retained variables that are not substantially correlated and that capture most of the environmental variation across the conterminous US. The latter was determined via PCA and some manual adjustments afterwards (*1_2_variable_selection.R*). We ended up with ten climate and four land use related variables.

#### Misc:

The script *1_3_match_BBS_to_env_data.R* extracts values of the environmental variables at the routes' centroids. These will be used to fit the models.
The script *1_3_env_data_df_contUS.R* extracts values of the variables at each grid cell's centroids across the conterminous US. These will be used to predict to the whole conterminous US using the fitted models.
The script *1_3_explore_occ_data.R* contains preliminary code for some exploratory plots of the data.

### Modelling:

Initial occupancy probability was modelled using linear and quadratic terms of the selected climatic and land use covariates summarising the three years preceeding the first survey considered.
Colonisation and extinction probability were modelled using the same covariates, but each summarising only the 12 months preceeding each survey.
Detection probability was modelled with route section (1-5) as covariate, which we consider a proxy for time of the day.
For each species, we considered only the routes that are located within a 750 km buffer around all routes where the species was recorded at least once.

To fit the models, we initially used JAGS, but then switched to the R package *flocker*, which was developed to fit DOMs via Stan, as it promised better performance / faster model fitting.
With JAGS we had some issues with chain convergence, we explored whether using different subsets of covariates improves chain convergence. We included calculating a goodness-of-fit estimate and used lasso priors.
With flocker we experimented with different priors, since lasso priors where not supported and not recommended.
We also experiment regarding routes that are identified as highly influential for the model by PSIS-LOO CV.
All *2_1_....R* scripts refer to this process.

### Model fit explorations:

- *2_2_DOM_explore_results.R*: explore and compare DOMs fitted with JAGS and flocker
- *2_2_DOM_flocker_results.R*: extract and plot outputs from flocker fit, including PSIS-LOO CV
- *2_2_DOM_predictions_USA.R*: plot maps of predicted occupancy, colonisation and extinction probability across the conterminous US
- *2_2_evaluate_occ_models.R*: some very preliminary code snippets for model evaluation