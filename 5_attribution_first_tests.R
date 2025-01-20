# first tests to attribute occupancy dynamics to change in environmental variables
# idea: detrend environmental variables, refit the model (or: predict to detrended data using the original model)
# see how prediction accuracy changes
# do we need trend in environmental variables to predict occupancy dynamics?

# example species for which model seems comparatively good:
# Black Vulture (spatial C index 0.89, temporal trend (increasing) roughly captured)
# occurs in Southern US and South America, expands its range in the US, one reason could be climate change (warmer temprature)

# first: use only variables that show a significant trend, based on a linear model, and
# for which parameter estimates don't overlap zero in the original model:
# bio1, bio3, pr_summer, pr_winter, sum_annual_crops, secdf, urban

# packages: ----
library(dplyr)
library(sf)
library(doParallel)
library(flocker)
library(cmdstanr)
set_cmdstan_path(path = NULL)#set_cmdstan_path("C:/Users/schifferle1/Documents/cmdstan-2.34.1") # xx


# functions: ----

source("0_functions.R")

# directories: ----

print(tempdir())
dir <- file.path("/import", "ecoc9z", "data-zurell", "schifferle", "BBS_occupancy_models_2023")
#dir <- getwd()


# load data: ----

# selected routes and focal years matched to environmental data:
# merged route, year, environment data:
load(file = file.path("data", "route_year_env_data.RData"))
route_sel_env_dt_final <- route_sel_env_dt_final %>% 
  # only selected variables:
  select(c("RTENO", "Year", "bio1", "bio2", "bio3", "bio7", "bio14", "bio15", 
           "pr_spring", "pr_summer","pr_autumn", "pr_winter",
           "bio1_3yrs", "bio2_3yrs", "bio3_3yrs", "bio7_3yrs", "bio14_3yrs", "bio15_3yrs",
           "pr_spring_3yrs", "pr_summer_3yrs", "pr_autumn_3yrs", "pr_winter_3yrs", 
           "sum_annual_crops", "secdf","pastr", "urban",
           "sum_annual_crops_3yrs", "secdf_3yrs", "pastr_3yrs", "urban_3yrs"))


# routes-years:
load(file = file.path("data", "BBS_for_occ_selection.RData")) # route_sel_dt; output of 1_3_match_BBS_to_env_data.R 

# selected routes spatial data (to buffer presences):
routes_sel_sf <- st_read(file.path("data", "route_selection_1991_2015_surv_beg_end_max_5y_miss_v2_spat_thin_100km_max_30_r_per_BCR_centroids.shp")) # output of 1_1_route_selection.R

# route-year-species information (only surveyed)
load(file = file.path("data", "BBS_for_occ_spec_records.RData")) # bbs_dt_occ; output of 1_0_reformat_BBS_data.R

# selected species:
load(file = file.path("data", "final_species_selection_eco_sorted.RData")) # final_species_eco_sorted; output of 1_2_species_selection.R
species_discard <- c("Golden Eagle", "Prairie Falcon", "Sharp-shinned Hawk", "Broad-winged Hawk",
                     "Cooper's Hawk", "Osprey", "Evening Grosbeak", "Golden-crowned Kinglet", "Olive-sided Flycatcher",
                     "Wilson's Warbler", "Bullock's Oriole", "Western Wood-Pewee", "Black-throated Gray Warbler",
                     "MacGillivray's Warbler", "Ruffed Grouse", "Northern Bobwhite",
                     "Common Raven", "Common Yellowthroat", "Hairy Woodpecker", "Indigo Bunting",
                     "Louisiana Waterthrush", "Northern Mockingbird", "Northern Parula", "Prairie Warbler",
                     "Bewick's Wren", "Brewer's Blackbird")
species_set <- final_species_eco_sorted[!final_species_eco_sorted %in% species_discard]


# detrend environmental data: ----

spec <- "Black Vulture"

# relevant routes for the species, within distance of 750 km of presences:
rel_routes <- training_routes(species = spec, buffer_km = 750, output = "RTENOs")

env_dt_spec <- route_sel_env_dt_final %>% 
  filter(RTENO %in% rel_routes)

# load file with trends in environmental data: (output of 1_3_explore_occ_data.R)
lm_env_trend_df <- read.csv(file = file.path("data", "env_vars_trends_species.csv"))
# reformat:
lm_env_trend_df_reformat <- lm_env_trend_df %>% 
  tidyr::pivot_longer(cols = !species, names_to = "variable", values_to = "values") %>%
  mutate(value = ifelse(grepl(pattern = "trend", x = variable), "slope", "pvalue")) %>%
  mutate(variable = gsub("(_trend)|(_p_value)", "", x = variable)) %>% 
  tidyr::pivot_wider(names_from = value, values_from = values)


# iterate over variables:

# detrend one variable at a time:

#var <- "bio1"

# register cores for parallel computation:
ncores <- 28 # models * 4 chains?
cl <- makeCluster(ncores, setup_timeout = 0.5)
registerDoParallel(cl)

foreach(var = c("bio1", "bio3", "pr_summer", "pr_winter", "sum_annual_crops", "secdf", "urban"),
        .packages = c("dplyr", "collapse", "flocker", "cmdstanr", "brms", "sf"), # xx
        .errorhandling = "pass", #"remove",
        .verbose = TRUE) %dopar% {

#for(var in c("bio1", "bio3", "pr_summer", "pr_winter", "sum_annual_crops", "secdf", "urban")){
  
  prog_log_file <- file(file.path("logfiles", "attribution", paste0(spec, "_", var, "_progress.txt")), open = "wt") # write console output here
  sink(prog_log_file, type = "message")
  sink(prog_log_file, type = "output")
  
  print(spec)
  print(var)
  
  print("detrending")
  
  # get linear trend:
  lm_trend_var <- lm_env_trend_df_reformat %>% 
    filter(species == spec) %>% 
    filter(variable == var) %>% 
    pull(slope)
  
  # remove linear trend from data:
  env_dt_spec_detr <- env_dt_spec %>% 
    mutate("{var}" := get(var) - lm_trend_var*(Year - 1991))
  
  # plots to test:
  env_dt_spec %>%
    ggplot(aes(x = Year, y = pr_winter)) +
    geom_line(aes(group = RTENO)) +
    geom_smooth(method = 'lm')
  env_dt_spec_detr %>%
    ggplot(aes(x = Year, y = pr_winter)) +
    geom_line(aes(group = RTENO)) +
    geom_smooth(method = 'lm')
  
  # scale resulting environmental variables:
  env_dt_spec_detr_scaled <- env_dt_spec_detr %>% 
    mutate(across(!c(RTENO, Year), ~ (scale(.)) %>% as.vector()))
  
  # refit model: ----
  
  print("format data for model")
  
  buffer_km <- 750 # 250
  nyears <- length(unique(env_dt_spec_detr_scaled$Year)) # 25
  nsurveys <- 5
  
  # species presences-absences:
  occ_dt_spec <- BBS_pres_abs_spec(species = spec)
  
  # relevant routes, within distance of 750 km of species records:
  rel_routes <- training_routes(species = spec, buffer_km = buffer_km, output = "RTENOs")
  
  # reformat obs. as array sites x surveys x years:
  
  years <- seq(min(occ_dt_spec$Year), max(occ_dt_spec$Year))
  nsites <- length(rel_routes)
  
  y_array <- array(NA, dim = c(nsites, nsurveys, nyears))
  for (t in 1:nyears){
    y_array[1:nsites, 1:nsurveys, t] <- as.matrix(occ_dt_spec[which(occ_dt_spec$Year == years[t] & occ_dt_spec$RTENO %in% rel_routes), 
                                                              c(paste0("Count", seq(10, 50, 10)))])
  }
  
  # reformat environmental covariates:
  env_cov <- vector("list", length = nyears)
  for (t in 1:nyears){
    env_cov[[t]] <- env_dt_spec_detr_scaled[which(env_dt_spec_detr_scaled$Year == years[t]), 
                                            c("bio1", "bio2", "bio3", "bio7", "bio14", "bio15", 
                                              "pr_spring", "pr_summer","pr_autumn", "pr_winter",
                                              "bio1_3yrs", "bio2_3yrs", "bio3_3yrs", "bio7_3yrs", "bio14_3yrs", "bio15_3yrs",
                                              "pr_spring_3yrs", "pr_summer_3yrs", "pr_autumn_3yrs", "pr_winter_3yrs", 
                                              "sum_annual_crops", "secdf","pastr", "urban",
                                              "sum_annual_crops_3yrs", "secdf_3yrs", "pastr_3yrs", "urban_3yrs")]
  }
  
  # covariate for detection probability:
  det_cov <- vector("list", length = 1)
  names(det_cov) <- "route_section"
  det_cov$route_section <- array(NA, dim = c(nsites, nsurveys, nyears))
  det_cov$route_section[ , , 1:nyears] <- matrix(rep(c("Sect1", "Sect2", "Sect3", "Sect4", "Sect5"), nsites), nsites, byrow = TRUE)
  
  # make flocker data:
  
  fd <- make_flocker_data_dynamic(
    obs = y_array,
    unit_covs = env_cov, 
    event_covs = det_cov, 
    quiet = TRUE
  )
  
  # fit model: 
  start.time <- Sys.time()
  
  print("fit model")
  
  out <- flock(
    f_occ = ~ bio1_3yrs + bio2_3yrs + bio3_3yrs + bio7_3yrs + bio14_3yrs + bio15_3yrs +
      pr_spring_3yrs + pr_summer_3yrs + pr_autumn_3yrs + pr_winter_3yrs  + 
      I(bio1_3yrs^2) + I(bio2_3yrs^2) + I(bio3_3yrs^2) + I(bio7_3yrs^2) + I(bio14_3yrs^2) + I(bio15_3yrs^2) + 
      I(pr_spring_3yrs^2) + I(pr_summer_3yrs^2) + I(pr_autumn_3yrs^2) + I(pr_winter_3yrs^2)  +
      sum_annual_crops_3yrs + secdf_3yrs + pastr_3yrs + urban_3yrs +
      I(sum_annual_crops_3yrs^2) + I(secdf_3yrs^2) + I(pastr_3yrs^2) + I(urban_3yrs^2),
    f_det = ~ route_section,
    f_col = ~ bio1 + bio2 + bio3 + bio7 + bio14 + bio15 +
      pr_spring + pr_summer + pr_autumn + pr_winter +
      I(bio1^2) + I(bio2^2) + I(bio3^2) + I(bio7^2)+ I(bio14^2)+ I(bio15^2) +
      I(pr_spring^2) + I(pr_summer^2) + I(pr_autumn^2) + I(pr_winter^2)  + 
      sum_annual_crops + secdf + pastr + urban +
      I(sum_annual_crops^2) + I(secdf^2) + I(pastr^2) + I(urban^2),
    f_ex = ~ bio1 + bio2 + bio3 + bio7 + bio14 + bio15 +
      pr_spring + pr_summer + pr_autumn + pr_winter +
      I(bio1^2) + I(bio2^2) + I(bio3^2) + I(bio7^2)+ I(bio14^2)+ I(bio15^2) +
      I(pr_spring^2) + I(pr_summer^2) + I(pr_autumn^2) + I(pr_winter^2)  + 
      sum_annual_crops + secdf + pastr + urban +
      I(sum_annual_crops^2) + I(secdf^2) + I(pastr^2) + I(urban^2),
    flocker_data = fd,
    prior = c(brms::set_prior("logistic(0,1)", class = "Intercept") + # flat on probability scale (https://cran.r-project.org/web/packages/flocker/vignettes/flocker_tutorial.html)
                brms::set_prior("logistic(0,1)", class = "Intercept", dpar = "occ"),
              brms::set_prior("logistic(0,1)", class = "Intercept", dpar = "colo"),
              brms::set_prior("logistic(0,1)", class = "Intercept", dpar = "ex"),
              brms::set_prior("normal(0,2)", class = "b"),
              brms::set_prior("normal(0,2)", dpar = "occ", class = "b"),
              brms::set_prior("normal(0,2)", dpar = "colo", class = "b"),
              brms::set_prior("normal(0,2)", dpar = "ex", class = "b")),
    multiseason = "colex",
    multi_init = "explicit",
    backend = "cmdstanr",
    cores = 4,
    chains = 4,
    warmup = 1000, # first round: 500
    iter = 1000 + 1000 # first round: 500 + 1000
  )
  
  print(out)
  
  save(out, file = file.path(dir, "results", "attribution", paste0("out_", spec, "_", var, ".RData")))
  
  end.time <- Sys.time()
  print(round(end.time - start.time, 2))
  
  print("postprocessing")
  
  start.time <- Sys.time()
  
  # calculate further results:
  fitted_initocc_col_ex <- fitted_flocker(out) 
  occ_posterior <- get_Z(out, history_condition = FALSE) # default
  y_predictions <- predict_flocker(out, history_condition = FALSE) # default, necessary for validation?
  loo_cv <- loo_flocker(out, thin = NULL)
  
  res_list <- list("fitted" = fitted_initocc_col_ex,
                   "occ_posterior" = occ_posterior,
                   "y_preds" = y_predictions,
                   "loo_cv" = loo_cv)
  
  print(res_list$loo_cv)
  
  high_k_routes <- loo::pareto_k_ids(loo_cv)
  high_k_RTENO <- rel_routes[high_k_routes]
  print(paste("RTENO with too high pareto k:", paste(high_k_RTENO, collapse = ", ")))
  
  rm(fitted_initocc_col_ex, occ_posterior, y_predictions, loo_cv)
  
  save(res_list, file = file.path(dir, "results", "attribution", paste0("postproc_", spec, "_", var, ".RData")))
  
  end.time <- Sys.time()
  print(round(end.time - start.time, 2))
  
  sink(type="message")
  sink(type="output")
}


#---

load(file.path(results_dir, "postproc_Black Vulture_bio1.RData"))
res_bio1 <- res_list

load(file.path("M:", "Documents", "DEBTs", "analysis", "Schifferle_BBS_occupancy_models_2023", "results", "full_model", "final_run_750km", "postproc_Black Vulture_fm_buffer750.RData"))
res_fm <- res_list

res_bio1$loo_cv
res_fm$loo_cv
