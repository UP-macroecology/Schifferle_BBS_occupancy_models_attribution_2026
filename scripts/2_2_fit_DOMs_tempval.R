# Model fitting for temporal validation:
# fit model with data of first 15 years
# predict y (observations) for all 25 years, only last ten years used later for model evaluation

# executed once (round 1), then MCMC checked with 2_3a_fit_DOMs_check_fit.R, 
# then executed again (round 2) for species with issues with larger number of iterations
# -> is this round 1 or 2:

# round <- 1
round <- 2

# packages: --------------------------------------------------------------------

library(dplyr)
library(sf)
library(doParallel)
library(flocker)
library(cmdstanr) 


# directories: -----------------------------------------------------------------

set_cmdstan_path(path = NULL) # for HPC; local: set_cmdstan_path("C:/Users/schifferle1/Documents/cmdstan-2.34.1")

print(tempdir())

# project directory:
dir <- file.path("//NAS-2-P-SN-01.ibb.uni-potsdam.de", "daten$", "AG26", "Transfer", "Schifferle_BBS_occupancy_models_2023")
#dir <- file.path("/mnt", "ibb_share", "zurell_transfer", "Schifferle_BBS_occupancy_models_2023")


# functions: -------------------------------------------------------------------

source(file.path("scripts", "0_functions.R"))


# load data: -------------------------------------------------------------------

# route-year-species information (only surveyed)
load(file = file.path("data", "BBS_for_occ_spec_records.RData")) # bbs_dt_occ; output of 1_0_dataprep_BBS_bird_data.R

# selected routes spatial data (to buffer presences):
routes_sel_sf <- st_read(file.path("data", "route_selection_1995_2019_surv_beg_end_max_5y_miss_v2_spat_thin_100km_max_30_r_per_BCR_centroids.shp")) # output of 1_1_dataprep_BBS_route_selection.R

# selected species:
load(file = file.path("data", "final_species_selection.RData")) # species_selection_final; output of 1_2_dataprep_BBS_species_selection.R

# selected variables:
load(file = file.path("data", "selected_variables.RData")) # selvar_final; output of 1_2a_dataprep_env_variable_selection.R

# routes-years:
load(file = file.path("data", "BBS_for_occ_selection.RData")) # route_sel_dt; output of 1_3_dataprep_match_BBS_routes_env_data.R

# selected routes and focal years matched to environmental data:
# merged route, year, environment data:
load(file = file.path("data", "route_year_env_data.RData")) # route_sel_env_dt_final; output 1_3_dataprep_match_BBS_routes_env_data.R


# settings: --------------------------------------------------------------------

if(round == 1){
  
  # directory for logfiles:
  log_dir <- file.path("logfiles", "temp_val_buffer_750_10yrs")
  
  # directory for results:
  res_dir <- file.path(dir, "results", "temp_val_buffer_750_10yrs")
  
  # species to fit models for:
  species_set <- species_selection_final
  
  # fitting iterations:
  iterations <- 1000
}


if(round == 2){
  
  # directory for logfiles:
  log_dir <- file.path("logfiles", "temp_val_buffer_750_10yrs", "refit_2000_2000")
  
  # directory for results:
  res_dir <- file.path(dir, "results", "temp_val_buffer_750_10yrs", "refit_2000_2000")
  
  # load species for which MCMC with 1000 + 1000 iterations failed (specs_MCMC_failed; output of 2_3a_fit_DOMs_check_fit.R)
  load(file.path(dir, "results", "temp_val_buffer_750_10yrs", "check_output", "specs_MCMC_failed.RData"))
  
  # species to fit models for:
  species_set <- specs_MCMC_failed
  
  # fitting iterations:
  iterations <- 2000
}


# assemble overall data: -------------------------------------------------------

# size of buffer area around presences considered for fitting DOMs:
buffer_km <- 750 # 250

# how many years in training data:
n_train_years <- 15

# scale covariates:

route_sel_env_dt_scaled <- route_sel_env_dt_final %>% 
  select(-c(Latitude, Longitude, BCR, ObsN, doy)) %>% 
  mutate(across(!c(RTENO, Year, Surveyed), ~ as.numeric(scale(., center=mean(.), scale = sd(.)))))

rm(route_sel_env_dt_final)


# fit models: ------------------------------------------------------------------

# register cores for parallel computation:
ncores <- 40 # species * 4 chains? 
cl <- makeCluster(ncores, setup_timeout = 0.5)
registerDoParallel(cl)


# iterate over species:

foreach(spec = species_set,
        .packages = c("dplyr", "collapse", "flocker", "cmdstanr", "brms", "sf"), 
        .errorhandling = "pass", #"remove",
        .verbose = TRUE) %dopar% {
          
          # log file:
          spec_log_file <- file(file.path(log_dir, paste0("out_", spec, "_temp_val_10yrs_buffer", buffer_km, ".txt")), open = "wt") # write console output here
          sink(spec_log_file, type = "message")
          sink(spec_log_file, type = "output")
          
          print(spec)
          
          # check whether species has run already:
          model_run <- file.exists(file.path(res_dir, paste0("out_", spec, "_temp_val_10yrs_buffer_", buffer_km, ".RData")))
          post_proc_run <- file.exists(file.path(res_dir,  paste0("preds_", spec, "_temp_val_10yrs_buffer_", buffer_km, ".RData")))
          
          if(model_run & post_proc_run) {
            print(paste(spec, "ran already."))
            next
          }
          
          # assemble data:
          # train with first 15 years, validate with last 10 years
          
          # species presence - absence:
          occ_dt_spec <- BBS_pres_abs_spec(species = spec)
          
          # relevant routes, within distance of 750 km of species records:
          rel_routes <- training_routes(species = spec, buffer_km = buffer_km, output = "RTENOs")
          
          print(paste("training routes:", length(rel_routes)))
          
          # subset data:
          occ_dt_spec_subset_train <- occ_dt_spec %>% 
            filter(RTENO %in% rel_routes) %>% 
            filter(Year <= 2019-10)
          
          # reformat observations as array sites x surveys x years:
          
          years <- seq(min(occ_dt_spec_subset_train$Year), max(occ_dt_spec_subset_train$Year))
          nyears <- length(years)
          nsites <- length(rel_routes)
          nsurveys <- 5
          
          y_array <- array(NA, dim = c(nsites, nsurveys, nyears))
          for (t in 1:nyears){
            y_array[1:nsites, 1:nsurveys, t] <- as.matrix(occ_dt_spec_subset_train[which(occ_dt_spec_subset_train$Year == years[t]), 
                                                                      c(paste0("Count", seq(10, 50, 10)))])
          }
          
          # reformat environmental covariates:
          
          env_cov <- vector("list", length = nyears)
          for (t in 1:nyears){
            env_cov[[t]] <- route_sel_env_dt_scaled[which(route_sel_env_dt_scaled$Year == years[t] & route_sel_env_dt_scaled$RTENO %in% rel_routes), 
                                                    c(selvar_final, paste0(selvar_final, "_3yrs"))]
          }
          
          # covariate for detection probability:
          
          det_cov <- vector("list", length = 1)
          names(det_cov) <- "route_section"
          det_cov$route_section <- array(NA, dim = c(nsites, nsurveys, nyears))
          det_cov$route_section[ , , 1:nyears] <- matrix(rep(c("Sect1", "Sect2", "Sect3", "Sect4", "Sect5"), nsites), nsites, byrow = TRUE)
          
          # make flocker data:
          
          fd_train_yrs <- make_flocker_data_dynamic(
            obs = y_array,
            unit_covs = env_cov, 
            event_covs = det_cov, 
            quiet = TRUE
          )

          # fit model:
          
          print("model fitting")
          
          start.time <- Sys.time()
          
          out <- flock(
            f_occ = ~ bio1_3yrs + bio2_3yrs + bio3_3yrs + bio7_3yrs + bio14_3yrs + bio15_3yrs +
              pr_mean_spring_3yrs + pr_mean_summer_3yrs + pr_mean_autumn_3yrs + pr_mean_winter_3yrs  + 
              I(bio1_3yrs^2) + I(bio2_3yrs^2) + I(bio3_3yrs^2) + I(bio7_3yrs^2) + I(bio14_3yrs^2) + I(bio15_3yrs^2) + 
              I(pr_mean_spring_3yrs^2) + I(pr_mean_summer_3yrs^2) + I(pr_mean_autumn_3yrs^2) + I(pr_mean_winter_3yrs^2)  +
              urbanareas_3yrs + managed_pastures_3yrs + primary_nonforests_3yrs + secondary_nonforests_3yrs + sum_annual_crops_3yrs +
              I(urbanareas_3yrs^2) + I(managed_pastures_3yrs^2) + I(primary_nonforests_3yrs^2) + I(secondary_nonforests_3yrs^2) + I(sum_annual_crops_3yrs^2),
            f_det = ~ route_section,
            f_col = ~ bio1 + bio2 + bio3 + bio7 + bio14 + bio15 +
              pr_mean_spring + pr_mean_summer + pr_mean_autumn + pr_mean_winter +
              I(bio1^2) + I(bio2^2) + I(bio3^2) + I(bio7^2)+ I(bio14^2)+ I(bio15^2) +
              I(pr_mean_spring^2) + I(pr_mean_summer^2) + I(pr_mean_autumn^2) + I(pr_mean_winter^2)  + 
              urbanareas + managed_pastures + primary_nonforests + secondary_nonforests + sum_annual_crops +
              I(urbanareas^2) + I(managed_pastures^2) + I(primary_nonforests^2) + I(secondary_nonforests^2) + I(sum_annual_crops^2),
            f_ex = ~ bio1 + bio2 + bio3 + bio7 + bio14 + bio15 +
              pr_mean_spring + pr_mean_summer + pr_mean_autumn + pr_mean_winter +
              I(bio1^2) + I(bio2^2) + I(bio3^2) + I(bio7^2)+ I(bio14^2)+ I(bio15^2) +
              I(pr_mean_spring^2) + I(pr_mean_summer^2) + I(pr_mean_autumn^2) + I(pr_mean_winter^2)  + 
              urbanareas + managed_pastures + primary_nonforests + secondary_nonforests + sum_annual_crops +
              I(urbanareas^2) + I(managed_pastures^2) + I(primary_nonforests^2) + I(secondary_nonforests^2) + I(sum_annual_crops^2),
            flocker_data = fd_train_yrs,
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
            warmup = iterations,
            iter = iterations + iterations
          )
          
          print(out)
          
          # save fitted model:
          save(out, file = file.path(res_dir, paste0("out_", spec, "_temp_val_10yrs_buffer_", buffer_km, ".RData")))
          
          end.time <- Sys.time()
          print(round(end.time - start.time, 2))

          
          # predict to all env. data (data used in training + data of last 10 years): ----
          
          print("start predictions")
          
          start.time <- Sys.time()

          # assemble data to predict to:
          
          # only relevant routes within buffer:
          occ_dt_spec_subset <- occ_dt_spec %>% 
            filter(RTENO %in% rel_routes)
          
          years <- seq(min(route_sel_env_dt_scaled$Year), max(route_sel_env_dt_scaled$Year))
          nyears <- length(years)

          # dummy observations as array sites x surveys x years:
          
          # https://cran.r-project.org/web/packages/flocker/vignettes/flocker_tutorial.html:
          # "Note that if predictions are desired at sites without observations, it is acceptable 
          # to pass an array of dummy observations (e.g. all zeros) to make_flocker_data() and 
          # then to set history_condition = FALSE in the call to get_Z()."
          
          y_array <- array(NA, dim = c(nsites, nsurveys, nyears))
          for (t in 1:nyears){
            y_array[1:nsites, 1:nsurveys, t] <- as.matrix(0)
          }
          
          print(paste("dim y_array:", paste(dim(y_array), collapse = ", ")))
          
          # reformat environmental covariates:
          
          env_cov <- vector("list", length = nyears)
          for (t in 1:nyears){
            env_cov[[t]] <- route_sel_env_dt_scaled[which(route_sel_env_dt_scaled$Year == years[t] & route_sel_env_dt_scaled$RTENO %in% rel_routes), 
                                                    c(selvar_final, paste0(selvar_final, "_3yrs"))]
          }
          
          # covariate for detection probability:
          # ("please add a dummy event covariate. You do not need to use this covariate in your model formula"):
          
          det_cov <- vector("list", length = 1)
          names(det_cov) <- "route_section"
          det_cov$route_section <- array(NA, dim = c(nsites, nsurveys, nyears))
          det_cov$route_section[ , , 1:nyears] <- matrix(rep(c("Sect1", "Sect2", "Sect3", "Sect4", "Sect5"), nsites), nsites, byrow = TRUE)
          
          
          # make flocker data:
          # (using the `new_data` argument within get_Z() for a multiseason model requires passing a `flocker_data` object)
          
          fd <- make_flocker_data_dynamic(
            obs = y_array, # no NAs allowed
            unit_covs = env_cov, 
            event_covs = det_cov, 
            quiet = TRUE
          )

          # predicted observations (y: 0 or 1) for all routes and years:
          
          y_predictions <- predict_flocker(out, 
                                           history_condition = FALSE,
                                           new_data = fd) # for each site, route section and year 4000 draws
          
          # discarded (saved, but not used for model evaluation):
          
          # # summarise detections across route sections since we look at route level 
          # y_predictions_sections_sum <- apply(y_predictions, MAR = c(1,3,4), FUN = max) 
          # 
          # # in which proportion of draws is species detected on route:
          # y_predictions_mean <- apply(y_predictions_sections_sum, MAR = c(1,2), FUN = function(x) as.integer(round(mean(x) * 100)))
          #
          # # fitted values for all routes and years:
          # 
          # fitted_initocc_col_ex_det <- fitted_flocker(out, components = c("occ", "col", "ex", "det"), new_data = fd)
          # 
          # # mean and median for each site and year (as integer to reduce storage space):
          # fitted_initocc_col_ex_det_mean <- list(
          #   "linpred_occ" = apply(fitted_initocc_col_ex_det$linpred_occ[,1,1,], MAR = 2, FUN = function(x) as.integer(round(mean(x)*100))), # same across route sections and years
          #   "linpred_col" = apply(fitted_initocc_col_ex_det$linpred_col[,1,,], MAR = c(1,2), FUN = function(x) as.integer(round(mean(x)*100))), # same across route sections
          #   "linpred_ex" = apply(fitted_initocc_col_ex_det$linpred_ex[,1,,], MAR = c(1,2), FUN = function(x) as.integer(round(mean(x)*100))), # same across route sections
          #   "linpred_det" = apply(fitted_initocc_col_ex_det$linpred_det[1,,1,], MAR = 1, FUN = function(x) as.integer(round(mean(x)*100))) # same across years and sites
          # )
          # fitted_initocc_col_ex_det_median <- list(
          #   "linpred_occ" = apply(fitted_initocc_col_ex_det$linpred_occ[,1,1,], MAR = 2, FUN = function(x) as.integer(round(median(x)*100))), 
          #   "linpred_col" = apply(fitted_initocc_col_ex_det$linpred_col[,1,,], MAR = c(1,2), FUN = function(x) as.integer(round(median(x)*100))), 
          #   "linpred_ex" = apply(fitted_initocc_col_ex_det$linpred_ex[,1,,], MAR = c(1,2), FUN = function(x) as.integer(round(median(x)*100))), 
          #   "linpred_det" = apply(fitted_initocc_col_ex_det$linpred_det[1,,1,], MAR = 1, FUN = function(x) as.integer(round(median(x)*100))) 
          # )
          # 
          # # occupancy probability for all routes and years:
          # 
          # occ_posterior <- get_Z(out, history_condition = FALSE, new_data = fd, sample = FALSE)
          # 
          # # mean and median for each site and year (as integer to reduce storage space):
          # occ_posterior_median <- apply(occ_posterior, MARGIN = c(1,2), FUN = function(x) as.integer(round(median(x) * 100)))
          # occ_posterior_mean <- apply(occ_posterior, MARGIN = c(1,2), FUN = function(x) as.integer(round(mean(x) * 100)))


          # save predictions:
          
          res_list <- list(
            "y_preds" = y_predictions,
            #"y_preds_mean" = y_predictions_mean,
            #"fitted_mean" = fitted_initocc_col_ex_det_mean,
            #"fitted_median" = fitted_initocc_col_ex_det_median,
            #"occ_posterior_mean" = occ_posterior_mean,
            #"occ_posterior_median" = occ_posterior_median,
            #"occ_posterior" = occ_posterior
            )
          
          save(res_list, file = file.path(res_dir,  paste0("preds_", spec, "_temp_val_10yrs_buffer_", buffer_km, ".RData")))

          rm(res_list)
          
          end.time <- Sys.time()
          print(paste("finished", spec, "; time taken", round(end.time - start.time, 2)))
          
          sink(type="message")
          sink(type="output")
      }

stopCluster(cl)

rm(list=ls())
gc()