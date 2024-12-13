# in cross validation we blocked in space, not in time
# here: "block" in time:
# fit model with data of first 15 years, predict to last 10 years (all routes)

# packages: ----

library(dplyr)
library(doParallel)
library(flocker)
library(cmdstanr)
set_cmdstan_path(path = NULL)
#set_cmdstan_path("C:/Users/schifferle1/Documents/cmdstan-2.34.1") # xx
library(sf)

# register cores for parallel computation:
ncores <- 40 # species * 4 chains? 
cl <- makeCluster(ncores, setup_timeout = 0.5)
registerDoParallel(cl)


# directories: ----

print(tempdir())
dir <- file.path("/import", "ecoc9z", "data-zurell", "schifferle", "BBS_occupancy_models_2023")
#dir <- getwd()

# directory for logfiles:
log_dir <- file.path("logfiles", "temp_val_buffer750km")
# directory for results:
res_dir <- file.path(dir, "results", "temp_val_buffer750km")

# functions: ----

source("0_functions.R")


# load data: ----

# how many years in training data:
n_train_years <- 15
# buffer around presences within which routes are considered:
buffer_km <- 750 # 250


# selected routes and focal years matched to environmental data:
# merged route, year, environment data:

# selected variables:
load(file = file.path("data", "selected_variables.RData")) # selvar_final; output of 1_2_variable_selection.R

load(file = file.path("data", "route_year_env_data.RData")) # route_sel_env_dt_final; output of 1_3_match_BBS_to_env_data.R

# routes-years:
load(file = file.path("data", "BBS_for_occ_selection.RData")) # route_sel_dt; output of 1_3_match_BBS_to_env_data.R 

# selected routes spatial data (to buffer presences):
routes_sel_sf <- st_read(file.path("data", "route_selection_1995_2019_surv_beg_end_max_5y_miss_v2_spat_thin_100km_max_30_r_per_BCR_centroids.shp")) # output of 1_1_route_selection.R

# route-year-species information (only surveyed)
load(file = file.path("data", "BBS_for_occ_spec_records.RData")) # bbs_dt_occ; output of 1_0_reformat_BBS_data.R

# selected species:
# sorted by ecoregion:
load(file = file.path("data", "final_species_selection_eco_sorted.RData")) # final_species_eco_sorted; output of 1_2_species_selection.R


# assemble overall data: ----

# scale covariates:

route_sel_env_dt_scaled <- route_sel_env_dt_final %>% 
  select(-c(Latitude, Longitude, BCR, ObsN, doy)) %>% 
  mutate(across(!c(RTENO, Year, Surveyed), ~ as.numeric(scale(., center=mean(.), scale = sd(.)))))

rm(route_sel_env_dt_final)


# make blocks:
# 10 species per block
spec_blocks_list <- split(final_species_eco_sorted, c(rep(1:(floor(length(final_species_eco_sorted)/10)), each = 10), 
                                                      rep(ceiling(length(final_species_eco_sorted)/10), length(final_species_eco_sorted) %% 10))) # xx change once it runs
names(spec_blocks_list) <- NULL

# log progress:
prog_log_file <- file(file.path(log_dir, paste0("temp_val_10yrs_buffer_", buffer_km, "_progress.txt")), open = "wt") # write console output here
sink(prog_log_file, type = "message")
sink(prog_log_file, type = "output")

for(i in 1:3){#length(spec_blocks_list)){

  print(paste("block", i, "of", length(spec_blocks_list)))
  
  # iterate over species:
  foreach(spec = spec_blocks_list[[i]],
          .packages = c("dplyr", "collapse", "flocker", "cmdstanr", "brms", "sf"), # xx
          .errorhandling = "pass", #"remove",
          .verbose = TRUE) %dopar% {
            
            
            # check whether species has run already:
            model_run <- file.exists(file.path(res_dir, paste0("out_", spec, "_temp_val_10yrs_buffer_", buffer_km, ".RData")))
            post_proc_run <- file.exists(file.path(res_dir,  paste0("preds_", spec, "_temp_val_10yrs_buffer_", buffer_km, ".RData")))
            
            if(model_run & post_proc_run) {
              print(paste(spec, "ran already."))
              next
            }
            
            spec_log_file <- file(file.path(log_dir, paste0("out_", spec, "_temp_val_10yrs_buffer", buffer_km, ".txt")), open = "wt") # write console output here
            sink(spec_log_file, type = "message")
            sink(spec_log_file, type = "output")

            print(spec)

            
            # 1) train with first 15 years, validate with last 10 years: ----
            
            print("training")
            
            # assemble data:

            # species presence - absence:
            occ_dt_spec <- BBS_pres_abs_spec(species = spec)
            
            # relevant routes, within distance of 750 km of species records:
            rel_routes <- training_routes(species = spec, buffer_km = buffer_km, output = "RTENOs")
            
            print(paste("routes:", length(rel_routes)))
            
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
            
            print(paste("dim y_array:", paste(dim(y_array), collapse = ", ")))
            
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
            
            print("start fitting model")
            
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
              warmup = 1000,
              iter = 1000 + 1000
            )
            
            print(out)
            
            # save fitted model:
            save(out, file = file.path(res_dir, paste0("out_", spec, "_temp_val_10yrs_buffer_", buffer_km, ".RData")))
            
            end.time <- Sys.time()
            time.taken <- round(end.time - start.time, 2)
            print(time.taken)

            
            # predict to all env. data (data used in training + data of last 10 years): ----
            
            print("start predictions")
            
            start.time <- Sys.time()
  
            # prepare data to predict to:
            
            # get_Z(): using the `new_data` argument for a multiseason model requires passing a `flocker_data` object
            
            # Note that if predictions are desired at sites without observations, it is acceptable 
            # to pass an array of dummy observations (e.g. all zeros) to make_flocker_data() and 
            # then to set history_condition = FALSE in the call to get_Z().
            # (https://cran.r-project.org/web/packages/flocker/vignettes/flocker_tutorial.html)
            # "please add a dummy event covariate. You do not need to use this covariate in your model formula":
            
            # subset data:
            occ_dt_spec_subset <- occ_dt_spec %>% 
              filter(RTENO %in% rel_routes)
            
            # reformat obs. as array sites x surveys x years:
            years <- seq(min(route_sel_env_dt_scaled$Year), max(route_sel_env_dt_scaled$Year))
            nyears <- length(years)

            # dummy obs. as array sites x surveys x years:
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
            det_cov <- vector("list", length = 1)
            names(det_cov) <- "route_section"
            det_cov$route_section <- array(NA, dim = c(nsites, nsurveys, nyears))
            det_cov$route_section[ , , 1:nyears] <- matrix(rep(c("Sect1", "Sect2", "Sect3", "Sect4", "Sect5"), nsites), nsites, byrow = TRUE)
            
            # make flocker data:
            fd <- make_flocker_data_dynamic(
              obs = y_array,# if NAs in here I don't get predictions!
              unit_covs = env_cov, 
              event_covs = det_cov, 
              quiet = TRUE
            )

            
            # fitted values for all routes and years (not necessary):
            
            fitted_initocc_col_ex_det <- fitted_flocker(out,
                                                        components = c("occ", "col", "ex", "det"),
                                                        new_data = fd)
            
            # mean and median for each site and year (as integer to reduce storage space):
            fitted_initocc_col_ex_det_mean <- list(
              "linpred_occ" = apply(fitted_initocc_col_ex_det$linpred_occ[,1,1,], MAR = 2, FUN = function(x) as.integer(round(mean(x)*100))), # same across route sections and years
              "linpred_col" = apply(fitted_initocc_col_ex_det$linpred_col[,1,,], MAR = c(1,2), FUN = function(x) as.integer(round(mean(x)*100))), # same across route sections
              "linpred_ex" = apply(fitted_initocc_col_ex_det$linpred_ex[,1,,], MAR = c(1,2), FUN = function(x) as.integer(round(mean(x)*100))), # same across route sections
              "linpred_det" = apply(fitted_initocc_col_ex_det$linpred_det[1,,1,], MAR = 1, FUN = function(x) as.integer(round(mean(x)*100))) # same across years and sites
            )
            fitted_initocc_col_ex_det_median <- list(
              "linpred_occ" = apply(fitted_initocc_col_ex_det$linpred_occ[,1,1,], MAR = 2, FUN = function(x) as.integer(round(median(x)*100))), # same across route sections and years
              "linpred_col" = apply(fitted_initocc_col_ex_det$linpred_col[,1,,], MAR = c(1,2), FUN = function(x) as.integer(round(median(x)*100))), # same across route sections
              "linpred_ex" = apply(fitted_initocc_col_ex_det$linpred_ex[,1,,], MAR = c(1,2), FUN = function(x) as.integer(round(median(x)*100))), # same across route sections
              "linpred_det" = apply(fitted_initocc_col_ex_det$linpred_det[1,,1,], MAR = 1, FUN = function(x) as.integer(round(median(x)*100))) # same across years and sites
            )
            
            
            # occupancy probability for all routes and years:
            
            occ_posterior <- get_Z(out, history_condition = FALSE, new_data = fd, sample = FALSE)
            print(paste("dim occ. posterior", paste(dim(occ_posterior), collapse = ", ")))
            
            # mean and median for each site and year (as integer to reduce storage space):
            occ_posterior_median <- apply(occ_posterior, MARGIN = c(1,2), FUN = function(x) as.integer(round(median(x) * 100)))
            occ_posterior_mean <- apply(occ_posterior, MARGIN = c(1,2), FUN = function(x) as.integer(round(mean(x) * 100)))
            
            
            # predicted observations (y: 0 or 1) for all routes and years:
            
            y_predictions <- predict_flocker(out, 
                                             history_condition = FALSE,
                                             new_data = fd) # for each site, route section and year 4000 draws
            
            y_predictions_sections_sum <- apply(y_predictions, MAR = c(1,3,4), FUN = max) # summarise detections across route sections since we look at route level
            
            # in which proportion of draws is species detected on route:
            y_predictions_mean <- apply(y_predictions_sections_sum, MAR = c(1,2), FUN = function(x) as.integer(round(mean(x) * 100)))
            
            
            # save results:
            res_list <- list("fitted_mean" = fitted_initocc_col_ex_det_mean,
                             "fitted_median" = fitted_initocc_col_ex_det_median,
                             "occ_posterior_mean" = occ_posterior_mean,
                             "occ_posterior_median" = occ_posterior_median,
                             "occ_posterior" = occ_posterior,
                             "y_preds" = y_predictions,
                             "y_preds_mean" = y_predictions_mean)
            
            save(res_list, file = file.path(res_dir,  paste0("preds_", spec, "_temp_val_10yrs_buffer_", buffer_km, ".RData")))

            rm(res_list)
            
            end.time <- Sys.time()
            print(paste("finished", spec, "; time taken", round(end.time - start.time, 2)))
            
            sink(type="message")
            sink(type="output")
        }

sink(type="message")
sink(type="output")

}

stopCluster(cl)

rm(list=ls())
gc()