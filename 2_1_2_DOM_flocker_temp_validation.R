# problem: in cross validation we blocked in space, not in time
# here: "block" in time:

# runs:
# 1) fit model with data of first 20 years, predict to last 5 years (all routes)
# 2) fit model with data of first 15 years, predict to last 10 years (all routes)
# both first for 750 km buffer, then for 250 km buffer!


# how many years in training data:
n_train_years <- 20 # 15
# buffer around presences within which routes are considered:
buffer_km <- 750 # 250


# packages: ----

library(dplyr)
library(doParallel)
library(flocker)
library(cmdstanr)
set_cmdstan_path(path = NULL)
#set_cmdstan_path("C:/Users/schifferle1/Documents/cmdstan-2.34.1") # xx
library(sf)

# directories: ----

print(tempdir())
dir <- file.path("/import", "ecoc9z", "data-zurell", "schifferle", "BBS_occupancy_models_2023")
#dir <- getwd()

# functions: ----

source("0_functions.R")


# load data: ----

# selected routes and focal years matched to environmental data:
# merged route, year, environment data:

load(file = file.path("data", "route_year_env_data.RData")) # route_sel_env_dt_final; output of 1_3_match_BBS_to_env_data.R
# scale covariates:
route_sel_env_dt_scaled <- route_sel_env_dt_final %>% 
  mutate(across(bio2:pr_winter_3yrs, ~ (scale(.)) %>% as.vector())) %>% 
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
#load(file = file.path("data", "final_species_selection.RData")) # species_selection_final; output of 1_2_species_selection.R
# sorted by ecoregion:
load(file = file.path("data", "final_species_selection_eco_sorted.RData")) # final_species_eco_sorted; output of 1_2_species_selection.R


# register cores for parallel computation:
ncores <- 37 # species * 4 chains? 
cl <- makeCluster(ncores, setup_timeout = 0.5)
registerDoParallel(cl)



# make blocks:
# 9 species per block:
spec_blocks_list <- split(final_species_eco_sorted, c(rep(1:19, each = 9), rep(20, 3)))
names(spec_blocks_list) <- NULL

# log progress:
prog_log_file <- file(paste0("temp_val/", "temp_val_", ifelse(n_train_years == 20, "5yrs", "10yrs"), "_buffer_", buffer_km, "_progress.txt"), open = "wt") # write console output here
sink(prog_log_file, type = "message")
sink(prog_log_file, type = "output")

for(i in 1:length(spec_blocks_list)){

  print(paste("block", i, "of", length(spec_blocks_list)))
  
  # iterate over species:
  foreach(spec = spec_blocks_list[[i]],
          .packages = c("dplyr", "collapse", "flocker", "cmdstanr", "brms", "sf"), # xx
          .errorhandling = "pass", #"remove",
          .verbose = TRUE) %dopar% {
            
            log_file_spec <- file(paste0("temp_val/", "temp_val_", ifelse(n_train_years == 20, "5yrs", "10yrs"), "_buffer_", buffer_km, "_", spec, ".txt"), open = "wt") # write console output here
            sink(log_file_spec, type = "message")
            sink(log_file_spec, type = "output")

            print(spec)

            
            # 1) train with first 20 (15) years, validate with last 5 (10) years: ----
            
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
              filter(Year <= if_else(n_train_years == 20, 2010, 2005))
            
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
              iter = 1000 + 1000 # 250 + 1000 # including warmup
            )
            
            print(out)
            
            # save fitted model:
            save(out, file = file.path(dir, "data", "temp_val", paste0("out_", spec, "_temp_val_", ifelse(n_train_years == 20, "5yrs", "10yrs"), ".RData")))
            
            end.time <- Sys.time()
            time.taken <- round(end.time - start.time, 2)
            print(time.taken)

            
            # predict to all env. data (data used in training + data of last 5 years): ----
            
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
                             "y_preds" = y_predictions,
                             "y_preds_mean" = y_predictions_mean)
            
            save(res_list, file = file.path(dir, "data", "temp_val",  paste0("test_preds_", spec, "_temp_val_", ifelse(n_train_years == 20, "5yrs", "10yrs"), ".RData")))

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
