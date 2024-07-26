# 5-fold spatial cross validation to assess predictive ability of dynamic occupancy models:

# folds created with package blockCV: 2_2_DOM_spatialCV.R

# 1.) for species load fold assignment
# 2.) refit model 5 times (5-fold cross validation)
# 3.) predict to test folds: occupancy probability per year and predicted y (0/1) (this is comparable to observed data)

# first for 750 km buffer, later repeat for 250 km buffer!

# packages: ----

library(dplyr)
library(doParallel)
library(flocker)
#install.packages("cmdstanr", repos = c("https://mc-stan.org/r-packages/", getOption("repos")))
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
load(file = file.path("data", "final_species_selection.RData")) # species_selection_final; output of 1_2_species_selection.R


# assemble data: ----

nyears <- length(unique(route_sel_env_dt_final$Year)) # 25
nsurveys <- 5
nsites <- length(unique(route_sel_env_dt_final$RTENO)) # 476


# register cores for parallel computation:
ncores <- 25 # fold * 4 chains? 
cl <- makeCluster(ncores, setup_timeout = 0.5)
registerDoParallel(cl)

# for species block...

# make blocks:
# 3 species per block:
#spec_blocks_list <- split(species_selection_final, rep(1:(length(species_selection_final)/2), each = 2))
#names(spec_blocks_list) <- NULL

#for(i in 1:length(spec_blocks_list)){
for(i in 1:20){#length(species_selection_final)){
 
  print(paste(i, "of", length(species_selection_final)))
  
  spec <- species_selection_final[i]
  
  sink(paste0(spec, "_CV_fitting.txt")) # write console output here
  sink(type = "message")

  print(spec)
  
  # relevant routes, within distance of 750 km of species records:
  rel_routes <- training_routes(species = spec, buffer_km = 750, output = "RTENOs")
  
  # load fold assignment:
  load(file.path("data", "CV_route_block_allocation", "block_size_500km", paste0(spec, ".RData")))
  
  
  # assemble data: ----
  
  occ_dt_spec <- BBS_pres_abs_spec(species = spec)
  
  
  foreach(fold = 1:5,
          .packages = c("dplyr", "collapse", "flocker", "cmdstanr", "brms"), # xx
          .errorhandling = "pass", #"remove",
          .verbose = TRUE) %dopar% {
            
            sink(paste0(spec, "_CV_fitting_fold", fold, ".txt")) # write console output here
            sink(type = "message")
            
            print(spec)
            print(paste("fold", fold))
            
            # data in current fold:
            training_RTENOs <- rel_routes[sb_US$folds_list[[fold]][[1]]] # training data fold 
            test_RTENOs <- rel_routes[sb_US$folds_list[[fold]][[2]]] # test data fold 
  
            occ_dt_spec_subset <- occ_dt_spec %>% 
              filter(RTENO %in% training_RTENOs)
            
            # reformat obs. as array sites x surveys x years:
            
            years <- seq(min(occ_dt_spec_subset$Year), max(occ_dt_spec_subset$Year))
            nsites <- length(unique(occ_dt_spec_subset$RTENO))
            
            y_array <- array(NA, dim = c(nsites, nsurveys, nyears))
            for (t in 1:nyears){
              y_array[1:nsites, 1:nsurveys, t] <- as.matrix(occ_dt_spec_subset[which(occ_dt_spec_subset$Year == years[t]), 
                                                                        c(paste0("Count", seq(10, 50, 10)))])
            }
            
            # reformat environmental covariates:
            
            env_cov <- vector("list", length = nyears)
            for (t in 1:nyears){
              env_cov[[t]] <- route_sel_env_dt_scaled[which(route_sel_env_dt_scaled$Year == years[t] & route_sel_env_dt_scaled$RTENO %in% occ_dt_spec_subset$RTENO), 
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
            
            
            # fit model: ----
  
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
              warmup = 250,
              iter = 250 + 1000
            )
            
            
            # save fitted model for each fold:
            save(out, file = file.path(dir, "data", paste0("out_", spec, "_CV_fold", fold, ".RData")))
            
            end.time <- Sys.time()
            time.taken <- round(end.time - start.time,2)
            print(time.taken)
            
            
            # calculate further results / predictions to test routes:
            
            print("start predictions")
            
            start.time <- Sys.time()
            
            # data for test folds:
            test_env_data <- route_sel_env_dt_scaled %>%
              filter(RTENO %in% test_RTENOs)
            
            # fitted values for test routes (not necessary):
            fitted_initocc_col_ex <- fitted_flocker(out,
                                                    components = c("occ", "col", "ex"),
                                                    new_data = as.data.frame(test_env_data))
            
            # occupancy probability for test data:
            
            # prepare data for get_Z():
            
            # get_Z(): using the `new_data` argument for a multiseason model requires passing a `flocker_data` object
            
            # Note that if predictions are desired at sites without observations, it is acceptable 
            # to pass an array of dummy observations (e.g. all zeros) to make_flocker_data() and 
            # then to set history_condition = FALSE in the call to get_Z().
            # (https://cran.r-project.org/web/packages/flocker/vignettes/flocker_tutorial.html)
            # "please add a dummy event covariate. You do not need to use this covariate in your model formula":
            
            years_test <- seq(min(test_env_data$Year), max(test_env_data$Year))
            nsites_test <- length(unique(test_env_data$RTENO))
            
            # reformat obs. as array sites x surveys x years:
            y_array_dummy_test <- array(NA, dim = c(nsites_test, nsurveys, nyears))
            for (t in 1:nyears){
              y_array_dummy_test[1:nsites_test, 1:nsurveys, t] <- as.matrix(0)
            }
  
            # covariate for detection probability:
            det_cov_test <- vector("list", length = 1)
            names(det_cov_test) <- "route_section"
            det_cov_test$route_section <- array(NA, dim = c(nsites_test, nsurveys, nyears))
            det_cov_test$route_section[ , , 1:nyears] <- matrix(rep(c("Sect1", "Sect2", "Sect3", "Sect4", "Sect5"), nsites_test), nsites_test, byrow = TRUE)
            
            
            # reformat environmental covariates:
            env_cov_test <- vector("list", length = nyears)
            for (t in 1:nyears){
              env_cov_test[[t]] <- route_sel_env_dt_scaled[which(route_sel_env_dt_scaled$Year == years[t] & route_sel_env_dt_scaled$RTENO %in% test_RTENOs), 
                                                      c("bio1", "bio2", "bio3", "bio7", "bio14", "bio15", 
                                                        "pr_spring", "pr_summer","pr_autumn", "pr_winter",
                                                        "bio1_3yrs", "bio2_3yrs", "bio3_3yrs", "bio7_3yrs", "bio14_3yrs", "bio15_3yrs",
                                                        "pr_spring_3yrs", "pr_summer_3yrs", "pr_autumn_3yrs", "pr_winter_3yrs", 
                                                        "sum_annual_crops", "secdf","pastr", "urban",
                                                        "sum_annual_crops_3yrs", "secdf_3yrs", "pastr_3yrs", "urban_3yrs")]
            }
            
            # make flocker data:
            fd_test <- make_flocker_data_dynamic(
              obs = y_array_dummy_test, 
              unit_covs = env_cov_test, 
              event_covs = det_cov_test, 
              quiet = TRUE
            )
        
            occ_posterior <- get_Z(out, history_condition = FALSE, new_data = fd_test, sample = FALSE) # for each site and year 4000 draws
            
            # predicted observations (y, 0 or 1) for test data:
            y_predictions <- predict_flocker(out, 
                                             history_condition = FALSE,
                                             new_data = fd_test) # for each site, route section and year 4000 draws
  
            # save results:
            res_list <- list("fitted" = fitted_initocc_col_ex,
                             "occ_posterior" = occ_posterior,
                             "y_preds" = y_predictions)
            
            save(res_list, file = file.path(dir, "data", paste0("test_preds_", spec, "_CV_fold", fold, ".RData")))
           
            end.time <- Sys.time()
            print(paste("finished", spec, "; time taken", round(end.time - start.time, 2)))
            
            sink(file = NULL)
          }
  
  sink(file = NULL)

}

stopCluster(cl)

rm(list=ls())
gc()
