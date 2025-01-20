# 5-fold spatial cross validation to assess predictive ability of dynamic occupancy models:

# folds created with package blockCV: 2_1_1_assign_spatial_blocks_CV_folds.R

# 1.) for species load fold assignment
# 2.) refit model 5 times (5-fold cross validation)
# 3.) predict to test folds: occupancy probability per year and predicted y (0/1) (this is comparable to observed data)

# first for 750 km buffer

# packages: ----

library(dplyr)
library(doParallel)
library(flocker)
#install.packages("cmdstanr", repos = c("https://mc-stan.org/r-packages/", getOption("repos")))
library(cmdstanr)
set_cmdstan_path(path = NULL)
#set_cmdstan_path("C:/Users/schifferle1/Documents/cmdstan-2.34.1") # xx
library(sf)

# register cores for parallel computation:
ncores <- 20 # 5 fold * 4 chains * species
cl <- makeCluster(ncores, setup_timeout = 0.5)
registerDoParallel(cl)


# directories: ----

print(tempdir())
dir <- file.path("/import", "ecoc9z", "data-zurell", "schifferle", "BBS_occupancy_models_2023")
#dir <- getwd()

# directory for logfiles:
log_dir <- file.path("logfiles", "CV_buffer750km")
# directory for results:
res_dir <- file.path(dir, "results", "CV_buffer750km")


# functions: ----

source("0_functions.R")


# load data: ----

# selected routes and focal years matched to environmental data:

# merged route, year, environment data:
load(file = file.path("data", "route_year_env_data.RData")) # route_sel_env_dt_final; output of 1_3_match_BBS_to_env_data.R

# selected variables:
load(file = file.path("data", "selected_variables.RData")) # selvar_final; output of 1_2_variable_selection.R

# routes-years:
load(file = file.path("data", "BBS_for_occ_selection.RData")) # route_sel_dt; output of 1_3_match_BBS_to_env_data.R 

# selected routes spatial data (to buffer presences):
routes_sel_sf <- st_read(file.path("data", "route_selection_1995_2019_surv_beg_end_max_5y_miss_v2_spat_thin_100km_max_30_r_per_BCR_centroids.shp")) # output of 1_1_route_selection.R

# route-year-species information (only surveyed)
load(file = file.path("data", "BBS_for_occ_spec_records.RData")) # bbs_dt_occ; output of 1_0_reformat_BBS_data.R

# selected species:
# sorted by ecoregion:
load(file = file.path("data", "final_species_selection_eco_sorted.RData")) # final_species_eco_sorted; output of 1_2_species_selection.R

# so far discarded species:
load(file = file.path("M:", "Documents", "DEBTs", "analysis", "Schifferle_BBS_occupancy_models_2023",
                      "results", "fm_buffer750km", "refit_2000_2000", "check_output", "specs_discard.RData")) # specs_discard_fm; output of 2_2_DOM_check_flocker_fit.R

species_set <- final_species_eco_sorted[-which(final_species_eco_sorted %in% specs_discard_fm)]

# assemble overall data: ----

# scale covariates:

route_sel_env_dt_scaled <- route_sel_env_dt_final %>% 
  select(-c(Latitude, Longitude, BCR, ObsN, doy)) %>% 
  mutate(across(!c(RTENO, Year, Surveyed), ~ as.numeric(scale(., center=mean(.), scale = sd(.)))))

rm(route_sel_env_dt_final)


nyears <- length(unique(route_sel_env_dt_scaled$Year)) # 25
nsurveys <- 5
#nsites <- length(unique(route_sel_env_dt_scaled$RTENO)) # 476


# make blocks of species:
# 3 species per block:
spec_blocks_list <- split(species_set, c(rep(1:(floor(length(species_set)/3)), each = 3), 
                                                      rep(ceiling(length(species_set)/3), length(species_set) %% 3))) # xx change once it runs
names(spec_blocks_list) <- NULL

prog_log_file <- file(file.path(log_dir, "CV_buffer_750_progress.txt"), open = "wt") # write console output here
sink(prog_log_file, type = "message")
sink(prog_log_file, type = "output")

#for(i in 1:3){
  
#  print(paste("block", i, "of", length(spec_blocks_list)))
  
  foreach(spec = spec_blocks_list[[1]][2], #spec_blocks_list[[i]],
          .packages = c("dplyr", "collapse", "flocker", "cmdstanr", "brms", "sf"), # xx
          .errorhandling = "pass", #"remove",
          .verbose = TRUE) %:% 
            
            foreach(fold = 1:5,
                    .packages = c("dplyr", "collapse", "flocker", "cmdstanr", "brms"), # xx
                    .errorhandling = "pass", #"remove",
                    .verbose = TRUE) %dopar% {
                      
                      spec_log_file <- file(file.path(log_dir, paste0(spec, "_CV_fitting.txt")), open = "wt") # write console output here
                      sink(spec_log_file, type = "message")
                      sink(spec_log_file, type = "output")
                      
                      # check whether species has run already:
                      CV_run <- length(list.files(path = file.path(res_dir), pattern = paste0(spec, "_CV_fold"))) == 10
                      if(CV_run) {
                        print(paste(spec, "ran already."))
                        next
                      }

                      print(spec)
                      
                      # relevant routes, within distance of 750 km of species records:
                      rel_routes <- training_routes(species = spec, buffer_km = 750, output = "RTENOs")
                      
                      print(length(rel_routes))
                      
                      # load fold assignment:
                      load(file.path("data", "CV_route_block_allocation", "block_size_500km", paste0(spec, ".RData")))
                      
                      # assemble data:
                      occ_dt_spec <- BBS_pres_abs_spec(species = spec)

                      log_file_spec_fold <- file(file.path(log_dir, paste0(spec, "_CV_fitting_fold", fold, ".txt")), open = "wt") # write console output here
                      sink(log_file_spec_fold, type = "message")
                      sink(log_file_spec_fold, type = "output")
                      
                      print(spec)
                      print(paste("fold", fold))
                      
                      # data in current fold:
                      training_RTENOs <- rel_routes[sb_US$folds_list[[fold]][[1]]] # training data fold 
                      test_RTENOs <- rel_routes[sb_US$folds_list[[fold]][[2]]] # test data fold 
            
                      print(paste("num rel_routes", length(rel_routes)))
                      print(c(length(training_RTENOs), length(test_RTENOs)))
                      
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
                      
                      print(paste("dim y_array:", paste(dim(y_array), collapse = ", ")))
                      
                      # reformat environmental covariates:
                      
                      env_cov <- vector("list", length = nyears)
                      for (t in 1:nyears){
                        env_cov[[t]] <- route_sel_env_dt_scaled[which(route_sel_env_dt_scaled$Year == years[t] & route_sel_env_dt_scaled$RTENO %in% occ_dt_spec_subset$RTENO), 
                                                                c(selvar_final, paste0(selvar_final, "_3yrs"))]
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
                        warmup = 1000,
                        iter = 1000 + 1000
                      )
                      
                      print(out)
                      
                      # save fitted model for each fold:
                      save(out, file = file.path(res_dir, paste0("out_", spec, "_CV_fold", fold, ".RData")))
                      
                      end.time <- Sys.time()
                      time.taken <- round(end.time - start.time, 2)
                      print(time.taken)
                      
                      
                      # calculate further results / predictions to test routes:
                      
                      print("start predictions")
                      
                      start.time <- Sys.time()
                      
                      # data for test folds:
                      test_env_data <- route_sel_env_dt_scaled %>%
                        filter(RTENO %in% test_RTENOs)
                      
                      print(paste("test_en_data N routes:", length(unique(test_env_data$RTENO))))
                      
              
                      # prepare test data:
          
                      
                      # get_Z(): using the `new_data` argument for a multiseason model requires passing a `flocker_data` object
                      
                      # Note that if predictions are desired at sites without observations, it is acceptable 
                      # to pass an array of dummy observations (e.g. all zeros) to make_flocker_data() and 
                      # then to set history_condition = FALSE in the call to get_Z().
                      # (https://cran.r-project.org/web/packages/flocker/vignettes/flocker_tutorial.html)
                      # "please add a dummy event covariate. You do not need to use this covariate in your model formula":
                      
                      #years_test <- seq(min(test_env_data$Year), max(test_env_data$Year))
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
                                                                c(selvar_final, paste0(selvar_final, "_3yrs"))]
                      }
                      
                      
                      # make flocker data:
                      fd_test <- make_flocker_data_dynamic(
                        obs = y_array_dummy_test, 
                        unit_covs = env_cov_test, 
                        event_covs = det_cov_test, 
                        quiet = TRUE
                      )
                      
                      
                      # fitted values for test routes (not necessary):
                      
                      fitted_initocc_col_ex_det <- fitted_flocker(out,
                                                                  components = c("occ", "col", "ex", "det"),
                                                                  new_data = fd_test)
                      
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
                      
          
                      # occupancy probability:
                      
                      occ_posterior <- get_Z(out, 
                                             history_condition = FALSE, 
                                             new_data = fd_test, 
                                             sample = FALSE) # for each site and year 4000 draws
                      
                      print(paste("dim occ. posterior", dim(occ_posterior)))
                      
                      # mean and median for each site and year (as integer to reduce storage space):
                      occ_posterior_median <- apply(occ_posterior, MARGIN = c(1,2), FUN = function(x) as.integer(round(median(x) * 100)))
                      occ_posterior_mean <- apply(occ_posterior, MARGIN = c(1,2), FUN = function(x) as.integer(round(mean(x) * 100)))
                      
          
                      # predicted observations (y: 0 or 1) for test data:
                      y_predictions <- predict_flocker(out, 
                                                       history_condition = FALSE,
                                                       new_data = fd_test) # for each site, route section and year 4000 draws
                      
                      y_predictions_sections_sum <- apply(y_predictions, MAR = c(1,3,4), FUN = max) # summarise detections across route sections since we look at route level
                      # if detection prob. is 0.4, probability that species is missed on whole route is 0.6*0.6*0.6*0.6*0.6 = 8% prob.
                      # -> we expect that prob. that species is detected on a route is 8% smaller than prob. that species is present on a route?
          
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
                      
                      save(res_list, file = file.path(res_dir, paste0("test_preds_", spec, "_CV_fold", fold, ".RData")))
          
                      end.time <- Sys.time()
                      print(paste("finished", spec, "; time taken", round(end.time - start.time, 2)))
                      
                      sink(type="message")
                      sink(type="output")
                      sink(type="message")
                      sink(type="output")
                    }

#}
sink(type="message")
sink(type="output")

stopCluster(cl)

rm(list=ls())
gc()
