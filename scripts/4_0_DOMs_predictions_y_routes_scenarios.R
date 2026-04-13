# use fitted dynamic occupancy models to simulate occupancy dynamics of species
# for which models show acceptable performance in space and time for the scenarios:
# - no climate change since 1995, but land use change
# - no land use change since 1995, but climate change
# - no climate and no land use change since 1995

# (simulate observations (y), not Z, to compare to real observations)

# counterfactual environmental data to simulate colonisation and extinction probability,
# factual data for initial occupancy


# packages: --------------------------------------------------------------------

library(dplyr)
library(sf)
library(doParallel)
library(flocker)
library(cmdstanr)


# directories: -----------------------------------------------------------------

set_cmdstan_path(path = NULL) # for HPC; local: set_cmdstan_path("C:/Users/schifferle1/Documents/cmdstan-2.34.1")

# project directory:
#dir <- file.path("//NAS-2-P-SN-01.ibb.uni-potsdam.de", "daten$", "AG26", "Transfer", "Schifferle_BBS_occupancy_models_2023")
dir <- file.path("/mnt", "ibb_share", "zurell_transfer", "Schifferle_BBS_occupancy_models_2023")

# logfiles:
log_dir <- file.path("logfiles", "attribution", "y_preds_routes_cf_1995")
if(!dir.exists(log_dir)){dir.create(log_dir, recursive = TRUE)}

# directory with fitted models:
res_dir <- file.path(dir, "results", "fm_buffer750km")

# directory to store predictions:
preds_dir <- file.path(dir, "results", "attribution", "fm_y_preds_routes_cf_1995_all")
if(!dir.exists(preds_dir)){dir.create(preds_dir, recursive = TRUE)}


# functions: -------------------------------------------------------------------

source(file.path("scripts", "0_functions.R"))


# load data: -------------------------------------------------------------------

# route-year-species information (only surveyed)
load(file = file.path("data", "BBS_for_occ_spec_records.RData")) # bbs_dt_occ; output of 1_0_dataprep_BBS_bird_data.R

# selected routes spatial data (to buffer presences):
routes_sel_sf <- st_read(file.path("data", "route_selection_1995_2019_surv_beg_end_max_5y_miss_v2_spat_thin_100km_max_30_r_per_BCR_centroids.shp")) # output of 1_1_dataprep_BBS_route_selection.R

# species:
load(file = file.path("data", "species_set_analysis.RData")) # output of 3_1_DOM_CV_evaluation_metrics.R
final_species

# selected variables:
load(file = file.path("data", "selected_variables.RData")) # selvar_final; output of 1_2a_dataprep_env_variable_selection.R

# routes-years:
load(file = file.path("data", "BBS_for_occ_selection.RData")) # route_sel_dt; output of 1_3_dataprep_match_BBS_routes_env_data.R

# selected routes and focal years matched to environmental data:
load(file = file.path("data", "route_year_env_data.RData")) # route_sel_env_dt_final; output 1_3_dataprep_match_BBS_routes_env_data.R

# environmental data:

# factual data:
load(file.path("data", "route_year_env_data.RData")) # route_sel_env_dt_final; output of 1_3_dataprep_match_BBS_routes_env_data.R
f_env_data <- route_sel_env_dt_final %>% 
  select(-c(Latitude, Longitude, BCR, ObsN, doy, Surveyed))

# counterfactual data:
load(file.path("data", "route_year_env_data_cf.RData")) # route_sel_env_dt_final; output of 1_3_dataprep_match_BBS_routes_env_data.R
cf_env_data <- route_sel_env_dt_final %>% 
  select(-c(Latitude, Longitude, BCR, ObsN, doy, Surveyed))


# scale environmental data: ----------------------------------------------------

# load mean and sd with which training data were scaled:
load(file.path("data", "route_env_dt_scale_pars.RData")) # output of 2_1_fit_DOMs_full_model.R; env_scale_pars

# factual data:
f_env_data_scaled <- f_env_data
# iterate over columns:
for(c in colnames(f_env_data[, 3:ncol(f_env_data)])){
  print(c)
  f_env_data_scaled[, c] <- as.numeric(scale(f_env_data[, c],
                                              center = as.numeric(env_scale_pars$center[c]),
                                              scale = as.numeric(env_scale_pars$scale[c])))
}

# counterfactual data:
cf_env_data_scaled <- cf_env_data
# iterate over columns:
for(c in colnames(cf_env_data[, 3:ncol(cf_env_data)])){
  print(c)
  cf_env_data_scaled[, c] <- as.numeric(scale(cf_env_data[, c],
                                              center = as.numeric(env_scale_pars$center[c]),
                                              scale = as.numeric(env_scale_pars$scale[c])))
}


# species for which DOMs show acceptable predictive performance: ---------------

# okay in time:
load(file.path(dir, "results", "temp_val_buffer_750_10yrs", "temp_eval", "10_years", "spec_set_temp_val_ok.RData")) # output of 3_2_eval_DOMs_temp.R
spec_temp_okay <- specs_thresh

# okay in space:
CV_eval_summary <- read.csv(file = file.path(dir, "results", "CV_buffer750km", "CV_eval", "CV_eval_summary.csv")) # output of 3_1_eval_DOMs_CV.R
spec_spat_okay <- CV_eval_summary %>% 
  filter(y_spat_auc_mean >= 0.7) %>% 
  pull(species)

# okay in both:
spec_attr <- intersect(spec_temp_okay, spec_spat_okay) # 80
#save(spec_attr, file = file.path("data", "species_DOM_val_okay.RData"))


# load model for species and predict: ------------------------------------------

# register cores for parallel computation:
ncores <- 17
cl <- makeCluster(ncores, setup_timeout = 0.5)
registerDoParallel(cl)


# iterate over species:

foreach(spec = spec_attr,
        .packages = c("dplyr", "collapse", "flocker", "cmdstanr", "brms", "sf", "terra"),
        .errorhandling = "pass", #"remove",
        .verbose = TRUE) %dopar% {
          
          # log progress:
          prog_log_file <- file(file.path(log_dir, paste0(spec, "_ypreds_cf_progress.txt")), open = "wt") # write console output here
          sink(prog_log_file, type = "message")
          sink(prog_log_file, type = "output")
          
          print(Sys.time())
          
          print(spec)
          
          # check where to look for model output:
          if(file.exists(file.path(res_dir, "refit_2000_2000", paste0("out_", spec, "_fm_buffer750.RData")))){
            output_dir <- file.path(res_dir, "refit_2000_2000")
          } else {
            output_dir <- res_dir
          }
          
          print(output_dir)
          
          # load fitted model:
          skip_to_next <- FALSE
          tryCatch(print(load(file.path(output_dir, paste0("out_", spec, "_fm_buffer750.RData")))), # output 2_1_fit_DOMs_full_model.R
                   error = function(e) {skip_to_next <<- TRUE})
          if(skip_to_next) {
            print("model output not found")
            next 
          }
          
          # route locations within distance of 750 km of species records:
          rel_routes <- training_routes(species = spec, buffer_km = 750, output = "RTENOs")

          # factual environmental data for relevant routes:
          f_env_data_scaled_spec <- f_env_data_scaled %>%
            filter(RTENO %in% rel_routes)

          # counterfactual environmental data for relevant routes:
          cf_env_data_scaled_spec <- cf_env_data_scaled %>%
            filter(RTENO %in% rel_routes)
          
          # iterate over counterfactual scenarios:
          
          # counterclim: counterfactual data for all climate variables, factual for land use variables
          # 1995soc: factual data for climate variables, land use variables kept constant at 1995
          # counterclim_1995soc: counterfactual data for all climate variables, land use variables kept constant at 1995

          scenarios <- c("counterclim", "1995soc", "counterclim_1995soc")

          for(v in scenarios){
            
            print(v)
            
            # gather environmental data:
            
            if (v == "counterclim"){
              
              new_data_lu <- f_env_data_scaled_spec %>%
                select(RTENO, Year, !matches("(bio)|(pr_mean)"))
              
              new_data_clim <- cf_env_data_scaled_spec %>%
                select(RTENO, Year, matches("(bio)|(pr_mean)"))

              # merge both:
              new_data <- new_data_clim %>% 
                left_join(new_data_lu)
              
            } else if (v == "1995soc") {
              
              new_data_lu <- cf_env_data_scaled_spec %>%
                select(RTENO, Year, !matches("(bio)|(pr_mean)"))
              
              new_data_clim <- f_env_data_scaled_spec %>%
                select(RTENO, Year, matches("(bio)|(pr_mean)"))
              
              # merge both:
              new_data <- new_data_clim %>% 
                left_join(new_data_lu)
              
            } else {
              
              new_data <- cf_env_data_scaled_spec
            
            }
            
            # predict to counterfactual data:

            years <- seq(min(new_data$Year), max(new_data$Year))
            nyears <- length(years)
            nsites <- length(unique(new_data$RTENO))
            nsurveys <- 5
            
            # dummy observations as array sites x surveys x years:
            
            # https://cran.r-project.org/web/packages/flocker/vignettes/flocker_tutorial.html:
            # "Note that if predictions are desired at sites without observations, it is acceptable 
            # to pass an array of dummy observations (e.g. all zeros) to make_flocker_data() and 
            # then to set history_condition = FALSE in the call to get_Z()."
            
            y_array_dummy <- array(NA, dim = c(nsites, nsurveys, nyears))
            for (t in 1:nyears){
              y_array_dummy[1:nsites, 1:nsurveys, t] <- as.matrix(0)
            }
            
            # reformat environmental covariates:
            env_cov <- vector("list", length = nyears)
            for (t in 1:nyears){
              env_cov[[t]] <- new_data[which(new_data$Year == years[t]), c(selvar_final, paste0(selvar_final, "_3yrs"))]
            }
            
            # covariate for detection probability:
            # ("please add a dummy event covariate. You do not need to use this covariate in your model formula"):
            
            det_cov <- vector("list", length = 1)
            names(det_cov) <- "route_section"
            det_cov$route_section <- array(NA, dim = c(nsites, nsurveys, nyears))
            det_cov$route_section[ , , 1:nyears] <- matrix(rep(c("Sect1", "Sect2", "Sect3", "Sect4", "Sect5"), nsites), nsites, byrow = TRUE)
            
            # make flocker data:
            # (using the `new_data` argument within get_Z() for a multiseason model requires passing a `flocker_data` object)
            
            fd_new <- make_flocker_data_dynamic(
              obs = y_array_dummy, 
              unit_covs = env_cov, 
              event_covs = det_cov, 
              quiet = TRUE
            )
            
            print(Sys.time())
            
            y_predictions <- predict_flocker(flocker_fit = out, 
                                             history_condition = FALSE,
                                             new_data = fd_new,
                                             draw_ids = seq(1, 4000, 4)) # 1000 draws 
            
            print(Sys.time())
            
            # summarise detections across route sections since we look at route level
            y_preds_route_cf <- apply(y_predictions, MAR = c(1,3,4), FUN = max)
            
            # save predictions:
            save(y_preds_route_cf, file = file.path(preds_dir,  paste0(spec, "_y_preds_cf_", v, ".RData")))
            
          }
          
          sink(type="message")
          sink(type="output")
        }

stopCluster(cl)
rm(list=ls())
gc()