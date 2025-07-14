# use DOMs to predict detections (y) for ISIMIP:
# obsclim + histsoc and 
# counterclim + histsoc
# 1901 - 2019
# for selected BBS routes

# obsclim + histsoc or counterclim + histsoc:
obsclim <- TRUE


# packages: --------------------------------------------------------------------

library(dplyr)
library(doParallel)
library(flocker)
library(cmdstanr)
set_cmdstan_path(path = NULL)
#set_cmdstan_path("C:/Users/schifferle1/Documents/cmdstan-2.34.1")
library(sf)
library(ggplot2)


# functions: -------------------------------------------------------------------

source("0_functions.R")

# register cores for parallel computation:
ncores <- 5
cl <- makeCluster(ncores, setup_timeout = 0.5)
registerDoParallel(cl)


# directories: -----------------------------------------------------------------

# logfiles:
log_dir <- file.path("logfiles", "ISIMIP", ifelse(obsclim, "obsclim_histsoc", "counterclim_histsoc"))
if(!dir.exists(log_dir)){dir.create(log_dir, recursive = TRUE)}

# directory with results, fitted models:

# res_dir <- file.path("//NAS-2-P-SN-01.ibb.uni-potsdam.de/daten$", "AG26", "Transfer", "Schifferle_BBS_occupancy_models_2023", "results")
res_dir <- file.path("/mnt", "ibb_share", "zurell_transfer", "Schifferle_BBS_occupancy_models_2023", "results")

# directory to store predictions:
preds_dir <- file.path(res_dir, "fm_preds_ISIMIP", ifelse(obsclim, "obsclim_histsoc", "counterclim_histsoc")) 
if(!dir.exists(preds_dir)){dir.create(preds_dir, recursive = TRUE)}


# load data: -------------------------------------------------------------------

# routes-years:
load(file = file.path("data", "BBS_for_occ_selection.RData")) # route_sel_dt; output of 1_3_match_BBS_to_env_data.R 

if(obsclim == TRUE){
  load(file = file.path("data", "route_sel_env_dt_ISIMIP_obsclim.RData")) 
} else{
  load(file = file.path("data", "route_sel_env_dt_ISIMIP_counterclim.RData"))
}  # route_sel_env_dt_ISIMIP; output of 1_3_match_routes_env_data_ISIMIP.R 

# route-year-species information (only surveyed)
load(file = file.path("data", "BBS_for_occ_spec_records.RData")) # bbs_dt_occ; output of 1_0_reformat_BBS_data.R

# selected variables:
load(file = file.path("data", "selected_variables.RData")) # selvar_final; output of 1_2_variable_selection.R
selvar_final

# selected routes spatial data (to buffer presences):
routes_sel_sf <- st_read(file.path("data", "route_selection_1995_2019_surv_beg_end_max_5y_miss_v2_spat_thin_100km_max_30_r_per_BCR_centroids.shp")) # output of 1_1_route_selection.R

# selected species:
load(file = file.path("data", "species_set_analysis.RData")) # output of 3_1_DOM_CV_evaluation_metrics.R
final_species


# scale env. data: ----

# load mean and sd with which training data were scaled:
load(file.path("data", "route_env_dt_scale_pars.RData")) # output of 2_1_DOM_flocker_fit_fm.R; env_scale_pars


# factual data:
route_sel_env_dt_ISIMIP_scaled <- route_sel_env_dt_ISIMIP
# iterate over columns:
for(c in colnames(route_sel_env_dt_ISIMIP[, 3:ncol(route_sel_env_dt_ISIMIP)])){
  print(c)
  route_sel_env_dt_ISIMIP_scaled[, c] <- as.numeric(scale(route_sel_env_dt_ISIMIP[, c], 
                                              center = as.numeric(env_scale_pars$center[c]),
                                              scale = as.numeric(env_scale_pars$scale[c])))
}


# # check:
# round(apply(route_sel_env_dt_ISIMIP_scaled, 2, mean), 1)
# round(apply(route_sel_env_dt_ISIMIP_scaled, 2, sd), 1)


# load model for species and predict: ------------------------------------------


# species:
foreach(spec = final_species[1:5],#length(final_species)], # xx
        .packages = c("dplyr", "collapse", "flocker", "cmdstanr", "brms", "sf", "terra"),
        .errorhandling = "pass", #"remove",
        .verbose = TRUE) %dopar% {
          
          # log progress:
          prog_log_file <- file(file.path(log_dir, paste0(spec, "_ypreds_progress.txt")), open = "wt") # write console output here
          sink(prog_log_file, type = "message")
          sink(prog_log_file, type = "output")
          
          print(spec)
          
          # check where to look for model output (did MCMC fitting work with less or only with more iterations?)
          if(file.exists(file.path(res_dir, "fm_buffer750km", "refit_2000_2000", paste0("out_", spec, "_fm_buffer750.RData")))){
            output_dir <- file.path(res_dir, "fm_buffer750km", "refit_2000_2000")
          } else {
            output_dir <- file.path(res_dir, "fm_buffer750km")
          }
          
          print(output_dir)
          
          # load fitted model:
          
          skip_to_next <- FALSE
          tryCatch(print(load(file.path(output_dir, paste0("out_", spec, "_fm_buffer750.RData")))), # output of 2_1_DOM_flocker_single_model.R
                   error = function(e) { skip_to_next <<- TRUE})
          if(skip_to_next) {
            print("model output not found")
            next 
          }
          
          # route locations within distance of 750 km of species records:
          rel_routes <- training_routes(species = spec, buffer_km = 750, output = "RTENOs")
          
          routes_in_buffer_sf <- routes_sel_sf %>% 
            filter(RTENO_BBS %in% rel_routes)
          #plot(st_geometry(routes_in_buffer_sf))

          # env. data for these cells:
          env_df_routes <- route_sel_env_dt_ISIMIP_scaled %>% 
            filter(RTENO_BBS %in% routes_in_buffer_sf$RTENO_BBS)

          # # test plot:
          # env_df_routes %>%
          #   ggplot() +
          #   geom_line(aes(x = Year, y = sum_annual_crops, group = RTENO_BBS, colour = as.factor(RTENO_BBS))) +
          #   theme_bw() +
          #   theme(legend.position = "none")
          

          # predict:
          # for predict_flocker new data must be formatted as a flocker data object:
          
          years <- seq(min(env_df_routes$Year), max(env_df_routes$Year))
          nyears <- length(years)
          nsites <- length(unique(env_df_routes$RTENO_BBS))
          nsurveys <- 5
          
          # reformat obs. as array sites x surveys x years:
          y_array_dummy <- array(NA, dim = c(nsites, nsurveys, nyears))
          for (t in 1:nyears){
            y_array_dummy[1:nsites, 1:nsurveys, t] <- as.matrix(0)
          }
          
          # covariate for detection probability:
          det_cov <- vector("list", length = 1)
          names(det_cov) <- "route_section"
          det_cov$route_section <- array(NA, dim = c(nsites, nsurveys, nyears))
          det_cov$route_section[ , , 1:nyears] <- matrix(rep(c("Sect1", "Sect2", "Sect3", "Sect4", "Sect5"), nsites), nsites, byrow = TRUE)
          
          # reformat environmental covariates:
          env_cov <- vector("list", length = nyears)
          for (t in 1:nyears){
            env_cov[[t]] <- env_df_routes[which(env_df_routes$Year == years[t]), c(selvar_final, paste0(selvar_final, "_3yrs"))]
          }
          
          print("start creating flocker data")
          print(Sys.time())
          
          # make flocker data:
          fd_new <- make_flocker_data_dynamic(
            obs = y_array_dummy, 
            unit_covs = env_cov, 
            event_covs = det_cov, 
            quiet = TRUE
          )
          
          print("created flocker data, start predictions")
          print(Sys.time())
          
          y_predictions <- predict_flocker(flocker_fit = out, 
                                           history_condition = FALSE,
                                           new_data = fd_new,
                                           draw_ids = seq(1, 4000, 4)) # 1000
          
          # summarise detections across route sections since we look at route level
          y_preds_route <- apply(y_predictions, MAR = c(1,3,4), FUN = max)
          
          # save predictions:
          
          # use counterfactual data also for predictions of initial occupancy:
          save(y_preds_route, file = file.path(preds_dir,  paste0(spec, "_y_preds_", ifelse(obsclim, "obsclim_histsoc", "counterclim_histsoc"), ".RData")))
          
          sink(type="message")
          sink(type="output")
          
          rm(y_predictions, y_preds_route, fd_new)
        }

stopCluster(cl)
rm(list=ls())
gc()