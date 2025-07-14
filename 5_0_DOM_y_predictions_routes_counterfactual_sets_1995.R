# use DOMs to predict detections (y) for counterfactual data on BBS routes
# to compare observations (thus y and not Z needed) to counterfactual expectations based on the DOMs
# use climate data from ISIMIP, but detrended with ATTRICI from 1995 onwards
# use land use data from ISIMIP, fix at 1995 level
# counterfactual data for colonisation and extinction probability
# factual data for initial occupancy
# predict once for all climate covariates counterfactual, once for all land use covariates counterfactual
# once for all counterfactual for colonisation and extinction probability


# packages: --------------------------------------------------------------------

library(doParallel)
library(dplyr)
library(flocker)
library(cmdstanr)
set_cmdstan_path(path = NULL)
#set_cmdstan_path("C:/Users/schifferle1/Documents/cmdstan-2.34.1")
library(sf)
library(ggplot2)
library(terra)


# register cores for parallel computation:
ncores <- 17
cl <- makeCluster(ncores, setup_timeout = 0.5)
registerDoParallel(cl)

# functions: -------------------------------------------------------------------

source("0_functions.R")

# directories: -----------------------------------------------------------------

# logfiles:
log_dir <- file.path("logfiles", "attribution", "y_preds_routes_cf_1995")
if(!dir.exists(log_dir)){dir.create(log_dir, recursive = TRUE)}

# directory with results, fitted models:
# res_dir <- file.path("T:", "Schifferle_BBS_occupancy_models_2023", "results")
res_dir <- file.path("/mnt", "ibb_share", "zurell_transfer", "Schifferle_BBS_occupancy_models_2023", "results")

# directory to store predictions:
preds_dir <- file.path(res_dir, "attribution", "fm_y_preds_routes_cf_1995_all")
if(!dir.exists(preds_dir)){dir.create(preds_dir, recursive = TRUE)}


# load data: -------------------------------------------------------------------

# routes-years:
load(file = file.path("data", "BBS_for_occ_selection.RData")) # route_sel_dt; output of 1_3_match_BBS_to_env_data.R 

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

# US cells:
US_cells_sf <- read_sf(file.path("data", "cell_centroids_US_ESRI102003.shp")) # output of 1_3_nev_data_df_contUS.R
US_cells_sf # 4149 cells

# factual environmental data for each US grid cell:
load(file = file.path("data", "US_grid_env_data.RData")) # output of 1_3_env_data_df_contUS.R; clim_lu_cells_sf
fact_env_df <- sf::st_drop_geometry(clim_lu_cells_sf)

# counterfactual environmental data for each US grid cell:
load(file = file.path("data", "US_grid_env_data_cf_1995_all.RData")) # output of 1_3_env_data_df_contUS.R; clim_lu_cells_cf_1995
clim_lu_cells_cf_1995


# scale env. data: -------------------------------------------------------------

# load mean and sd with which training data were scaled:
load(file.path("data", "route_env_dt_scale_pars.RData")) # output of 2_1_DOM_flocker_fit_fm.R; env_scale_pars

# factual data:
fact_env_df_scaled <- fact_env_df
# iterate over columns:
for(c in colnames(fact_env_df[, 3:ncol(fact_env_df)])){
  print(c)
  fact_env_df_scaled[, c] <- as.numeric(scale(fact_env_df[, c], 
                                              center = as.numeric(env_scale_pars$center[c]),
                                              scale = as.numeric(env_scale_pars$scale[c])))
}

# counterfactual data:
counterfact_env_df_scaled <- clim_lu_cells_cf_1995
# iterate over columns:
for(c in colnames(clim_lu_cells_cf_1995[, 3:ncol(clim_lu_cells_cf_1995)])){
  print(c)
  counterfact_env_df_scaled[, c] <- as.numeric(scale(clim_lu_cells_cf_1995[, c], 
                                                     center = as.numeric(env_scale_pars$center[c]),
                                                     scale = as.numeric(env_scale_pars$scale[c])))
}

# # check:
# round(apply(fact_env_df_scaled, 2, mean), 1)
# round(apply(fact_env_df_scaled, 2, sd), 1)
# round(apply(counterfact_env_df_scaled, 2, mean), 1)
# round(apply(counterfact_env_df_scaled, 2, sd), 1)


# load model for species and predict: ------------------------------------------


# species:
foreach(#spec = c("Cassin's Sparrow", "Brown Creeper", "Bobolink", "Brewer's Sparrow", "Olive-sided Flycatcher"),#final_species[1:5],#length(final_species)], # xx
        spec = subset(final_species, !final_species %in% c("Cassin's Sparrow", "Brown Creeper", "Bobolink", "Brewer's Sparrow", "Olive-sided Flycatcher")),
        .packages = c("dplyr", "collapse", "flocker", "cmdstanr", "brms", "sf", "terra"),
        .errorhandling = "pass", #"remove",
        .verbose = TRUE) %dopar% {
          
          # log progress:
          prog_log_file <- file(file.path(log_dir, paste0(spec, "_ypreds_cf_progress.txt")), open = "wt") # write console output here
          sink(prog_log_file, type = "message")
          sink(prog_log_file, type = "output")
          
          print(Sys.time())
          
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
          
          # find US cell IDs where routes are located:
          
          # route locations within distance of 750 km of species records:
          rel_routes <- training_routes(species = spec, buffer_km = 750, output = "RTENOs")
          
          routes_in_buffer_sf <- routes_sel_sf %>% 
            filter(RTENO_BBS %in% rel_routes)
          #plot(st_geometry(routes_in_buffer_sf))
          
          # convert US cells to raster to then find cell ID for each route:
          US_cells_rast <- US_cells_sf %>%
            mutate(x = st_coordinates(.)[,1],
                   y = st_coordinates(.)[,2]) %>%
            st_drop_geometry() %>%
            select(x, y, everything()) %>%
            rast(., type='xyz', crs=crs(US_cells_sf))
          
          # extract cellID for each route location:
          obs_cells <- extract(x = US_cells_rast, y = routes_in_buffer_sf, cells = FALSE, ID = FALSE)
          #plot(US_cells_rast, main = "cellID")
          #plot(st_geometry(routes_in_buffer_sf), add = TRUE)
          
          # env. data for these cells:
          fact_env_df_routes <- fact_env_df_scaled %>% 
            filter(cellID %in% obs_cells$cellID)
          
          counterfact_env_df_routes <- counterfact_env_df_scaled %>% 
            filter(cellID %in% obs_cells$cellID)
          
          # # test plot:
          # fact_env_df_routes %>%
          #   mutate(scen = "factual") %>%
          #   rbind(counterfact_env_df_routes %>%  mutate(scen = "cf")) %>%
          #   ggplot() +
          #   facet_wrap(~scen) +
          #   geom_line(aes(x = year, y = sum_annual_crops, group = cellID, colour = as.factor(cellID))) +
          #   theme_bw() +
          #   theme(legend.position = "none")
          
          
          # iterate over counterfactual scenarios:
          
          # either counterfactual values for all climate, all land use or all variables at a time:

          #scenarios <- c("counterclim", "1995soc", "counterclim_1995soc")
          #scenarios <- "1995soc"
          scenarios <- c("counterclim", "counterclim_1995soc")
          
          for(v in scenarios){
            
            print(v)
            
            # check whether scenario ran already:
            # if(file.exists(file.path(preds_dir,  paste0(spec, "_y_preds_cf_", v, ".RData")))){
            #   print(paste(spec, v, "ran already."))
            #   next
            # }
            
            # gather new data:
            
            if (v == "counterclim"){
              
              new_data_lu <- fact_env_df_routes %>%
                select(cellID, year, !matches("(bio)|(pr_mean)"))
              
              new_data_clim <- counterfact_env_df_routes %>%
                select(cellID, year, matches("(bio)|(pr_mean)"))

              # merge both:
              new_data <- new_data_clim %>% 
                left_join(new_data_lu)
              
            } else if (v == "1995soc") {
              
              new_data_lu <- counterfact_env_df_routes %>%
                select(cellID, year, !matches("(bio)|(pr_mean)"))
              
              new_data_clim <- fact_env_df_routes %>%
                select(cellID, year, matches("(bio)|(pr_mean)"))
              
              # merge both:
              new_data <- new_data_clim %>% 
                left_join(new_data_lu)
              
            } else {
              
              new_data <- counterfact_env_df_routes
            
            }
            
            
            # predict:
            # for predict_flocker new data must be formatted as a flocker data object:
            
            years <- seq(min(new_data$year), max(new_data$year))
            nyears <- length(years)
            nsites <- length(unique(new_data$cellID))
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
              env_cov[[t]] <- new_data[which(new_data$year == years[t]), c(selvar_final, paste0(selvar_final, "_3yrs"))]
            }
            
            # make flocker data:
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
                                             draw_ids = seq(1, 4000, 4)) # 1000
            
            print(Sys.time())
            
            # summarise detections across route sections since we look at route level
            y_preds_route_cf <- apply(y_predictions, MAR = c(1,3,4), FUN = max)
            
            # save predictions:
            
            # use counterfactual data also for predictions of initial occupancy:
            save(y_preds_route_cf, file = file.path(preds_dir,  paste0(spec, "_y_preds_cf_", v, ".RData")))
            
          }
          
          sink(type="message")
          sink(type="output")
        }

stopCluster(cl)
rm(list=ls())
gc()
