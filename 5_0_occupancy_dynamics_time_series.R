# calculate time series of sum of routes with observed or predicted species presence:
# median per year, CI per year, 100 draws of posterior per year

# packages: --------------------------------------------------------------------

library(dplyr)
library(sf)
library(ggplot2)
library(doParallel)
library(bayestestR)

# register cores for parallel computation:
ncores <- 2
cl <- makeCluster(ncores, setup_timeout = 0.5)
registerDoParallel(cl)

# functions: -------------------------------------------------------------------

source("0_functions.R")


# directories: -----------------------------------------------------------------

# main directory:
# main_dir <- file.path("//NAS-2-P-SN-01.ibb.uni-potsdam.de", "daten$", "AG26", "Transfer", 
#                       "Schifferle_BBS_occupancy_models_2023")

main_dir <- file.path("/mnt", "ibb_share", "zurell_transfer", "Schifferle_BBS_occupancy_models_2023")


# save observations time series:
obs_dir <- file.path(main_dir, "data", "observed_time_series_1995_2019")

# save predicted time series for factual data:
fact_dir <- file.path(main_dir, "results", "fm_buffer750km", "fact_pred_time_series_1995_2019")

# save predicted time series for counterfactual data:
cfact_dir <- file.path(main_dir, "results", "attribution", "cfact_pred_time_series_1995_2019")

if(!dir.exists(obs_dir)){dir.create(obs_dir)}
if(!dir.exists(fact_dir)){dir.create(fact_dir)}
if(!dir.exists(cfact_dir)){dir.create(cfact_dir)}


# load data: -------------------------------------------------------------------

# selected species:
load(file = file.path("data", "species_set_analysis.RData")) # output of 3_1_DOM_CV_evaluation_metrics.R
final_species

# routes-years:
load(file = file.path("data", "BBS_for_occ_selection.RData")) # route_sel_dt; output of 1_3_match_BBS_to_env_data.R 

# selected routes spatial data (to buffer presences):
routes_sel_sf <- st_read(file.path("data", "route_selection_1995_2019_surv_beg_end_max_5y_miss_v2_spat_thin_100km_max_30_r_per_BCR_centroids.shp")) # output of 1_1_route_selection.R

# route-year-species information (only surveyed)
load(file = file.path("data", "BBS_for_occ_spec_records.RData")) # bbs_dt_occ; output of 1_0_reformat_BBS_data.R


# calculate time series: -------------------------------------------------------

# iterate over species:

foreach(s = 1:length(final_species),
  .packages = c("dplyr", "collapse", "sf", "bayestestR"),
  .errorhandling = "pass", #"remove",
  .verbose = TRUE) %dopar% {
  
  spec <- final_species[s]
  print(paste(s, spec))
  
  # observations: -----------------------------
  
  print("observations")
  
  rel_routes <- training_routes(species = spec, buffer_km = 750, output = "RTENOs")
  occ_dt_spec <- BBS_pres_abs_spec(species = spec) %>%
    filter(RTENO %in% rel_routes)
  
  # route-level presence:
  # sum all routes for each year (temporal trend)
  ts_obs <- occ_dt_spec %>%
    rename("year" = Year) %>% 
    group_by(year) %>%
    summarise(Npres = sum(presence, na.rm = TRUE))
  
  save(ts_obs, file = file.path(obs_dir, paste0(spec, "_obs_ts_sum_occ_routes.RData")))
  
  
  # predictions for factual data incl. CI: ---------
  
  print("factual predictions")
  
  # check where to look for model output (did MCMC fitting work with less or only with more iterations?)
  if(file.exists(file.path(main_dir, "results", "fm_buffer750km", "refit_2000_2000", paste0("out_", spec, "_fm_buffer_750.RData")))){
    output_dir <- file.path(main_dir, "results", "fm_buffer750km", "refit_2000_2000")
  } else {
    output_dir <- file.path(main_dir, "results", "fm_buffer750km")
  }
  
  #load(file.path(output_dir, paste0("postproc_", spec, "_fm_buffer750.RData")))
  # sum across route sections:
  #preds_routes <- apply(res_list$y_preds, MAR = c(1,3,4), FUN = max) # as soon as one predicted presence -> 1
  # save this:
  #save(preds_routes, file = file.path(res_dir, "fm_buffer750km", "y_preds_route_level_section_sum", paste0(spec, "_y_preds_route_level_section_sum.RData")))
  
  load(file.path(main_dir, "results", "fm_buffer750km", "y_preds_route_level_section_sum", 
                 paste0(spec, "_y_preds_route_level_section_sum.RData")))
  preds_routes # sites, years, draws
  
  # sum across routes for each year:
  preds_years <- apply(preds_routes, MAR = c(2,3), FUN = sum, na.rm = TRUE) 
  
  # median and 90% credible interval:
  ts_median <- apply(preds_years, MAR = 1, FUN = median)
  
  ts_ci90 <- apply(preds_years, MAR = 1, FUN = bayestestR::ci, ci = 0.9, method = "ETI")
  ts_ci90_low <- unlist(lapply(ts_ci90, FUN = function(x) x$CI_low))
  ts_ci90_high <- unlist(lapply(ts_ci90, FUN = function(x) x$CI_high))
  
  # add 100 draws of posterior distribution for each year:
  n_draws <- 100
  draws <- t(apply(preds_years, MAR = 1, FUN = function(x) sample(x = x, size = n_draws, replace = FALSE))) 
  colnames(draws) <- paste0("draw", 1:n_draws)
  draws <- draws %>% 
    as_tibble() %>% 
    mutate(year = 1995:2019) %>% 
    select(year, everything())
    
  # draws %>% 
  #   tidyr::pivot_longer(starts_with("draw"), names_to = "draw", values_to = "value") %>% 
  #   ggplot(aes(x = year, y = value)) +
  #   geom_point() +
  #   geom_smooth(method = "lm")

  # assemble df:
  ts_preds_fact <- tibble(year = 1995:2019, median_Nocc_f = ts_median, CI90low_f = ts_ci90_low, CI90high_f = ts_ci90_high) %>% 
    left_join(draws)
  
  # # plot:
  # ggplot(ts_preds_fact) +
  #   geom_line(aes(x = year, y = median_Nocc_f)) +
  #   geom_ribbon(aes(x = year, ymax = CI90high_f, ymin = CI90low_f),
  #               alpha = 0.2, fill = "cornflowerblue") +
  #   ggtitle(spec) +
  #   theme_bw()
  
  save(ts_preds_fact, file = file.path(fact_dir, paste0(spec, "_ts_sum_occ_routes_f_preds.RData")))
  
  
  # predictions for counterfactual data incl. CI: -----------
  
  print("counterfactual predictions")
  
  # directory with predictions for counterfactual scenarios:
  cf_dir <- file.path(main_dir, "results", "attribution", "fm_y_preds_routes_cf_1995_all") # output of 5_0_DOM_y_predictions_routes_counterfactual_sets_1995.R
  
  # counterfactual prediction for this species:
  cf_files <- list.files(cf_dir, pattern = spec)
  
  # iterate over counterfactual scenarios:
  
  for(v in 1:length(cf_files)){
    
    print(v)
    
    # extract scenario from file name:
    scen <- gsub(cf_files[v], pattern = paste0("(", paste0(spec, "_y_preds_cf_"), ")|(.RData)"), replacement = "")
    
    load(file.path(cf_dir, cf_files[v]))
    y_preds_route_cf
    
    # sum across routes for each year:
    preds_years <- apply(y_preds_route_cf, MAR = c(2,3), FUN = sum, na.rm = TRUE)
    
    # median and 90% credible interval:
    ts_median <- apply(preds_years, MAR = 1, FUN = median)
    
    ts_ci90 <- apply(preds_years, MAR = 1, FUN = bayestestR::ci, ci = 0.9, method = "ETI")
    ts_ci90_low <- unlist(lapply(ts_ci90, FUN = function(x) x$CI_low))
    ts_ci90_high <- unlist(lapply(ts_ci90, FUN = function(x) x$CI_high))
    
    # add 100 draws of posterior distribution for each year:
    n_draws <- 100
    draws <- t(apply(preds_years, MAR = 1, FUN = function(x) sample(x = x, size = n_draws, replace = FALSE))) 
    colnames(draws) <- paste0("draw", 1:n_draws)
    draws <- draws %>% 
      as_tibble() %>% 
      mutate(year = 1995:2019) %>% 
      select(year, everything())
    
    if(scen == "counterclim"){

      ts_preds_cfact_clim <- tibble(year = 1995:2019, median_Nocc_cf_clim = ts_median, 
                                    CI90low_cf_clim = ts_ci90_low, CI90high_cf_clim = ts_ci90_high) %>% 
        left_join(draws)
      
    } else if(scen == "1995soc"){
      
      ts_preds_cfact_1995soc <- tibble(year = 1995:2019, median_Nocc_cf_1995soc = ts_median, 
                               CI90low_cf_1995soc = ts_ci90_low, CI90high_cf_1995soc = ts_ci90_high) %>% 
        left_join(draws)
      
      
    } else {
      
      ts_preds_cfact_clim_1995soc <- tibble(year = 1995:2019, median_Nocc_cf_clim_1995soc = ts_median, 
                                            CI90low_cf_clim_1995soc = ts_ci90_low, CI90high_cf_clim_1995soc = ts_ci90_high) %>% 
        left_join(draws)
      
    }
  }
  
  # add to list:
  ts_preds_cfact <- list(ts_preds_cfact_clim, ts_preds_cfact_1995soc, ts_preds_cfact_clim_1995soc)
  names(ts_preds_cfact) <- c("cf_clim", "cf_1995soc", "cf_clim_1995soc")
  
  save(ts_preds_cfact, file = file.path(cfact_dir, paste0(spec, "_ts_sum_occ_routes_cf_preds.RData")))
  
}
