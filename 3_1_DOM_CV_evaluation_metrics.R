# calculate evaluation measures for cross validation 
# (model predictions generated with 2_1_DOM_flocker_CV.R)

#  calculate Harrel's C indices
# - spatially, temporally


# packages: ----

library(dplyr)
library(flocker)
library(sf)
library(ggplot2)
library(cmdstanr)
#set_cmdstan_path(path = NULL)#set_cmdstan_path("C:/Users/schifferle1/Documents/cmdstan-2.34.1") # xx
set_cmdstan_path("C:/Users/schifferle1/Documents/cmdstan-2.34.1")

# functions: ----

source("0_functions.R")

# directories: ----

print(tempdir())
#dir <- file.path("/import", "ecoc9z", "data-zurell", "schifferle", "BBS_occupancy_models_2023")
dir <- getwd()
results_dir <- file.path("M:", "Documents", "DEBTs", "analysis", "Schifferle_BBS_occupancy_models_2023", "results", "CV_cluster") 
#results_dir <- file.path("results", "CV_cluster")


# load data: ----

# selected species, sorted by ecoregion:
load(file = file.path("data", "final_species_selection_eco_sorted.RData")) # final_species_eco_sorted; output of 1_2_species_selection.R

# selected routes spatial data (to buffer presences):
routes_sel_sf <- st_read(file.path("data", "route_selection_1991_2015_surv_beg_end_max_5y_miss_v2_spat_thin_100km_max_30_r_per_BCR_centroids.shp")) # output of 1_1_route_selection.R

# routes-years:
load(file = file.path("data", "BBS_for_occ_selection.RData")) # route_sel_dt; output of 1_3_match_BBS_to_env_data.R 

# route-year-species information (only surveyed)
load(file = file.path("data", "BBS_for_occ_spec_records.RData")) # bbs_dt_occ; output of 1_0_reformat_BBS_data.R

years <- 1991:2015


# data prep.: ----


for(i in 1:length(final_species_eco_sorted)){
  
  spec <- final_species_eco_sorted[i]
  
  print(paste(i, spec))
  
  # load model predictions of each test fold:
  res_lists_folds <- vector(mode = "list", length = 5)
  
  # test whether species data are there:
  skip_to_next <- FALSE
  tryCatch(print(load(file = file.path(results_dir, paste0("test_preds_", spec, "_CV_fold5.RData")))),
           error = function(e) { skip_to_next <<- TRUE})
  if(skip_to_next) { next }
  
  for(fold in 1:5){
    
    print(fold)
    
    # assemble fold results:
    load(file = file.path(results_dir, paste0("test_preds_", spec, "_CV_fold", fold, ".RData")))
    res_lists_folds[[fold]] <- res_list
  }
  
  
  # find route IDs matching the predictions:
  
  # relevant routes for the species, within distance of 750 km of presences:
  rel_routes <- training_routes(species = spec, buffer_km = 750, output = "RTENOs")
  
  # load fold assignment (sorting important!):
  load(file.path("data", "CV_route_block_allocation", "block_size_500km", paste0(spec, ".RData")))
  
  test_RTENOs <- c(sort(rel_routes[sb_US$folds_list[[1]][[2]]]), # test data fold 1
                   sort(rel_routes[sb_US$folds_list[[2]][[2]]]),
                   sort(rel_routes[sb_US$folds_list[[3]][[2]]]),
                   sort(rel_routes[sb_US$folds_list[[4]][[2]]]),
                   sort(rel_routes[sb_US$folds_list[[5]][[2]]]))
  
  
  # observations can be compared either to the predicted occupancy probability or to predictions of y, in the latter 
  # case, imperfect detection is taken into account
  
  # observations:
  occ_dt_spec <- BBS_pres_abs_spec(species = spec)
  
  
  # spatial C index: ---
  
  
  # do we get differences between sites right?
  
  # for each year separately, calculate Harrel's C index
  # C index near 1 -> good
  # C index = 0.5 -> model as good as random guessing of which of two routes has a higher probability of being occupied
  
  
  # 1) based on predictions of y (= occupancy (0/1) * detection prob.):
  
  # proportion of draws where species is detected on a route:
  y_preds_all_routes <- rbind(res_lists_folds[[1]]$y_preds_mean,
                              res_lists_folds[[2]]$y_preds_mean,
                              res_lists_folds[[3]]$y_preds_mean,
                              res_lists_folds[[4]]$y_preds_mean,
                              res_lists_folds[[5]]$y_preds_mean)
  
  y_preds_obs_df <- as.data.frame(y_preds_all_routes) %>% 
    cbind(test_RTENOs)
  colnames(y_preds_obs_df) <- c(years, "RTENO")
  
  y_preds_obs_df <- y_preds_obs_df %>% 
    tidyr::pivot_longer(cols = !RTENO, names_to = "Year", values_to = "pred_y_mean") %>% 
    mutate(Year = as.integer(Year)) %>% 
    left_join(occ_dt_spec, by = c("RTENO", "Year")) %>% 
    select(c(RTENO, Year, pred_y_mean, Count10:Count50, presence))
  
  # C-index per year:
  
  for(y in unique(y_preds_obs_df$Year)){
    
    dt <- y_preds_obs_df %>% 
      filter(Year == y)
    
    y_rank_corr <- c("Year" = y, Hmisc::rcorr.cens(x = dt$presence, S = dt$pred_y_mean, outx=FALSE)) # xx
    
    if(y == unique(y_preds_obs_df$Year)[1]) {
      y_rank_corr_spat <- y_rank_corr} else{
        y_rank_corr_spat <- rbind(y_rank_corr_spat, y_rank_corr)
      }
  }
  rownames(y_rank_corr_spat) <- NULL
  y_rank_corr_spat
  
  
  # 2) based on predictions of occupancy probability:
  
  # proportion of draws where species is detected on a route:
  occ_preds_all_routes <- rbind(res_lists_folds[[1]]$occ_posterior_mean,
                                res_lists_folds[[2]]$occ_posterior_mean,
                                res_lists_folds[[3]]$occ_posterior_mean,
                                res_lists_folds[[4]]$occ_posterior_mean,
                                res_lists_folds[[5]]$occ_posterior_mean)
  
  # occ_preds_all_routes <- rbind(res_lists_folds[[1]]$occ_posterior_median,
  #                               res_lists_folds[[2]]$occ_posterior_median,
  #                               res_lists_folds[[3]]$occ_posterior_median,
  #                               res_lists_folds[[4]]$occ_posterior_median,
  #                               res_lists_folds[[5]]$occ_posterior_median)
  
  occ_preds_obs_df <- as.data.frame(occ_preds_all_routes) %>% 
    cbind(test_RTENOs)
  colnames(occ_preds_obs_df) <- c(years, "RTENO")
  
  occ_preds_obs_df <- occ_preds_obs_df %>% 
    tidyr::pivot_longer(cols = !RTENO, names_to = "Year", values_to = "pred_occ_mean") %>% 
    mutate(Year = as.integer(Year)) %>% 
    left_join(occ_dt_spec, by = c("RTENO", "Year")) %>% 
    select(c(RTENO, Year, pred_occ_mean, Count10:Count50, presence))
  
  # C-index per year:
  
  for(y in unique(occ_preds_obs_df$Year)){
    
    dt <- occ_preds_obs_df %>% 
      filter(Year == y)
    
    occ_rank_corr <- c("Year" = y, Hmisc::rcorr.cens(x = dt$presence, S = dt$pred_occ_mean, outx=FALSE)) # xx
    
    if(y == unique(occ_preds_obs_df$Year)[1]) {
      occ_rank_corr_spat <- occ_rank_corr} else{
        occ_rank_corr_spat <- rbind(occ_rank_corr_spat, occ_rank_corr)
      }
  }
  rownames(occ_rank_corr_spat) <- NULL
  occ_rank_corr_spat
  
  
  # temporal C indices: ---
  
  
  # do we get trends right?
  
  # sum predicted detections across all routes for each year
  # vs.
  # sum observations across all routes for each year
  
  y_preds_obs_df_temp <- y_preds_obs_df %>% 
    group_by(Year) %>% 
    summarise(sum_obs = sum(presence, na.rm = TRUE),
              sum_preds = sum(pred_y_mean))
  
  C_Ind_y_temp <- Hmisc::rcorr.cens(x = y_preds_obs_df_temp$sum_preds,
                                    S = y_preds_obs_df_temp$sum_obs, outx=FALSE)
  
  
  # same for comparing observations with occupancy probability:
  
  occ_preds_obs_df_temp <- occ_preds_obs_df %>% 
    group_by(Year) %>% 
    summarise(sum_obs = sum(presence, na.rm = TRUE),
              sum_preds = sum(pred_occ_mean))
  
  
  C_Ind_occ_temp <- Hmisc::rcorr.cens(x = occ_preds_obs_df_temp$sum_preds,
                                      S = occ_preds_obs_df_temp$sum_obs, outx=FALSE)
  
  
  # save evaluation outputs:
  CV_eval <- list(y_rank_corr_spat, occ_rank_corr_spat, C_Ind_y_temp, C_Ind_occ_temp)
  names(CV_eval) <- c("C_spat_y", "C_spat_occ", "C_temp_y", "C_temp_occ")
  
  save(CV_eval, file = file.path(dir, "results", "CV_cluster", "CV_eval", paste0("CV_eval_", spec, ".RData")))
  
}


# assemble C-values for all species: ----

CV_eval_summary <- data.frame("species" = final_species_eco_sorted,
                              "y_spatial_C_median" = NA,
                              "occ_spatial_C_median" = NA,
                              "y_temp_C" = NA,
                              "occ_temp_C" = NA)


for(i in 1:length(final_species_eco_sorted)){
  
  spec <- final_species_eco_sorted[i]
  
  print(paste(i, spec))
  
  skip_to_next <- FALSE
  tryCatch(print(load(file = file.path(dir, "results", "CV_cluster", "CV_eval", paste0("CV_eval_", spec, ".RData")))),
           error = function(e) { skip_to_next <<- TRUE})
           if(skip_to_next) { next }
  
  CV_eval_summary$y_spatial_C_median[which(CV_eval_summary$species == spec)] <- median(CV_eval$C_spat_y[, "C Index"])
  CV_eval_summary$occ_spatial_C_median[which(CV_eval_summary$species == spec)] <- median(CV_eval$C_spat_occ[, "C Index"])        
  CV_eval_summary$y_temp_C[which(CV_eval_summary$species == spec)] <- CV_eval$C_temp_y["C Index"]
  CV_eval_summary$occ_temp_C[which(CV_eval_summary$species == spec)] <- CV_eval$C_temp_occ["C Index"]
  
}

write.csv(CV_eval_summary, file = file.path(dir, "results", "CV_cluster", "CV_eval", "CV_eval_summary.csv"), row.names = FALSE)

CV_eval_summary





# Briscoe metrics: ---
# xx

# plots: ----

# boxplots:
y_preds_obs_df %>% 
  filter(!is.na(presence)) %>% 
  ggplot() + 
  facet_wrap(~Year) +
  geom_boxplot(aes(group = as.factor(presence), y = pred_y_mean, x = as.factor(presence)))

# plot map comparing mean predictions and observations for single year: ---

obs_preds_sf <- y_preds_obs_df %>% 
  left_join(occ_preds_obs_df[c("RTENO", "Year", "pred_occ_mean")]) %>% 
  left_join(routes_sel_sf, by = c(RTENO = "RTENO_BBS")) %>%
  st_as_sf()

# the more similar the colours of observation and prediction dots, the better:
ggplot(obs_preds_sf) +
  geom_sf(aes(fill = as.factor(presence)), pch = 21, size = 3) +
  scale_fill_viridis_d(name = "observation") +
  geom_sf(aes(color = pred_y_mean), size = 1) +
  scale_color_viridis_c(name = "mean y. prediction") +
  theme_bw() +
  ggtitle(spec)

ggplot(obs_preds_sf) +
  geom_sf(aes(fill = as.factor(presence)), pch = 21, size = 3) +
  scale_fill_viridis_d(name = "observation") +
  geom_sf(aes(color = pred_occ_mean), size = 1) +
  scale_color_viridis_c(name = "mean occ. prediction") +
  theme_bw() +
  ggtitle(spec)


ggplot(occ_preds_obs_df_temp, aes(x = Year)) +
  geom_line(aes(y = scale(sum_preds))) +
  geom_line(aes(y = scale(sum_obs)), color = "blue")
