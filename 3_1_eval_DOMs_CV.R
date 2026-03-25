# model evaluation: spatial predictive performance based on cross validation 

# measures:
# (spatio-temporal: AUC)
# mean yearly AUC

# based on predictions of y (occ. prob. * det. prob.) for all sites (+ only for sites with change)


# packages: --------------------------------------------------------------------

library(dplyr)
library(sf)
library(ggplot2)


# directories: -----------------------------------------------------------------

print(tempdir())

# project directory:
dir <- file.path("//NAS-2-P-SN-01.ibb.uni-potsdam.de", "daten$", "AG26", "Transfer", "Schifferle_BBS_occupancy_models_2023")

results_dir <- file.path(dir, "results", "CV_buffer750km")

# directory to store CV evaluation metrics:
if(!dir.exists(file.path(results_dir, "CV_eval"))){dir.create(file.path(results_dir, "CV_eval"))}


# functions: -------------------------------------------------------------------

source("0_functions.R")


# load data: -------------------------------------------------------------------

# selected species:
load(file = file.path("data", "final_species_selection.RData")) # species_selection_final; output of 1_2_dataprep_BBS_species_selection.R

# route-year-species information (only surveyed)
load(file = file.path("data", "BBS_for_occ_spec_records.RData")) # bbs_dt_occ; output of 1_0_dataprep_BBS_bird_data.R

# routes-years:
load(file = file.path("data", "BBS_for_occ_selection.RData")) # route_sel_dt; output of 1_3_dataprep_match_BBS_routes_env_data.R

# selected routes spatial data (to buffer presences):
routes_sel_sf <- st_read(file.path("data", "route_selection_1995_2019_surv_beg_end_max_5y_miss_v2_spat_thin_100km_max_30_r_per_BCR_centroids.shp")) # output of 1_1_dataprep_BBS_route_selection.R


# species set for which MCMC fitting worked for full model, cross validation and temporal validation: ----

# species for which fitting full models failed:
load(file = file.path(dir, "results", "fm_buffer750km", "refit_2000_2000", "check_output", "specs_MCMC_failed.RData")) # output of 2_3a_fit_DOMs_check_fit.R
specs_discard_fm <- specs_MCMC_failed

# species for which fitting first 15 years failed:
load(file = file.path(dir, "results", "temp_val_buffer_750_10yrs", "refit_2000_2000", "check_output", "specs_MCMC_failed.RData")) # output of 2_3a_fit_DOMs_check_fit.R
specs_discard_tv <- specs_MCMC_failed

# species for which fitting cross validation failed for at least one fold:
load(file = file.path(dir, "results", "CV_buffer750km", "refit_2000_2000", "check_output", "specs_folds_MCMC_failed.RData"))
specs_discard_cv <- names(which(lengths(spec_folds_MCMC_fail) != 0))

final_species <- sort(subset(species_selection_final,
       !species_selection_final %in% c(specs_discard_fm, specs_discard_cv, specs_discard_tv))) # 159

#save(final_species, file = file.path("data", "species_set_analysis.RData"))

#load(file = file.path("data", "species_set_analysis.RData"))


# evaluation metrics: ----------------------------------------------------------

# iterate over species:

for(i in 1:length(final_species)){

  spec <- final_species[i]
  
  print(paste(i, spec))
  
  
  # assemble data frame with predicted mean value for each route and year and observations:
  
  # load model predictions of each test fold:
  res_lists_folds <- vector(mode = "list", length = 5)

  # test whether species data are there:
  skip_to_next <- FALSE
  tryCatch(print(load(file = file.path(results_dir, paste0("test_preds_", spec, "_CV_fold5.RData")))),
           error = function(e) { skip_to_next <<- TRUE})
  if(skip_to_next) { next }

  # test whether evaluation ran already for this species:
  if(file.exists(file.path(results_dir, "CV_eval", paste0("CV_eval_", spec, ".RData")))){
    print(paste(spec, "ran alrady."))
    next
  }

  # check where to look for model output (did MCMC fitting work with less or only with more iterations?)
  if(file.exists(file.path(results_dir, "refit_2000_2000", paste0("out_", spec, "_CV_fold1.RData")))){
    output_dir <- file.path(results_dir, "refit_2000_2000")
  } else {
    output_dir <- results_dir
  }

  print(output_dir)

  for(fold in 1:5){

    print(fold)

    # assemble fold results:
    load(file = file.path(output_dir, paste0("test_preds_", spec, "_CV_fold", fold, ".RData")))
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
  
  # observations are compared to predictions of y, where imperfect detection is taken into account

  # proportion of draws where species is predicted to be detected on a route in a year:
  y_preds_all_routes <- rbind(res_lists_folds[[1]]$y_preds_mean,
                              res_lists_folds[[2]]$y_preds_mean,
                              res_lists_folds[[3]]$y_preds_mean,
                              res_lists_folds[[4]]$y_preds_mean,
                              res_lists_folds[[5]]$y_preds_mean)
  
  y_preds_obs_df <- as.data.frame(y_preds_all_routes) %>%
    cbind(test_RTENOs)
  colnames(y_preds_obs_df) <- c(min(route_sel_dt$Year):max(route_sel_dt$Year) , "RTENO")

  # add observations:
  occ_dt_spec <- BBS_pres_abs_spec(species = spec)
  
  y_preds_obs_df <- y_preds_obs_df %>%
    tidyr::pivot_longer(cols = !RTENO, names_to = "Year", values_to = "pred_y_mean") %>%
    mutate(Year = as.integer(Year)) %>%
    left_join(occ_dt_spec, by = c("RTENO", "Year")) %>%
    select(c(RTENO, Year, pred_y_mean, Count10:Count50, presence))

  # save:
  #save(y_preds_obs_df, file = file.path(results_dir, "CV_eval", "obs_preds", paste0(spec, "_obs_ymean_preds.RData")))

  
  # spatio-temporal evaluation:

  # how well do we overall discriminate between occupied and non-occupied sites:

  # calculate overall AUC based on predictions of y (= occupancy (0/1) * detection prob.):

  load(file = file.path(results_dir, "CV_eval", "obs_preds", paste0(spec, "_obs_ymean_preds.RData")))
  y_preds_obs_df

  # all routes:

  y_auc_overall <- pROC::roc(response = y_preds_obs_df$presence, predictor = y_preds_obs_df$pred_y_mean)$auc

  # only routes with change in detection status:
  
  routes_occ_change <- y_preds_obs_df %>%
    filter(complete.cases(.)) %>%
    group_by(RTENO) %>%
    summarise(occ_change = length(unique(presence))-1) %>% # 1 = occupancy change on that route, 0 = no change
    filter(occ_change == 1) %>%
    pull(RTENO)

  y_preds_obs_df_co <- y_preds_obs_df %>%
    filter(RTENO %in% routes_occ_change)

  y_auc_overall_co <- pROC::roc(response = y_preds_obs_df_co$presence, predictor = y_preds_obs_df_co$pred_y_mean)$auc

  
  # spatial evaluation:

  # do we get differences between sites right?

  # AUC for each year, based on predictions of y (= occupancy (0/1) * detection prob.):

  # all routes:
  
  for(y in unique(y_preds_obs_df$Year)){

    print(y)
    dt <- y_preds_obs_df %>%
      filter(Year == y)

    # AUC:
    y_auc <- c("Year" = y,
               "auc" = pROC::roc(response = dt$presence, predictor = dt$pred_y_mean)$auc) 

    if(y == unique(y_preds_obs_df$Year)[1]) {
      y_auc_spat <- y_auc} else{
        y_auc_spat <- rbind(y_auc_spat, y_auc)
      }

  }

  rownames(y_auc_spat) <- NULL
  y_auc_spat

  
  # only routes with change in detection status:

  for(y in unique(y_preds_obs_df_co$Year)){

    print(y)
    dt <- y_preds_obs_df_co %>%
      filter(Year == y)

    # AUC:
    y_auc <- c("Year" = y,
               "auc" = pROC::roc(response = dt$presence, predictor = dt$pred_y_mean)$auc)

    if(y == unique(y_preds_obs_df_co$Year)[1]) {
      y_auc_spat_co <- y_auc} else{
        y_auc_spat_co <- rbind(y_auc_spat_co, y_auc)
      }
  }

  rownames(y_auc_spat_co) <- NULL
  y_auc_spat_co
  
  # save evaluation outputs:
  
  CV_eval <- list("auc_spattemp_y" = y_auc_overall, 
                  "auc_spattemp_y_c" = y_auc_overall_co, 
                  "auc_spat_y" = y_auc_spat, 
                  "auc_spat_y_c" = y_auc_spat_co)
  
  # save(CV_eval, file = file.path(results_dir, "CV_eval", paste0("CV_eval_", spec, ".RData")))
}


# assemble evaluation metrics for all species: ---- 

CV_eval_summary <- data.frame("species" = final_species,
                              "y_spattemp_auc" = NA,
                              "y_spattemp_auc_c" = NA,
                              "y_spat_auc_mean" = NA,
                              "y_spat_auc_mean_c" = NA
                              )

# iterate over species:

for(i in 1:length(final_species)){
  
  spec <- final_species[i]
  
  print(paste(i, spec))
  
  skip_to_next <- FALSE
  tryCatch(print(load(file = file.path(results_dir, "CV_eval", paste0("CV_eval_", spec, ".RData")))),
           error = function(e) { skip_to_next <<- TRUE})
  if(skip_to_next) { next }

  # overall AUC:
  CV_eval_summary$y_spattemp_auc[which(CV_eval_summary$species == spec)] <- CV_eval$auc_spattemp_y
  CV_eval_summary$y_spattemp_auc_c[which(CV_eval_summary$species == spec)] <- CV_eval$auc_spattemp_y_c

  # mean yearly AUC:
  CV_eval_summary$y_spat_auc_mean[which(CV_eval_summary$species == spec)] <- mean(CV_eval$auc_spat_y[ ,"auc"], na.rm = TRUE)
  CV_eval_summary$y_spat_auc_mean_c[which(CV_eval_summary$species == spec)] <- mean(CV_eval$auc_spat_y_c[ ,"auc"], na.rm = TRUE)

}
CV_eval_summary


#write.csv(CV_eval_summary, file = file.path(results_dir, "CV_eval", "CV_eval_summary.csv"), row.names = FALSE)
CV_eval_summary <- read.csv(file = file.path(results_dir, "CV_eval", "CV_eval_summary.csv"))
summary(CV_eval_summary)


# explorative plots: ----

# maps comparing mean CV predictions and observations:

obs_preds_sf <- y_preds_obs_df %>%
  left_join(routes_sel_sf, by = c(RTENO = "RTENO_BBS")) %>%
  st_as_sf()

# the more similar the colours of observation and prediction dots, the better:

# multiple years:
obs_preds_sf %>%
  filter(Year %in% seq(1995, 2019, length = 4)) %>%
  ggplot() +
  facet_wrap(~Year) +
  geom_sf(aes(fill = as.factor(presence)), pch = 21, size = 3) +
  scale_fill_viridis_d(name = "observation") +
  geom_sf(aes(color = pred_y_mean), size = 1) +
  scale_color_viridis_c(name = "mean y prediction") +
  theme_bw() +
  ggtitle(spec)