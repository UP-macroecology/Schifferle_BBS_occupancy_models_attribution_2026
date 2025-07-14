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
library(ggrepel)


# functions: ----

source("0_functions.R")


# directories: ----

print(tempdir())
#dir <- file.path("/import", "ecoc9z", "data-zurell", "schifferle", "BBS_occupancy_models_2023")
dir <- getwd()

# results_dir <- file.path("//NAS-2-P-SN-01.ibb.uni-potsdam.de", "users$", "schifferle1", "Documents", "DEBTs", "analysis",
#                          "Schifferle_BBS_occupancy_models_2023", "results", "CV_cluster")
results_dir <- file.path("M:", "Documents", "DEBTs", "analysis", "Schifferle_BBS_occupancy_models_2023", 
                         "results", "CV_buffer750km")

# directory to store CV evaluation metrics:
if(!dir.exists(file.path(results_dir, "CV_eval"))){
  dir.create(file.path(results_dir, "CV_eval"))
}

# load data: ----

# selected species, sorted by ecoregion:
load(file = file.path("data", "final_species_selection_eco_sorted.RData")) # final_species_eco_sorted; output of 1_2_species_selection.R

# selected routes spatial data (to buffer presences):
routes_sel_sf <- st_read(file.path("data", "route_selection_1995_2019_surv_beg_end_max_5y_miss_v2_spat_thin_100km_max_30_r_per_BCR_centroids.shp")) # output of 1_1_route_selection.R

# routes-years:
load(file = file.path("data", "BBS_for_occ_selection.RData")) # route_sel_dt; output of 1_3_match_BBS_to_env_data.R 

# route-year-species information (only surveyed)
load(file = file.path("data", "BBS_for_occ_spec_records.RData")) # bbs_dt_occ; output of 1_0_reformat_BBS_data.R

years <- min(route_sel_dt$Year):max(route_sel_dt$Year) 


# # determine final species set for which MCMC fitting worked for full model, cross validation and temporal validation: ---
# 
# # species for which fitting full models failed:
# load(file = file.path("M:", "Documents", "DEBTs", "analysis", "Schifferle_BBS_occupancy_models_2023",
#                       "results", "fm_buffer750km", "refit_2000_2000", "check_output", 
#                       "specs_MCMC_failed.RData"))
# specs_discard_fm <- specs_MCMC_failed # black-billed magpie and brown-headed nuthatch with divergent transitions
# 
# # species for which fitting first 15 years failed:
# load(file = file.path("M:", "Documents", "DEBTs", "analysis", "Schifferle_BBS_occupancy_models_2023",
#                       "results", "temp_val_buffer_750_10yrs", "refit_2000_2000", "check_output", 
#                       "specs_MCMC_failed.RData"))
# specs_discard_tv <- specs_MCMC_failed
# 
# # species for which fitting cross validation failed for at least one fold:
# load(file = file.path("M:", "Documents", "DEBTs", "analysis", "Schifferle_BBS_occupancy_models_2023",
#                       "results", "CV_buffer750km", "refit_2000_2000", "check_output", 
#                       "specs_folds_MCMC_failed.RData"))
# specs_discard_cv <- names(which(lengths(spec_folds_MCMC_fail) != 0))
# 
# final_species <- sort(subset(final_species_eco_sorted, 
#        !final_species_eco_sorted %in% c(specs_discard_fm, specs_discard_cv, specs_discard_tv))) # 159
# 
# save(final_species, file = file.path("data", "species_set_analysis.RData"))

# selected species:
load(file = file.path("data", "species_set_analysis.RData"))
final_species

# evaluation metrics calculation: ----

for(i in 1:length(final_species)){

  spec <- final_species[i]
  
  print(paste(i, spec))
  
  # # load model predictions of each test fold:
  # res_lists_folds <- vector(mode = "list", length = 5)
  # 
  # # test whether species data are there:
  # skip_to_next <- FALSE
  # tryCatch(print(load(file = file.path(results_dir, paste0("test_preds_", spec, "_CV_fold5.RData")))),
  #          error = function(e) { skip_to_next <<- TRUE})
  # if(skip_to_next) { next }
  # 
  # # # test whether script ran already for this species:
  # # if(file.exists(file.path(results_dir, "CV_eval", paste0("CV_eval_", spec, ".RData")))){
  # #   print(paste(spec, "ran alrady."))
  # #   next
  # # }
  # 
  # # check where to look for model output (did MCMC fitting work with less or only with more iterations?)
  # if(file.exists(file.path(results_dir, "refit_2000_2000", paste0("out_", spec, "_CV_fold1.RData")))){
  #   output_dir <- file.path(results_dir, "refit_2000_2000")
  # } else {
  #   output_dir <- results_dir
  # }
  # 
  # print(output_dir)
  
  # for(fold in 1:5){
  #   
  #   print(fold)
  #   
  #   # assemble fold results:
  #   load(file = file.path(output_dir, paste0("test_preds_", spec, "_CV_fold", fold, ".RData")))
  #   res_lists_folds[[fold]] <- res_list
  # }
  # 
  # # find route IDs matching the predictions:
  # 
  # # relevant routes for the species, within distance of 750 km of presences:
  # rel_routes <- training_routes(species = spec, buffer_km = 750, output = "RTENOs")
  # 
  # # load fold assignment (sorting important!):
  # load(file.path("data", "CV_route_block_allocation", "block_size_500km", paste0(spec, ".RData")))
  # 
  # test_RTENOs <- c(sort(rel_routes[sb_US$folds_list[[1]][[2]]]), # test data fold 1
  #                  sort(rel_routes[sb_US$folds_list[[2]][[2]]]),
  #                  sort(rel_routes[sb_US$folds_list[[3]][[2]]]),
  #                  sort(rel_routes[sb_US$folds_list[[4]][[2]]]),
  #                  sort(rel_routes[sb_US$folds_list[[5]][[2]]]))
  # 
  # 
  # # observations can be compared either to the predicted occupancy probability or to predictions of y, in the latter 
  # # case, imperfect detection is taken into account
  # 
  # # observations:
  # occ_dt_spec <- BBS_pres_abs_spec(species = spec)
  # 
  # 
  # # spatio-temporal C index: ---
  # 
  # # how well do we overall discriminate between occupied and non-occupied sites:
  # 
  # # calculate overall Harrel's C index:
  # 
  # # 1) based on predictions of y (= occupancy (0/1) * detection prob.):
  # 
  # # proportion of draws where species is predicted to be detected on a route:
  # y_preds_all_routes <- rbind(res_lists_folds[[1]]$y_preds_mean,
  #                             res_lists_folds[[2]]$y_preds_mean,
  #                             res_lists_folds[[3]]$y_preds_mean,
  #                             res_lists_folds[[4]]$y_preds_mean,
  #                             res_lists_folds[[5]]$y_preds_mean)
  # 
  # y_preds_obs_df <- as.data.frame(y_preds_all_routes) %>% 
  #   cbind(test_RTENOs)
  # colnames(y_preds_obs_df) <- c(years, "RTENO")
  # 
  # # add observations:
  # y_preds_obs_df <- y_preds_obs_df %>% 
  #   tidyr::pivot_longer(cols = !RTENO, names_to = "Year", values_to = "pred_y_mean") %>% 
  #   mutate(Year = as.integer(Year)) %>% 
  #   left_join(occ_dt_spec, by = c("RTENO", "Year")) %>% 
  #   select(c(RTENO, Year, pred_y_mean, Count10:Count50, presence))
  # 
  # # save:
  # save(y_preds_obs_df, file = file.path(results_dir, "CV_eval", "obs_preds", 
  #                                       paste0(spec, "_obs_ymean_preds.RData")))
  # 
  load(file = file.path(results_dir, "CV_eval", "obs_preds", paste0(spec, "_obs_ymean_preds.RData")))
  y_preds_obs_df
  
  # routes with occupancy change (evaluate also only with regard to those):
  routes_occ_change <- y_preds_obs_df %>%
    filter(complete.cases(.)) %>%
    group_by(RTENO) %>%
    summarise(occ_change = length(unique(presence))-1) %>% # 1 = occupancy change on that route, 0 = no change
    filter(occ_change == 1) %>%
    pull(RTENO)

  y_preds_obs_df_co <- y_preds_obs_df %>%
    filter(RTENO %in% routes_occ_change)

  # all routes:
  y_Cind <- Hmisc::rcorr.cens(x = y_preds_obs_df$pred_y_mean, S = y_preds_obs_df$presence, outx=FALSE) # xx
  y_auc_overall <- pROC::roc(response = y_preds_obs_df$presence, predictor = y_preds_obs_df$pred_y_mean)$auc

  # only routes with change in occupancy:
  y_Cind_co <- Hmisc::rcorr.cens(x = y_preds_obs_df_co$pred_y_mean, S = y_preds_obs_df_co$presence, outx=FALSE) # xx
  y_auc_overall_co <- pROC::roc(response = y_preds_obs_df_co$presence, predictor = y_preds_obs_df_co$pred_y_mean)$auc

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
  
  # add observations:
  occ_preds_obs_df <- occ_preds_obs_df %>% 
    tidyr::pivot_longer(cols = !RTENO, names_to = "Year", values_to = "pred_occ_mean") %>% 
    mutate(Year = as.integer(Year)) %>% 
    left_join(occ_dt_spec, by = c("RTENO", "Year")) %>% 
    select(c(RTENO, Year, pred_occ_mean, Count10:Count50, presence))
  
  # save:
  save(occ_preds_obs_df, file = file.path(results_dir, "CV_eval", "obs_preds", 
                                        paste0(spec, "_obs_occmean_preds.RData")))
  
  # # all routes:
  # occ_Cind <- Hmisc::rcorr.cens(x = occ_preds_obs_df$pred_occ_mean, S = occ_preds_obs_df$presence, outx=FALSE)
  # occ_auc_overall <- pROC::roc(response = occ_preds_obs_df$presence, predictor = occ_preds_obs_df$pred_occ_mean)$auc
  # 
  # 
  # # only routes with change in occupancy:
  # 
  # occ_preds_obs_df_co <- occ_preds_obs_df %>% 
  #   filter(RTENO %in% routes_occ_change)
  # 
  # occ_Cind_co <- Hmisc::rcorr.cens(x = occ_preds_obs_df_co$pred_occ_mean, S = occ_preds_obs_df_co$presence, outx=FALSE) # xx
  # occ_auc_overall_co <- pROC::roc(response = occ_preds_obs_df_co$presence, predictor = occ_preds_obs_df_co$pred_occ_mean)$auc
  # 
  # 
  # # spatial C index: ---
  # 
  # # do we get differences between sites right?
  # 
  # # for each year separately, calculate Harrel's C index
  # # C index near 1 -> good
  # # C index = 0.5 -> model as good as random guessing of which of two routes has a higher probability of being occupied
  # 
  # 
  # # 1) based on predictions of y (= occupancy (0/1) * detection prob.):
  # 
  # # C-index / AUC per year:
  # 
  # # all routes:
  # for(y in unique(y_preds_obs_df$Year)){
  #   
  #   print(y)
  #   dt <- y_preds_obs_df %>% 
  #     filter(Year == y)
  #   
  #   # C index:
  #   y_rank_corr <- c("Year" = y, Hmisc::rcorr.cens(x = dt$pred_y_mean, S = dt$presence , outx=FALSE)) # xx
  #   
  #   if(y == unique(y_preds_obs_df$Year)[1]) {
  #     y_rank_corr_spat <- y_rank_corr} else{
  #       y_rank_corr_spat <- rbind(y_rank_corr_spat, y_rank_corr)
  #     }
  #   
  #   # AUC:
  #   y_auc <- c("Year" = y, 
  #              "auc" = tryCatch(pROC::roc(response = dt$presence, predictor = dt$pred_y_mean)$auc,
  #                               error = function(e) {NA})) # xx
  #   
  #   if(y == unique(y_preds_obs_df$Year)[1]) {
  #     y_auc_spat <- y_auc} else{
  #       y_auc_spat <- rbind(y_auc_spat, y_auc)
  #     }
  #   
  # }
  # rownames(y_rank_corr_spat) <- NULL
  # y_rank_corr_spat
  # rownames(y_auc_spat) <- NULL
  # y_auc_spat
  # 
  # # only routes with occupancy change:
  # 
  # for(y in unique(y_preds_obs_df_co$Year)){
  #   
  #   print(y)
  #   dt <- y_preds_obs_df_co %>% 
  #     filter(Year == y)
  #   
  #   # C index:
  #   y_rank_corr <- c("Year" = y, Hmisc::rcorr.cens(x = dt$pred_y_mean, S = dt$presence , outx=FALSE)) # xx
  #   
  #   if(y == unique(y_preds_obs_df_co$Year)[1]) {
  #     y_rank_corr_spat_co <- y_rank_corr} else{
  #       y_rank_corr_spat_co <- rbind(y_rank_corr_spat_co, y_rank_corr)
  #     }
  #   
  #   # AUC:
  #   y_auc <- c("Year" = y, 
  #              "auc" = tryCatch(pROC::roc(response = dt$presence, predictor = dt$pred_y_mean)$auc,
  #                               error = function(e) {NA})) # xx
  #   
  #   if(y == unique(y_preds_obs_df_co$Year)[1]) {
  #     y_auc_spat_co <- y_auc} else{
  #       y_auc_spat_co <- rbind(y_auc_spat_co, y_auc)
  #     }
  #   
  # }
  # rownames(y_rank_corr_spat_co) <- NULL
  # y_rank_corr_spat_co
  # rownames(y_auc_spat_co) <- NULL
  # y_auc_spat_co
  # 
  # 
  # # 2) based on predictions of occupancy probability:
  # 
  # # C-index / AUC per year:
  # 
  # # all routes:
  # 
  # for(y in unique(occ_preds_obs_df$Year)){
  #   
  #   dt <- occ_preds_obs_df %>% 
  #     filter(Year == y)
  #   
  #   occ_rank_corr <- c("Year" = y, Hmisc::rcorr.cens(x = dt$pred_occ_mean, S = dt$presence, outx=FALSE)) # xx
  #   
  #   if(y == unique(occ_preds_obs_df$Year)[1]) {
  #     occ_rank_corr_spat <- occ_rank_corr} else{
  #       occ_rank_corr_spat <- rbind(occ_rank_corr_spat, occ_rank_corr)
  #     }
  #   
  #   occ_auc <- c("Year" = y, 
  #                "auc" = tryCatch(pROC::roc(response = dt$presence, predictor = dt$pred_occ_mean)$auc,
  #                                             error = function(e) {NA})) # xx
  #   
  #   if(y == unique(y_preds_obs_df$Year)[1]) {
  #     occ_auc_spat <- occ_auc} else{
  #       occ_auc_spat <- rbind(occ_auc_spat, occ_auc)
  #     }
  #   
  #   
  # }
  # rownames(occ_rank_corr_spat) <- NULL
  # occ_rank_corr_spat
  # rownames(occ_auc_spat) <- NULL
  # occ_auc_spat
  # 
  # # only routes with occupancy change:
  # 
  # for(y in unique(occ_preds_obs_df_co$Year)){
  #   
  #   dt <- occ_preds_obs_df_co %>% 
  #     filter(Year == y)
  #   
  #   occ_rank_corr <- c("Year" = y, Hmisc::rcorr.cens(x = dt$pred_occ_mean, S = dt$presence, outx=FALSE)) # xx
  #   
  #   if(y == unique(occ_preds_obs_df_co$Year)[1]) {
  #     occ_rank_corr_spat_co <- occ_rank_corr} else{
  #       occ_rank_corr_spat_co <- rbind(occ_rank_corr_spat_co, occ_rank_corr)
  #     }
  #   
  #   occ_auc <- c("Year" = y, 
  #                "auc" = tryCatch(pROC::roc(response = dt$presence, predictor = dt$pred_occ_mean)$auc,
  #                                 error = function(e) {NA})) # xx
  #   
  #   if(y == unique(occ_preds_obs_df_co$Year)[1]) {
  #     occ_auc_spat_co <- occ_auc} else{
  #       occ_auc_spat_co <- rbind(occ_auc_spat_co, occ_auc)
  #     }
  #   
  #   
  # }
  # rownames(occ_rank_corr_spat_co) <- NULL
  # occ_rank_corr_spat_co
  # rownames(occ_auc_spat_co) <- NULL
  # occ_auc_spat_co
  # 
  # 
  # 
  # temporal C indices: ---


  # do we get trends right?

  # sum predicted detections across all routes for each year
  # vs.
  # sum observations across all routes for each year

  # all routes:

  y_preds_obs_df_temp <- y_preds_obs_df %>%
    group_by(Year) %>%
    summarise(sum_obs = sum(presence, na.rm = TRUE),
              sum_preds = sum(pred_y_mean)) # xx



  C_Ind_y_temp <- Hmisc::rcorr.cens(x = y_preds_obs_df_temp$sum_preds,
                                    S = y_preds_obs_df_temp$sum_obs, outx=FALSE)

  # only routes with occupancy change:

  y_preds_obs_df_temp_co <- y_preds_obs_df_co %>%
    group_by(Year) %>%
    summarise(sum_obs = sum(presence, na.rm = TRUE),
              sum_preds = sum(pred_y_mean)) # xx

  C_Ind_y_temp_co <- Hmisc::rcorr.cens(x = y_preds_obs_df_temp_co$sum_preds,
                                    S = y_preds_obs_df_temp_co$sum_obs, outx=FALSE)


  # # plot time series:
  # # plot predicted y against observations summed over years:
  #
  # # # directory to save plots:
  # # if(!dir.exists(file.path(results_dir, "CV_eval", "plots_temp_performance"))){
  # #   dir.create(file.path(results_dir, "CV_eval", "plots_temp_performance"), recursive = TRUE)
  # # }
  # #
  # # # jpeg(file = file.path(results_dir, "CV_eval", "plots_temp_performance", paste0("CV_temp_C_", spec, "2.jpg")),
  # # #      width = 1000, height = 700, quality = 100)
  # print(
  #   ggplot(y_preds_obs_df_temp, aes(x = Year)) +
  #     geom_line(aes(y = sum_obs)) +
  #     geom_point(aes(y = sum_obs), size = 3) +
  #     geom_line(aes(y = sum_preds/100), color = "cornflowerblue") + # / 100 because proportion was saved as percent to be integer
  #     geom_point(aes(y = sum_preds/100), color = "cornflowerblue", size = 3) +
  #     #geom_smooth(aes(x = Year, y = sum_obs), method = "lm", color = "black") +
  #     #geom_smooth(aes(x = Year, y = sum_preds/100), method = "lm", color = "cornflowerblue") +
  #     ylab("N routes with presence") +
  #     theme_bw() +
  #     theme(text = element_text(size = 20)) +
  #     #ggtitle(paste(spec, "C index temp. CV", round(C_Ind_y_temp["C Index"],2)))
  #   ggtitle(paste0(spec, ", C-index ", round(C_Ind_y_temp["C Index"],2)))
  # )
  # # # dev.off()
  # 
  # 
  # # same for comparing observations with occupancy probability:
  # 
  # # all routes:
  # 
  # occ_preds_obs_df_temp <- occ_preds_obs_df %>% 
  #   group_by(Year) %>% 
  #   summarise(sum_obs = sum(presence, na.rm = TRUE),
  #             sum_preds = sum(pred_occ_mean))
  # 
  # C_Ind_occ_temp <- Hmisc::rcorr.cens(x = occ_preds_obs_df_temp$sum_preds,
  #                                     S = occ_preds_obs_df_temp$sum_obs, outx=FALSE)
  # 
  # # only routes with occupancy change:
  # occ_preds_obs_df_temp_co <- occ_preds_obs_df_co %>% 
  #   group_by(Year) %>% 
  #   summarise(sum_obs = sum(presence, na.rm = TRUE),
  #             sum_preds = sum(pred_occ_mean))
  # 
  # C_Ind_occ_temp_co <- Hmisc::rcorr.cens(x = occ_preds_obs_df_temp_co$sum_preds,
  #                                     S = occ_preds_obs_df_temp_co$sum_obs, outx=FALSE)
  # 
  # # overall trend captured?
  # 
  # # all routes:
  # 
  # lm_obs <- lm(sum_obs  ~ Year, data = y_preds_obs_df_temp)
  # lm_preds <- lm(sum_preds  ~ Year, data = y_preds_obs_df_temp)
  # lm_obs_preds <- predict(lm_obs, y_preds_obs_df_temp)
  # lm_preds_preds <- predict(lm_preds, y_preds_obs_df_temp)
  # trend_corr_CV <- cor(lm_obs_preds, lm_preds_preds) # 1 = trend captured, -1 = not captured
  # 
  # # only routes with occupancy change:
  # lm_obs <- lm(sum_obs  ~ Year, data = y_preds_obs_df_temp_co)
  # lm_preds <- lm(sum_preds  ~ Year, data = y_preds_obs_df_temp_co)
  # lm_obs_preds <- predict(lm_obs, y_preds_obs_df_temp_co)
  # lm_preds_preds <- predict(lm_preds, y_preds_obs_df_temp_co)
  # trend_corr_CV_c <- cor(lm_obs_preds, lm_preds_preds) # 1 = trend captured, -1 = not captured
  # 
  # # quantify absolute error, find whether there is a trend in the error: 
  # 
  # y_preds_obs_df_temp <- y_preds_obs_df_temp %>% 
  #   mutate(abs_error = abs((sum_preds/100) - sum_obs)) %>% 
  #   mutate(diff_preds_obs = (sum_preds/100) - sum_obs) %>% 
  #   mutate(Year_Ind = 1:25)
  # 
  # # linear trend in absolute error:
  # lm_error <- lm(abs_error  ~ Year_Ind, data = y_preds_obs_df_temp)
  # slope_abs_error <- lm_error$coefficients[2] # slope
  # p_slope_abs_error <- summary(lm_error)$coefficients[2,4] # slope significance
  # intercept_abs_error <- lm_error$coefficients[1] # intercept
  # 
  # # linear trend in (predictions - observations):
  # lm_error <- lm(diff_preds_obs  ~ Year_Ind, data = y_preds_obs_df_temp)
  # slope_diff_p_o <- lm_error$coefficients[2] # slope
  # p_slope_diff_p_o <- summary(lm_error)$coefficients[2,4] # slope significance
  # intercept_diff_p_o <- lm_error$coefficients[1] # intercept
  # 
  # # same for only routes with change:
  # 
  # y_preds_obs_df_temp_co <- y_preds_obs_df_temp_co %>% 
  #   mutate(abs_error = abs((sum_preds/100) - sum_obs)) %>% 
  #   mutate(diff_preds_obs = (sum_preds/100) - sum_obs) %>% 
  #   mutate(Year_Ind = 1:25)
  # 
  # # linear trend in absolute error:
  # lm_error <- lm(abs_error  ~ Year_Ind, data = y_preds_obs_df_temp_co)
  # slope_abs_error_co <- lm_error$coefficients[2] # slope
  # p_slope_abs_error_co <- summary(lm_error)$coefficients[2,4] # slope significance
  # intercept_abs_error_co <- lm_error$coefficients[1] # intercept
  # 
  # # linear trend in (predictions - observations):
  # lm_error <- lm(diff_preds_obs  ~ Year_Ind, data = y_preds_obs_df_temp_co)
  # slope_diff_p_o_co <- lm_error$coefficients[2] # slope
  # p_slope_diff_p_o_co <- summary(lm_error)$coefficients[2,4] # slope significance
  # intercept_diff_p_o_co <- lm_error$coefficients[1] # intercept
  # 
  # 
  # 
  # # save evaluation outputs:
  # CV_eval <- list(y_Cind, y_Cind_co, y_auc_overall, y_auc_overall_co, 
  #                 occ_Cind, occ_Cind_co, occ_auc_overall, occ_auc_overall_co,
  #                 y_rank_corr_spat, y_rank_corr_spat_co, y_auc_spat, y_auc_spat_co, 
  #                 occ_rank_corr_spat, occ_rank_corr_spat_co, occ_auc_spat, occ_auc_spat_co,
  #                 C_Ind_y_temp, C_Ind_y_temp_co, C_Ind_occ_temp, C_Ind_occ_temp_co,
  #                 trend_corr_CV, trend_corr_CV_c,
  #                 slope_abs_error, p_slope_abs_error, intercept_abs_error,
  #                 slope_diff_p_o, p_slope_diff_p_o, intercept_diff_p_o,
  #                 slope_abs_error_co, p_slope_abs_error_co, intercept_abs_error_co,
  #                 slope_diff_p_o_co, p_slope_diff_p_o_co, intercept_diff_p_o_co)
  # 
  # names(CV_eval) <- c("C_spattemp_y", "C_spattemp_y_c", "auc_spattemp_y", "auc_spattemp_y_c",
  #                     "C_spattemp_occ", "C_spattemp_occ_c", "auc_spattemp_occ", "auc_spattemp_occ_c",
  #                     "C_spat_y", "C_spat_y_c", "auc_spat_y", "auc_spat_y_c", 
  #                     "C_spat_occ", "C_spat_occ_c", "auc_spat_occ", "auc_spat_occ_c", 
  #                     "C_temp_y", "C_temp_y_c", "C_temp_occ", "C_temp_occ_c",
  #                     "trend_corr_CV", "trend_corr_CV_c",
  #                     "slope_abs_err", "p_slope_abs_err", "int_abs_err",
  #                     "slope_diffpo", "p_slope_diffpo", "int_diffpo",
  #                     "slope_abs_err_c", "p_slope_abs_err_c", "int_abs_err_c",
  #                     "slope_diffpo_c", "p_slope_diffpo_c", "int_diff_p_o_c")
  # 
  # save(CV_eval, file = file.path(results_dir, "CV_eval", paste0("CV_eval_", spec, "_2.RData")))
  # 
}


# assemble C-values for all species: ---- xx

CV_eval_summary <- data.frame("species" = final_species,
                              "y_spattemp_C" = NA,
                              "y_spattemp_C_c" = NA, # only for routes with occupancy change
                              "occ_spattemp_C" = NA,
                              "occ_spattemp_C_c" = NA,
                              "y_spattemp_auc" = NA,
                              "y_spattemp_auc_c" = NA,
                              "occ_spattemp_auc" = NA,
                              "occ_spattemp_auc_c" = NA,
                              "y_spat_C_mean" = NA,
                              "y_spat_C_mean_c" = NA,
                              "occ_spat_C_mean" = NA,
                              "occ_spat_C_mean_c" = NA,
                              "y_spat_auc_mean" = NA,
                              "y_spat_auc_mean_c" = NA,
                              "occ_spat_auc_mean" = NA,
                              "occ_spat_auc_mean_c" = NA,
                              "y_temp_C" = NA,
                              "y_temp_C_c" = NA,
                              "occ_temp_C" = NA,
                              "occ_temp_C_c" = NA,
                              "trend_corr_CV"= NA,
                              "trend_corr_CV_c" = NA)


for(i in 1:length(final_species)){
  
  spec <- final_species[i]
  
  print(paste(i, spec))
  
  skip_to_next <- FALSE
  tryCatch(print(load(file = file.path(results_dir, "CV_eval", paste0("CV_eval_", spec, ".RData")))),
           error = function(e) { skip_to_next <<- TRUE})
           if(skip_to_next) { next }
  
  CV_eval_summary$y_spattemp_C[which(CV_eval_summary$species == spec)] <- CV_eval$C_spattemp_y["C Index"]
  CV_eval_summary$y_spattemp_C_c[which(CV_eval_summary$species == spec)] <- CV_eval$C_spattemp_y_c["C Index"]
  CV_eval_summary$occ_spattemp_C[which(CV_eval_summary$species == spec)] <- CV_eval$C_spattemp_occ["C Index"]
  CV_eval_summary$occ_spattemp_C_c[which(CV_eval_summary$species == spec)] <- CV_eval$C_spattemp_occ_c["C Index"]
  
  CV_eval_summary$y_spattemp_auc[which(CV_eval_summary$species == spec)] <- CV_eval$auc_spattemp_y
  CV_eval_summary$y_spattemp_auc_c[which(CV_eval_summary$species == spec)] <- CV_eval$auc_spattemp_y_c
  CV_eval_summary$occ_spattemp_auc[which(CV_eval_summary$species == spec)] <- CV_eval$auc_spattemp_occ
  CV_eval_summary$occ_spattemp_auc_c[which(CV_eval_summary$species == spec)] <- CV_eval$auc_spattemp_occ_c
  
  CV_eval_summary$y_spat_C_mean[which(CV_eval_summary$species == spec)] <- mean(CV_eval$C_spat_y[, "C Index"], na.rm = TRUE)
  CV_eval_summary$y_spat_C_mean_c[which(CV_eval_summary$species == spec)] <- mean(CV_eval$C_spat_y_c[, "C Index"], na.rm = TRUE)
  CV_eval_summary$occ_spat_C_mean[which(CV_eval_summary$species == spec)] <- mean(CV_eval$C_spat_occ[, "C Index"], na.rm = TRUE)  
  CV_eval_summary$occ_spat_C_mean_c[which(CV_eval_summary$species == spec)] <- mean(CV_eval$C_spat_occ_c[, "C Index"], na.rm = TRUE)
  
  CV_eval_summary$y_spat_auc_mean[which(CV_eval_summary$species == spec)] <- mean(CV_eval$auc_spat_y[ ,"auc"], na.rm = TRUE)
  CV_eval_summary$y_spat_auc_mean_c[which(CV_eval_summary$species == spec)] <- mean(CV_eval$auc_spat_y_c[ ,"auc"], na.rm = TRUE)
  CV_eval_summary$occ_spat_auc_mean[which(CV_eval_summary$species == spec)] <- mean(CV_eval$auc_spat_occ[ ,"auc"], na.rm = TRUE)
  CV_eval_summary$occ_spat_auc_mean_c[which(CV_eval_summary$species == spec)] <- mean(CV_eval$auc_spat_occ_c[ ,"auc"], na.rm = TRUE)
  
  CV_eval_summary$y_temp_C[which(CV_eval_summary$species == spec)] <- CV_eval$C_temp_y["C Index"]
  CV_eval_summary$y_temp_C_c[which(CV_eval_summary$species == spec)] <- CV_eval$C_temp_y_c["C Index"]
  CV_eval_summary$occ_temp_C[which(CV_eval_summary$species == spec)] <- CV_eval$C_temp_occ["C Index"]
  CV_eval_summary$occ_temp_C_c[which(CV_eval_summary$species == spec)] <- CV_eval$C_temp_occ_c["C Index"]
  
  CV_eval_summary$trend_corr_CV[which(CV_eval_summary$species == spec)] <- CV_eval$trend_corr_CV
  CV_eval_summary$trend_corr_CV_c[which(CV_eval_summary$species == spec)] <- CV_eval$trend_corr_CV_c
  
}
CV_eval_summary

#write.csv(CV_eval_summary, file = file.path(results_dir, "CV_eval", "CV_eval_summary.csv"), row.names = FALSE)
CV_eval_summary <- read.csv(file = file.path(results_dir, "CV_eval", "CV_eval_summary.csv"))
summary(CV_eval_summary)

CV_eval_summary %>% 
  filter(y_temp_C_c > 0.9) %>% 
  pull(species)



load(file.path(results_dir, "spec_set_temp_val_ok1.RData")) # xx changed

CV_eval_summary %>% 
  select(species, y_spat_C_mean, y_spat_C_mean_c, y_spat_auc_mean, y_spat_auc_mean_c, 
         occ_spat_C_mean, occ_spat_C_mean_c, occ_spat_auc_mean, occ_spat_auc_mean_c) %>% 
  View

CV_eval_summary %>% 
  select(species, y_spat_C_mean, y_spat_C_mean_c, y_spat_auc_mean, y_spat_auc_mean_c, 
         occ_spat_C_mean, occ_spat_C_mean_c, occ_spat_auc_mean, occ_spat_auc_mean_c) %>% 
  filter(species %in% specs_thresh4) %>% 
  View


# xx calculate correlation and error metrics: ----

CV_temp_metrics <- data.frame("species" = final_species, 
                               "cor_p" = NA,
                               "cor_s" = NA,
                               "mean_error" = NA,
                               "mean_square_error" = NA,
                               "root_mean_square_error" = NA,
                               "mean_abs_error" = NA,
                               "mean_n_routes_obs" = NA,
                              "slope_abs_error" = NA, 
                              "p_slope_abs_error" = NA,
                              "intercept_abs_error" = NA,
                              "slope_diffpo" = NA,
                              "p_slope_diffpo" = NA,
                              "intercept_diffpo" = NA)
# load intermediate file observations - predictions:

for(i in 1:nrow(CV_temp_metrics)){
  
  spec <- final_species[i]
  
  print(paste(i, spec))
  
  load(file = file.path(results_dir, "CV_eval", "obs_preds", paste0(spec, "_obs_ymean_preds.RData"))) # xx occ
  
  # sum across years: 
  y_preds_obs_df_temp <- y_preds_obs_df %>%
    group_by(Year) %>%
    summarise(sum_obs = sum(presence, na.rm = TRUE),
              sum_preds = sum(pred_y_mean)/100)
  
  # # plot to check:
  # print(
  #   ggplot(y_preds_obs_df_temp, aes(x = Year)) +
  #     geom_line(aes(y = sum_obs)) +
  #     geom_point(aes(y = sum_obs), size = 3) +
  #     geom_line(aes(y = sum_preds), color = "cornflowerblue") + # / 100 because proportion was saved as percent to be integer
  #     geom_point(aes(y = sum_preds), color = "cornflowerblue", size = 3) +
  #     geom_smooth(aes(x = Year, y = sum_obs), method = "lm", color = "black") +
  #     geom_smooth(aes(x = Year, y = sum_preds), method = "lm", color = "cornflowerblue") +
  #     ylab("N routes with presence") +
  #     theme_bw() +
  #     theme(text = element_text(size = 20)) +
  #     ggtitle(spec)
  # )
  
  # pearson correlation:
  CV_temp_metrics$cor_p[i] <- cor(y_preds_obs_df_temp$sum_obs, y_preds_obs_df_temp$sum_preds, method = "p")
  
  # rank correlation:
  CV_temp_metrics$cor_s[i] <- cor(y_preds_obs_df_temp$sum_obs, y_preds_obs_df_temp$sum_preds, method = "s")
  
  # mean error:
  CV_temp_metrics$mean_error[i] <- mean(y_preds_obs_df_temp$sum_preds - y_preds_obs_df_temp$sum_obs)
  
  # mean square error:
  CV_temp_metrics$mean_square_error[i] <- mean((y_preds_obs_df_temp$sum_preds - y_preds_obs_df_temp$sum_obs)^2)
  
  # root mean square error:
  CV_temp_metrics$root_mean_square_error[i] <- sqrt(mean((y_preds_obs_df_temp$sum_preds - y_preds_obs_df_temp$sum_obs)^2))
  
  # mean absolute error:
  CV_temp_metrics$mean_abs_error[i] <- mean(abs(y_preds_obs_df_temp$sum_preds - y_preds_obs_df_temp$sum_obs))
  
  # mean number of routes with presences observed (for scaling?):
  CV_temp_metrics$mean_n_routes_obs[i] <- mean(y_preds_obs_df_temp$sum_obs)
  
  # linear trend in absolute error:
  trend_ae_df <- data.frame("Year_Ind" = 1:25, "abs_error" = abs(y_preds_obs_df_temp$sum_preds - y_preds_obs_df_temp$sum_obs))
  lm_error <- lm(abs_error  ~ Year_Ind, data = trend_ae_df)
  CV_temp_metrics$slope_abs_error[i] <- lm_error$coefficients[2] # slope
  CV_temp_metrics$p_slope_abs_error[i] <- summary(lm_error)$coefficients[2,4] # slope significance
  CV_temp_metrics$intercept_abs_error[i] <- lm_error$coefficients[1] # intercept
  
  # linear trend in diff. preds. - obs.:
  trend_diffpo_df <- data.frame("Year_Ind" = 1:25, "diffpo" = y_preds_obs_df_temp$sum_preds - y_preds_obs_df_temp$sum_obs)
  lm_error2 <- lm(diffpo  ~ Year_Ind, data = trend_diffpo_df)
  CV_temp_metrics$slope_diffpo[i] <- lm_error2$coefficients[2] # slope
  CV_temp_metrics$p_slope_diffpo[i] <- summary(lm_error2)$coefficients[2,4] # slope significance
  CV_temp_metrics$intercept_diffpo[i] <- lm_error2$coefficients[1] # intercept
  
  print(CV_temp_metrics[i,])
}

#save(CV_temp_metrics, file = file.path(results_dir, "CV_eval", "CV_temp_metrics.RData"))
load(file = file.path(results_dir, "CV_eval", "CV_temp_metrics.RData"))
CV_temp_metrics
summary(CV_temp_metrics)

# plots: ----

## maps comparing mean CV predictions and observations: ----

obs_preds_sf <- y_preds_obs_df %>%
  left_join(occ_preds_obs_df[c("RTENO", "Year", "pred_occ_mean")]) %>%
  left_join(routes_sel_sf, by = c(RTENO = "RTENO_BBS")) %>%
  st_as_sf()

# the more similar the colours of observation and prediction dots, the better:

# for single years:
# for predicted y:
obs_preds_sf %>% 
  filter(Year == 2019) %>% 
  ggplot() +
  geom_sf(aes(fill = as.factor(presence)), pch = 21, size = 3) +
  scale_fill_viridis_d(name = "observation") +
  geom_sf(aes(color = pred_y_mean), size = 1) +
  scale_color_viridis_c(name = "mean y. prediction") +
  theme_bw() +
  ggtitle(spec)

# for predicted occ.:
obs_preds_sf %>% 
  filter(Year == 2019) %>% 
  ggplot() +
  geom_sf(aes(fill = as.factor(presence)), pch = 21, size = 3) +
  scale_fill_viridis_d(name = "observation") +
  geom_sf(aes(color = pred_occ_mean), size = 1) +
  scale_color_viridis_c(name = "mean occ. prediction") +
  theme_bw() +
  ggtitle(spec)

# more years:
obs_preds_sf %>%
  filter(Year %in% seq(1995, 2019, length = 4)) %>%
ggplot() +
  facet_wrap(~Year) +
  geom_sf(aes(fill = as.factor(presence)), pch = 21, size = 3) +
  scale_fill_viridis_d(name = "observation") +
  geom_sf(aes(color = pred_occ_mean), size = 1) +
  scale_color_viridis_c(name = "mean occ. prediction") +
  theme_bw() +
  ggtitle(spec)
