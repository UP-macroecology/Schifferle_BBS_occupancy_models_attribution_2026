# model evaluation: temporal predictive performance of DOMs

# model fitted to 15 years of training data and evaluated based on the following 10 years


# packages: --------------------------------------------------------------------

library(sf)
library(dplyr)
library(ggplot2)
library(gridExtra)


# directories: -----------------------------------------------------------------

# project directory:
dir <- file.path("//NAS-2-P-SN-01.ibb.uni-potsdam.de", "daten$", "AG26", "Transfer", "Schifferle_BBS_occupancy_models_2023")

results_dir <- file.path(dir, "results", "temp_val_buffer_750_10yrs")

# directory to save time series calculated from predictions:

intermediate_time_series_dir <- file.path(results_dir, "temp_eval", "10_years", "pred_occ_sum")
if(!dir.exists(intermediate_time_series_dir)){dir.create(intermediate_time_series_dir)}


# functions: -------------------------------------------------------------------

source("0_functions.R")


# load data: -------------------------------------------------------------------

# selected species:
load(file = file.path("data", "species_set_analysis.RData")) # output of 3_1_eval_DOMs_CV.R
final_species

# route-year-species information (only surveyed)
load(file = file.path("data", "BBS_for_occ_spec_records.RData")) # bbs_dt_occ; output of 1_0_dataprep_BBS_bird_data.R

# routes-years:
load(file = file.path("data", "BBS_for_occ_selection.RData")) # route_sel_dt; output of 1_3_dataprep_match_BBS_routes_env_data.R

# selected routes spatial data (to buffer presences):
routes_sel_sf <- st_read(file.path("data", "route_selection_1995_2019_surv_beg_end_max_5y_miss_v2_spat_thin_100km_max_30_r_per_BCR_centroids.shp")) # output of 1_1_dataprep_BBS_route_selection.R

n_train_years <- 15 


# extract time series of summed observations and summed y predictions: ---------

for(i in 1:length(final_species)){

  spec <- final_species[i]
  print(paste(i, spec))

  # check where to look for model output (did MCMC fitting work with less or only with more iterations?)
  if(file.exists(file.path(results_dir, "refit_2000_2000", paste0("out_", spec, "_temp_val_10yrs_buffer_750.RData")))){
    output_dir <- file.path(results_dir, "refit_2000_2000")
  } else {
    output_dir <- results_dir
  }

  print(output_dir)

  # observations:

  # relevant routes for the species, within distance of 750 km of presences:
  rel_routes <- training_routes(species = spec, buffer_km = 750, output = "RTENOs")
  occ_dt_spec <- BBS_pres_abs_spec(species = spec) %>%
    filter(RTENO %in% rel_routes)

  # route-level presence:
  # sum all routes for each year (temporal trend)
  obs_temp_trend <- occ_dt_spec %>%
    group_by(Year) %>%
    summarise(pres_sum = sum(presence, na.rm = TRUE))

  # model predictions:

  load(file.path(output_dir, paste0("preds_", spec, "_temp_val_10yrs_buffer_750.RData")))

  # sum across route sections:
  dim(res_list$y_preds) # routes - sections - years - draws
  preds_routes <- apply(res_list$y_preds, MAR = c(1,3,4), FUN = max) # as soon as one predicted presence -> 1
  # sum across routes for each year:
  preds_years <- apply(preds_routes, MAR = c(2,3), FUN = sum, na.rm = TRUE)
  # mean:
  preds_years_mean <- apply(preds_years, MAR = 1, FUN = mean)
  # sd:
  preds_years_sd <- apply(preds_years, MAR = 1, FUN = sd)
  # 95% range:
  preds_years_0025_0975 <- apply(preds_years, MAR = 1, FUN = function(x) quantile(x, probs = c(0.025, 0.975)))

  obs_temp_trend$preds_years_mean <- preds_years_mean
  obs_temp_trend$preds_years_sd <- preds_years_sd
  obs_temp_trend$preds_years_0025 <- preds_years_0025_0975[1,]
  obs_temp_trend$preds_years_0975 <- preds_years_0025_0975[2,]

  # plot observations, predictions:

  # ggplot(obs_temp_trend, aes(x = Year)) +
  #   geom_line(aes(y = pres_sum)) +
  #   geom_point(aes(y = pres_sum), size = 3) +
  #   geom_line(aes(y = preds_years_mean), color = "cornflowerblue") +
  #   geom_point(aes(y = preds_years_mean), color = "cornflowerblue", size = 3) +
  #   ylab("N routes with presence") +
  #   theme_bw() +
  #   theme(text = element_text(size = 20)) +
  #   geom_vline(xintercept = years[1]+(n_train_years-1), linetype = "dashed") +
  #   ggtitle(spec) +
  #   geom_smooth(aes(y = pres_sum), method = "loess", color = "black")+
  #   geom_smooth(aes(y = preds_years_mean), method = "loess", color = "cornflowerblue")

  # save:
  save(obs_temp_trend, file = file.path(intermediate_time_series_dir,
                                        paste0(spec, "_occ_sum_series.RData")))
  }


# validation metrics: ----------------------------------------------------------

temp_val_metrics <- data.frame("species" = final_species, 
                               "mape" = NA, # mean absolute percentage error
                               "obs_in_CI95" = NA, # observation within 95 % credible interval of predictions
                               "cor_p" = NA, # correlation predictions - observations
                               "mae" = NA, # mean absolute error
                               "trend_mae" = NA, # trend in absolute error (slope of lm mae ~ Year)
                               "p_trend_mae" = NA, # p value of trend
                               "p_cor_p_neg" = NA # p value Pearson correlation significantly < 0
                               )


# load intermediate file observations - predictions:

for(i in 1:nrow(temp_val_metrics)){
  
  spec <- final_species[i]
  
  print(paste(i, spec))
  
  load(file = file.path(results_dir, "temp_eval", "10_years", "pred_occ_sum", paste0(spec, "_occ_sum_series.RData"))) 
  
  # test data:
  obs_temp_trend_test <- obs_temp_trend[16:25,] %>% # only test years
    mutate(aerr = abs(preds_years_mean - pres_sum)) # absolute error for each year
           
  # mean absolute percentage error:
  temp_val_metrics$mape[i] <- Metrics::mape(actual = obs_temp_trend_test$pres_sum, 
                                            predicted = obs_temp_trend_test$preds_years_mean)
    
  # observation within 95 % credible interval of predictions for every year:
  temp_val_metrics$obs_in_CI95[i] <- all(obs_temp_trend_test$pres_sum < obs_temp_trend_test$preds_years_0975 & 
                                           obs_temp_trend_test$pres_sum > obs_temp_trend_test$preds_years_0025)
  
  # Pearson correlation:
  temp_val_metrics$cor_p[i] <- cor(obs_temp_trend_test$pres_sum, obs_temp_trend_test$preds_years_mean, method = "p")
  
  # mean absolute error:
  temp_val_metrics$mae[i] <- mean(obs_temp_trend_test$aerr)
  
  # trend in absolute error:
  lm_error <- lm(aerr ~ c(1:10), data = obs_temp_trend_test) # years 1:10
  temp_val_metrics$trend_mae[i] <- lm_error$coefficients[2] # slope of lm
  temp_val_metrics$p_trend_mae[i] <- summary(lm_error)$coefficients[2,4] # slope significance
  
  # significant negative correlation test:
  temp_val_metrics$p_cor_p_neg[i] <- cor.test(obs_temp_trend_test$pres_sum, obs_temp_trend_test$preds_years_mean, method = "p",
           alternative = "less")$p.value 
  
  print(temp_val_metrics[i,])
}

#save(temp_val_metrics, file = file.path(results_dir, "temp_eval", "10_years", "temp_val_metrics.RData"))


# filter species based on validation criteria: ---------------------------------

specs_thresh <- temp_val_metrics %>% 
  # either error is comparatively small or overall trend is captured:
  # mean absolute error below 10 % or observations in 95 % prediction credible interval or correlation between obs and preds > 0.5:
  filter((mape < 0.1 | obs_in_CI95) | (cor_p > 0.5)) %>% 
  # and no large deviations:
  # no sign. positive trend in mean absolute error and correlation is not significantly negative 
  filter(!(trend_mae > 0 & p_trend_mae <= 0.05) & !(p_cor_p_neg <= 0.05)) %>% 
  pull(species) # 81 species left

#save(specs_thresh, file = file.path(results_dir, "temp_eval", "10_years", "spec_set_temp_val_ok.RData"))


# # (plot time series with metrics:) -----
# 
# # predicted y against observations summed over years:
# 
# buffer_km <- 750
# 
# # directory to save plots:
# if(!dir.exists(file.path(results_dir, "temp_eval", "10_years", "time_series_aggr_across_routes"))){
#   dir.create(file.path(results_dir, "temp_eval", "10_years", "time_series_aggr_across_routes"), recursive = TRUE)
# }
# 
# for(i in 1:length(final_species)){
#   
#   spec <- final_species[i]
#   print(paste(i, spec)) 
#   
#   load(file = file.path(results_dir, "temp_eval", "10_years", "pred_occ_sum", paste0(spec, "_occ_sum_series.RData")))
#   
#   metrics_table <- temp_val_metrics %>% 
#     filter(species == spec) %>% 
#     tidyr::pivot_longer(cols = !species, 
#                         names_to = "metric",
#                         values_to = "value") %>% 
#     select(-species) %>% 
#     mutate(value = round(value, 2))
# 
# 
#   jpeg(file = file.path(results_dir, "temp_eval", "10_years", "time_series_aggr_across_routes", paste0("temp_val_10yrs_", spec,"_", buffer_km, "km.jpg")), 
#        width = 1000, height = 700, quality = 100)
#   
#   print(
#     ggplot(obs_temp_trend, aes(x = Year)) +
#       geom_line(aes(y = pres_sum)) +
#       geom_point(aes(y = pres_sum), size = 3) +
#       geom_ribbon(aes(y = preds_years_mean, ymax = preds_years_0975, ymin = preds_years_0025), 
#                   alpha = 0.2, fill = "cornflowerblue") +
#       geom_line(aes(y = preds_years_mean), color = "cornflowerblue") +
#       geom_point(aes(y = preds_years_mean), color = "cornflowerblue", size = 3) +
#       ylab("N routes with presence") +
#       theme_bw() +
#       theme(text = element_text(size = 20)) +
#       geom_vline(xintercept = 2009.5, linetype = "dashed") +
#       ggtitle(spec) +
#       expand_limits(x = 2026) + 
#       # add metrics table:
#       annotation_custom(tableGrob(metrics_table,
#                                   rows = NULL,
#                                   theme = ttheme_minimal(core=list(fg_params=list(hjust=0, x=0.1)),
#                                                          colhead=list(fg_params=list(hjust=0, x=0.1)),
#                                                          base_size = 18)
#                                   ),
#                         xmin=2021, xmax=2025)
#   )
#   dev.off()
# }
