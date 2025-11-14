# quantify temporal predictive performance of DOMs with different metrics
# based on model fitted to 15 years of training data and predictions for following 10 years


# packages: --------------------------------------------------------------------

library(sf)
library(dplyr)
library(ggplot2)
library(gridExtra)

# directories: ----

results_dir <- file.path("//NAS-2-P-SN-01.ibb.uni-potsdam.de/users$/schifferle1", "Documents", "DEBTs", "analysis", "Schifferle_BBS_occupancy_models_2023",
                         "results", "temp_val_buffer_750_10yrs")

# directory to save time series calculated from predictions:

intermediate_time_series_dir <- file.path(results_dir, "temp_eval", "10_years", "pred_occ_sum")
if(!dir.exists(intermediate_time_series_dir)){dir.create(intermediate_time_series_dir)}


# functions: ----

source("0_functions.R")


# load data: ----

# selected species:
load(file = file.path("data", "species_set_analysis.RData"))
final_species

# selected routes spatial data (to buffer presences):
routes_sel_sf <- st_read(file.path("data", "route_selection_1995_2019_surv_beg_end_max_5y_miss_v2_spat_thin_100km_max_30_r_per_BCR_centroids.shp")) # output of 1_1_route_selection.R

# routes-years:
load(file = file.path("data", "BBS_for_occ_selection.RData")) # route_sel_dt; output of 1_3_match_BBS_to_env_data.R 

# route-year-species information (only surveyed)
load(file = file.path("data", "BBS_for_occ_spec_records.RData")) # bbs_dt_occ; output of 1_0_reformat_BBS_data.R

n_train_years <- 15
years <- min(route_sel_dt$Year):max(route_sel_dt$Year) 


# save time series of summed observations and mean summed y predictions: ----

# for(i in 1:length(final_species)){
#   
#   spec <- final_species[i]
#   print(paste(i, spec))
# 
#   # check where to look for model output (did MCMC fitting work with less or only with more iterations?)
#   if(file.exists(file.path(results_dir, "refit_2000_2000", paste0("out_", spec, "_temp_val_10yrs_buffer_750.RData")))){
#     output_dir <- file.path(results_dir, "refit_2000_2000")
#   } else {
#     output_dir <- results_dir
#   }
# 
#   print(output_dir)
# 
#   # observations: ---
# 
#   # relevant routes for the species, within distance of 750 km of presences:
#   rel_routes <- training_routes(species = spec, buffer_km = 750, output = "RTENOs")
#   occ_dt_spec <- BBS_pres_abs_spec(species = spec) %>%
#     filter(RTENO %in% rel_routes)
#   
#   # route-level presence:
#   # sum all routes for each year (temporal trend)
#   obs_temp_trend <- occ_dt_spec %>%
#     group_by(Year) %>%
#     summarise(pres_sum = sum(presence, na.rm = TRUE))
#   
#   # model predictions: ---
#   
#   load(file.path(output_dir, paste0("preds_", spec, "_temp_val_10yrs_buffer_750.RData")))
#   
#   # sum across route sections:
#   dim(res_list$y_preds) # routes - sections - years - draws
#   preds_routes <- apply(res_list$y_preds, MAR = c(1,3,4), FUN = max) # as soon as one predicted presence -> 1
#   # sum across routes for each year:
#   preds_years <- apply(preds_routes, MAR = c(2,3), FUN = sum, na.rm = TRUE) 
#   # mean:
#   preds_years_mean <- apply(preds_years, MAR = 1, FUN = mean)
#   # sd:
#   preds_years_sd <- apply(preds_years, MAR = 1, FUN = sd)
#   # 95% range: 
#   preds_years_0025_0975 <- apply(preds_years, MAR = 1, FUN = function(x) quantile(x, probs = c(0.025, 0.975)))
# 
#   obs_temp_trend$preds_years_mean <- preds_years_mean
#   obs_temp_trend$preds_years_sd <- preds_years_sd
#   obs_temp_trend$preds_years_0025 <- preds_years_0025_0975[1,]
#   obs_temp_trend$preds_years_0975 <- preds_years_0025_0975[2,]
# 
#   # plot observations, predictions: ---
#   
#   ggplot(obs_temp_trend, aes(x = Year)) +
#     geom_line(aes(y = pres_sum)) +
#     geom_point(aes(y = pres_sum), size = 3) +
#     geom_line(aes(y = preds_years_mean), color = "cornflowerblue") +
#     geom_point(aes(y = preds_years_mean), color = "cornflowerblue", size = 3) +
#     ylab("N routes with presence") +
#     theme_bw() +
#     theme(text = element_text(size = 20)) +
#     geom_vline(xintercept = years[1]+(n_train_years-1), linetype = "dashed") +
#     ggtitle(spec) +
#     geom_smooth(aes(y = pres_sum), method = "loess", color = "black")+
#     geom_smooth(aes(y = preds_years_mean), method = "loess", color = "cornflowerblue")#+
# 
#   # save:
#   save(obs_temp_trend, file = file.path(intermediate_time_series_dir,
#                                         paste0(spec, "_occ_sum_series.RData")))
# 
#   }



# metrics based on time series: ----

temp_val_metrics <- data.frame("species" = final_species, 
                               "mae" = NA, # mean absolute error
                               "trend_mae" = NA, # trend in absolute error (slope of lm mae ~ Year)
                               "int_mae" = NA, # intercept of lm mae ~ Year
                               "p_trend_mae" = NA, # p value of trend
                               "me" = NA, # mean error predictions minus observations
                               "trend_me" = NA,  # trend in error (slope of lm me ~ Year)
                               "int_me" = NA, # intercept of lm me ~ Year
                               "p_trend_me" = NA, # p value of trend
                               "mse" = NA, # mean squared error
                               "trend_mse" = NA,  # trend in squared error (slope of lm mse ~ Year)
                               "int_mse" = NA, # intercept of lm mse ~ Year
                               "p_trend_mse" = NA, # p value of trend
                               "rmse" = NA, # root mean squared error
                               "mape" = NA, # mean absolute percentage error
                               "mase" = NA, # mean absolute scaled error
                               "cor_p" = NA, # correlation predictions - observations
                               "p_cor_p_pos" = NA, # p value Pearson correlation significantly > 0
                               "p_cor_p_neg" = NA, # p value Pearson correlation significantly < 0
                               "cor_s" = NA, # rank correlation predictions - observations 
                               "p_cor_s_pos" = NA, # p value rank correlation significantly > 0
                               "p_cor_s_neg" = NA, # p value rank correlation significantly < 0
                               "C_temp_val" = NA, # Harrel's C index
                               "obs_in_CI95" = NA # is observation within 95 % credible interval (?) of predictions
                               )


# load intermediate file observations - predictions:

for(i in 1:nrow(temp_val_metrics)){
  
  spec <- final_species[i]
  
  print(paste(i, spec))
  
  load(file = file.path(results_dir, "temp_eval", "10_years", "pred_occ_sum", paste0(spec, "_occ_sum_series.RData"))) 
  
  # only test data:
  obs_temp_trend_test <- obs_temp_trend[16:25,] %>% 
    mutate(Year_Ind = 1:10, # index for years, to fit lms
           aerr = abs(preds_years_mean - pres_sum), # absolute error
           err_pred_obs = preds_years_mean - pres_sum, # error predictions minus observations
           se = (preds_years_mean - pres_sum)^2, # squared error
           aerr_scaled = aerr/ pres_sum,
           err_pred_obs_scaled = err_pred_obs / pres_sum)
           
    
  # mean absolute error:
  temp_val_metrics$mae[i] <- mean(obs_temp_trend_test$aerr)
  
  # # plot absolute error:
  # ggplot(obs_temp_trend_test, aes(x = Year)) +
  #   geom_point(aes(y = aerr)) +
  #   theme_bw() +
  #   theme(text = element_text(size = 20)) +
  #   ggtitle(spec) +
  #   geom_smooth(aes(x = Year, y = aerr), method = "lm", color = "black")
  
  
  # trend in absolute error:
  lm_error <- lm(aerr ~ Year_Ind, data = obs_temp_trend_test)
  temp_val_metrics$trend_mae[i] <- lm_error$coefficients[2] # slope of lm
  temp_val_metrics$int_mae[i] <- lm_error$coefficients[1] # intercept of lm
  temp_val_metrics$p_trend_mae[i] <- summary(lm_error)$coefficients[2,4] # slope significance
  
  # mean error:
  temp_val_metrics$me[i] <- mean(obs_temp_trend_test$err_pred_obs)
  
  # trend in error:
  lm_error <- lm(err_pred_obs ~ Year_Ind, data = obs_temp_trend_test)
  temp_val_metrics$trend_me[i] <- lm_error$coefficients[2] # slope of lm
  temp_val_metrics$int_me[i] <- lm_error$coefficients[1] # intercept of lm
  temp_val_metrics$p_trend_me[i] <- summary(lm_error)$coefficients[2,4] # slope significance
  
  # mean squared error:
  temp_val_metrics$mse[i] <- mean(obs_temp_trend_test$se)
  
  # trend in squared error:
  lm_error <- lm(se ~ Year_Ind, data = obs_temp_trend_test)
  temp_val_metrics$trend_mse[i] <- lm_error$coefficients[2] # slope of lm
  temp_val_metrics$int_mse[i] <- lm_error$coefficients[1] # intercept of lm
  temp_val_metrics$p_trend_mse[i] <- summary(lm_error)$coefficients[2,4] # slope significance
  
  # root mean squared error
  temp_val_metrics$rmse[i] <- sqrt(mean(obs_temp_trend_test$se))
  
  # mean absolute percent error:
  temp_val_metrics$mape[i] <- Metrics::mape(actual = obs_temp_trend_test$pres_sum, 
                                            predicted = obs_temp_trend_test$preds_years_mean)
  

  # mean absolute scaled error:
  temp_val_metrics$mase[i] <- Metrics::mase(actual = obs_temp_trend_test$pres_sum, 
                                            predicted = obs_temp_trend_test$preds_years_mean)
  
  # correlations: 
  
  # Pearson correlation:
  temp_val_metrics$cor_p[i] <- cor(obs_temp_trend_test$pres_sum, obs_temp_trend_test$preds_years_mean, method = "p")
  
  # tests:
  temp_val_metrics$p_cor_p_pos[i] <- cor.test(obs_temp_trend_test$pres_sum, obs_temp_trend_test$preds_years_mean, method = "p",
           alternative = "greater")$p.value # sign. positive correlation

  temp_val_metrics$p_cor_p_neg[i] <- cor.test(obs_temp_trend_test$pres_sum, obs_temp_trend_test$preds_years_mean, method = "p",
           alternative = "less")$p.value # sign. negative correlation

  # rank correlation:
  temp_val_metrics$cor_s[i] <- cor(obs_temp_trend_test$pres_sum, obs_temp_trend_test$preds_years_mean, method = "s")
  
  # tests:
  temp_val_metrics$p_cor_s_pos[i] <- cor.test(obs_temp_trend_test$pres_sum, obs_temp_trend_test$preds_years_mean, method = "s",
                                              alternative = "greater")$p.value # sign. negative correlation
  
  temp_val_metrics$p_cor_s_neg[i] <- cor.test(obs_temp_trend_test$pres_sum, obs_temp_trend_test$preds_years_mean, method = "s",
                                              alternative = "less")$p.value # sign. positive correlation
  
  # C-index of last 10 years:
  C_ind <- Hmisc::rcorr.cens(x = obs_temp_trend_test$preds_years_mean, 
                             S = obs_temp_trend_test$pres_sum)
  temp_val_metrics$C_temp_val[i] <- C_ind["C Index"]
  
  # is observation within 95 % credible interval (?) of predictions for every year:
  temp_val_metrics$obs_in_CI95[i] <- all(obs_temp_trend_test$pres_sum < obs_temp_trend_test$preds_years_0975 & 
                                              obs_temp_trend_test$pres_sum > obs_temp_trend_test$preds_years_0025)
  
  print(temp_val_metrics[i,])
}

#save(temp_val_metrics, file = file.path(results_dir, "temp_eval", "10_years", "temp_val_metrics_final.RData"))

# note: only positive trends in error important, a negative trend (= better predictions further into the future) is fine


# plot time series with metrics: -----

# predicted y against observations summed over years:

buffer_km <- 750

# directory to save plots:
if(!dir.exists(file.path(results_dir, "temp_eval", "10_years", "time_series_plots_sum_presence_across_routes"))){
  dir.create(file.path(results_dir, "temp_eval", "10_years", "time_series_plots_sum_presence_across_routes"), recursive = TRUE)
}

for(i in 1:length(final_species)){
  
  spec <- final_species[i]
  print(paste(i, spec)) 
  
  load(file = file.path(results_dir, "temp_eval", "10_years", "pred_occ_sum", paste0(spec, "_occ_sum_series.RData")))
  
  metrics_table <- temp_val_metrics %>% 
    filter(species == spec) %>% 
    tidyr::pivot_longer(cols = !species, 
                        names_to = "metric",
                        values_to = "value") %>% 
    select(-species) %>% 
    filter(metric %in% c("mae", "me", "mape", "mase",
                         "trend_mae",  "trend_me",
                         "p_trend_mae", "p_trend_me",
                         "cor_p", "p_cor_p_pos", "p_cor_p_neg", "C_temp_val")) %>% 
    mutate(value = round(value, 2)) %>% 
    arrange(factor(metric, levels = c("mae", "me", "mape", "mase",
                                      "trend_mae",  "trend_me",
                                      "p_trend_mae", "p_trend_me",
                                      "cor_p", "p_cor_p_pos", "p_cor_p_neg", "C_temp_val")))
  
  
  jpeg(file = file.path(results_dir, "temp_eval", "10_years", "time_series_plots_sum_presence_across_routes", paste0("temp_val_10yrs_", spec,"_", buffer_km, "km.jpg")), 
       width = 1000, height = 700, quality = 100)
  
  print(
    ggplot(obs_temp_trend, aes(x = Year)) +
      geom_line(aes(y = pres_sum)) +
      geom_point(aes(y = pres_sum), size = 3) +
       #geom_ribbon(aes(y = preds_years_mean, ymax = preds_years_mean + 1*preds_years_sd, ymin = preds_years_mean - 1*preds_years_sd), 
      #             alpha = 0.2, fill = "cornflowerblue") +
      geom_ribbon(aes(y = preds_years_mean, ymax = preds_years_0975, ymin = preds_years_0025), 
                  alpha = 0.2, fill = "cornflowerblue") +
      geom_line(aes(y = preds_years_mean), color = "cornflowerblue") +
      geom_point(aes(y = preds_years_mean), color = "cornflowerblue", size = 3) +
      #ylab("sum of mean obs. per section across routes") +
      ylab("N routes with presence") +
      theme_bw() +
      theme(text = element_text(size = 20)) +
      geom_vline(xintercept = 2009.5, linetype = "dashed") +
      ggtitle(spec) +
      expand_limits(x = 2026) + 
      # add metrics table:
      annotation_custom(tableGrob(metrics_table,
                                  rows = NULL,
                                  theme = ttheme_minimal(core=list(fg_params=list(hjust=0, x=0.1)),
                                                         colhead=list(fg_params=list(hjust=0, x=0.1)),
                                                         base_size = 18)
      ), 
      xmin=2021, xmax=2025)
  )
  dev.off()
}


# explorations: ----

load(file.path(results_dir, "temp_eval", "10_years", "temp_val_metrics_final.RData"))


# correlations among evaluation metrics:
M <- cor(temp_val_metrics %>% 
           select(-species) %>% 
           select(mae, me, mse, rmse, mape, mase, cor_p, cor_s, C_temp_val), 
         method = "p")
corrplot::corrplot(M, method = "square", order = "hclust",
                   addCoef.col = "black",
                   diag = FALSE,
                   tl.cex = 1,#1
                   number.cex = 0.8, # 0.8
                   number.digits= 2)

# highly correlated are:
# - cor_p and cor_s and C_temp_val
# - mae, mse, rmse, (mase)

plot(mae ~ mase, data = temp_val_metrics)
plot(mae ~ mape, data = temp_val_metrics)
plot(mase ~ C_temp_val, data = temp_val_metrics)
plot(mase ~ rmse, data = temp_val_metrics)
plot(mape ~ rmse, data = temp_val_metrics)
plot(mape ~ mse, data = temp_val_metrics)

# find suitable thresholds ----

# that keep species for which model works okay and flags species
# for which model doesn't predict okay in time:

summary(temp_val_metrics)

# proposal:

specs_thresh <- temp_val_metrics %>% 
  # either error is comparatively small or overall trend is captured:
  # mean absolute error below 10 % or observations in 95 % prediction credible interval or correlation between obs and preds > 0.5:
  filter((mape < 0.1 | obs_in_CI95) | (cor_p > 0.5)) %>% 
  # and no large deviations:
  # no sign. positive trend in mean absolute error and correlation is not significantly negative 
  filter(!(trend_mae > 0 & p_trend_mae <= 0.05) & !(p_cor_p_neg <= 0.05)) %>% 
  pull(species) # 81 species left

save(specs_thresh, file = file.path(results_dir, "temp_eval", "10_years", "spec_set_temp_val_ok1.RData"))
