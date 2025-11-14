# assemble metrics used to quantify the contribution of climate and land use change to estimated occupancy dynamics
# based on counterfactual and factual predictions

# packages: --------------------------------------------------------------------

library(dplyr)
library(ggplot2)

# directories: -----------------------------------------------------------------
main_dir <- file.path("//NAS-2-P-SN-01.ibb.uni-potsdam.de", "daten$", "AG26", "Transfer",
                      "Schifferle_BBS_occupancy_models_2023")


# load data: -------------------------------------------------------------------

# species for which fitting worked:
load(file = file.path("data", "species_set_analysis.RData")) # output of 3_1_DOM_CV_evaluation_metrics.R
final_species

# species for which models worked fine:

# okay in time:
load(file.path("//NAS-2-P-SN-01.ibb.uni-potsdam.de/users$/schifferle1", "Documents", "DEBTs", "analysis", "Schifferle_BBS_occupancy_models_2023",
               "results", "temp_val_buffer_750_10yrs", "temp_eval", "10_years", "spec_set_temp_val_ok1.RData")) # output of 3_1_DOM_temp_evaluation_metrics.R
spec_temp_okay <- specs_thresh
# okay in space:
CV_eval_summary <- read.csv(file = file.path("//NAS-2-P-SN-01.ibb.uni-potsdam.de/users$/schifferle1", "Documents", "DEBTs", "analysis", "Schifferle_BBS_occupancy_models_2023",
                                             "results", "CV_buffer750km", "CV_eval", "CV_eval_summary.csv")) # 3_1_DOM_CV_evaluation_metrics.R
spec_spat_okay <- CV_eval_summary %>% 
  filter(y_spat_auc_mean >= 0.7) %>% 
  pull(species)

# okay in both:
spec_okay <- intersect(spec_temp_okay, spec_spat_okay) # 80
#save(spec_okay, file = file.path("data", "species_DOM_val_okay.RData"))

# directories: -----------------------------------------------------------------

obs_dir <- file.path(main_dir, "data", "observed_time_series_1995_2019")
fact_dir <- file.path(main_dir, "results", "fm_buffer750km", "fact_pred_time_series_1995_2019")
cfact_dir <- file.path(main_dir, "results", "attribution", "cfact_pred_time_series_1995_2019") # output of 5_0_occupancy_dynamics_time_series.R

# assemble df: -----------------------------------------------------------------

# df to store results:
attr_metr_df <- tibble(species = final_species,
                      slope_obs = NA,
                      p_obs = NA,
                      slope_fact = NA,
                      slope_CIlow_fact = NA,
                      slope_CIhigh_fact = NA,
                      p_fact = NA,
                      slope_cfclim = NA,
                      slope_CIlow_cfclim = NA,
                      slope_CIhigh_cfclim = NA,
                      p_cfclim = NA,
                      slope_cflu = NA,
                      slope_CIlow_cflu = NA,
                      slope_CIhigh_cflu = NA,
                      p_cflu = NA,
                      slope_cfclimlu = NA,
                      slope_CIlow_cfclimlu = NA,
                      slope_CIhigh_cfclimlu = NA,
                      p_cfclimlu = NA,
                      mape_fact = NA,
                      mape_cfclim = NA,
                      mape_cflu = NA,
                      mape_cfclimlu = NA,
                      mae_fact = NA,
                      mae_cfclim = NA,
                      mae_cflu = NA,
                      mae_cfclimlu = NA)
# species:

for(s in 1:nrow(attr_metr_df)){
  
  spec <- final_species[s]
  
  print(paste(s, spec))
  
  # observations time series:
  load(file.path(obs_dir, paste0(spec, "_obs_ts_sum_occ_routes.RData")))
  ts_obs
  
  # time series predictions factual:
  load(file.path(fact_dir, paste0(spec, "_ts_sum_occ_routes_f_preds.RData")))
  ts_preds_fact
  
  # time series predictions counterfactual:
  load(file.path(cfact_dir, paste0(spec, "_ts_sum_occ_routes_cf_preds.RData")))
  ts_preds_cfact
  
  # determine overall linear trend:
  
  # ## observed: not used further?
  # trend_obs <- summary(lm(Npres ~ year, data = ts_obs))
  # attr_metr_df$slope_obs[s] <- trend_obs$coefficients["year", "Estimate"]
  # attr_metr_df$p_obs[s] <- trend_obs$coefficients["year", "Pr(>|t|)"]
  
  
  ## predictions:
  # 100 posterior draws per year;
  # fit linear model without intercept
  # scale data by subtracting median of first year
  
  # data:
  
  ## factual:
  lm_df_f <- ts_preds_fact %>% 
    select(year, starts_with("draw")) %>% 
    tidyr::pivot_longer(starts_with("draw"), names_to = "draw", values_to = "value") %>% 
    mutate(value_scaled = value - ts_preds_fact$median_Nocc[1],
           year_scaled = year - 1995,
           scenario = "fact")
  
  ## climate counterfactual:
  lm_df_clim <- ts_preds_cfact$cf_clim %>% 
    select(year, starts_with("draw")) %>% 
    tidyr::pivot_longer(starts_with("draw"), names_to = "draw", values_to = "value") %>% 
    mutate(value_scaled = value - ts_preds_cfact$cf_clim$median_Nocc[1],
           year_scaled = year - 1995,
           scenario = "clim")
  
  ## land use counterfactual:
  lm_df_lu <- ts_preds_cfact$cf_1995soc %>% 
    select(year, starts_with("draw")) %>% 
    tidyr::pivot_longer(starts_with("draw"), names_to = "draw", values_to = "value") %>% 
    mutate(value_scaled = value - ts_preds_cfact$cf_1995soc$median_Nocc[1],
           year_scaled = year - 1995,
           scenario = "lu")
  
  ## climate + land use counterfactual:
  lm_df_climlu <- ts_preds_cfact$cf_clim_1995soc %>% 
    select(year, starts_with("draw")) %>% 
    tidyr::pivot_longer(starts_with("draw"), names_to = "draw", values_to = "value") %>% 
    mutate(value_scaled = value - ts_preds_cfact$cf_clim_1995soc$median_Nocc[1],
           year_scaled = year - 1995,
           scenario = "climlu")
  
  lm_df_all <- rbind(lm_df_f, lm_df_clim, lm_df_lu, lm_df_climlu) %>% 
    mutate(scenario = factor(scenario))
  
  # test plot:
  # xx improve for manuscript!
  plot_dir <- file.path("plots", "attribution", "explorations", "lm_posterior_draws_scaled_no_int") # xx
  if(!dir.exists(plot_dir)) dir.create(plot_dir, recursive = TRUE)
  jpeg(file = file.path(plot_dir, paste0(spec, "_lm_posterior_draws_opt3.jpg")), # xx
       width = 1200, height = 800, quality = 100)
  print(lm_df_all %>%
          mutate(scenario = factor(scenario, levels = c("fact", "clim", "lu", "climlu")),
                 scenario = recode(scenario, fact = "factual", clim = "counterfactual\nclimate",
                                   lu = "counterfactual\nland use", climlu = "counterfactual\nclimate + land use")) %>%
          ggplot(aes(x = year_scaled, y = value_scaled, colour = scenario, fill = scenario)) +
          geom_point(size = 1.5, position = position_jitter(width = 0.3), alpha = 0.2, shape = 19) +
          geom_smooth(method = "lm", formula = y ~ x + 0, se = TRUE) +
          labs(y = "scaled number of occupied routes",
               title = paste0(spec, ", 100 draws of posterior distribution"),
               x = "scaled year") +
          scale_colour_manual(values = c("counterfactual\nclimate + land use" = "#046865",
                                         "counterfactual\nland use" = "#B7410E",
                                         "counterfactual\nclimate" = "#0D98BA",
                                         "factual" = "#85CB33")) +
          scale_fill_manual(values = c("counterfactual\nclimate + land use" = "#046865",
                                       "counterfactual\nland use" = "#B7410E",
                                       "counterfactual\nclimate" = "#0D98BA",
                                       "factual" = "#85CB33")) +
          guides(fill = "none") +
          theme_bw() +
          theme(text = element_text(size = 20))
  )
  dev.off()
  
  # linear models to quantify trend:
  
  ## factual
  fm <- lm(value_scaled ~ 0 + year_scaled, data = lm_df_f)
  attr_metr_df$slope_fact[s] <- summary(fm)$coefficients["year_scaled", "Estimate"]
  attr_metr_df$p_fact[s] <- summary(fm)$coefficients["year_scaled", "Pr(>|t|)"]
  # confidence interval around slope
  CI_f <- confint(fm, "year_scaled", level = 0.95)
  attr_metr_df$slope_CIlow_fact[s] <- CI_f[1]
  attr_metr_df$slope_CIhigh_fact[s] <- CI_f[2]
  
  ## climate counterfactual:
  fm <- lm(value_scaled ~ 0 + year_scaled, data = lm_df_clim)
  attr_metr_df$slope_cfclim[s] <- summary(fm)$coefficients["year_scaled", "Estimate"]
  attr_metr_df$p_cfclim[s] <- summary(fm)$coefficients["year_scaled", "Pr(>|t|)"]
  # confidence interval around slope
  CI_f <- confint(fm, "year_scaled", level = 0.95)
  attr_metr_df$slope_CIlow_cfclim[s] <- CI_f[1]
  attr_metr_df$slope_CIhigh_cfclim[s] <- CI_f[2]
  
  ## land use counterfactual:
  fm <- lm(value_scaled ~ 0 + year_scaled, data = lm_df_lu)
  attr_metr_df$slope_cflu[s] <- summary(fm)$coefficients["year_scaled", "Estimate"]
  attr_metr_df$p_cflu[s] <- summary(fm)$coefficients["year_scaled", "Pr(>|t|)"]
  # confidence interval around slope
  CI_f <- confint(fm, "year_scaled", level = 0.95)
  attr_metr_df$slope_CIlow_cflu[s] <- CI_f[1]
  attr_metr_df$slope_CIhigh_cflu[s] <- CI_f[2]
  
  ## climate + land use counterfactual:
  fm <- lm(value_scaled ~ 0 + year_scaled, data = lm_df_climlu)
  attr_metr_df$slope_cfclimlu[s] <- summary(fm)$coefficients["year_scaled", "Estimate"]
  attr_metr_df$p_cfclimlu[s] <- summary(fm)$coefficients["year_scaled", "Pr(>|t|)"]
  # confidence interval around slope
  CI_f <- confint(fm, "year_scaled", level = 0.95)
  attr_metr_df$slope_CIlow_cfclimlu[s] <- CI_f[1]
  attr_metr_df$slope_CIhigh_cfclimlu[s] <- CI_f[2]
  
  # mean absolute percent error:
  # xx based on mean? or rather with uncertainty? 
  attr_metr_df$mape_fact[s] <- Metrics::mape(actual = ts_obs$Npres, 
                                            predicted = ts_preds_fact$median_Nocc)
  
  attr_metr_df$mape_cfclim[s] <- Metrics::mape(actual = ts_obs$Npres, 
                                              predicted = ts_preds_cfact$cf_clim$median_Nocc)
  
  attr_metr_df$mape_cflu[s] <- Metrics::mape(actual = ts_obs$Npres, 
                                            predicted = ts_preds_cfact$cf_1995soc$median_Nocc)
  
  attr_metr_df$mape_cfclimlu[s] <- Metrics::mape(actual = ts_obs$Npres, 
                                                predicted = ts_preds_cfact$cf_clim_1995soc$median_Nocc)
  
  # mean absolute error:
  attr_metr_df$mae_fact[s] <- Metrics::mae(actual = ts_obs$Npres, 
                                          predicted = ts_preds_fact$median_Nocc)
  
  attr_metr_df$mae_cfclim[s] <- Metrics::mae(actual = ts_obs$Npres, 
                                            predicted = ts_preds_cfact$cf_clim$median_Nocc)
  
  attr_metr_df$mae_cflu[s] <- Metrics::mae(actual = ts_obs$Npres, 
                                          predicted = ts_preds_cfact$cf_1995soc$median_Nocc)
  
  attr_metr_df$mae_cfclimlu[s] <- Metrics::mae(actual = ts_obs$Npres, 
                                              predicted = ts_preds_cfact$cf_clim_1995soc$median_Nocc)
  
  
}

#save(attr_metr_df, file = file.path(main_dir, "results", "attribution", "attribution_metrics_final.RData"))

