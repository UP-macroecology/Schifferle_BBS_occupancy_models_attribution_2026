# Script:   5_1_attribution_metrics.R
# Purpose:  Assemble metrics to compare occupancy time series for different scenarios to quantify the contribution of climate and land use change
# Inputs:   results/species_DOM_val_okay.RData
#           data/observed_time_series_1995_2019/<species>_obs_ts_sum_occ_routes.RData
#           results/fm_buffer750km/fact_pred_time_series_1995_2019/<species>_ts_sum_occ_routes_f_preds.RData (one file per species)
#           results/attribution/cfact_pred_time_series_1995_2019/<species>_ts_sum_occ_routes_cf_preds.RData (one file per species)
# Outputs:  results/attribution/attribution_metrics_final.RData
#           plots/attribution/trend_calculation_example.svg (Fig. S2)
# Runs on:  Local


source(file.path("scripts", "0_paths.R"))


# packages: --------------------------------------------------------------------

library(dplyr)
library(ggplot2)


# directories: -----------------------------------------------------------------

# observed time series:
obs_dir <- file.path(dir, "data", "observed_time_series_1995_2019") # output of 4_1_DOMs_predictions_time_series.R
# predicted time series for factual data:
fact_dir <- file.path(dir, "results", "fm_buffer750km", "fact_pred_time_series_1995_2019") # output of 4_1_DOMs_predictions_time_series.R
# predicted time series for counterfactual data:
cfact_dir <- file.path(dir, "results", "attribution", "cfact_pred_time_series_1995_2019") # output of 4_1_DOMs_predictions_time_series.R


# load data: -------------------------------------------------------------------

# species for attribution:
load(file = file.path(dir, "results", "species_DOM_val_okay.RData")) # output of 4_0_DOMs_predictions_y_routes_scenarios.R
spec_attr


# assemble df: -----------------------------------------------------------------

attr_metr_df <- tibble(species = spec_attr,
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
                      mape_cfclimlu = NA)

# iterate over species:

for(s in 1:nrow(attr_metr_df)){
  
  spec <- spec_attr[s]
  
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
  
  # # plot to check:
  # print(lm_df_all %>%
  #         mutate(scenario = factor(scenario, levels = c("fact", "clim", "lu", "climlu")),
  #                scenario = recode(scenario, fact = "factual", clim = "counterfactual\nclimate",
  #                                  lu = "counterfactual\nland use", climlu = "counterfactual\nclimate + land use")) %>%
  #         ggplot(aes(x = year_scaled, y = value_scaled, colour = scenario, fill = scenario)) +
  #         geom_point(size = 1.5, position = position_jitter(width = 0.3), alpha = 0.2, shape = 19) +
  #         geom_smooth(method = "lm", formula = y ~ x + 0, se = TRUE) +
  #         labs(y = "scaled number of occupied routes",
  #              title = spec,
  #              x = "scaled year") +
  #         scale_colour_manual(values = c("counterfactual\nclimate + land use" = "#046865",
  #                                        "counterfactual\nland use" = "#B7410E",
  #                                        "counterfactual\nclimate" = "#0D98BA",
  #                                        "factual" = "#85CB33")) +
  #         scale_fill_manual(values = c("counterfactual\nclimate + land use" = "#046865",
  #                                      "counterfactual\nland use" = "#B7410E",
  #                                      "counterfactual\nclimate" = "#0D98BA",
  #                                      "factual" = "#85CB33")) +
  #         guides(fill = "none") +
  #         theme_bw() +
  #         theme(text = element_text(size = 20))
  # )

  
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
  attr_metr_df$mape_fact[s] <- Metrics::mape(actual = ts_obs$Npres, 
                                            predicted = ts_preds_fact$median_Nocc)
  
  attr_metr_df$mape_cfclim[s] <- Metrics::mape(actual = ts_obs$Npres, 
                                              predicted = ts_preds_cfact$cf_clim$median_Nocc)
  
  attr_metr_df$mape_cflu[s] <- Metrics::mape(actual = ts_obs$Npres, 
                                            predicted = ts_preds_cfact$cf_1995soc$median_Nocc)
  
  attr_metr_df$mape_cfclimlu[s] <- Metrics::mape(actual = ts_obs$Npres, 
                                                predicted = ts_preds_cfact$cf_clim_1995soc$median_Nocc)
}

save(attr_metr_df, file = file.path(dir, "results", "attribution", "attribution_metrics_final.RData"))



# figure for manuscript to illustrate trend calculation (Fig. S2): -----------------------

# directory to save plot:
plot_dir <- file.path(dir, "plots", "attribution")

# example species:
spec <- "Dickcissel"

# observations time series:
load(file.path(obs_dir, paste0(spec, "_obs_ts_sum_occ_routes.RData")))
# time series predictions factual:
load(file.path(fact_dir, paste0(spec, "_ts_sum_occ_routes_f_preds.RData")))
# time series predictions counterfactual:
load(file.path(cfact_dir, paste0(spec, "_ts_sum_occ_routes_cf_preds.RData")))

# determine overall linear trend:

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

# plot:
trend_calc_plot <- lm_df_all %>%
        mutate(scenario = factor(scenario, levels = c("fact", "clim", "lu", "climlu")),
               scenario = recode(scenario, fact = "factual", clim = "counterfactual\nclimate",
                                 lu = "counterfactual\nland use", climlu = "counterfactual\nclimate + land use")) %>%
        ggplot(aes(x = year_scaled, y = value_scaled, colour = scenario, fill = scenario)) +
        geom_point(size = 1.5, position = position_jitter(width = 0.3), alpha = 0.2, shape = 19) +
        geom_smooth(method = "lm", formula = y ~ x + 0, se = TRUE) +
        labs(y = "scaled number of occupied routes",
             x = "scaled year") +
        scale_colour_manual(values = c("counterfactual\nclimate + land use" = "#046865",
                                       "counterfactual\nland use" = "#B7410E",
                                       "counterfactual\nclimate" = "#0D98BA",
                                       "factual" = "#85CB33")) +
        scale_fill_manual(values = c("counterfactual\nclimate + land use" = "#046865",
                                     "counterfactual\nland use" = "#B7410E",
                                     "counterfactual\nclimate" = "#0D98BA",
                                     "factual" = "#85CB33")) +
        guides(fill = guide_legend(title = "Scenario"),
               colour = guide_legend(title = "Scenario")) +
        theme_bw() +
        theme(text = element_text(size = 20),
              legend.key.height =  unit(1, units = "cm"),
              legend.key.width =  unit(1, units = "cm"),
              legend.key.spacing.y = unit(0.5, units = "cm"))
trend_calc_plot
ggsave(filename = file.path(plot_dir, "trend_calculation_example.svg"), 
       plot = trend_calc_plot,
       device = "svg",
       width = 29.7,
       height = 18, # A4
       units = "cm")

# session info:
writeLines(capture.output(sessionInfo()), file.path(dir, "results", "sessionInfo", "5_1_attribution_metrics.txt"))
