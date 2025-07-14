# quantify the contribution of factual climate and land use change to estimated occupancy dynamics
# based on counterfactual and factual predictions?

# packages: --------------------------------------------------------------------

library(dplyr)
library(sf)
library(ggplot2)
library(ggalluvial) # alluvial plots
library(ggnewscale) # to have multiple colour fill scales

# functions: -------------------------------------------------------------------

source("0_functions.R")

# directories: ----
#main_dir <- file.path("T:", "Schifferle_BBS_occupancy_models_2023")
main_dir <- file.path("//NAS-2-P-SN-01.ibb.uni-potsdam.de", "daten$", "AG26", "Transfer",
                      "Schifferle_BBS_occupancy_models_2023")
 

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

# for traits:
load(file.path("data", "BBS_data_merged.RData")) # bbs_dt

# species for which models worked fine:

# okay in time:
load(file.path("//NAS-2-P-SN-01.ibb.uni-potsdam.de/users$/schifferle1", "Documents", "DEBTs", "analysis", "Schifferle_BBS_occupancy_models_2023",
               "results", "temp_val_buffer_750_10yrs", "temp_eval", "10_years", "spec_set_temp_val_ok1.RData")) # 3_2_quant_temp_performance.R
spec_temp_okay <- specs_thresh
# okay in space:
CV_eval_summary <- read.csv(file = file.path("//NAS-2-P-SN-01.ibb.uni-potsdam.de/users$/schifferle1", "Documents", "DEBTs", "analysis", "Schifferle_BBS_occupancy_models_2023",
                                             "results", "CV_buffer750km", "CV_eval", "CV_eval_summary.csv")) # 3_1_DOM_CV_evaluation_metrics.R
spec_spat_okay <- CV_eval_summary %>% 
  filter(occ_spat_auc_mean >= 0.7) %>% 
  pull(species)
# okay in both:
spec_okay <- intersect(spec_temp_okay, spec_spat_okay) # 80


# directories: -----------------------------------------------------------------

obs_dir <- file.path(main_dir, "data", "observed_time_series_1995_2019")
fact_dir <- file.path(main_dir, "results", "fm_buffer750km", "fact_pred_time_series_1995_2019")
cfact_dir <- file.path(main_dir, "results", "attribution", "cfact_pred_time_series_1995_2019") # output of 5_0_occupancy_dynamics_time_series.R


# assemble df: -----------------------------------------------------------------

# df to store results:
rel_infl_df <- tibble(species = final_species,
                      slope_obs = NA,
                      p_obs = NA,
                      slope_fact = NA,
                      p_fact = NA,
                      slope_cfclim = NA,
                      p_cfclim = NA,
                      slope_cflu = NA,
                      p_cflu = NA,
                      slope_cfclimlu = NA,
                      p_cfclimlu = NA,
                      opt2_slope_fact = NA,
                      opt2_p_fact = NA,
                      opt2_slope_cfclim = NA,
                      opt2_p_cfclim = NA,
                      opt2_slope_cflu = NA,
                      opt2_p_cflu = NA,
                      opt2_slope_cfclimlu = NA,
                      opt2_p_cfclimlu = NA,
                      mape_fact = NA,
                      mape_cfclim = NA,
                      mape_cflu = NA,
                      mape_cfclimlu = NA,
                      mae_fact = NA,
                      mae_cfclim = NA,
                      mae_cflu = NA,
                      mae_cfclimlu = NA)
# species:

for(s in 1:nrow(rel_infl_df)){
  
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
  
  ## option 1: linear trend based on median number of occupied routes per year: ----
  
  ## observed:
  trend_obs <- summary(lm(Npres ~ year, data = ts_obs))
  rel_infl_df$slope_obs[s] <- trend_obs$coefficients["year", "Estimate"]
  rel_infl_df$p_obs[s] <- trend_obs$coefficients["year", "Pr(>|t|)"]
  ## factual:
  trend_fact <- summary(lm(median_Nocc_f ~ year, data = ts_preds_fact))
  rel_infl_df$slope_fact[s] <- trend_fact$coefficients["year", "Estimate"]
  rel_infl_df$p_fact[s] <- trend_fact$coefficients["year", "Pr(>|t|)"]
  ## climate counterfactual:
  trend_cfclim <- summary(lm(median_Nocc_cf_clim ~ year, data = ts_preds_cfact$cf_clim))
  rel_infl_df$slope_cfclim[s] <- trend_cfclim$coefficients["year", "Estimate"]
  rel_infl_df$p_cfclim[s] <- trend_cfclim$coefficients["year", "Pr(>|t|)"]
  ## land use counterfactual:
  trend_cflu <- summary(lm(median_Nocc_cf_1995soc ~ year, data = ts_preds_cfact$cf_1995soc))
  rel_infl_df$slope_cflu[s] <- trend_cflu$coefficients["year", "Estimate"]
  rel_infl_df$p_cflu[s] <- trend_cflu$coefficients["year", "Pr(>|t|)"]
  ## climate + land use counterfactual:
  trend_cfclimlu <- summary(lm(median_Nocc_cf_clim_1995soc ~ year, data = ts_preds_cfact$cf_clim_1995soc))
  rel_infl_df$slope_cfclimlu[s] <- trend_cfclimlu$coefficients["year", "Estimate"]
  rel_infl_df$p_cfclimlu[s] <- trend_cfclimlu$coefficients["year", "Pr(>|t|)"]
  
  # option 2: linear trend based on 100 draws of posterior for number of occupied routes per year: ----
  
  ## factual:
  lm_df <- ts_preds_fact %>% 
    select(year, starts_with("draw")) %>% 
    tidyr::pivot_longer(starts_with("draw"), names_to = "draw", values_to = "value")
  
  trend <- summary(lm(value ~ year, data = lm_df))
  
  rel_infl_df$opt2_slope_fact[s] <- trend$coefficients["year", "Estimate"]
  rel_infl_df$opt2_p_fact[s] <- trend$coefficients["year", "Pr(>|t|)"]
  
  ## climate counterfactual:
  lm_df <- ts_preds_cfact$cf_clim %>% 
    select(year, starts_with("draw")) %>% 
    tidyr::pivot_longer(starts_with("draw"), names_to = "draw", values_to = "value")

  trend <- summary(lm(value ~ year, data = lm_df))
  
  rel_infl_df$opt2_slope_cfclim[s] <- trend$coefficients["year", "Estimate"]
  rel_infl_df$opt2_p_cfclim[s] <- trend$coefficients["year", "Pr(>|t|)"]
  
  ## land use counterfactual:
  lm_df <- ts_preds_cfact$cf_1995soc %>% 
    select(year, starts_with("draw")) %>% 
    tidyr::pivot_longer(starts_with("draw"), names_to = "draw", values_to = "value")
  
  trend <- summary(lm(value ~ year, data = lm_df))
  
  rel_infl_df$opt2_slope_cflu[s] <- trend$coefficients["year", "Estimate"]
  rel_infl_df$opt2_p_cflu[s] <- trend$coefficients["year", "Pr(>|t|)"]
  
  ## climate + land use counterfactual:
  lm_df <- ts_preds_cfact$cf_clim_1995soc %>% 
    select(year, starts_with("draw")) %>% 
    tidyr::pivot_longer(starts_with("draw"), names_to = "draw", values_to = "value")
  
  trend <- summary(lm(value ~ year, data = lm_df))
  
  rel_infl_df$opt2_slope_cfclimlu[s] <- trend$coefficients["year", "Estimate"]
  rel_infl_df$opt2_p_cfclimlu[s] <- trend$coefficients["year", "Pr(>|t|)"]

  
  # mean absolute percent error:
  rel_infl_df$mape_fact[s] <- Metrics::mape(actual = ts_obs$Npres, 
                                            predicted = ts_preds_fact$median_Nocc_f)
  
  rel_infl_df$mape_cfclim[s] <- Metrics::mape(actual = ts_obs$Npres, 
                                              predicted = ts_preds_cfact$cf_clim$median_Nocc_cf_clim)
  
  rel_infl_df$mape_cflu[s] <- Metrics::mape(actual = ts_obs$Npres, 
                                            predicted = ts_preds_cfact$cf_1995soc$median_Nocc_cf_1995soc)
  
  rel_infl_df$mape_cfclimlu[s] <- Metrics::mape(actual = ts_obs$Npres, 
                                                predicted = ts_preds_cfact$cf_clim_1995soc$median_Nocc_cf_clim_1995soc)
  
  # mean absolute error:
  rel_infl_df$mae_fact[s] <- Metrics::mae(actual = ts_obs$Npres, 
                                               predicted = ts_preds_fact$median_Nocc_f)
  
  rel_infl_df$mae_cfclim[s] <- Metrics::mae(actual = ts_obs$Npres, 
                                                 predicted = ts_preds_cfact$cf_clim$median_Nocc_cf_clim)
  
  rel_infl_df$mae_cflu[s] <- Metrics::mae(actual = ts_obs$Npres, 
                                               predicted = ts_preds_cfact$cf_1995soc$median_Nocc_cf_1995soc)
  
  rel_infl_df$mae_cfclimlu[s] <- Metrics::mae(actual = ts_obs$Npres, 
                                                   predicted = ts_preds_cfact$cf_clim_1995soc$median_Nocc_cf_clim_1995soc)
  
  
}

#save(rel_infl_df, file = file.path(main_dir, "results", "attribution", "rel_influence_010725.RData"))
load(file = file.path(main_dir, "results", "attribution", "rel_influence_010725.RData"))


# occ. dyn. trend for scenarios across species: --------------------------------


# for how many species does occupancy increase, decrease, remain stable for each scenario:
# (based on slope and p-value)

sum_dyn_scen <- rel_infl_df %>% 
  #filter(species %in% spec_okay) %>% 
  select(c(species, starts_with("slope"), starts_with("p_"))) %>% 
  tidyr::pivot_longer(cols = matches("(slope_)|(p_)"), names_to = c("metric", "scenario"), 
                      values_to = "value", names_pattern = "(.*)_(.*)") %>% 
  tidyr::pivot_wider(names_from = metric, values_from = value) %>% 
  mutate(dynamics = ifelse(slope > 0 & p < 0.05, "increase",
                           ifelse(slope < 0 & p < 0.05, "decrease", "stable"))) %>% 
  group_by(scenario, dynamics) %>% 
  summarise(n = n()) %>% 
  arrange(scenario, dynamics)

## lollipop plots: ----

sum_dyn_scen %>%
  mutate(scenario = factor(scenario, levels = c("cfclimlu", "cflu", "cfclim", "fact"))) %>% 
  filter(scenario != "obs") %>% 
  ggplot(aes(x = n, y = dynamics)) +
  geom_errorbarh(aes(xmin = 0, xmax = n, colour = scenario),
                 height = 0,
                 position = position_dodge(width = .7)) +
  geom_point(aes(colour = scenario), size=4, position = position_dodge(width = .7)) +
  theme_bw() +
  xlab("N species") +
  ylab("") +
  scale_colour_manual(values = c("#046865", "#B7410E", "#0D98BA", "#85CB33"),
                      labels = c("counterfactual climate + land use",
                                 "counterfactual land use", "counterfactual climate", "factual"))+
  guides(colour = guide_legend(reverse = TRUE)) +
  ggtitle("Occupancy dynamics 1995 - 2019") +
  theme(text = element_text(size = 15))

# plot other way around:
sum_dyn_scen %>% 
  mutate(scenario = factor(scenario, levels = c("cfclimlu", "cflu", "cfclim", "fact"))) %>% 
  mutate(scenario = recode(scenario, fact = "factual", cfclim = "counterfactual\nclimate",
                           cflu = "counterfactual\nland use", cfclimlu = "counterfactual\nclimate + land use")) %>% 
  filter(scenario != "obs") %>% 
  mutate(n_perc = n/length(final_species)*100) %>% 
  ggplot(aes(x = scenario, y = n_perc)) +
  geom_errorbar(aes(ymin = 0, ymax = n_perc, colour = dynamics),
                position = position_dodge(width = .6), width = 0, linewidth = 1) +
  geom_point(aes(colour = dynamics), size=4, position = position_dodge(width = .6)) +
  ylab("N species [%]") +
  xlab("") +
  ggtitle("Occupancy dynamics 1995 - 2019") +
  theme_bw() +
  theme(text = element_text(size = 15), legend.title=element_blank()) +
  coord_flip()

# plot number of selected routes over time (do they increase?) ---
occ_dt_spec <- BBS_pres_abs_spec(species = final_species[1])
routes_year <- occ_dt_spec %>% 
  filter(Surveyed == 1) %>% 
  select(c(RTENO, Year)) %>%
  filter(!Year %in% c(1995, 2019)) %>% # remove first and last year for which all routes were surveyed
  group_by(Year) %>% 
  summarise(n_routes = n()) 
ggplot(data = routes_year, aes(x = Year, y = n_routes)) +
  geom_point() +
  geom_smooth(method = "lm") +
  theme_bw() +
  ylab("Number of surveyed routes") +
  theme(text = element_text(size = 18)) +
  labs(title = paste(" slope =", signif(lm(n_routes ~ Year, data = routes_year)$coef[[2]], 2),
                     " p =", signif(summary(lm(n_routes ~ Year, data = routes_year))$coef[2,4], 2)))
# significant increase (when not considering first and last year) ---


## alluvial plots: ----

# https://r-charts.com/flow/ggalluvial/
# change categories based on linear trend

### option 1: linear trend based on median of occupied routes: ----

flow_df <- rel_infl_df %>% 
  #filter(species %in% spec_okay) %>% 
  select(c(species, starts_with("slope"), starts_with("p_"))) %>% 
  tidyr::pivot_longer(cols = matches("(slope_)|(p_)"), names_to = c("metric", "scenario"), 
                      values_to = "value", names_pattern = "(.*)_(.*)") %>% 
  tidyr::pivot_wider(names_from = metric, values_from = value) %>% 
  mutate(dynamics = ifelse(slope > 0 & p < 0.05, "increase",
                           ifelse(slope < 0 & p < 0.05, "decrease", "stable"))) %>% 
  mutate(dynamics = factor(dynamics, levels = c("increase", "stable", "decrease"))) %>% 
  select(-c(slope, p)) %>% 
  tidyr::pivot_wider(names_from = scenario, values_from = dynamics)

# alluvial plot change in AOO factual vs. counterfactual climate:
jpeg(file = file.path("plots", "attribution", "occ_dyn_flow_climate.jpg"), 
     width = 1000, height = 800, quality = 100)
flow_df %>%
  filter(species %in% spec_okay) %>% 
  group_by(fact, cfclim) %>% 
  summarise(n = n()) %>% 
  ggplot(aes(axis1 = cfclim, axis2 = fact, y = n)) +
  geom_alluvium(aes(fill = cfclim), show.legend = FALSE, width = 1/5, colour = "grey50") +
  geom_stratum(width = 1/5, aes(fill = cfclim), show.legend = FALSE) + 
  scale_fill_manual(values = c("increase" = "#EFCA08", "decrease" =  "#A17C6B", "stable" =  "#93C6D6"), na.value = NA) +
  geom_stratum(width = 1/5, aes(fill = fact), show.legend = FALSE) +
  scale_fill_manual(values = c("increase" = "#EFCA08", "decrease" =  "#A17C6B", "stable" =  "#93C6D6"), na.value = NA) +
  geom_text(stat = "stratum", aes(label = after_stat(stratum)), size = 8) +
  geom_text(stat = "flow",
            #aes(label = n),
            aes(label = ifelse(test = after_stat(x) == "1",
                               yes = flow_df %>% 
                                 filter(species %in% spec_okay) %>% 
                                 group_by(fact, cfclim) %>% summarise(n = n()) %>% arrange(cfclim) %>% pull(n) %>% rev,
                               no = "")),
            size = 8, nudge_x = 0.15) +
  scale_x_discrete(limits = c("cfclim", "fact"),
                   expand = c(0.15, 0.05),
                   labels = c("counterfactual\nclimate", "factual\nclimate")) +
  xlab("") +
  #scale_fill_viridis_d(direction = -1) +
  ylab("Number of species") +
  theme_bw() +
  ggtitle("Change in area of occupancy, considering only slope + significance") +
  theme(axis.ticks.x = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text.y = element_blank(),
        axis.text.x = element_blank(),
        text = element_text(size = 24),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank()) +
  annotate("text", x = 1, y = 85, label = "Counterfactual\nclimate", hjust = 0.5, size = 8) + # y = 85 for spec_okay
  annotate("text", x = 2, y = 85, label = "Factual\nclimate", hjust = 0.5, size = 8) 
dev.off()

# same alluvial plot for land use factual vs. counterfactual:
jpeg(file = file.path("plots", "attribution", "occ_dyn_flow_landuse_spec_okay.jpg"), 
     width = 1000, height = 800, quality = 100)
flow_df %>%
  filter(species %in% spec_okay) %>% 
  group_by(fact, cflu) %>% 
  summarise(n = n()) %>% 
  ggplot(aes(axis1 = cflu, axis2 = fact, y = n)) +
  geom_alluvium(aes(fill = cflu), show.legend = FALSE, width = 1/5, colour = "grey50") +
  geom_stratum(width = 1/5, aes(fill = cflu), show.legend = FALSE) + 
  scale_fill_manual(values = c("increase" = "#EFCA08", "decrease" =  "#A17C6B", "stable" =  "#93C6D6"), na.value = NA) +
  geom_stratum(width = 1/5, aes(fill = fact), show.legend = FALSE) +
  scale_fill_manual(values = c("increase" = "#EFCA08", "decrease" =  "#A17C6B", "stable" =  "#93C6D6"), na.value = NA) +
  geom_text(stat = "stratum", aes(label = after_stat(stratum)), size = 8) +
  geom_text(stat = "flow",
            aes(label = ifelse(test = after_stat(x) == "1",
                               yes = flow_df %>% 
                                 filter(species %in% spec_okay) %>% 
                                 group_by(fact, cflu) %>% summarise(n = n()) %>% arrange(cflu) %>% pull(n) %>% rev,
                               no = "")),
            size = 8, nudge_x = 0.15) +
  scale_x_discrete(limits = c("cflu", "fact"), expand = c(0.15, 0.05)) +
  xlab("") +
  #scale_fill_viridis_d(direction = -1) +
  ylab("Number of species") +
  theme_bw() +
  ggtitle("Change in area of occupancy") +
  theme(axis.ticks.x = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text.y = element_blank(),
        axis.text.x = element_blank(),
        text = element_text(size = 24),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank()) +
  annotate("text", x = 1, y = 85, label = "Counterfactual\nland use", hjust = 0.5, size = 8) +
  annotate("text", x = 2, y = 85, label = "Factual\nland use", hjust = 0.5, size = 8) 
dev.off()

### option 2: linear trend based on draws of posterior: ----


flow_df <- rel_infl_df %>% 
  select(c(species, starts_with("opt2_slope"), starts_with("opt2_p_"))) %>% 
  tidyr::pivot_longer(cols = matches("(slope_)|(p_)"), names_to = c("metric", "scenario"), 
                      values_to = "value", names_pattern = "(.*)_(.*)") %>% 
  tidyr::pivot_wider(names_from = metric, values_from = value) %>% 
  mutate(dynamics = ifelse(opt2_slope > 0 & opt2_p < 0.05, "increase",
                           ifelse(opt2_slope < 0 & opt2_p < 0.05, "decrease", "stable"))) %>% 
  mutate(dynamics = factor(dynamics, levels = c("increase", "stable", "decrease"))) %>% 
  select(-c(opt2_slope, opt2_p)) %>% 
  tidyr::pivot_wider(names_from = scenario, values_from = dynamics)

# alluvial plot change in AOO factual vs. counterfactual climate:

plot_df <- flow_df %>%
  #filter(species %in% spec_okay) %>% 
  group_by(fact, cfclim) %>% 
  summarise(n = n()) %>% 
  mutate(rel_change = case_when(
    cfclim == "increase" & fact == "stable" ~ "relative decrease",
    cfclim == "increase" & fact == "decrease" ~ "absolute decrease",
    cfclim == "stable" & fact == "increase" ~ "absolute increase",
    cfclim == "stable" & fact == "decrease" ~ "absolute decrease",
    cfclim == "decrease" & fact == "increase" ~ "absolute increase",
    cfclim == "decrease" & fact == "stable" ~ "relative increase",
    .default = "no change")) %>% 
  mutate(rel_change = factor(rel_change, levels = c("absolute increase", "relative increase", "no change", "relative decrease", "absolute decrease"))) %>% 
  arrange(forcats::fct_rev(rel_change)) %>% # change plotting order
  arrange(cfclim, fact)

jpeg(file = file.path("plots", "attribution", "occ_dyn_flow_climate_linear_trend_incl_uncertainty.jpg"),
    width = 900, height = 1200, quality = 100)

plot_df %>% 
  ggplot(aes(axis1 = cfclim, axis2 = fact, y = n)) +
  geom_alluvium(aes(fill = rel_change), show.legend = TRUE, width = 1/5, colour = "grey50") +
  #scale_fill_brewer(palette = "RdBu", drop = FALSE, name = "relative change") +
  # scale_fill_manual(values = c("relative decrease" = "#E4572E", 
  #                              "absolute decrease" =  "#714955", 
  #                              "absolute increase" = "#53B3CB", #"#17BEBB",# "#00b2ca", #"#006DAA", #,
  #                              "relative increase" = "#7dcfb6",
  #                              "no change" = "gray80"), drop = FALSE, na.value = NA, name = "relative change") +
  scale_fill_manual(values = c("relative decrease" = "#abd9e9",
                               "absolute decrease" =  "#2c7bb6",
                               "absolute increase" = "#d7191c", #"#17BEBB",# "#00b2ca", #"#006DAA", #,
                               "relative increase" = "#fdae61",
                               "no change" = "gray90"), drop = FALSE, na.value = NA, name = "relative change") +
  
  new_scale_fill() +
  geom_stratum(width = 1/5, aes(fill = cfclim), show.legend = FALSE, colour = "grey20") + 
  scale_fill_manual(values = c("increase" = "#EFCA08", "decrease" =  "#7DC4C9", "stable" =  "#EDE6F2"), na.value = NA) +
  #scale_fill_manual(values = c("increase" = "#EFCA08", "decrease" =  "#A17C6B", "stable" =  "#A5C4D4"), na.value = NA) +
  geom_stratum(width = 1/5, aes(fill = fact), show.legend = FALSE, colour = "grey20") +
  geom_text(stat = "stratum", 
            fontface = "italic",
            aes(label = ifelse(test = after_stat(x) == "1",
                               yes = paste(stratum, "\n",
                                 flow_df %>%
                                 #filter(species %in% spec_okay) %>% 
                                 group_by(cfclim) %>% 
                                 summarise(n = n()) %>% 
                                 pull(n) %>% 
                                 rev),
                               no = paste(stratum, "\n",
                                 flow_df %>%
                                 #filter(species %in% spec_okay) %>% 
                                 group_by(fact) %>% 
                                 summarise(n = n()) %>% 
                                 pull(n) %>% 
                                 rev))),
            size = 8) +
  geom_text(stat = "flow", # species numbers in flow
            fontface = "italic",
            aes(label = ifelse(test = after_stat(x) == "1",
                               yes = plot_df %>% pull(n) %>% rev,
                               no = "")),
            size = 7, nudge_x = 0.15) +
  scale_x_discrete(limits = c("cfclim", "fact"), expand = c(0.13, 0.0)) +
  xlab("") + ylab("") +
  theme_bw() +
  ggtitle("Change in area of occupancy,\nsign. linear trend in 100 posterior draws per year") +
  labs(title = "Change in area of occupancy",
       subtitle = "linear trend in 100 posterior draws of N. occ. routes per year\nnumbers = number of species") +
  theme(axis.ticks.x = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text.y = element_blank(),
        axis.text.x = element_blank(),
        text = element_text(size = 28),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        legend.text = element_text(size = 26),
        panel.border = element_blank(),
        plot.subtitle = element_text(margin=margin(0,0,30,0)),
        legend.position = "bottom",
        legend.direction="vertical",
        legend.box.spacing = unit(-50, "pt")) +
  annotate("text", x = 1, y = 170, label = "Counterfactual\nclimate", hjust = 0.5, size = 9) + # y = 85 for spec_okay
  annotate("text", x = 2, y = 170, label = "Factual\nclimate", hjust = 0.5, size = 9) +
  guides(fill = "none")
dev.off()

#BE6E46
#9C7866

# # https://matthewdharris.com/2017/11/11/a-brief-diversion-into-static-alluvial-sankey-diagrams-in-r/
# 
# library(ggforce)
# plot_df %>% 
#   select(cfclim, fact, n, rel_change) %>% 
#   mutate(cfclim = as.character(cfclim),
#           fact = as.character(fact)
#           ) %>% 
#   gather_set_data(1:2) %>% 
#   arrange(x, fact, desc(cfclim)) %>% # xx
# ggplot(aes(x = x, id = id, split = y, value = n)) +
#   geom_parallel_sets(aes(fill = rel_change),
#                          alpha = 0.5, axis.width = 0.2,
#                      n=100, strength = 0.5) +
#   geom_parallel_sets_axes(axis.width = 0.25, fill = "gray95", color = "gray80", size = 0.15) +
#   geom_parallel_sets_labels(colour = 'gray35', size = 4.5, angle = 0, fontface="bold") +
#   #scale_fill_manual(values  = c(A_col, B_col, C_col)) +
#   #scale_color_manual(values = c(A_col, B_col, C_col)) +
#   scale_x_continuous(breaks = 1:2, labels = c("counterfactual\nclimate", "factual climate")) +
#   theme_minimal() +
#   theme(
#     legend.position = "none",
#     panel.grid.major = element_blank(),
#     panel.grid.minor = element_blank(),
#     axis.text.y = element_blank(),
#     axis.text.x = element_text(size = 20, face = "bold"),
#     axis.title.x  = element_blank()
#   )


# same for land use:
plot_df <- flow_df %>%
  #filter(species %in% spec_okay) %>% 
  group_by(fact, cflu) %>% 
  summarise(n = n()) %>% 
  mutate(rel_change = case_when(
    cflu == "increase" & fact == "stable" ~ "relative decrease",
    cflu == "increase" & fact == "decrease" ~ "absolute decrease",
    cflu == "stable" & fact == "increase" ~ "absolute increase",
    cflu == "stable" & fact == "decrease" ~ "absolute decrease",
    cflu == "decrease" & fact == "increase" ~ "absolute increase",
    cflu == "decrease" & fact == "stable" ~ "relative increase",
    .default = "no change")) %>% 
  mutate(rel_change = factor(rel_change, levels = c("absolute increase", "relative increase", "no change", "relative decrease", "absolute decrease"))) %>% 
  arrange(forcats::fct_rev(rel_change)) %>% # change plotting order
  arrange(cflu, fact)

jpeg(file = file.path("plots", "attribution", "occ_dyn_flow_landuse_linear_trend_incl_uncertainty.jpg"),
     width = 900, height = 1200, quality = 100)

plot_df %>% 
  ggplot(aes(axis1 = cflu, axis2 = fact, y = n)) +
  geom_alluvium(aes(fill = rel_change), show.legend = TRUE, width = 1/5, colour = "grey50") +
  scale_fill_manual(values = c("relative decrease" = "#abd9e9",
                               "absolute decrease" =  "#2c7bb6",
                               "absolute increase" = "#d7191c", #"#17BEBB",# "#00b2ca", #"#006DAA", #,
                               "relative increase" = "#fdae61",
                               "no change" = "gray90"), drop = FALSE, na.value = NA, name = "relative change") +
  new_scale_fill() +
  geom_stratum(width = 1/5, aes(fill = cflu), show.legend = FALSE, colour = "grey20") + 
  scale_fill_manual(values = c("increase" = "#EFCA08", "decrease" =  "#7DC4C9", "stable" =  "#EDE6F2"), na.value = NA) +
  geom_stratum(width = 1/5, aes(fill = fact), show.legend = FALSE, colour = "grey20") +
  geom_text(stat = "stratum", 
            fontface = "italic",
            aes(label = ifelse(test = after_stat(x) == "1",
                               yes = paste(stratum, "\n",
                                           flow_df %>%
                                             #filter(species %in% spec_okay) %>% 
                                             group_by(cflu) %>% 
                                             summarise(n = n()) %>% 
                                             pull(n) %>% 
                                             rev),
                               no = paste(stratum, "\n",
                                          flow_df %>%
                                            #filter(species %in% spec_okay) %>% 
                                            group_by(fact) %>% 
                                            summarise(n = n()) %>% 
                                            pull(n) %>% 
                                            rev))),
            size = 8) +
  geom_text(stat = "flow", # species numbers in flow
            fontface = "italic",
            aes(label = ifelse(test = after_stat(x) == "1",
                               yes = plot_df %>% pull(n) %>% rev,
                               no = "")),
            size = 7, nudge_x = 0.15) +
  scale_x_discrete(limits = c("cflu", "fact"), expand = c(0.13, 0.0)) +
  xlab("") + ylab("") +
  theme_bw() +
  ggtitle("Change in area of occupancy,\nsign. linear trend in 100 posterior draws per year") +
  labs(title = "Change in area of occupancy",
       subtitle = "linear trend in 100 posterior draws of N. occ. routes per year\nnumbers = number of species") +
  theme(axis.ticks.x = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text.y = element_blank(),
        axis.text.x = element_blank(),
        text = element_text(size = 28),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        legend.text = element_text(size = 26),
        panel.border = element_blank(),
        plot.subtitle = element_text(margin=margin(0,0,30,0)),
        legend.position = "bottom",
        legend.direction="vertical",
        legend.box.spacing = unit(-50, "pt")) +
  annotate("text", x = 1, y = 170, label = "Counterfactual\nland use", hjust = 0.5, size = 9) + # y = 85 for spec_okay
  annotate("text", x = 2, y = 170, label = "Factual\nland use", hjust = 0.5, size = 9) +
  guides(fill = "none")
dev.off()




## occupancy trend and traits: -------------------------------------------------

# occupancy trend:
occ_dyn_df <- rel_infl_df %>% 
  #filter(species %in% spec_okay) %>% 
  select(c(species, starts_with("slope"), starts_with("p_"))) %>% 
  tidyr::pivot_longer(cols = matches("(slope_)|(p_)"), names_to = c("metric", "scenario"), 
                      values_to = "value", names_pattern = "(.*)_(.*)") %>% 
  tidyr::pivot_wider(names_from = metric, values_from = value) %>% 
  mutate(dynamics = ifelse(slope > 0 & p < 0.05, "increase",
                           ifelse(slope < 0 & p < 0.05, "decrease", "stable"))) %>% 
  mutate(scenario = factor(scenario, levels = c("obs", "fact", "cfclim", "cflu", "cfclimlu"))) %>% 
  mutate(scenario = recode(scenario, 
                           obs = "observed",
                           fact = "factual predictions", 
                           cfclim = "counterfactual climate",
                           cflu = "counterfactual land use",
                           cfclimlu = "counterfactual climate & land use")) 
# add traits:
occ_dyn_traits <- bbs_dt %>% 
  select(English_Common_Name:Range.Size) %>% 
  distinct() %>% 
  right_join(occ_dyn_df, by = c(English_Common_Name = "species")) %>% 
  rename("species" = English_Common_Name)

plot_dir <- file.path("plots", "attribution", "explorations")
if(!dir.exists(plot_dir)) dir.create(plot_dir, recursive = TRUE)

### lollipop plots: ----
category <- c("ORDER", "Habitat", "Migration", "Trophic.Level",
              "Trophic.Niche", "Primary.Lifestyle")

for(c in 1:length(category)){

  print(category[c])
  
  jpeg(file = file.path(plot_dir, paste0("occ_trend_", category[c], ".jpg")), 
  width = 1600, height = 1200, quality = 100)
  
  print(occ_dyn_traits %>% 
    filter(scenario != "observed") %>% 
    mutate(dynamics = as.factor(dynamics)) %>% 
    group_by(.data[[category[c]]], dynamics, scenario) %>% 
    summarise(n = n()) %>% 
    ggplot() +
    facet_wrap(~scenario) +
    geom_errorbar(aes(xmin = 0, xmax = n, y = reorder(.data[[category[c]]], n), colour = dynamics),
                  position = position_dodge(width = .6), width = 0, linewidth = 1.1) +
    geom_point(aes(y = .data[[category[c]]], x = n, colour = dynamics), size = 6, position = position_dodge(width = .6)) +
    ggtitle(paste("Occupancy trend &", category[c])) +
    xlab("Number of species") +
    ylab("") +
    theme_bw() +
    guides(colour = guide_legend(reverse = TRUE)) +
    theme(text = element_text(size = 35), legend.title = element_blank())
  )
  
  dev.off()
  
  }


# other way around: dynamics in facets:
for(c in 1:length(category)){
  
  print(category[c])
  
  jpeg(file = file.path(plot_dir, paste0("occ_trend_", category[c], "2.jpg")), 
       width = 1600, height = 1200, quality = 100)
  
  # print(occ_dyn_traits %>% 
  #         filter(scenario != "observed") %>% 
  #         mutate(dynamics = as.factor(dynamics)) %>% 
  #         group_by(.data[[category[c]]], dynamics, scenario) %>% 
  #         summarise(n = n()) %>% 
  #         ggplot() +
  #         facet_wrap(~dynamics) +
  #         geom_errorbar(aes(xmin = 0, xmax = n, y = reorder(.data[[category[c]]], n), colour = scenario),
  #                       position = position_dodge(width = .6), width = 0, linewidth = 1.1) +
  #         geom_point(aes(y = .data[[category[c]]], x = n, colour = scenario), size = 6, position = position_dodge(width = .6)) +
  #         ggtitle(paste("Occupancy trend &", category[c])) +
  #         xlab("Number of species") +
  #         ylab("") +
  #         theme_bw() +
  #         scale_colour_manual(values = c("#046865", "#B7410E", "#0D98BA", "#85CB33"),
  #                             labels = c("counterfactual climate + land use",
  #                                        "counterfactual land use", "counterfactual climate", "factual")) +
  #         guides(colour = guide_legend(reverse = TRUE)) +
  #         theme(text = element_text(size = 18), legend.title = element_blank())
  # )
  
  print(occ_dyn_traits %>% 
          filter(scenario != "observed") %>% 
          mutate(dynamics = as.factor(dynamics)) %>% 
          group_by(.data[[category[c]]], dynamics, scenario) %>% 
          summarise(n = n()) %>% 
          ggplot() +
          facet_wrap(~reorder(.data[[category[c]]], -n)) +
          geom_errorbar(aes(xmin = 0, xmax = n, y = dynamics, colour = forcats::fct_rev(scenario)),
                        position = position_dodge(width = .6), width = 0, linewidth = 1.1) +
          geom_point(aes(y = dynamics, x = n, colour = forcats::fct_rev(scenario)), size = 6, position = position_dodge(width = .6)) +
          ggtitle(paste("Occupancy trend &", category[c])) +
          xlab("Number of species") +
          ylab("") +
          theme_bw() +
          scale_colour_manual(values = c("#046865", "#B7410E","#0D98BA", "#85CB33")) +
          guides(colour = guide_legend(reverse = TRUE)) +
          theme(text = element_text(size = 35), legend.title = element_blank())
  )
  dev.off()
  
}


### alluvial plots: ----

category <- c("Habitat", "Migration", "Trophic.Level",
              "Trophic.Niche", "Primary.Lifestyle")

# merge occ. dyn. and traits:
dyn_traits_df <- flow_df %>% # from above!
  left_join(bbs_dt %>% select(English_Common_Name:Range.Size), by = c(species = "English_Common_Name")) %>% 
  distinct()

for(c in 1:length(category)){
  
  print(category[c])
  
  filter_cat <- dyn_traits_df %>% 
    #filter(species %in% spec_okay) %>% 
    group_by(.data[[category[c]]]) %>% 
    summarise(n = n()) %>% 
    filter(n >= 5) %>% 
    pull(.data[[category[c]]])
  
  # climate
  trait_df <- dyn_traits_df %>% 
    filter(species %in% spec_okay) %>% 
    group_by(fact, cfclim, .data[[category[c]]]) %>% 
    summarise(n = n())
  
  jpeg(file = file.path(plot_dir, paste0("occ_dyn_flow_climate_", category[c], "_prelim_spec_okay.jpg")), 
       width = 1400, height = 1000, quality = 100)
  
  print(ggplot(data = trait_df %>% filter(.data[[category[c]]] %in% filter_cat), 
         aes(axis1 = cfclim, axis2 = fact, y = n, label = n)) +
    geom_alluvium(aes(fill = cfclim), show.legend = FALSE, width = 1/5) +
    geom_stratum(width = 1/5) +
    geom_text(stat = "stratum", aes(label = after_stat(stratum)), size = 6) +
    scale_x_discrete(limits = c("cfclim", "fact"),
                     expand = c(0.15, 0.05),
                     labels = c("Counterfactual\nclimate", "Factual\nclimate"),
                     #position = "top"
    ) +
    facet_wrap(~.data[[category[c]]], scales = 'free_y', axes = "all_x") +
    geom_text(stat = "flow", size = 6, nudge_x = 0.15) +
    xlab("") +
    scale_fill_viridis_d(direction = -1) +
    ylab("Number of species") +
    theme_bw() +
    ggtitle("Change in area of occupancy") +
    theme(axis.ticks.x = element_blank(),
          axis.ticks.y = element_blank(),
          axis.text.y = element_blank(),
          text = element_text(size = 22),
          panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank())
  )
  dev.off()
  
  
  # land use:
  trait_df <- dyn_traits_df %>% 
    filter(species %in% spec_okay) %>% 
    group_by(fact, cflu, .data[[category[c]]]) %>% 
    summarise(n = n())
  
  jpeg(file = file.path(plot_dir, paste0("occ_dyn_flow_landuse_", category[c], "_prelim_spec_okay.jpg")), 
       width = 1400, height = 1000, quality = 100)
  
  print(ggplot(data = trait_df %>% filter(.data[[category[c]]] %in% filter_cat),
         aes(axis1 = cflu, axis2 = fact, y = n, label = n)) +
    geom_alluvium(aes(fill = cflu), show.legend = FALSE, width = 1/5) +
    geom_stratum(width = 1/5) +
    geom_text(stat = "stratum", aes(label = after_stat(stratum)), size = 6) +
    scale_x_discrete(limits = c("cflu", "fact"),
                     expand = c(0.15, 0.05),
                     labels = c("Counterfactual\nland use", "Factual\nland use"),
                     #position = "top"
    ) +
    facet_wrap(~.data[[category[c]]], scales = 'free_y', axes = "all_x") +
    geom_text(stat = "flow", size = 6, nudge_x = 0.15) +
    xlab("") +
    scale_fill_viridis_d(direction = -1) +
    ylab("Number of species") +
    theme_bw() +
    ggtitle("Change in area of occupancy") +
    theme(axis.ticks.x = element_blank(),
          axis.ticks.y = element_blank(),
          axis.text.y = element_blank(),
          text = element_text(size = 22),
          panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank())
  )
  dev.off()
  
  }



# body mass: rather for relative categories? xx

# # to add number of cases above boxplots:
# get_box_stats <- function(y, upper_limit = 1.15) {
#   return(data.frame(
#     y = 0.95 * upper_limit,
#     label = length(y)
#   ))
# }
# 

# jpeg(file = file.path(plot_dir, "occ_dyn_body_mass.jpg"),
#      width = 1200, height = 900, quality = 100)
# occ_dyn_traits %>% 
#   select(species, scenario, dynamics, Mass) %>% 
#   filter(scenario != "observed") %>% 
#   mutate(scenario = factor(scenario, levels = c("counterfactual climate & land use",
#                                                   "counterfactual land use", "counterfactual climate",
#                                                   "factual predictions"))) %>%
#   mutate(dynamics = factor(dynamics, levels = c("increase", "stable", "decrease"))) %>% 
#   ggplot(aes(y = log(Mass), x = dynamics, group = dynamics:scenario, colour = scenario)) +
#     geom_boxplot() +
#     ggtitle("Occupancy trend and body mass") +
#     ylab("log(body mass) [g]") +
#     xlab("") +
#     #stat_summary(fun.data = get_box_stats, geom = "text", hjust = 0.5, vjust = 0.9, size = 8) +
#   guides(colour = guide_legend(reverse = TRUE)) +  
#   scale_colour_manual(values = c("#046865", "#B7410E", "#0D98BA", "#85CB33"),
#                       labels = c("counterfactual climate + land use",
#                                  "counterfactual land use", "counterfactual climate", "factual")) +
#   theme_bw() +
#   theme(text = element_text(size = 15))
# dev.off()

# alternative: compare number of occupied routes at beginning and end of time series:  ----
# summarise first and last three years, respectively:

rel_inf_change_df <- tibble(species = final_species,
                            Nocc_f_mean_first_3yrs = NA,
                            Nocc_f_90CIlow_first_3yrs = NA,
                            Nocc_f_90CIhigh_first_3yrs = NA,
                            Nocc_f_mean_last_3yrs = NA,
                            Nocc_f_90CIlow_last_3yrs = NA,
                            Nocc_f_90CIhigh_last_3yrs = NA,
                            Nocc_cfclim_mean_first_3yrs = NA,
                            Nocc_cfclim_90CIlow_first_3yrs = NA,
                            Nocc_cfclim_90CIhigh_first_3yrs = NA,
                            Nocc_cfclim_mean_last_3yrs = NA,
                            Nocc_cfclim_90CIlow_last_3yrs = NA,
                            Nocc_cfclim_90CIhigh_last_3yrs = NA,
                            Nocc_cflu_mean_first_3yrs = NA,
                            Nocc_cflu_90CIlow_first_3yrs = NA,
                            Nocc_cflu_90CIhigh_first_3yrs = NA,
                            Nocc_cflu_mean_last_3yrs = NA,
                            Nocc_cflu_90CIlow_last_3yrs = NA,
                            Nocc_cflu_90CIhigh_last_3yrs = NA,
                            Nocc_cfclimlu_mean_first_3yrs = NA,
                            Nocc_cfclimlu_90CIlow_first_3yrs = NA,
                            Nocc_cfclimlu_90CIhigh_first_3yrs = NA,
                            Nocc_cfclimlu_mean_last_3yrs = NA,
                            Nocc_cfclimlu_90CIlow_last_3yrs = NA,
                            Nocc_cfclimlu_90CIhigh_last_3yrs = NA)

for(s in 1:length(final_species)){
  
  spec <- final_species[s]
  
  print(paste(s, spec))
  
  # factual predictions:
  load(file.path(main_dir, "results", "fm_buffer750km", "y_preds_route_level_section_sum", 
                 paste0(spec, "_y_preds_route_level_section_sum.RData"))) # output of xx
  # sum across routes for each year:
  preds_years <- apply(preds_routes, MAR = c(2,3), FUN = sum, na.rm = TRUE) 
  
  # first 3 years:
  rel_inf_change_df$Nocc_f_mean_first_3yrs[s] <- mean(preds_years[1:3,])
  ts_ci90 <- bayestestR::ci(preds_years[1:3,], ci = 0.7, method = "ETI") # 0.9 XXXX
  rel_inf_change_df$Nocc_f_90CIlow_first_3yrs[s] <- ts_ci90$CI_low
  rel_inf_change_df$Nocc_f_90CIhigh_first_3yrs[s] <- ts_ci90$CI_high
  
  # last 3 years:
  rel_inf_change_df$Nocc_f_mean_last_3yrs[s] <- mean(preds_years[23:25,])
  ts_ci90 <- bayestestR::ci(preds_years[23:25,], ci = 0.7, method = "ETI") 
  rel_inf_change_df$Nocc_f_90CIlow_last_3yrs[s] <- ts_ci90$CI_low
  rel_inf_change_df$Nocc_f_90CIhigh_last_3yrs[s] <- ts_ci90$CI_high
  
  
  # counterfactual climate:
  
  load(file.path(main_dir, "results", "attribution", "fm_y_preds_routes_cf_1995_all",
                 paste0(spec, "_y_preds_cf_counterclim.RData")))
  # sum across routes for each year:
  preds_years <- apply(y_preds_route_cf, MAR = c(2,3), FUN = sum, na.rm = TRUE)
  
  # first 3 years:
  rel_inf_change_df$Nocc_cfclim_mean_first_3yrs[s] <- mean(preds_years[1:3,])
  ts_ci90 <- bayestestR::ci(preds_years[1:3,], ci = 0.7, method = "ETI") 
  rel_inf_change_df$Nocc_cfclim_90CIlow_first_3yrs[s] <- ts_ci90$CI_low
  rel_inf_change_df$Nocc_cfclim_90CIhigh_first_3yrs[s] <- ts_ci90$CI_high
  
  # last 3 years:
  rel_inf_change_df$Nocc_cfclim_mean_last_3yrs[s] <- mean(preds_years[23:25,])
  ts_ci90 <- bayestestR::ci(preds_years[23:25,], ci = 0.7, method = "ETI") 
  rel_inf_change_df$Nocc_cfclim_90CIlow_last_3yrs[s] <- ts_ci90$CI_low
  rel_inf_change_df$Nocc_cfclim_90CIhigh_last_3yrs[s] <- ts_ci90$CI_high
  
  # counterfactual land use:
  
  load(file.path(main_dir, "results", "attribution", "fm_y_preds_routes_cf_1995_all",
                 paste0(spec, "_y_preds_cf_1995soc.RData")))
  # sum across routes for each year:
  preds_years <- apply(y_preds_route_cf, MAR = c(2,3), FUN = sum, na.rm = TRUE)
  
  # first 3 years:
  rel_inf_change_df$Nocc_cflu_mean_first_3yrs[s] <- mean(preds_years[1:3,])
  ts_ci90 <- bayestestR::ci(preds_years[1:3,], ci = 0.7, method = "ETI") 
  rel_inf_change_df$Nocc_cflu_90CIlow_first_3yrs[s] <- ts_ci90$CI_low
  rel_inf_change_df$Nocc_cflu_90CIhigh_first_3yrs[s] <- ts_ci90$CI_high
  
  # last 3 years:
  rel_inf_change_df$Nocc_cflu_mean_last_3yrs[s] <- mean(preds_years[23:25,])
  ts_ci90 <- bayestestR::ci(preds_years[23:25,], ci = 0.7, method = "ETI") 
  rel_inf_change_df$Nocc_cflu_90CIlow_last_3yrs[s] <- ts_ci90$CI_low
  rel_inf_change_df$Nocc_cflu_90CIhigh_last_3yrs[s] <- ts_ci90$CI_high
  
  # counterfactual climate & land use:
  
  load(file.path(main_dir, "results", "attribution", "fm_y_preds_routes_cf_1995_all",
                 paste0(spec, "_y_preds_cf_counterclim_1995soc.RData")))
  # sum across routes for each year:
  preds_years <- apply(y_preds_route_cf, MAR = c(2,3), FUN = sum, na.rm = TRUE)
  
  # first 3 years:
  rel_inf_change_df$Nocc_cfclimlu_mean_first_3yrs[s] <- mean(preds_years[1:3,])
  ts_ci90 <- bayestestR::ci(preds_years[1:3,], ci = 0.7, method = "ETI") 
  rel_inf_change_df$Nocc_cfclimlu_90CIlow_first_3yrs[s] <- ts_ci90$CI_low
  rel_inf_change_df$Nocc_cfclimlu_90CIhigh_first_3yrs[s] <- ts_ci90$CI_high
  
  # last 3 years:
  rel_inf_change_df$Nocc_cfclimlu_mean_last_3yrs[s] <- mean(preds_years[23:25,])
  ts_ci90 <- bayestestR::ci(preds_years[23:25,], ci = 0.7, method = "ETI") # 0.9
  rel_inf_change_df$Nocc_cfclimlu_90CIlow_last_3yrs[s] <- ts_ci90$CI_low
  rel_inf_change_df$Nocc_cfclimlu_90CIhigh_last_3yrs[s] <- ts_ci90$CI_high
  
}

save(rel_inf_change_df, file = file.path(main_dir, "results", "attribution", "rel_inf_change_df_CI70.RData"))
#load(file = file.path(main_dir, "results", "attribution", "rel_inf_change_df.RData"))


## plot mean number of occupied routes + CI in first 3 years and last 3 years per species: ----

rel_inf_change_df_ref <- rel_inf_change_df %>% 
  filter(species %in% spec_okay) %>% 
  tidyr::pivot_longer(cols = -species, values_to = "value", names_to = "measure") %>%
  mutate(timeperiod = ifelse(grepl("first", measure), "1995 - 1997", "2017 - 2019"),
         timeperiod = factor(timeperiod, levels = c("2017 - 2019", "1995 - 1997"))) %>%
  mutate(measure = gsub("_first_3yrs", "", measure)) %>%  
  mutate(measure = gsub("_last_3yrs", "", measure)) %>% 
  tidyr::pivot_wider(names_from = measure, values_from = value) %>% 
  mutate(species = forcats::fct_reorder(species, Nocc_f_mean))

pdf(file = file.path("plots", "attribution", "change_Nocc_routes_first_last3yrs.pdf"), 
    width = 8.27, height = 11.69)
ggplot(data = rel_inf_change_df_ref, aes(y = species)) + 
  geom_point(aes(x = Nocc_f_mean, colour = timeperiod),
             position = position_dodge(0.8), size = 1.8) +
  geom_errorbar(aes(xmin = Nocc_f_90CIlow, xmax=Nocc_f_90CIhigh, colour = timeperiod), 
                width = 0.4, linewidth = 0.4, position = position_dodge(0.8)) +
  scale_colour_manual(values = c("#119DA4", "#084B83"), breaks = c("1995 - 1997", "2017 - 2019"), name = "time period") +
  theme_bw() +
  xlab("Number of occupied routes & 90 % CI") +
  ylab("") +
  theme(text = element_text(size = 14))
dev.off()


# (same for counterfactual climate):
ggplot(data = rel_inf_change_df_ref, aes(y = species)) + 
  geom_point(aes(x = Nocc_cfclim_mean, colour = timeperiod),
             position = position_dodge(0.8), size = 1.8) +
  geom_errorbar(aes(xmin = Nocc_cfclim_90CIlow, xmax=Nocc_cfclim_90CIhigh, colour = timeperiod), 
                width = 0.4, linewidth = 0.4, position = position_dodge(0.8)) +
  scale_colour_manual(values = c("#119DA4", "#084B83"), breaks = c("1995 - 1997", "2017 - 2019"), name = "time period") +
  theme_bw() +
  xlab("Number of occupied routes & 90 % CI") +
  ylab("") +
  theme(text = element_text(size = 14))

## alluvial plots: ----
# for category "stable" we need to define what's stable -> CIs

rel_inf_change_df <- rel_inf_change_df %>% 
  mutate(change_fact = ifelse(Nocc_f_mean_last_3yrs > Nocc_f_90CIhigh_first_3yrs,
                              "increase", 
                              ifelse(Nocc_f_mean_last_3yrs < Nocc_f_90CIlow_first_3yrs,
                                     "decrease", "stable"))) %>% 
  mutate(change_cfclim = ifelse(Nocc_cfclim_mean_last_3yrs > Nocc_cfclim_90CIhigh_first_3yrs,
                                "increase", 
                                ifelse(Nocc_cfclim_mean_last_3yrs < Nocc_cfclim_90CIlow_first_3yrs,
                                       "decrease", "stable"))) %>%
  mutate(change_cflu = ifelse(Nocc_cflu_mean_last_3yrs > Nocc_cflu_90CIhigh_first_3yrs,
                              "increase", 
                              ifelse(Nocc_cflu_mean_last_3yrs < Nocc_cflu_90CIlow_first_3yrs,
                                     "decrease", "stable"))) %>%  
  mutate(change_cfclimlu = ifelse(Nocc_cfclimlu_mean_last_3yrs > Nocc_cfclimlu_90CIhigh_first_3yrs,
                                  "increase", 
                                  ifelse(Nocc_cfclimlu_mean_last_3yrs < Nocc_cfclimlu_90CIlow_first_3yrs,
                                         "decrease", "stable"))) %>% 
  mutate(change_fact = factor(change_fact, levels = c("increase", "stable", "decrease")),
         change_cfclim = factor(change_cfclim, levels = c("increase", "stable", "decrease")),
         change_cflu = factor(change_cflu, levels = c("increase", "stable", "decrease")),
         change_cfclimlu = factor(change_cfclimlu, levels = c("increase", "stable", "decrease")))


jpeg(file = file.path("plots", "attribution", "change_Nocc_routes_first_last3yrs_climate_flow_val_okay_CI70.jpg"), 
     width = 1000, height = 800, quality = 100)
rel_inf_change_df %>%
  filter(species %in% spec_okay) %>% 
  # alluvial plot:
  group_by(change_fact, change_cfclim) %>% 
  summarise(n = n()) %>% 
  ggplot(aes(axis1 = change_cfclim, axis2 = change_fact, y = n)) +
  geom_alluvium(aes(fill = change_cfclim), show.legend = FALSE, width = 1/5, colour = "grey50") +
  geom_stratum(width = 1/5, aes(fill = change_cfclim), show.legend = FALSE) + 
  scale_fill_manual(values = c("increase" = "#EFCA08", "decrease" =  "#A17C6B", "stable" =  "#93C6D6"), na.value = NA) +
  geom_stratum(width = 1/5, aes(fill = change_fact), show.legend = FALSE) +
  scale_fill_manual(values = c("increase" = "#EFCA08", "decrease" =  "#A17C6B", "stable" =  "#93C6D6"), na.value = NA) +
  geom_text(stat = "stratum", aes(label = after_stat(stratum)), size = 8) +
  geom_text(stat = "flow",
            aes(label = ifelse(test = after_stat(x) == "1",
                               yes = rel_inf_change_df %>% 
                                 filter(species %in% spec_okay) %>%
                                 group_by(change_fact, change_cfclim) %>% 
                                 summarise(n = n()) %>% arrange(change_cfclim) %>% pull(n) %>% rev,
                               no = "")), size = 8, nudge_x = 0.15) +
  scale_x_discrete(limits = c("change_cfclim", "change_fact"), expand = c(0.15, 0.05)) +
  xlab("") + ylab("Number of species") +
  #scale_fill_viridis_d(direction = -1, na.value = NA) +
  theme_bw() +
  ggtitle("Change in area of occupancy, first 3 years vs. last 3 years, cons. 70% CI") +
  theme(axis.ticks.x = element_blank(), axis.ticks.y = element_blank(),
        axis.text.y = element_blank(), axis.text.x = element_blank(),
        text = element_text(size = 24),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  annotate("text", x = 1, y = 85, label = "Counterfactual\nclimate", hjust = 0.5, size = 8) + #170 y = 85 for spec_okay
  annotate("text", x = 2, y = 85, label = "Factual\nclimate", hjust = 0.5, size = 8) # y = 170
dev.off()


# same for land use:
jpeg(file = file.path("plots", "attribution", "change_Nocc_routes_first_last3yrs_land_use_flow.jpg"), 
     width = 1000, height = 800, quality = 100)
rel_inf_change_df %>%
  #filter(species %in% spec_okay) %>% 
  # alluvial plot:
  group_by(change_fact, change_cflu) %>% 
  summarise(n = n()) %>% 
  ggplot(aes(axis1 = change_cflu, axis2 = change_fact, y = n)) +
  geom_alluvium(aes(fill = change_cflu), show.legend = FALSE, width = 1/5, colour = "grey50") +
  geom_stratum(width = 1/5, aes(fill = change_cflu), show.legend = FALSE) + 
  scale_fill_manual(values = c("increase" = "#EFCA08", "decrease" =  "#A17C6B", "stable" =  "#93C6D6"), na.value = NA) +
  geom_stratum(width = 1/5, aes(fill = change_fact), show.legend = FALSE) +
  scale_fill_manual(values = c("increase" = "#EFCA08", "decrease" =  "#A17C6B", "stable" =  "#93C6D6"), na.value = NA) +
  geom_text(stat = "stratum", aes(label = after_stat(stratum)), size = 8) +
  geom_text(stat = "flow",
            aes(label = ifelse(test = after_stat(x) == "1",
                               yes = rel_inf_change_df %>% 
                                 #filter(species %in% spec_okay) %>%
                                 group_by(change_fact, change_cflu) %>% 
                                 summarise(n = n()) %>% arrange(change_cflu) %>% pull(n) %>% rev,
                               no = "")), size = 8, nudge_x = 0.15) +
  scale_x_discrete(limits = c("change_cflu", "change_fact"), expand = c(0.15, 0.05)) +
  xlab("") + ylab("Number of species") +
  #scale_fill_viridis_d(direction = -1, na.value = NA) +
  theme_bw() +
  ggtitle("Change in area of occupancy, first 3 years vs. last 3 years, cons. 90% CI") +
  theme(axis.ticks.x = element_blank(), axis.ticks.y = element_blank(),
        axis.text.y = element_blank(), axis.text.x = element_blank(),
        text = element_text(size = 24),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  annotate("text", x = 1, y = 170, label = "Counterfactual\nland use", hjust = 0.5, size = 8) + #170 y = 85 for spec_okay
  annotate("text", x = 2, y = 170, label = "Factual\nland use", hjust = 0.5, size = 8) # y = 170
dev.off()

rel_inf_change_df %>% 
  filter(species %in% spec_okay) %>% View


## traits: ----

category <- c("Habitat", "Migration", "Trophic.Level",
              "Trophic.Niche", "Primary.Lifestyle")

# merge occ. dyn. and traits:
dyn_traits_df <- rel_inf_change_df %>% 
  left_join(bbs_dt %>% select(English_Common_Name:Range.Size), by = c(species = "English_Common_Name")) %>% 
  distinct() 

for(c in 1:length(category)){
  
  print(category[c])
  
  
  filter_cat <- dyn_traits_df %>% 
    filter(species %in% spec_okay) %>% 
    group_by(.data[[category[c]]]) %>% 
    summarise(n = n()) %>% 
    filter(n >= 5) %>% 
    pull(.data[[category[c]]])
  
  # climate
  trait_df <- dyn_traits_df %>% 
    filter(species %in% spec_okay) %>% 
    group_by(change_fact, change_cfclim, .data[[category[c]]]) %>% 
    summarise(n = n())
  
  jpeg(file = file.path(plot_dir, paste0("change_Nocc_routes_climate_", category[c], "_prelim_spec_okay.jpg")), 
       width = 1400, height = 1000, quality = 100)
  # xx change to https://www.jbarg.net/posts/2024-08-05-how-the-heck-does-ggalluvial-work/
  print(ggplot(data = trait_df %>% 
                 filter(.data[[category[c]]] %in% filter_cat), 
               aes(axis1 = change_cfclim, axis2 = change_fact, y = n, label = n)) +
          geom_alluvium(aes(fill = change_cfclim), show.legend = FALSE, width = 1/5, colour = "grey50") +
          geom_stratum(width = 1/5, aes(fill = change_fact), show.legend = FALSE) +
          geom_stratum(width = 1/5, aes(fill = change_cfclim), show.legend = FALSE) + 
          scale_fill_manual(values = c("increase" = "#EFCA08", "decrease" =  "#A17C6B", "stable" =  "#93C6D6"), 
                           na.value = NA) +
          geom_text(stat = "stratum", aes(label = after_stat(stratum)), size = 6) +
          scale_x_discrete(limits = c("change_cfclim", "change_fact"),
                           expand = c(0.15, 0.05),
                           labels = c("Counterfactual\nclimate", "Factual\nclimate"),
                           #position = "top"
          ) +
          facet_wrap(~.data[[category[c]]], scales = 'free_y', axes = "all_x") +
          geom_text(stat = "flow", size = 6, nudge_x = 0.15) +
          ylab("Number of species") +  xlab("") +
          theme_bw() +
          ggtitle("Change in area of occupancy") +
          theme(axis.ticks.x = element_blank(),
                axis.ticks.y = element_blank(),
                axis.text.y = element_blank(),
                text = element_text(size = 22),
                panel.grid.major = element_blank(), 
                panel.grid.minor = element_blank())
  )
  dev.off()
  
  
  # land use:
  trait_df <- dyn_traits_df %>% 
    filter(species %in% spec_okay) %>% 
    group_by(change_fact, change_cflu, .data[[category[c]]]) %>% 
    summarise(n = n())
  
  jpeg(file = file.path(plot_dir, paste0("change_Nocc_routes_landuse_", category[c], "_prelim_spec_okay.jpg")), 
       width = 1400, height = 1000, quality = 100)
  
  print(ggplot(data = trait_df %>% filter(.data[[category[c]]] %in% filter_cat),
               aes(axis1 = change_cflu, axis2 = change_fact, y = n, label = n)) +
          geom_alluvium(aes(fill = change_cflu), show.legend = FALSE, width = 1/5) +
          geom_stratum(width = 1/5, aes(fill = change_fact), show.legend = FALSE) +
          geom_stratum(width = 1/5, aes(fill = change_cflu), show.legend = FALSE) + 
          scale_fill_manual(values = c("increase" = "#EFCA08", "decrease" =  "#A17C6B", "stable" =  "#93C6D6"), 
                            na.value = NA) +
          geom_text(stat = "stratum", aes(label = after_stat(stratum)), size = 6) +
          scale_x_discrete(limits = c("change_cflu", "change_fact"),
                           expand = c(0.15, 0.05),
                           labels = c("Counterfactual\nland use", "Factual\nland use"),
                           #position = "top"
          ) +
          facet_wrap(~.data[[category[c]]], scales = 'free_y', axes = "all_x") +
          geom_text(stat = "flow", size = 6, nudge_x = 0.15) +
          xlab("") + ylab("Number of species") +
          theme_bw() +
          ggtitle("Change in area of occupancy") +
          theme(axis.ticks.x = element_blank(),
                axis.ticks.y = element_blank(),
                axis.text.y = element_blank(),
                text = element_text(size = 22),
                panel.grid.major = element_blank(), 
                panel.grid.minor = element_blank())
  )
  dev.off()
  
}



# relative occupancy dynamics; factual vs. counterfactual: ---------------------
# inspired by Langhammer et al. categories

## option 1: linear trend based on median number of occupied routes per year ------

# some visual assessment:


#   # climate change: absolute increase = 1, relative increase = 2, relative decrease = 3, absolute decrease = 4
#   # comparable to Langhammer et al. impact categories
# 
# rel_infl_df$impact_cat[rel_infl_df$species == "Acadian Flycatcher"] <- 1
# rel_infl_df$impact_cat[rel_infl_df$species == "American Kestrel"] <- 4
# rel_infl_df$impact_cat[rel_infl_df$species == "American Redstart"] <- 2
# rel_infl_df$impact_cat[rel_infl_df$species == "Baltimore Oriole"] <- 4
# rel_infl_df$impact_cat[rel_infl_df$species == "Bell's Vireo"] <- 1
# rel_infl_df$impact_cat[rel_infl_df$species == "Black Vulture"] <- 1
# rel_infl_df$impact_cat[rel_infl_df$species == "Black-billed Cuckoo"] <- 4
# rel_infl_df$impact_cat[rel_infl_df$species == "Blackburnian Warbler"] <- 2
# rel_infl_df$impact_cat[rel_infl_df$species == "Black-capped Chickadee"] <- 2
# rel_infl_df$impact_cat[rel_infl_df$species == "Black-headed Grosbeak"] <- 2
# rel_infl_df$impact_cat[rel_infl_df$species == "Black-throated Green Warbler"] <- 2
# rel_infl_df$impact_cat[rel_infl_df$species == "Black-throated Sparrow"] <- 1
# rel_infl_df$impact_cat[rel_infl_df$species == "Blue-gray Gnatcatcher"] <- 1
# rel_infl_df$impact_cat[rel_infl_df$species == "Bobolink"] <- 4
# rel_infl_df$impact_cat[rel_infl_df$species == "Brown Creeper"] <- 2
# rel_infl_df$impact_cat[rel_infl_df$species == "Bullock's Oriole"] <- 2
# rel_infl_df$impact_cat[rel_infl_df$species == "Carolina Chickadee"] <- 1
# rel_infl_df$impact_cat[rel_infl_df$species == "Carolina Wren"] <- 1
# rel_infl_df$impact_cat[rel_infl_df$species == "Cassin's Finch"] <- 2
# rel_infl_df$impact_cat[rel_infl_df$species == "Clark's Nutcracker"] <- 4
# rel_infl_df$impact_cat[rel_infl_df$species == "Clay-colored Sparrow"] <- 2
# rel_infl_df$impact_cat[rel_infl_df$species == "Cliff Swallow"] <- 2
# rel_infl_df$impact_cat[rel_infl_df$species == "Common Raven"] <- 2
# rel_infl_df$impact_cat[rel_infl_df$species == "Common Yellowthroat"] <- 4
# rel_infl_df$impact_cat[rel_infl_df$species == "Dark-eyed Junco"] <- 4
# rel_infl_df$impact_cat[rel_infl_df$species == "Dickcissel"] <- 1
# rel_infl_df$impact_cat[rel_infl_df$species == "Eastern Bluebird"] <- 2
# rel_infl_df$impact_cat[rel_infl_df$species == "Eastern Kingbird"] <- 4
# rel_infl_df$impact_cat[rel_infl_df$species == "Eastern Meadowlark"] <- 3 # xx
# rel_infl_df$impact_cat[rel_infl_df$species == "Eastern Phoebe"] <- 2
# rel_infl_df$impact_cat[rel_infl_df$species == "Fish Crow"] <- 1
# rel_infl_df$impact_cat[rel_infl_df$species == "Golden-crowned Kinglet"] <- 4
# rel_infl_df$impact_cat[rel_infl_df$species == "Grasshopper Sparrow"] <- 4
# rel_infl_df$impact_cat[rel_infl_df$species == "Great-tailed Grackle"] <- 1
# rel_infl_df$impact_cat[rel_infl_df$species == "Hammond's Flycatcher"] <- 2
# rel_infl_df$impact_cat[rel_infl_df$species == "Hermit Thrush"] <- 4
# rel_infl_df$impact_cat[rel_infl_df$species == "House Wren"] <- 2 # xx
# rel_infl_df$impact_cat[rel_infl_df$species == "Ladder-backed Woodpecker"] <- 1
# rel_infl_df$impact_cat[rel_infl_df$species == "Lark Sparrow"] <- 2
# rel_infl_df$impact_cat[rel_infl_df$species == "Lazuli Bunting"] <- 2
# rel_infl_df$impact_cat[rel_infl_df$species == "Least Flycatcher"] <- 4
# rel_infl_df$impact_cat[rel_infl_df$species == "Lesser Goldfinch"] <- 2
# rel_infl_df$impact_cat[rel_infl_df$species == "Marsh Wren"] <- 2
# rel_infl_df$impact_cat[rel_infl_df$species == "Mississippi Kite"] <- 1
# rel_infl_df$impact_cat[rel_infl_df$species == "Mountain Bluebird"] <- 4
# rel_infl_df$impact_cat[rel_infl_df$species == "Mountain Chickadee"] <- 4
# rel_infl_df$impact_cat[rel_infl_df$species == "Mourning Warbler"] <- 4
# rel_infl_df$impact_cat[rel_infl_df$species == "Northern Bobwhite"] <- 3 # xx
# rel_infl_df$impact_cat[rel_infl_df$species == "Northern Flicker"] <- 4
# rel_infl_df$impact_cat[rel_infl_df$species == "Northern Harrier"] <- 4
# rel_infl_df$impact_cat[rel_infl_df$species == "Northern Mockingbird"] <- 1
# rel_infl_df$impact_cat[rel_infl_df$species == "Olive-sided Flycatcher"] <- 4
# rel_infl_df$impact_cat[rel_infl_df$species == "Pine Siskin"] <- 4
# rel_infl_df$impact_cat[rel_infl_df$species == "Pine Warbler"] <- 2
# rel_infl_df$impact_cat[rel_infl_df$species == "Prothonotary Warbler"] <- 3 # xx
# rel_infl_df$impact_cat[rel_infl_df$species == "Red-eyed Vireo"] <- 2
# rel_infl_df$impact_cat[rel_infl_df$species == "Red-shouldered Hawk"] <- 1
# rel_infl_df$impact_cat[rel_infl_df$species == "Ring-necked Pheasant"] <- 4
# rel_infl_df$impact_cat[rel_infl_df$species == "Rock Pigeon"] <- 4
# rel_infl_df$impact_cat[rel_infl_df$species == "Savannah Sparrow"] <- 4
# rel_infl_df$impact_cat[rel_infl_df$species == "Scissor-tailed Flycatcher"] <- 1
# rel_infl_df$impact_cat[rel_infl_df$species == "Sedge Wren"] <- 4
# rel_infl_df$impact_cat[rel_infl_df$species == "Summer Tanager"] <- 1
# rel_infl_df$impact_cat[rel_infl_df$species == "Tree Swallow"] <- 2
# rel_infl_df$impact_cat[rel_infl_df$species == "Vesper Swallow"] <- 4
# rel_infl_df$impact_cat[rel_infl_df$species == "Western Meadowlark"] <- 4
# rel_infl_df$impact_cat[rel_infl_df$species == "White-breasted Nuthatch"] <- 2
# rel_infl_df$impact_cat[rel_infl_df$species == "White-eyed Vireo"] <- 1
# rel_infl_df$impact_cat[rel_infl_df$species == "White-throated Sparrow"] <- 4
# rel_infl_df$impact_cat[rel_infl_df$species == "White-winged Dove"] <- 1
# rel_infl_df$impact_cat[rel_infl_df$species == "Yellow-bellied Sapsucker"] <- 2
# rel_infl_df$impact_cat[rel_infl_df$species == "Yellow-headed Blackbird"] <- 4
# table(rel_infl_df$impact_cat, useNA = "always")

# -> less subjectively:

rel_infl_df_lf <- rel_infl_df %>% 
  #filter(species %in% spec_okay) %>% 
  select(c(species, starts_with("slope"), starts_with("p_"))) %>% 
  tidyr::pivot_longer(cols = matches("(slope_)|(p_)"), names_to = c("metric", "scenario"), 
                      values_to = "value", names_pattern = "(.*)_(.*)") %>% 
  tidyr::pivot_wider(names_from = metric, values_from = value) %>% 
  mutate(dynamics = ifelse(slope > 0 & p < 0.05, "increase",
                           ifelse(slope < 0 & p < 0.05, "decrease", "stable")))

# factual - cf

# increase and difference positive = absolute increase
# decrease and difference negative = absolute decrease
rel_infl_df$cat_relimp_clim <- NA
rel_infl_df$cat_relimp_lu <- NA
rel_infl_df$cat_relimp_climlu <- NA

for(s in 1:length(final_species)){
  
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
  
  dyn <- rel_infl_df_lf %>% 
    filter(species == spec & scenario == "fact") %>% 
    pull(dynamics)
  
  # compare sign of sum of difference between factual and counterfactual:
  direct <- sign(sum(ts_preds_fact$median_Nocc_f - ts_preds_cfact$cf_clim$median_Nocc_cf_clim))
  
  if(dyn == "increase" & direct == 1) rel_infl_df$cat_relimp_clim[s] <- 1 # absolute increase
  if(dyn == "decrease" & direct == -1) rel_infl_df$cat_relimp_clim[s] <- 4 # absolute decrease
  if(dyn %in% c("stable", "decrease") & direct == 1) rel_infl_df$cat_relimp_clim[s] <- 2 # relative increase (occ. is not increasing, but higher than without change, increase only relative to cf.)
  if(dyn %in% c("stable", "increase") & direct == -1) rel_infl_df$cat_relimp_clim[s] <- 3 # relative decrease (occ. is not decreasing, but lower than without change)
  
  direct_lu <- sign(sum(ts_preds_fact$median_Nocc_f - ts_preds_cfact$cf_1995soc$median_Nocc_cf_1995soc))
  
  if(dyn == "increase" & direct_lu == 1) rel_infl_df$cat_relimp_lu[s] <- 1 # absolute increase
  if(dyn == "decrease" & direct_lu == -1) rel_infl_df$cat_relimp_lu[s] <- 4 # absolute decrease
  if(dyn %in% c("stable", "decrease") & direct_lu == 1) rel_infl_df$cat_relimp_lu[s] <- 2 # relative increase (occ. is not increasing, but higher than without change, increase only relative to cf.)
  if(dyn %in% c("stable", "increase") & direct_lu == -1) rel_infl_df$cat_relimp_lu[s] <- 3 # relative decrease (occ. is not decreasing, but lower than without change)
  
  direct_climlu <- sign(sum(ts_preds_fact$median_Nocc_f - ts_preds_cfact$cf_clim_1995soc$median_Nocc_cf_clim_1995soc))
  
  if(dyn == "increase" & direct_climlu == 1) rel_infl_df$cat_relimp_climlu[s] <- 1 # absolute increase
  if(dyn == "decrease" & direct_climlu == -1) rel_infl_df$cat_relimp_climlu[s] <- 4 # absolute decrease
  if(dyn %in% c("stable", "decrease") & direct_climlu == 1) rel_infl_df$cat_relimp_climlu[s] <- 2 # relative increase (occ. is not increasing, but higher than without change, increase only relative to cf.)
  if(dyn %in% c("stable", "increase") & direct_climlu == -1) rel_infl_df$cat_relimp_climlu[s] <- 3 # relative decrease (occ. is not decreasing, but lower than without change)
  
  
}
rel_infl_df %>% 
  filter(species %in% spec_okay) %>%  View


# climate change "winners" (for which model performance is okay)
cc_winners <- rel_infl_df %>% 
  filter(cat_relimp_clim %in% c(1,2)) %>% # absolute (1) and relative increase (2)
  filter(species %in% spec_okay) %>%
  select(-matches("(slope)|(^p_)")) %>% 
  mutate(infl_clim = round(infl_clim, 2),
         infl_lu  = round(infl_lu , 2),
         infl_climlu  = round(infl_climlu , 2))
summary(cc_winners)

cc_losers <- rel_infl_df %>% 
  filter(cat_relimp_clim %in% c(3,4)) %>% # absolute (4) and relative decrease (3)
  filter(species %in% spec_okay) %>%
  select(-matches("(slope)|(^p_)")) %>% 
  mutate(infl_clim = round(infl_clim, 2),
         infl_lu  = round(infl_lu , 2),
         infl_climlu  = round(infl_climlu , 2))
summary(cc_losers)


### lollipop plots: ----
rel_infl_df %>% 
  #filter(species %in% spec_okay) %>%
  select(species, cat_relimp_clim, cat_relimp_lu, cat_relimp_climlu) %>% 
  tidyr::pivot_longer(cols = starts_with("cat"), names_to = "scenario", values_to = "rel_impact_cat") %>% 
  group_by(scenario, rel_impact_cat) %>% 
  summarise(n = n()) %>% 
  mutate(n_perc = n/length(final_species)*100) %>% 
  mutate(scenario = factor(scenario, levels = c("cat_relimp_climlu", "cat_relimp_lu", "cat_relimp_clim"))) %>% 
  mutate(scenario = recode(scenario, 
                           cat_relimp_clim = "factual vs.\ncounterfactual climate",
                           cat_relimp_lu = "factual vs.\ncounterfactual land use", 
                           cat_relimp_climlu = "factual vs.\ncounterfactual climate + land use")) %>% 
  mutate(rel_impact_cat = factor(rel_impact_cat, levels = c("4","3","2","1"))) %>% 
  mutate(rel_impact_cat = recode(rel_impact_cat, 
                           "1" = "absolute increase",
                           "2" = "relative increase", 
                           "3" = "relative decrease",
                           "4" = "absolute decrease")) %>% 
  ggplot(aes(x = scenario, y = n_perc)) +
  geom_errorbar(aes(ymin = 0, ymax = n_perc, colour = rel_impact_cat),
                position = position_dodge(width = .6), width = 0, linewidth = 1) +
  geom_point(aes(colour = rel_impact_cat), size=4, position = position_dodge(width = .6)) +
  ylab("N species [%]") +
  xlab("") +
  ggtitle("Occupancy dynamics 1995 - 2019") +
  guides(colour = guide_legend(reverse = TRUE)) +
  theme_bw() +
  theme(text = element_text(size = 15), legend.title=element_blank()) +
  coord_flip()

#save(rel_infl_df, file = file.path(main_dir, "results", "attribution", "rel_influence1.RData"))
#load(file = file.path(main_dir, "results", "attribution", "rel_influence1.RData"))


### alternative: compare slopes to quantify relative change categories: ----
# 
# rel_infl_df <- rel_infl_df %>% 
#   # counterfactual climate:
#   mutate(cat_relimp_clim2 = case_when(
#     slope_fact > slope_cfclim ~ "1", # absolute increase
#     slope_fact < 0 & slope_fact > slope_cfclim ~ "2", # relative increase (occ. is not increasing, but higher than without change, increase only relative to cf.)
#     slope_fact > 0 & slope_fact < slope_cfclim ~ "3", # relative decrease (occ. is not decreasing, but lower than without change)
#     slope_fact < slope_cfclim ~ "4",  # absolute decrease
#     TRUE ~ NA_character_)) %>% 
#   # counterfactual land use:
#   mutate(cat_relimp_lu2 = case_when(
#     slope_fact > slope_cflu ~ "1", # absolute increase
#     slope_fact < 0 & slope_fact > slope_cflu ~ "2", # relative increase (occ. is not increasing, but higher than without change, increase only relative to cf.)
#     slope_fact > 0 & slope_fact < slope_cflu ~ "3", # relative decrease (occ. is not decreasing, but lower than without change)
#     slope_fact < slope_cflu ~ "4",  # absolute decrease
#     TRUE ~ NA_character_)) %>% 
#   # counterfactual climate and land use:
#   mutate(cat_relimp_climlu2 = case_when(
#     slope_fact > slope_cfclimlu ~ "1", # absolute increase
#     slope_fact < 0 & slope_fact > slope_cfclimlu ~ "2", # relative increase (occ. is not increasing, but higher than without change, increase only relative to cf.)
#     slope_fact > 0 & slope_fact < slope_cfclimlu ~ "3", # relative decrease (occ. is not decreasing, but lower than without change)
#     slope_fact < slope_cfclimlu ~ "4",  # absolute decrease
#     TRUE ~ NA_character_)) %>% 
#   mutate(cat_relimp_clim2 = as.numeric(cat_relimp_clim2),
#          cat_relimp_lu2 = as.numeric(cat_relimp_lu2),
#          cat_relimp_climlu2 = as.numeric(cat_relimp_climlu2))
# 
# rel_infl_df %>% 
#   mutate(test = ifelse(cat_relimp_clim2 == cat_relimp_clim & 
#                          cat_relimp_lu2 == cat_relimp_lu &
#                          cat_relimp_climlu2 == cat_relimp_climlu, 0,1)) %>%  View
# 
# # same categories only for 29 species, only comparing slopes not sufficient

##  option 2: linear trend based on 100 draws of posterior for number of occupied routes per year: ----


rel_infl_df_lf
rel_infl_df_lf_opt2 <- rel_infl_df %>% 
  #filter(species %in% spec_okay) %>% 
  select(c(species, starts_with("opt2_slope"), starts_with("opt2_p_"))) %>% 
  tidyr::pivot_longer(cols = matches("(slope_)|(p_)"), names_to = c("metric", "scenario"), 
                      values_to = "value", names_pattern = "(.*)_(.*)") %>% 
  tidyr::pivot_wider(names_from = metric, values_from = value) %>% 
  mutate(dynamics = ifelse(opt2_slope > 0 & opt2_p < 0.05, "increase",
                           ifelse(opt2_slope < 0 & opt2_p < 0.05, "decrease", "stable")))

# factual - cf

# increase and difference positive = absolute increase
# decrease and difference negative = absolute decrease
rel_infl_df$opt2_cat_relimp_clim <- NA
rel_infl_df$opt2_cat_relimp_lu <- NA
rel_infl_df$opt2_cat_relimp_climlu <- NA

for(s in 1:length(final_species)){
  
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
  
  dyn <- rel_infl_df_lf_opt2 %>% 
    filter(species == spec & scenario == "fact") %>% 
    pull(dynamics)
  
  # compare sign of sum of difference between factual and counterfactual:
  #direct <- sign(sum(ts_preds_fact$median_Nocc_f - ts_preds_cfact$cf_clim$median_Nocc_cf_clim))
  # (since I look at sign of sum of differences, I think it should not make a difference
  # whether I quantify them based on the median for each year or 100 draws of the posterior for each year)
  direct <- sign(sum(ts_preds_fact %>% select(starts_with("draw")) - ts_preds_cfact$cf_clim %>% select(starts_with("draw"))))

  if(dyn == "increase" & direct == 1) rel_infl_df$opt2_cat_relimp_clim[s] <- 1 # absolute increase
  if(dyn == "decrease" & direct == -1) rel_infl_df$opt2_cat_relimp_clim[s] <- 4 # absolute decrease
  if(dyn %in% c("stable", "decrease") & direct == 1) rel_infl_df$opt2_cat_relimp_clim[s] <- 2 # relative increase (occ. is not increasing, but higher than without change, increase only relative to cf.)
  if(dyn %in% c("stable", "increase") & direct == -1) rel_infl_df$opt2_cat_relimp_clim[s] <- 3 # relative decrease (occ. is not decreasing, but lower than without change)
  
  direct_lu <- sign(sum(ts_preds_fact %>% select(starts_with("draw")) - ts_preds_cfact$cf_1995soc %>% select(starts_with("draw"))))
  
  if(dyn == "increase" & direct_lu == 1) rel_infl_df$opt2_cat_relimp_lu[s] <- 1 # absolute increase
  if(dyn == "decrease" & direct_lu == -1) rel_infl_df$opt2_cat_relimp_lu[s] <- 4 # absolute decrease
  if(dyn %in% c("stable", "decrease") & direct_lu == 1) rel_infl_df$opt2_cat_relimp_lu[s] <- 2 # relative increase (occ. is not increasing, but higher than without change, increase only relative to cf.)
  if(dyn %in% c("stable", "increase") & direct_lu == -1) rel_infl_df$opt2_cat_relimp_lu[s] <- 3 # relative decrease (occ. is not decreasing, but lower than without change)
  
  direct_climlu <- sign(sum(ts_preds_fact %>% select(starts_with("draw")) - ts_preds_cfact$cf_clim_1995soc %>% select(starts_with("draw"))))
  
  if(dyn == "increase" & direct_climlu == 1) rel_infl_df$opt2_cat_relimp_climlu[s] <- 1 # absolute increase
  if(dyn == "decrease" & direct_climlu == -1) rel_infl_df$opt2_cat_relimp_climlu[s] <- 4 # absolute decrease
  if(dyn %in% c("stable", "decrease") & direct_climlu == 1) rel_infl_df$opt2_cat_relimp_climlu[s] <- 2 # relative increase (occ. is not increasing, but higher than without change, increase only relative to cf.)
  if(dyn %in% c("stable", "increase") & direct_climlu == -1) rel_infl_df$opt2_cat_relimp_climlu[s] <- 3 # relative decrease (occ. is not decreasing, but lower than without change)
  
  
}
#save(rel_infl_df, file = file.path(main_dir, "results", "attribution", "rel_influence_080725.RData"))

# rel_infl_df %>% 
#   mutate(diff1 = cat_relimp_clim - opt2_cat_relimp_clim,
#          diff2 = cat_relimp_lu - opt2_cat_relimp_lu,
#          diff3 = cat_relimp_climlu - opt2_cat_relimp_climlu) %>%  
#   filter(species %in% spec_okay) %>% View



# quantify relative influence of climate and land use change for single species: ----

## option 1: difference mean absolute percentage error: ----

# how much does mean absolute percentage error differ between predictions
# for factual scenario and predictions for counterfactual scenario

# influence = difference in mean absolute percentage error:
rel_infl_df2 <- rel_infl_df %>% 
  filter(species %in% spec_okay) %>%
  select(-matches("(slope)|(p_)")) %>% 
  mutate(infl_clim = round(mape_cfclim - mape_fact, 2),
         infl_lu = round(mape_cflu - mape_fact, 2),
         infl_climlu = round(mape_cfclimlu - mape_fact, 2))

## option 2: difference mean absolute error: ----


# ratio mean absolute error for counterfactual vs. mean absolute error for factual
# rel_infl_df2 <- rel_infl_df %>% 
#   filter(species %in% spec_okay) %>%
#   select(-matches("(slope)|(p_)")) %>% 
#   mutate(infl_clim = mae_cfclim/mae_fact,
#          infl_lu = mae_cflu/mae_fact,
#          infl_climlu = mae_cfclimlu/mae_fact)

# difference between mean absolute error (MAE_cf - MAE_f):
# xx

## plot relative influence: ----

### sort by factual overall occ. dyn. (opt2) ----
spec_fact_dyn_df <- rel_infl_df %>% 
  filter(species %in% spec_okay) %>%
  select(c(species, starts_with("opt2_slope"), starts_with("opt2_p_"))) %>% 
  mutate(fact = ifelse(opt2_slope_fact > 0 & opt2_p_fact  < 0.05, "increase",
                       ifelse(opt2_slope_fact < 0 & opt2_p_fact  < 0.05, "decrease", "stable"))) %>% 
  mutate(fact = factor(fact, levels = c("increase", "stable", "decrease"))) %>% 
  select(species, fact) %>% 
  arrange(fact, species) %>% 
  mutate(colour = ifelse(fact == "increase", "#EFCA08",
                         ifelse(fact == "decrease", "#7DC4C9", "#EDE6F2"))) %>%
  # add species orders etc. as sorting options:
  left_join(bbs_dt %>% select(English_Common_Name, Scientific_Name, ORDER, Family, Habitat, Migration, Trophic.Level, Trophic.Niche, Primary.Lifestyle) %>% distinct,
            by = c(species = "English_Common_Name")) %>% 
  mutate(ORDER2 = factor(ORDER, levels = names(sort(table(ORDER), decreasing = TRUE))))

plot_df <- rel_infl_df2 %>%
  select(c(species, starts_with("infl"))) %>% 
  left_join(spec_fact_dyn_df) %>% 
  tidyr::pivot_longer(cols = starts_with("infl"), names_to = "scenario", values_to = "value") %>% 
  mutate(scenario = factor(scenario, levels = c("infl_clim", "infl_lu", "infl_climlu"))) %>% 
  mutate(scenario = recode(scenario, infl_clim = "climate", infl_lu = "land use", infl_climlu = "climate + land use")) %>% 
  #arrange(fact, ORDER2, Family) %>%
  #arrange(Habitat, fact, ORDER2) %>%
  arrange(Trophic.Level, fact, ORDER2) %>%
  mutate(plot_order = row_number())

p <- plot_df %>% 
  mutate(species = as.factor(species)) %>% 
  ggplot(aes(x = scenario, y = forcats::fct_reorder(species, desc(plot_order)), fill = value)) + 
  geom_tile() +
  #facet_grid(rows = vars(fact), scales = "free", space = "free") +
  facet_grid(rows = vars(Trophic.Level), scales = "free", space = "free") +
  viridis::scale_fill_viridis(discrete=FALSE#, 
                                #trans = scales::transform_log()
                                #direction = -1
                                ) +
  theme_bw() +
  labs(title = "Relative influence of climate and land use change on occupancy dynamics",
       subtitle = "Species ordered by factual occupancy dynamics, order and family") +
  xlab("") + ylab("") +
  labs(fill = expression(MAPE [counterfactual] - MAPE [factual])) +
  scale_x_discrete(expand = expansion(mult = c(0, 0.3))) + # 0 space on the left, 30 % of space between plot and facet strips
  theme(plot.title = element_text(hjust = 0),
        plot.subtitle = element_text(hjust = 0),
        plot.title.position = "plot",
        text = element_text(size = 12),
        axis.text.x = element_text(size = 10),
        axis.text.y = element_text(family = "serif", vjust = 0.3),
        strip.text.y = element_text(angle = 0),
        strip.background = element_rect(colour=NA),
        panel.border =  element_rect(colour=NA),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()#,
        #axis.ticks = element_blank()
        )
 p

#plot_df %>% arrange(desc(Habitat), desc(fact), desc(ORDER2), desc(species)) %>% filter(scenario == "climate") %>% View 
plot_df %>% filter(scenario == "climate") %>% View

# add different colours for facet strips:
g <- ggplot_gtable(ggplot_build(p))
stripr <- which(grepl('strip-r', g$layout$name))
fills <- c("#EFCA08","#EDE6F2","#7DC4C9")
k <- 1
for (i in stripr) {
  j <- which(grepl('rect', g$grobs[[i]]$grobs[[1]]$childrenOrder))
  g$grobs[[i]]$grobs[[1]]$children[[j]]$gp$fill <- fills[k]
  k <- k+1
}
grid::grid.draw(g)

# other attempts to get occ. dyn. categories next to y axis:
# p2 <- p + theme(plot.margin = unit(c(1,1,1,5), "lines"),
#           panel.background = element_rect(fill='transparent')) +
#   annotation_custom(
#     grob = grid::rectGrob(gp = grid::gpar(fill = "#7DC4C9", col = "#7DC4C9", lty = "solid")),
#     xmin = -1, xmax = -Inf, ymin = -Inf, ymax = 32.5) +
#   annotation_custom(
#     grob = grid::rectGrob(gp = grid::gpar(fill = "#EDE6F2", col = "#EDE6F2", lty = "solid")),
#     xmin = -1, xmax = -Inf, ymin = 32.5, ymax = 35.5) +
#   annotation_custom(
#     grob = grid::rectGrob(gp = grid::gpar(fill = "#EFCA08", col = "#EFCA08", lty = "solid")),
#     xmin = -1, xmax = -Inf, ymin = 35.5, ymax = Inf)
# 
# gt <- ggplot_gtable(ggplot_build(p2))
# gt$layout$clip[gt$layout$name == "panel"] <- "off"
# grid::grid.draw(gt)

# # patchwork package based:
# library(patchwork)
# p2 <- p + 
#   annotate("rect", ymin = c(-Inf, 32, 35), xmin = c(-Inf, -Inf, -Inf),
#            ymax = c(32, 32+3, Inf), xmax = c(Inf, Inf, Inf),
#            fill = c("#7DC4C9", "#EDE6F2", "#EFCA08")) +
#   annotate("text", label = c("decrease", "stable", "increase"), size = 6,
#            x = rep(2, 3),
#            y = c(32/2, 32+(3/2), 32+3+(45/2))) +
#   theme_void()
# p2 | p + plot_layout(widths = 1, byrow = TRUE) # both plots side-by-side



### sort by rel. change in overall occ. dyn. (~ Langhammer def.)----

# for sorting categories and colours:

spec_fact_dyn_df_opt2 <- rel_infl_df %>% 
  filter(species %in% spec_okay) %>%
  select(c(species, starts_with("opt2_cat_relimp"))) %>% # xx one colour column for each scenario xx
  mutate(across(opt2_cat_relimp_clim:opt2_cat_relimp_climlu, ~ 
                  recode(.x, "1" = "absolute\nincrease",
                             "2" = "relative\nincrease", 
                             "3" = "relative\ndecrease",
                             "4" = "absolute\ndecrease"))) %>% 
  mutate(across(opt2_cat_relimp_clim:opt2_cat_relimp_climlu, ~ 
                  factor(.x, levels = c("absolute\nincrease",
                                        "relative\nincrease",
                                        "relative\ndecrease",
                                        "absolute\ndecrease")))) %>% 
  # choose scenario:
  #arrange(opt2_cat_relimp_clim, species) %>% 
  #arrange(opt2_cat_relimp_lu, species) %>% 
  arrange(opt2_cat_relimp_climlu, species) %>% 
  
  # not needed if colouring facet strips:
  # mutate(colour = case_when(
  #   opt2_cat_relimp_clim == "relative decrease" ~ "#abd9e9",
  #   opt2_cat_relimp_clim == "absolute decrease" ~  "#2c7bb6",
  #   opt2_cat_relimp_clim == "absolute increase" ~ "#d7191c", #"#17BEBB",# "#00b2ca", #"#006DAA", #,
  #   opt2_cat_relimp_clim == "relative increase" ~ "#fdae61",
  #   .default = "gray90")) %>% 
  # add species orders etc. as sorting options:
  left_join(bbs_dt %>% select(English_Common_Name, Scientific_Name, ORDER, Family) %>% distinct,
            by = c(species = "English_Common_Name")) %>% 
  mutate(ORDER2 = factor(ORDER, levels = names(sort(table(ORDER), decreasing = TRUE))))


# relative influence values:
plot_df <- rel_infl_df2 %>%
  select(c(species, starts_with("infl"))) %>%
  left_join(spec_fact_dyn_df_opt2) %>% 
  tidyr::pivot_longer(cols = starts_with("infl"), names_to = "scenario", values_to = "value") %>% 
  mutate(scenario = factor(scenario, levels = c("infl_clim", "infl_lu", "infl_climlu"))) %>% 
  mutate(scenario = recode(scenario, infl_clim = "climate", infl_lu = "land use", infl_climlu = "climate + land use")) %>% 
  # choose scenario to plot: 
  #arrange(opt2_cat_relimp_clim, ORDER2, Family) %>%
  #arrange(opt2_cat_relimp_lu, ORDER2, Family) %>%
  arrange(opt2_cat_relimp_climlu , ORDER2, Family) %>%
  mutate(plot_order = row_number())


# plot 2: categorize species by relative change in occupancy dynamics
# since this depends on whether looking at counterfactual climate or land use
# I'll plot these separately:



p_rel <- plot_df %>% 
  # choose scenario to plot: 
  #filter(scenario == "climate") %>% 
  #filter(scenario == "land use") %>% 
  filter(scenario == "climate + land use") %>% 
  mutate(species = as.factor(species)) %>% 
  ggplot(aes(x = scenario, y = forcats::fct_reorder(species, desc(plot_order)), fill = value)) + 
  geom_tile() +
  # choose scenario to plot: 
  #facet_grid(rows = vars(opt2_cat_relimp_clim), scales = "free", space = "free") +
  #facet_grid(rows = vars(opt2_cat_relimp_lu), scales = "free", space = "free") +
  facet_grid(rows = vars(opt2_cat_relimp_climlu), scales = "free", space = "free") +
  viridis::scale_fill_viridis(discrete=FALSE#, 
                              #trans = scales::transform_log()
                              #direction = -1
  ) +
  theme_bw() +
  # choose scenario to plot:
  #labs(title = "Relative influence of climate change on occupancy dynamics") +
  #labs(title = "Relative influence of land use change on occupancy dynamics") +
  labs(title = "Relative influence of climate and land use change on occupancy dynamics") +
  labs(subtitle = "Species ordered by relative occupancy dynamics, order and family") +
  xlab("") + ylab("") +
  labs(fill = expression(MAPE [counterfactual] - MAPE [factual])) +
  scale_x_discrete(expand = expansion(mult = c(0, 0.55))) + # 0 space on the left, 30 % of space between plot and facet strips
  theme(plot.title = element_text(hjust = 0),
        plot.subtitle = element_text(hjust = 0),
        plot.title.position = "plot",
        text = element_text(size = 12),
        axis.text.x = element_text(size = 10),
        axis.text.y = element_text(family = "serif", vjust = 0.3),
        strip.text.y = element_text(angle = 0),
        strip.background = element_rect(colour=NA),
        panel.border =  element_rect(colour=NA),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()#,
        #axis.ticks = element_blank()
  )
#theme(axis.text.y = element_text(colour=rev(spec_fact_dyn_df$colour)))
p_rel
# add different colours for facet strips:
g <- ggplot_gtable(ggplot_build(p_rel))
stripr <- which(grepl('strip-r', g$layout$name))
fills <- c("#d7191c","#abd9e9","#2c7bb6")
k <- 1
for (i in stripr) {
  j <- which(grepl('rect', g$grobs[[i]]$grobs[[1]]$childrenOrder))
  g$grobs[[i]]$grobs[[1]]$children[[j]]$gp$fill <- fills[k]
  k <- k+1
}
grid::grid.draw(g)



### sort by rel. change in overall occ. dyn., compare categories based on 100 draws of posterior: ----

# for sorting categories:

spec_fact_dyn_df_opt3 <- flow_df %>%
  filter(species %in% spec_okay) %>% 
  mutate(rel_change_cfclim = case_when(
    cfclim == "increase" & fact == "stable" ~ "relative\ndecrease",
    cfclim == "increase" & fact == "decrease" ~ "absolute\ndecrease",
    cfclim == "stable" & fact == "increase" ~ "absolute\nincrease",
    cfclim == "stable" & fact == "decrease" ~ "absolute\ndecrease",
    cfclim == "decrease" & fact == "increase" ~ "absolute\nincrease",
    cfclim == "decrease" & fact == "stable" ~ "relative\nincrease",
    .default = "no change"),
    rel_change_cflu = case_when(
      cflu == "increase" & fact == "stable" ~ "relative\ndecrease",
      cflu == "increase" & fact == "decrease" ~ "absolute\ndecrease",
      cflu == "stable" & fact == "increase" ~ "absolute\nincrease",
      cflu == "stable" & fact == "decrease" ~ "absolute\ndecrease",
      cflu == "decrease" & fact == "increase" ~ "absolute\nincrease",
      cflu == "decrease" & fact == "stable" ~ "relative\nincrease",
      .default = "no change"),
    rel_change_cfclimlu = case_when(
      cfclimlu == "increase" & fact == "stable" ~ "relative\ndecrease",
      cfclimlu == "increase" & fact == "decrease" ~ "absolute\ndecrease",
      cfclimlu == "stable" & fact == "increase" ~ "absolute\nincrease",
      cfclimlu == "stable" & fact == "decrease" ~ "absolute\ndecrease",
      cfclimlu == "decrease" & fact == "increase" ~ "absolute\nincrease",
      cfclimlu == "decrease" & fact == "stable" ~ "relative\nincrease",
      .default = "no change")) %>%
  mutate(across(rel_change_cfclim:rel_change_cfclimlu,
                ~ factor(.x, levels = c("absolute\nincrease", "relative\nincrease", "no change", "relative\ndecrease", "absolute\ndecrease")))) %>% 
  # choose scenario:
  #arrange(rel_change_cfclim , species) %>% 
  arrange(rel_change_cflu, species) %>% 
  #arrange(rel_change_cfclimlu, species) %>%  
  # add species orders etc. as sorting options:
  left_join(bbs_dt %>% select(English_Common_Name, Scientific_Name, ORDER, Family) %>% distinct,
            by = c(species = "English_Common_Name")) %>% 
  mutate(ORDER2 = factor(ORDER, levels = names(sort(table(ORDER), decreasing = TRUE))))

# relative influence values:
plot_df <- rel_infl_df2 %>%
  select(c(species, starts_with("infl"))) %>%
  left_join(spec_fact_dyn_df_opt3) %>% 
  tidyr::pivot_longer(cols = starts_with("infl"), names_to = "scenario", values_to = "value") %>% 
  mutate(scenario = factor(scenario, levels = c("infl_clim", "infl_lu", "infl_climlu"))) %>% 
  mutate(scenario = recode(scenario, infl_clim = "climate", infl_lu = "land use", infl_climlu = "climate + land use")) %>% 
  # choose scenario to plot: 
  #arrange(rel_change_cfclim, ORDER2, Family) %>%
  arrange(rel_change_cflu, ORDER2, Family) %>%
  #arrange(rel_change_cfclimlu , ORDER2, Family) %>%
  mutate(plot_order = row_number())

# plot:
p_rel <- plot_df %>% 
  # choose scenario to plot: 
  #filter(scenario == "climate") %>% 
  filter(scenario == "land use") %>% 
  #filter(scenario == "climate + land use") %>% 
  mutate(species = as.factor(species)) %>% 
  ggplot(aes(x = scenario, y = forcats::fct_reorder(species, desc(plot_order)), fill = value)) + 
  geom_tile() +
  # choose scenario to plot: 
  #facet_grid(rows = vars(rel_change_cfclim), scales = "free", space = "free") +
  facet_grid(rows = vars(rel_change_cflu), scales = "free", space = "free") +
  #facet_grid(rows = vars(rel_change_cfclimlu), scales = "free", space = "free") +
  viridis::scale_fill_viridis(discrete=FALSE) +
  theme_bw() +
  # choose scenario to plot:
  #labs(title = "Relative influence of climate change on occupancy dynamics") +
  labs(title = "Relative influence of land use change on occupancy dynamics") +
  #labs(title = "Relative influence of climate and land use change on occupancy dynamics") +
  labs(subtitle = "Species ordered by relative occupancy dynamics, order and family") +
  xlab("") + ylab("") +
  labs(fill = expression(MAPE [counterfactual] - MAPE [factual])) +
  scale_x_discrete(expand = expansion(mult = c(0, 0.55))) + # 0 space on the left, 30 % of space between plot and facet strips
  theme(plot.title = element_text(hjust = 0),
        plot.subtitle = element_text(hjust = 0),
        plot.title.position = "plot",
        text = element_text(size = 12),
        axis.text.x = element_text(size = 10),
        axis.text.y = element_text(family = "serif", vjust = 0.3),
        strip.text.y = element_text(angle = 0),
        strip.background = element_rect(colour=NA),
        panel.border =  element_rect(colour=NA),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()
  )
p_rel
# add different colours for facet strips:
g <- ggplot_gtable(ggplot_build(p_rel))
stripr <- which(grepl('strip-r', g$layout$name))
fills <- c("#d7191c","gray90", "#abd9e9","#2c7bb6")
k <- 1
for (i in stripr) {
  j <- which(grepl('rect', g$grobs[[i]]$grobs[[1]]$childrenOrder))
  g$grobs[[i]]$grobs[[1]]$children[[j]]$gp$fill <- fills[k]
  k <- k+1
}
grid::grid.draw(g)


#"relative increase" = "#fdae61",









# MISC ---


#testspecs <- c("Brewer's Sparrow", "Cassin's Sparrow", "Brown Creeper", "Olive-sided Flycatcher", "Bobolink")
#spec <- testspecs[5]
spec <- final_species[1]
spec <- "Wilson's Warbler"
spec <- "White-winged Dove"
# observations time series:
load(file.path(obs_dir, paste0(spec, "_obs_ts_sum_occ_routes.RData")))
ts_obs

# time series predictions factual:
load(file.path(fact_dir, paste0(spec, "_ts_sum_occ_routes_f_preds.RData")))
ts_preds_fact

# time series predictions counterfactual:
load(file.path(cfact_dir, paste0(spec, "_ts_sum_occ_routes_cf_preds.RData")))
ts_preds_cfact

# plot:

ggplot() +
  geom_point(aes(x = year, y = Npres), data = ts_obs) +
  geom_line(aes(x = year, y = Npres, colour = "black"), data = ts_obs) +
  geom_line(aes(x = year, y = median_Nocc_f, colour = "#85CB33"), data = ts_preds_fact, linewidth = 0.8) +
  #geom_ribbon(aes(x = year, ymax = CI90high_f, ymin = CI90low_f ), data =  ts_preds_fact,
  #                           alpha = 0.2, fill = "#85CB33") +
  geom_line(aes(x = year, y = median_Nocc_cf_clim, colour = "#0D98BA"), data = ts_preds_cfact, linetype = "dashed", linewidth = 0.8) +
  #geom_ribbon(aes(x = year, ymax = CI90high_cf_clim, ymin = CI90low_cf_clim), data = ts_preds_cfact,
  #            alpha = 0.2, fill = "#0D98BA") +
  geom_line(aes(x = year, y = median_Nocc_cf_1995soc, colour = "#B7410E"), data = ts_preds_cfact, linetype = "dashed", linewidth = 0.8) +
  #geom_ribbon(aes(x = year, ymax = CI90high_cf_1995soc, ymin = CI90low_cf_1995soc), data = ts_preds_cfact,
  #            alpha = 0.2, fill = "#B7410E") +
  geom_line(aes(x = year, y = median_Nocc_cf_clim_1995soc, colour = "#046865"), data = ts_preds_cfact, linetype = "dashed", linewidth = 0.8) +
  #geom_ribbon(aes(x = year, ymax = CI90high_cf_clim_1995soc, ymin = CI90low_cf_clim_1995soc), data = ts_preds_cfact,
  #            alpha = 0.2, fill = "#046865") +
  #geom_smooth(aes(x = year, y = median_Nocc_f), data = ts_preds_fact, colour = "#85CB33", method = "lm") +
  theme_bw() +
  ggtitle(spec) +
  xlab("") +
  ylab("N routes with presence") +
  scale_color_identity(name = "",
                     breaks = c("black", "#85CB33","#0D98BA", "#B7410E", "#046865"),
                     labels = c("observations", "factual", "cfact climate", "cfact land use", "cfact climate + land use"),
                     guide = "legend") +
  guides(color = guide_legend(nrow = 2)) +
  theme(text = element_text(size = 18), legend.position = "bottom") +
  ylim(c(15,50))



# calculate different metrics, evaluate how well they match visual impression

# mean absolute percent error:
mape_fact <- Metrics::mape(actual = ts_obs$Npres, 
              predicted = ts_preds_fact$median_Nocc_f)

sum(abs(ts_obs$Npres - ts_preds_fact$median_Nocc_f)/abs(ts_obs$Npres)) / 25
(ts_obs$Npres - ts_preds_cfact$median_Nocc_cf_clim)/ts_obs$Npres

# test for all species whether for single years mape is > abs(100 %)
mape_fact
mape_cfclim <- Metrics::mape(actual = ts_obs$Npres, 
                           predicted = ts_preds_cfact$median_Nocc_cf_clim)
mape_cfclim
mape_cflu <- Metrics::mape(actual = ts_obs$Npres, 
                             predicted = ts_preds_cfact$median_Nocc_cf_1995soc)
mape_cflu
mape_cfclimlu <- Metrics::mape(actual = ts_obs$Npres, 
                             predicted = ts_preds_cfact$median_Nocc_cf_clim_1995soc)
mape_cfclimlu

# change in error compared to factual:
mape_cfclim - mape_fact # problem, wenn cf besser als fact ist, neg. wenn besser, bei 0 kappen und sagen kein Effekt?
mape_cflu - mape_fact
mape_cfclimlu - mape_fact

# einheit mape?
# percent error that comes on top when using counterfactual data?

# how to include direction of effect?
# other metrics

# check whether I have cases where error for a year is > abs(100%)
# where MAPE could be biased, as far as I understand:

spec_ape_greater100 <- c()

for(s in 1:length(final_species)){
  
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
  
  
  # absolute percentage error for each year:
  fact_ape <- abs(ts_obs$Npres - ts_preds_fact$median_Nocc_f)/abs(ts_obs$Npres)
  clim_ape <- abs(ts_obs$Npres - ts_preds_cfact$median_Nocc_cf_clim)/abs(ts_obs$Npres)
  lu_ape <- abs(ts_obs$Npres - ts_preds_cfact$median_Nocc_cf_1995soc)/abs(ts_obs$Npres)
  climlu_ape <- abs(ts_obs$Npres - ts_preds_cfact$median_Nocc_cf_clim_1995soc)/abs(ts_obs$Npres)
  
  if(any(fact_ape > 1) | any(clim_ape > 1) | any(lu_ape > 1) | any(climlu_ape > 1)){
    
    print(spec)
    spec_ape_greater100 <- c(spec_ape_greater100, spec)
    
  }
  
}
spec_ape_greater100
spec_ape_greater100[which(spec_ape_greater100 %in% spec_okay)]
