# Script:   5_2_attribution_plots_trend_categories.R
# Purpose:  Generate plots to compare factual and counterfactual occupancy trends
# Inputs:   results/species_DOM_val_okay.RData
#           results/attribution/attribution_metrics_final.RData
# Outputs:  results/attribution/trend_categories.RData
#           plots/attribution/barplot_trend_cats_stacked.svg (Fig. 2)
# Runs on:  Local

# Steps:
# define categories of absolute and relative winners and losers of change
# Chi-square test of categories
# barplots
# (alternative: alluvial plots)

source(file.path("scripts", "0_paths.R"))


# packages: --------------------------------------------------------------------

library(dplyr)
library(ggplot2)
library(ggalluvial) # alluvial plots
library(ggnewscale) # for multiple fill scales in one plot
library(cowplot)
library(rphylopic)
library(grid)


# load data: -------------------------------------------------------------------

# species for attribution:
load(file = file.path(dir, "results", "species_DOM_val_okay.RData")) # output of 4_0_DOMs_predictions_y_routes_scenarios.R
spec_attr

# attribution metrics:
load(file = file.path(dir, "results", "attribution", "attribution_metrics_final.RData")) # output of 5_1_attribution_metrics.R


# reformat data: ---------------------------------------------------------------

flow_df <- attr_metr_df %>% 
  # categorize relative impact (based on Langhammer et al. 2024):
  select(c(species, starts_with("slope"), starts_with("p_"))) %>% 
  ## factual vs. counterfactual climate:
  mutate(trend_change_clim = case_when(slope_fact > slope_cfclim & slope_fact > 0 & p_fact < 0.05 ~ "absolute climate change winner",
                                       slope_fact > slope_cfclim & slope_fact > 0 & p_fact >= 0.05 ~ "relative climate change winner",
                                       slope_fact > slope_cfclim & slope_fact < 0 ~ "relative climate change winner",
                                       slope_fact < slope_cfclim & slope_fact > 0 ~ "relative climate change loser",
                                       slope_fact < slope_cfclim & slope_fact < 0 & p_fact >= 0.05 ~ "relative climate change loser",
                                       slope_fact < slope_cfclim & slope_fact < 0 & p_fact < 0.05 ~ "absolute climate change loser",
                                       .default = NA)) %>% 
  ## factual vs. counterfactual land use:
  mutate(trend_change_lu = case_when(slope_fact > slope_cflu & slope_fact > 0 & p_fact < 0.05 ~ "absolute land use change winner",
                                     slope_fact > slope_cflu & slope_fact > 0 & p_fact >= 0.05 ~ "relative land use change winner",
                                     slope_fact > slope_cflu & slope_fact < 0 ~ "relative land use change winner",
                                     slope_fact < slope_cflu & slope_fact > 0 ~ "relative land use change loser",
                                     slope_fact < slope_cflu & slope_fact < 0 & p_fact >= 0.05 ~ "relative land use change loser",
                                     slope_fact < slope_cflu & slope_fact < 0 & p_fact < 0.05 ~ "absolute land use change loser",
                                     .default = NA)) %>% 
  ## factual vs. counterfactual climate + land use:
  mutate(trend_change_climlu = case_when(slope_fact > slope_cfclimlu & slope_fact > 0 & p_fact < 0.05 ~ "absolute global change winner",
                                         slope_fact > slope_cfclimlu & slope_fact > 0 & p_fact >= 0.05 ~ "relative global change winner",
                                         slope_fact > slope_cfclimlu & slope_fact < 0 ~ "relative global change winner",
                                         slope_fact < slope_cfclimlu & slope_fact > 0 ~ "relative global change loser",
                                         slope_fact < slope_cfclimlu & slope_fact < 0 & p_fact >= 0.05 ~ "relative global change loser",
                                         slope_fact < slope_cfclimlu & slope_fact < 0 & p_fact < 0.05 ~ "absolute global change loser",
                                         .default = NA)) %>% 
  ## add "no change"-category if confidence intervals of slopes overlap (= smaller max. is larger than larger min.):
  mutate(stable_cfclim = pmin(slope_CIhigh_fact, slope_CIhigh_cfclim) > pmax(slope_CIlow_fact, slope_CIlow_cfclim),
         stable_cflu = pmin(slope_CIhigh_fact, slope_CIhigh_cflu) > pmax(slope_CIlow_fact, slope_CIlow_cflu),
         stable_cfclimlu = pmin(slope_CIhigh_fact, slope_CIhigh_cfclimlu) > pmax(slope_CIlow_fact, slope_CIlow_cfclimlu)) %>%
  mutate(trend_change_clim = ifelse(stable_cfclim, "no change", trend_change_clim),
         trend_change_lu = ifelse(stable_cflu, "no change", trend_change_lu),
         trend_change_climlu = ifelse(stable_cfclimlu, "no change", trend_change_climlu)) %>%
  select(c(-starts_with("stable"))) %>% 
  # convert to factors:
  mutate(trend_change_clim = factor(trend_change_clim, levels = c("absolute climate change winner", "relative climate change winner", 
                                                                  "no change", "relative climate change loser", "absolute climate change loser")),
         trend_change_lu = factor(trend_change_lu, levels = c("absolute land use change winner", "relative land use change winner", 
                                                              "no change", "relative land use change loser", "absolute land use change loser")),
         trend_change_climlu = factor(trend_change_climlu, levels = c("absolute global change winner", "relative global change winner", 
                                                                      "no change", "relative global change loser", "absolute global change loser"))) %>% 
  # categorize trend within each scenario:
  select(-c(matches("_CI"))) %>% 
  tidyr::pivot_longer(cols = matches("(slope_)|(p_)"), names_to = c("metric", "scenario"), 
                      values_to = "value", names_pattern = "(.*)_(.*)") %>% 
  tidyr::pivot_wider(names_from = metric, values_from = value) %>% 
  mutate(dynamics = case_when(slope > 0 & p < 0.05 ~ "positive\ntrend",
                              slope < 0 & p < 0.05 ~ "negative\ntrend",
                              p >= 0.05 ~ "stable",
                              .default = NA),
         dynamics = factor(dynamics, levels = c("positive\ntrend", "stable", "negative\ntrend"))) %>% 
  select(-c(slope, p)) %>% 
  tidyr::pivot_wider(names_from = scenario, values_from = dynamics)

#save(flow_df, file =  file.path(dir, "results", "attribution", "trend_categories.RData"))


# Chi-square test against equal expected classes: ------------------------------

clim_test <- flow_df %>% 
  group_by(trend_change_clim) %>% 
  summarise(n = n())
clim_test$n
chisq.test(clim_test$n, p = rep(0.2, 5))

lu_test <- flow_df %>% 
  count(trend_change_lu, .drop = FALSE)
lu_test$n
chisq.test(lu_test$n, p = rep(0.2, 5))

climlu_test <- flow_df %>% 
  group_by(trend_change_climlu) %>% 
  summarise(n = n())
climlu_test$n
chisq.test(climlu_test$n, p = rep(0.2, 5))


# barplots of number of winner and loser species under different scenarios: ----

bp_dt <- flow_df %>% 
  mutate(fact = forcats::fct_recode(fact, "negative\nN = 45" = "negative\ntrend",
                                    "positive\nN = 31" = "positive\ntrend",
                                    "stable\nN = 4" = "stable")) %>% 
  select(species, fact, trend_change_clim, trend_change_lu, trend_change_climlu) %>% 
  # trend categories: absolute / relative winner / loser without driver:
  mutate(trend_change_clim = gsub(pattern = " climate change", replacement = "", x = trend_change_clim),
         trend_change_lu = gsub(pattern = " land use change", replacement = "", x = trend_change_lu),
         trend_change_climlu = gsub(pattern = " global change", replacement = "", x = trend_change_climlu)) %>% 
  # convert to factor:
  mutate(trend_change_clim = factor(trend_change_clim, levels = c("absolute winner", "relative winner", "no change", "relative loser", "absolute loser")),
         trend_change_lu = factor(trend_change_lu, levels = c("absolute winner", "relative winner", "no change", "relative loser", "absolute loser")),
         trend_change_climlu  = factor(trend_change_climlu , levels = c("absolute winner", "relative winner", "no change", "relative loser", "absolute loser"))) %>% 
  # convert to long format:
  tidyr::pivot_longer(cols = starts_with("trend"), names_to = "scenario", values_to = "category") %>% 
  mutate(scenario = factor(scenario, levels = c("trend_change_clim", "trend_change_lu", "trend_change_climlu")),
         scenario = recode(scenario, trend_change_clim = "A: Climate change", trend_change_lu = "B: Land use change", trend_change_climlu = "C: Climate & land use change"))

bp_dt2 <- bp_dt %>% 
  group_by(fact, scenario, category) %>% 
  summarise(n = n()) 


# stacked barplots:
bp_trends_stacked <- ggplot(bp_dt2, aes(x = fact, y = n, group = scenario)) +
  geom_bar(aes(fill = category), stat = "identity", position = "stack", width = 0.6) +
  geom_text(aes(label = n, group = fact), 
             size = 4, colour = "grey20", position = position_stack(vjust = 0.5)) +
  facet_wrap(~scenario) +
  scale_fill_manual(values = c("relative loser" = "#abd9e9", #,
                                 "absolute loser" =  "#2c7bb6",#",
                                 "no change" = "#B18985", #, #"gray80",
                                 "absolute winner" = "#d7191c",
                                 "relative winner" = "#fdae61"), 
                      drop = FALSE, na.value = NA, name = "trend change") +
  scale_y_continuous(expand = expansion(mult = c(0, 0.15))) +
  labs(y = "N species", x = "factual occupancy trend") +
  theme_bw() +
  theme(strip.background = element_blank(),
        text = element_text(size = 14),
        axis.title.x = element_text(margin = margin(t = 10), size = 13),
        axis.title.y = element_text(margin = margin(r = 10), size = 13),
        legend.position = "bottom",
        legend.title.position = "top",
        legend.margin = margin(t = -5, r = 0, b = 0, l = 0),
        strip.text = element_text(size = 13))


ggsave(filename = file.path(dir, "plots", "attribution", "barplot_trend_cats_stacked.svg"), 
       plot = bp_trends_stacked,
       device = "svg",
       width = 21,
       height = 11, # A4
       units = "cm")

# (xx maybe add bird icon)


# alternative: alluvial plots: -------------------------------------------------

# https://r-charts.com/flow/ggalluvial/
# change categories based on linear trend

## climate change: ----

plot_df_clim <- flow_df %>%
  group_by(fact, cfclim, trend_change_clim) %>% 
  summarise(n = n()) %>% 
  arrange(cfclim, trend_change_clim) # change plotting order

# version for presentation:

# jpeg(file = file.path(dir, "plots", "attribution", "alluvial_climate_presentation.jpg"),
#      width = 910, height = 1200, quality = 100)
# svg(file = file.path(dir, "plots", "attribution", "alluvial_climate_presentation.svg"),
#     width = 14, height = 18)

plot_df_clim %>% 
  ggplot(aes(axis1 = cfclim, axis2 = fact, y = n)) +
  geom_alluvium(aes(fill = trend_change_clim), show.legend = TRUE, width = 1/5, colour = "grey50", alpha = 0.6) +
  scale_fill_manual(values = c("relative climate change loser" = "#abd9e9",
                               "absolute climate change loser" =  "#2c7bb6",
                               "no change" = "gray90",
                               "absolute climate change winner" = "#d7191c",
                               "relative climate change winner" = "#fdae61"), drop = FALSE, na.value = NA, name = "impact of drivers") +
  new_scale_fill() +
  scale_fill_manual(values = c("positive\ntrend" = "#EFCA08", "negative\ntrend" =  "#7DC4C9", "stable" =  "#EDE6F2"), na.value = NA) +
  geom_stratum(width = 1/5, aes(fill = cfclim), show.legend = FALSE, colour = "grey20") + 
  geom_stratum(width = 1/5, aes(fill = fact), show.legend = FALSE, colour = "grey20") +
  geom_text(stat = "stratum", 
            fontface = "italic",
            aes(label = ifelse(test = after_stat(x) == "1",
                               yes = paste0(stratum, " (",
                                            flow_df %>%
                                              filter(species %in% spec_attr) %>% 
                                              group_by(cfclim) %>% 
                                              summarise(n = n()) %>% 
                                              arrange(cfclim) %>%  
                                              pull(n) %>% 
                                              rev, ")"),
                               no = paste0(stratum, " (",
                                           flow_df %>%
                                             filter(species %in% spec_attr) %>%
                                             group_by(fact) %>%
                                             summarise(n = n()) %>%
                                             arrange(as.character(fact)) %>%
                                             pull(n) %>%
                                             rev, ")"))),
            size = 9) +
  # species numbers in flow:
  geom_text(stat = "alluvium", 
            fontface = "italic",
            aes(label = ifelse(test = after_stat(x) == "1",
                               yes = plot_df_clim %>% pull(n),
                               no = "")),
            size = 8, nudge_x = 0.15) +
  scale_x_discrete(limits = c("cfclim", "fact"), expand = c(0.13, 0.0)) +
  labs(title = "Change in area of occupancy",
       subtitle = "numbers = number of species",
       x = "", y = "") +
  annotate("text", x = 1, y = 85, label = "Counterfactual\nclimate", hjust = 0.5, size = 10) + # y = 85 for spec_attr
  annotate("text", x = 2, y = 85, label = "Observed\nclimate change", hjust = 0.5, size = 10) + # size = 9 for paper version
  guides(fill = "none") +
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        axis.ticks = element_blank(), axis.text = element_blank(),
        text = element_text(size = 35), # size = 28 for paper version
        plot.subtitle = element_text(margin = margin(0,0,30,0), size = 30),
        legend.text = element_text(size = 30),
        legend.position = "bottom", legend.direction = "vertical", legend.box.spacing = unit(-50, "pt"),
        legend.key.spacing.y = unit(0.2, "cm"), 
        legend.title = element_text(margin=margin(0,0,10,0)))

#dev.off()

# version for manuscript:
# climate and land use change side-by-side, common legend

clim_plot <- plot_df_clim %>% 
  ggplot(aes(axis1 = cfclim, axis2 = fact, y = n)) +
  geom_alluvium(aes(fill = trend_change_clim), show.legend = TRUE, width = 1/5, colour = "grey50", alpha = 0.6) +
  scale_fill_manual(values = c("relative climate change loser" = "#abd9e9",
                               "absolute climate change loser" =  "#2c7bb6",
                               "no change" = "gray90",
                               "absolute climate change winner" = "#d7191c",
                               "relative climate change winner" = "#fdae61"), 
                    drop = FALSE, na.value = NA) +
  new_scale_fill() +
  scale_fill_manual(values = c("positive\ntrend" = "#EFCA08", "negative\ntrend" =  "#7DC4C9", "stable" =  "#EDE6F2"), na.value = NA) +
  geom_stratum(width = 1/4, aes(fill = cfclim), show.legend = FALSE, colour = "grey20") + 
  geom_stratum(width = 1/4, aes(fill = fact), show.legend = FALSE, colour = "grey20") +
  # species numbers in bars:
  geom_text(stat = "stratum", 
            fontface = "italic",
            aes(label = ifelse(test = after_stat(x) == "1",
                               yes = paste0(stratum, " (",
                                            flow_df %>%
                                              filter(species %in% spec_attr) %>% 
                                              group_by(cfclim) %>% 
                                              summarise(n = n()) %>% 
                                              arrange(cfclim) %>%  
                                              pull(n) %>% 
                                              rev, ")"),
                               no = paste0(stratum, " (",
                                           flow_df %>%
                                             filter(species %in% spec_attr) %>%
                                             group_by(fact) %>%
                                             summarise(n = n()) %>%
                                             arrange(as.character(fact)) %>%
                                             pull(n) %>%
                                             rev, ")"))),
            size = 8) +
  # species numbers in flow:
  geom_text(stat = "alluvium", 
            fontface = "italic",
            aes(label = ifelse(test = after_stat(x) == "1",
                               yes = plot_df_clim %>% pull(n),
                               no = "")),
            size = 7, nudge_x = 0.15) +
  scale_x_discrete(limits = c("cfclim", "fact"), expand = c(0.13, 0.0)) +
  labs(#title = "Change in area of occupancy",
       #subtitle = "numbers = number of species",
       x = "", y = "") +
  annotate("text", x = 1, y = 85, label = "Counterfactual\nclimate", hjust = 0.5, size = 9) + # y = 85 for spec_attr
  annotate("text", x = 2, y = 85, label = "Factual\nclimate", hjust = 0.5, size = 9) + # size = 9 for paper version
  coord_cartesian(clip = 'off') + # prevent annotations to get cut
  guides(fill = "none") +
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        plot.margin = margin(40,15,0,5), # trbl # top margin to later have label when arranging multiple plots
        panel.border = element_blank(),
        axis.ticks = element_blank(), axis.text = element_blank(),
        text = element_text(size = 24),
        plot.subtitle = element_text(margin = margin(0,0,30,0), size = 30),
        legend.text = element_text(size = 24),
        legend.position = "bottom", legend.direction = "vertical", 
        legend.box.margin = margin(-50, 0, 0, 0),
        legend.key.spacing.y = unit(0.2, "cm"), 
        legend.title=element_blank())

# ggsave(filename = file.path(dir, "plots", "attribution", "alluvial_climate_manuscript.svg"),
#        plot = clim_plot,
#        device = "svg",
#        width = 21,
#        height = 29.7, # A4
#        units = "cm")


## land use change: ----

plot_df_lu <- flow_df %>%
  group_by(fact, cflu, trend_change_lu) %>% 
  summarise(n = n()) %>% 
  arrange(cflu, fact, trend_change_lu) # change plotting order


# with lodes form, seems more stable:
plot_df_lodes <- to_lodes_form(plot_df_lu, axes = c(1:2), id = trend_change_lu ) %>% 
  mutate(x = factor(x, levels = c("cflu", "fact"))) %>% 
  # same order as in plot_df_lu:
  mutate(trend_change_lu2 = rep(plot_df_lu$trend_change_lu, 2))


# version for presentation:

# jpeg(file = file.path(dir, "plots", "attribution", "alluvial_lu_presentation.jpg"),
#      width = 920, height = 1200, quality = 100)
# svg(file = file.path(dir, "plots", "attribution", "alluvial_lu_presentation.svg"),
#     width = 14, height = 18)

ggplot(plot_df_lodes, aes(x = x, stratum = stratum, alluvium = trend_change_lu, y = n)) +
  geom_flow(stat = "alluvium", aes(fill = trend_change_lu2),
            color = "darkgray",width = 1/5, show.legend = TRUE, alpha = 0.6) +
  scale_fill_manual(values = c("relative land use change loser" = "#abd9e9",
                               "absolute land use change loser" =  "#2c7bb6",
                               "no change" = "gray90",
                               "absolute land use change winner" = "#d7191c",
                               "relative land use change winner" = "#fdae61"), 
                    drop = FALSE, na.value = NA, name = "impact of drivers") +
  new_scale_fill() +
  geom_stratum(aes(fill = stratum), width = 1/5, colour = "grey20") + 
  scale_fill_manual(values = c("positive\ntrend" = "#EFCA08", "negative\ntrend" =  "#7DC4C9", "stable" =  "#EDE6F2"), na.value = NA) +
  # text in bars:
  geom_text(stat = "stratum", 
            fontface = "italic",
            size = 9,
            aes(label = ifelse(test = after_stat(x) == "1",
                               yes = paste0(stratum, " (",
                                            plot_df_lodes %>% 
                                              filter(x == "cflu") %>% 
                                              group_by(stratum) %>% 
                                              summarise(n = sum(n)) %>% 
                                              pull(n) %>% 
                                              rev, ")"),
                               no = paste0(stratum, " (",
                                           plot_df_lodes %>% 
                                             filter(x == "fact") %>% 
                                             group_by(stratum) %>% 
                                             summarise(n = sum(n)) %>% 
                                             pull(n) %>% 
                                             rev, ")")))) +
  # species numbers in flow:
  geom_text(stat = "alluvium", 
            fontface = "italic",
            aes(label = ifelse(test = after_stat(x) == "1",
                               yes = plot_df_lodes %>% filter(x == "cflu") %>% pull(n),
                               no = "")),
            size = 8, nudge_x = 0.15) +
  scale_x_discrete(limits = c("cflu", "fact"), expand = c(0.13, 0.0)) +
  labs(title = "Change in area of occupancy",
       subtitle = "numbers = number of species",
       x = "", y = "") +
  annotate("text", x = 1, y = 85, label = "Counterfactual\nland use", hjust = 0.5, size = 10) + 
  annotate("text", x = 2, y = 85, label = "Observed\nland use change", hjust = 0.5, size = 10) +
  guides(fill = "none") +
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        axis.ticks = element_blank(), axis.text = element_blank(),
        text = element_text(size = 35), # size = 28 for paper version
        plot.subtitle = element_text(margin = margin(0,0,30,0), size = 30),
        legend.text = element_text(size = 30),
        legend.position = "bottom", legend.direction = "vertical", legend.box.spacing = unit(-50, "pt"),
        legend.key.spacing.y = unit(0.2, "cm"), 
        legend.title = element_text(margin=margin(0,0,10,0)))

#dev.off()


# version for manuscript:

lu_plot <- ggplot(plot_df_lodes, aes(x = x, stratum = stratum, alluvium = trend_change_lu, y = n)) +
  geom_flow(stat = "alluvium", aes(fill = trend_change_lu2),
            color = "darkgray",width = 1/4, show.legend = TRUE, alpha = 0.6) +
  scale_fill_manual(values = c("relative land use change loser" = "#abd9e9",
                               "absolute land use change loser" =  "#2c7bb6",
                               "no change" = "gray90",
                               "absolute land use change winner" = "#d7191c",
                               "relative land use change winner" = "#fdae61"), 
                    drop = FALSE, na.value = NA, name = "impact of drivers") +
  new_scale_fill() +
  geom_stratum(aes(fill = stratum), width = 1/4, colour = "grey20") + 
  scale_fill_manual(values = c("positive\ntrend" = "#EFCA08", "negative\ntrend" =  "#7DC4C9", "stable" =  "#EDE6F2"), na.value = NA) +
  # species numbers in bars:
  geom_text(stat = "stratum", 
            fontface = "italic",
            size = 8,
            aes(label = ifelse(test = after_stat(x) == "1",
                               yes = paste0(stratum, " (",
                                            plot_df_lodes %>% 
                                              filter(x == "cflu") %>% 
                                              group_by(stratum) %>% 
                                              summarise(n = sum(n)) %>% 
                                              pull(n) %>% 
                                              rev, ")"),
                               no = paste0(stratum, " (",
                                           plot_df_lodes %>% 
                                             filter(x == "fact") %>% 
                                             group_by(stratum) %>% 
                                             summarise(n = sum(n)) %>% 
                                             pull(n) %>% 
                                             rev, ")")))) +
  # species numbers in flow:
  geom_text(stat = "alluvium", 
            fontface = "italic",
            aes(label = ifelse(test = after_stat(x) == "1",
                               yes = plot_df_lu %>% pull(n),
                               no = "")),
            size = 7, nudge_x = 0.15) +
  scale_x_discrete(limits = c("cflu", "fact"), expand = c(0.13, 0.0)) +
  labs(#title = "Change in area of occupancy",
    #subtitle = "numbers = number of species",
    x = "", y = "") +
  annotate("text", x = 1, y = 85, label = "Counterfactual\nland use", hjust = 0.5, size = 9) + 
  annotate("text", x = 2, y = 85, label = "Factual\nland use", hjust = 0.5, size = 9) +
  coord_cartesian(clip = 'off') + # prevent annotations to get cut
  guides(fill = "none") +
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        plot.margin = margin(40,15,0,5), # trbl
        panel.border = element_blank(),
        axis.ticks = element_blank(), axis.text = element_blank(),
        text = element_text(size = 24),
        plot.subtitle = element_text(margin = margin(0,0,30,0), size = 30),
        legend.text = element_text(size = 24),
        legend.position = "bottom", legend.direction = "vertical", 
        legend.box.margin = margin(-50, 0, 0, 0),
        legend.key.spacing.y = unit(0.2, "cm"), 
        legend.title=element_blank())

# ggsave(filename = file.path(dir, "plots", "attribution", "alluvial_lu_manuscript.svg"),
#        plot = lu_plot,
#        device = "svg",
#        width = 21,
#        height = 29.7, # A4
#        units = "cm")


## climate + land use change: ----

plot_df_climlu <- flow_df %>%
  group_by(fact, cfclimlu, trend_change_climlu) %>% 
  summarise(n = n()) %>% 
  arrange(cfclimlu, fact, trend_change_climlu) # change plotting order


# version for presentation:

# jpeg(file = file.path(dir, "plots", "attribution", "alluvial_climlu_presentation.jpg"),
#      width = 900, height = 1200, quality = 100)

plot_df_climlu %>% 
  ggplot(aes(axis1 = cfclimlu, axis2 = fact, y = n)) +
  geom_alluvium(aes(fill = trend_change_climlu), show.legend = TRUE, width = 1/5, colour = "grey50", alpha = 0.6) +
  scale_fill_manual(values = c("relative global change loser" = "#abd9e9",
                               "absolute global change loser" =  "#2c7bb6",
                               "no change" = "gray90",
                               "absolute global change winner" = "#d7191c",
                               "relative global change winner" = "#fdae61"),
                    drop = FALSE, na.value = NA, name = "impact of drivers") +
  new_scale_fill() +
  scale_fill_manual(values = c("positive\ntrend" = "#EFCA08", "negative\ntrend" =  "#7DC4C9", "stable" =  "#EDE6F2"), na.value = NA) +
  geom_stratum(width = 1/5, aes(fill = fact), show.legend = FALSE, colour = "grey20") +
  geom_stratum(width = 1/5, aes(fill = cfclimlu), show.legend = FALSE, colour = "grey20") + # order of lines somehow important
  # text in bars:
  geom_text(stat = "stratum", 
            fontface = "italic",
            aes(label = ifelse(test = after_stat(x) == "1",
                               yes = paste0(stratum, " (",
                                            flow_df %>%
                                              filter(species %in% spec_attr) %>% 
                                              group_by(cfclimlu) %>% 
                                              summarise(n = n()) %>% 
                                              arrange(cfclimlu) %>%  
                                              pull(n) %>% 
                                              rev, ")"),
                               no = paste0(stratum, " (",
                                           flow_df %>%
                                             filter(species %in% spec_attr) %>%
                                             group_by(fact) %>%
                                             summarise(n = n()) %>%
                                             arrange(fact) %>%
                                             pull(n) %>% 
                                             rev, ")"))),
            size = 9) +
  # species numbers in flow:
  geom_text(stat = "alluvium",
            fontface = "italic",
            aes(label = ifelse(test = after_stat(x) == "1",
                               yes = plot_df_climlu %>% pull(n),
                               no = "")),
            size = 8, nudge_x = 0.15) +
  scale_x_discrete(limits = c("cfclimlu", "fact"), expand = c(0.13, 0.0)) +
  labs(title = "Change in area of occupancy",
       subtitle = "numbers = number of species",
       x = "", y = "") +
  annotate("text", x = 1, y = 85, label = "Counterfactual\nclimate & land use", hjust = 0.5, size = 10) +
  annotate("text", x = 2, y = 85, label = "Observed \nglobal change", hjust = 0.5, size = 10) + 
  guides(fill = "none") +
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        axis.ticks = element_blank(), axis.text = element_blank(),
        text = element_text(size = 35),
        plot.subtitle = element_text(margin = margin(0,0,30,0), size = 30),
        legend.text = element_text(size = 30),
        legend.position = "bottom", legend.direction = "vertical", legend.box.spacing = unit(-50, "pt"),
        legend.key.spacing.y = unit(0.2, "cm"), 
        legend.title = element_text(margin=margin(0,0,10,0)))

#dev.off()


# version for manuscript:

climlu_plot <- plot_df_climlu %>% 
  ggplot(aes(axis1 = cfclimlu, axis2 = fact, y = n)) +
  geom_alluvium(aes(fill = trend_change_climlu), show.legend = TRUE, width = 1/5, colour = "grey50", alpha = 0.6) +
  scale_fill_manual(values = c("relative global change loser" = "#abd9e9",
                               "absolute global change loser" =  "#2c7bb6",
                               "no change" = "gray90",
                               "absolute global change winner" = "#d7191c",
                               "relative global change winner" = "#fdae61"),
                    drop = FALSE, na.value = NA, name = "impact of drivers") +
  new_scale_fill() +
  scale_fill_manual(values = c("positive\ntrend" = "#EFCA08", "negative\ntrend" =  "#7DC4C9", "stable" =  "#EDE6F2"), na.value = NA) +
  geom_stratum(width = 1/4, aes(fill = fact), show.legend = FALSE, colour = "grey20") +
  geom_stratum(width = 1/4, aes(fill = cfclimlu), show.legend = FALSE, colour = "grey20") + # order of lines somehow important
  # text in bars:
  geom_text(stat = "stratum", 
            fontface = "italic",
            aes(label = ifelse(test = after_stat(x) == "1",
                               yes = paste0(stratum, " (",
                                            flow_df %>%
                                              filter(species %in% spec_attr) %>% 
                                              group_by(cfclimlu) %>% 
                                              summarise(n = n()) %>% 
                                              arrange(cfclimlu) %>%  
                                              pull(n) %>% 
                                              rev, ")"),
                               no = paste0(stratum, " (",
                                           flow_df %>%
                                             filter(species %in% spec_attr) %>%
                                             group_by(fact) %>%
                                             summarise(n = n()) %>%
                                             arrange(fact) %>%
                                             pull(n) %>% 
                                             rev, ")"))),
            size = 8) +
  # species numbers in flow:
  geom_text(stat = "alluvium", 
            fontface = "italic",
            aes(label = ifelse(test = after_stat(x) == "1",
                               yes = plot_df_climlu %>% pull(n),
                               no = "")),
            size = 7, nudge_x = 0.15) +
  scale_x_discrete(limits = c("cfclimlu", "fact"), expand = c(0.13, 0.0)) +
  labs(#title = "Change in area of occupancy",
    #subtitle = "numbers = number of species",
    x = "", y = "") +
  annotate("text", x = 1, y = 85, label = "Counterfactual\nclimate & land use", hjust = 0.5, size = 9) + 
  annotate("text", x = 2, y = 85, label = "Factual climate\n& land use", hjust = 0.5, size = 9) + 
  coord_cartesian(clip = 'off') + # prevent annotations to get cut
  guides(fill = "none") +
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        plot.margin = margin(0,25,10,25), # trbl
        panel.border = element_blank(),
        axis.ticks = element_blank(), axis.text = element_blank(),
        text = element_text(size = 24),
        plot.subtitle = element_text(margin = margin(0,0,30,0), size = 30),
        legend.text = element_text(size = 24),
        legend.position = "bottom", legend.direction = "vertical", legend.box.spacing = unit(-50, "pt"),
        legend.key.spacing.y = unit(0.2, "cm"), 
        legend.title=element_blank())

# ggsave(filename = file.path(dir, "plots", "attribution", "alluvial_climlu_manuscript.svg"), 
#        plot = climlu_plot,
#        device = "svg",
#        width = 21,
#        height = 29.7, # A4
#        units = "cm")


## arrange plots for manuscript with conceptual legend and bird icon: -----

# concept figure as legend:

dummy_df <- expand.grid(impact = factor(c("absolute winners", "relative winners", 
                                          "relative losers", "absolute losers"), 
                                        levels = c("absolute winners", "relative winners", 
                                                   "relative losers", "absolute losers")),
                        scenario = factor(c("             observed" ,"            counterfactual")), # added spaces to get text placed in the wider part of the shaded area
                        x = c(0, 1)) %>% 
  as_tibble() %>% 
  arrange(impact) %>% 
  mutate(ymin = c(0, 0, 2, -2, # abs. winners
                  0, 0, -1, -6,  # rel. winners
                  0, 0, 1, 6, # rel. losers
                  0, 0, -6, -2), # abs. losers
         ymax = c(0, 0, 5, 1, # abs. winners
                  0, 0, -4, -9, # rel. winners
                  0, 0, 4, 9, # rel. losers
                  0, 0, -3, 1)) %>% # abs. losers
  # for label placement:
  mutate(ymean = ymin + (ymax-ymin)/2)

# plot:
conc_plot_m <- ggplot(dummy_df) +
  geom_ribbon(aes(x = x, ymin = ymin, ymax = ymax, fill = scenario)) +
  facet_wrap(~impact, ncol = 4, # change for horizontal version
             axes = "all", axis.labels = "all") +
  geom_hline(yintercept = 0, linetype = "longdash") +
  
  geomtextpath::geom_textpath(aes(x = x, y = ymean, label = scenario, colour = scenario), 
                              size = 8, text_only = TRUE) + # requires some manual adjustments
  
  scale_colour_manual(values = c("white", "black"), guide = "none") +
  scale_x_continuous(expand = c(0, 0)) +
  scale_fill_manual(values = c(alpha("grey10", 0.6), alpha("grey60", 0.6)), guide = "none" ) +
  labs(x = "time", y = "occupied area") +
  theme_bw() +
  theme(axis.ticks = element_blank(),
        axis.text = element_blank(),
        panel.grid = element_blank(),
        strip.background = element_rect(fill = NA),# colour added later
        text = element_text(size = 24),
        strip.text.x = element_text(size = 24),
        plot.margin = margin(5,5,5,5)
  ) 
conc_plot_m

# add colours for facet strips:
g <- ggplot_gtable(ggplot_build(conc_plot_m))
stripr <- which(grepl('strip', g$layout$name))
fills <- c(
  alpha("#d7191c", alpha = 0.6),
  alpha("#fdae61", alpha = 0.6),
  alpha("#abd9e9", alpha = 0.6),
  alpha("#2c7bb6", alpha = 0.6)
)

k <- 1
for (i in stripr) {
  j <- which(grepl('rect', g$grobs[[i]]$grobs[[1]]$childrenOrder))
  g$grobs[[i]]$grobs[[1]]$children[[j]]$gp$fill <- fills[k]
  k <- k+1
}
grid::grid.draw(g)


# climate change:

clim_plot <- plot_df_clim %>% 
  ggplot(aes(axis1 = cfclim, axis2 = fact, y = n)) +
  geom_alluvium(aes(fill = trend_change_clim), show.legend = TRUE, width = 1/5, colour = "grey50", alpha = 0.6) +
  scale_fill_manual(values = c("relative climate change loser" = "#abd9e9",
                               "absolute climate change loser" =  "#2c7bb6",
                               "no change" = "gray90",
                               "absolute climate change winner" = "#d7191c",
                               "relative climate change winner" = "#fdae61"), 
                    drop = FALSE, na.value = NA, guide = "none") +
  new_scale_fill() +
  scale_fill_manual(values = c("positive\ntrend" = "#EFCA08", "negative\ntrend" =  "#7DC4C9", "stable" =  "#EDE6F2"), na.value = NA, guide = "none") +
  geom_stratum(width = 1/4, aes(fill = cfclim), show.legend = FALSE, colour = "grey20") + 
  geom_stratum(width = 1/4, aes(fill = fact), show.legend = FALSE, colour = "grey20") +
  # species numbers in bars:
  geom_text(stat = "stratum", 
            fontface = "italic",
            aes(label = ifelse(test = after_stat(x) == "1",
                               yes = paste0(stratum, " (",
                                            flow_df %>%
                                              filter(species %in% spec_attr) %>% 
                                              group_by(cfclim) %>% 
                                              summarise(n = n()) %>% 
                                              arrange(cfclim) %>%  
                                              pull(n) %>% 
                                              rev, ")"),
                               no = paste0(stratum, " (",
                                           flow_df %>%
                                             filter(species %in% spec_attr) %>%
                                             group_by(fact) %>%
                                             summarise(n = n()) %>%
                                             arrange(as.character(fact)) %>%
                                             pull(n) %>%
                                             rev, ")"))),
            size = 8) +
  # species numbers in flow:
  geom_text(stat = "alluvium", 
            fontface = "italic",
            aes(label = ifelse(test = after_stat(x) == "1",
                               yes = plot_df_clim %>% pull(n),
                               no = "")),
            size = 7, nudge_x = 0.15) +
  scale_x_discrete(limits = c("cfclim", "fact"), expand = c(0.13, 0.0)) +
  labs(x = "", y = "") +
  annotate("text", x = 1, y = 85, label = "Counterfactual\nclimate", hjust = 0.5, size = 9) + 
  annotate("text", x = 2, y = 85, label = "Factual\nclimate", hjust = 0.5, size = 9) + 
  coord_cartesian(clip = 'off') + # prevent annotations to get cut
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        plot.margin = margin(40,15,-50,5), # trbl # top margin for labels when arranging multiple plots
        panel.border = element_blank(),
        axis.ticks = element_blank(), axis.text = element_blank(),
        text = element_text(size = 24),
        legend.position = "none")

clim_plot

# land use change:

lu_plot <- ggplot(plot_df_lodes, aes(x = x, stratum = stratum, alluvium = trend_change_lu, y = n)) +
  geom_flow(stat = "alluvium", aes(fill = trend_change_lu2),
            color = "darkgray",width = 1/4, show.legend = TRUE, alpha = 0.6) +
  scale_fill_manual(values = c("relative land use change loser" = "#abd9e9",
                               "absolute land use change loser" =  "#2c7bb6",
                               "no change" = "gray90",
                               "absolute land use change winner" = "#d7191c",
                               "relative land use change winner" = "#fdae61"), 
                    drop = FALSE, na.value = NA, name = "impact of drivers") +
  new_scale_fill() +
  geom_stratum(aes(fill = stratum), width = 1/4, colour = "grey20") + 
  scale_fill_manual(values = c("positive\ntrend" = "#EFCA08", "negative\ntrend" =  "#7DC4C9", "stable" =  "#EDE6F2"), na.value = NA) +
  # species numbers in bars:
  geom_text(stat = "stratum", 
            fontface = "italic",
            size = 8,
            aes(label = ifelse(test = after_stat(x) == "1",
                               yes = paste0(stratum, " (",
                                            plot_df_lodes %>% 
                                              filter(x == "cflu") %>% 
                                              group_by(stratum) %>% 
                                              summarise(n = sum(n)) %>% 
                                              pull(n) %>% 
                                              rev, ")"),
                               no = paste0(stratum, " (",
                                           plot_df_lodes %>% 
                                             filter(x == "fact") %>% 
                                             group_by(stratum) %>% 
                                             summarise(n = sum(n)) %>% 
                                             pull(n) %>% 
                                             rev, ")")))) +
  # species numbers in flow:
  geom_text(stat = "alluvium", 
            fontface = "italic",
            aes(label = ifelse(test = after_stat(x) == "1",
                               yes = plot_df_lu %>% pull(n),
                               no = "")),
            size = 7, nudge_x = 0.15) +
  scale_x_discrete(limits = c("cflu", "fact"), expand = c(0.13, 0.0)) +
  labs(x = "", y = "") +
  annotate("text", x = 1, y = 85, label = "Counterfactual\nland use", hjust = 0.5, size = 9) + 
  annotate("text", x = 2, y = 85, label = "Factual\nland use", hjust = 0.5, size = 9) + 
  coord_cartesian(clip = 'off') + # prevent annotations to get cut
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        axis.ticks = element_blank(), axis.text = element_blank(),
        plot.margin = margin(40,15,-50,5), # trbl
        panel.border = element_blank(),
        text = element_text(size = 24),
        legend.position = "none")

lu_plot

# arrange side-by-side:

combined_plot <- plot_grid(clim_plot,
                           lu_plot,
                           align = 'vh',
                           labels = "AUTO",
                           label_size = 26)

# get bird icon:

uuid <- get_uuid(name = "Contopus virens", n = 1)
img <- get_phylopic(uuid = uuid)
icon <- ggplot() +
  add_phylopic(x = 1, y = 1, img = img, remove_background = TRUE) + # vjust for alignment with legend in final plot
  theme_nothing()

# arrange all parts:

svg(file = file.path(dir, "plots", "attribution", "alluvial_climate_lu_manuscript_conceptual_legend.svg"),
    width = 17, # in inches
    height = 14)

# clear plot area:
grid.newpage()

# define area for main plot:
vp1 <- viewport(x = 0.5, y = 0.25, 
                height = 0.75, width = 1,
                just = c("center", "bottom"),
                name = "combined_plots")

# enter vp1 
pushViewport(vp1)
# add plot:
print(combined_plot, newpage = FALSE)

# leave vp1 - up one level (into root viewport)
upViewport(1)

# define legend area
vp2 <- viewport(x = 0.493, y = 0.02, 
                height = 0.25, width = 0.8,
                just = c("center", "bottom"),
                name = "legend")
# enter vp2
pushViewport(vp2)
# add plot:
grid.draw(g)

# leave vp2 - up one level (into root viewport)
upViewport(1)

# define icon area
vp3 <- viewport(x = 0.01, y = 0.09, 
                height = 0.15, width = 0.07,
                just = c("left", "bottom"),
                name = "icon")
pushViewport(vp3)
print(icon, newpage = FALSE)

dev.off()


## # version without conceptual legend: ----
# 
# # arrange plots:
# combined_plot <- plot_grid(clim_plot,
#                            lu_plot,
#                            align = 'vh',
#                            labels = "AUTO",
#                            label_size = 26)
# combined_plot
# 
# # add common legend:
# # plot to get legend from:
# leg_plot <- plot_df_clim %>% 
#   ggplot(aes(axis1 = cfclim, axis2 = fact, y = n)) +
#   geom_alluvium(aes(fill = trend_change_clim), show.legend = TRUE, width = 1/5, colour = "grey50", alpha = 0.6) +
#   scale_fill_manual(values = c("relative climate change loser" = "#abd9e9",
#                                "absolute climate change loser" =  "#2c7bb6",
#                                "no change" = "gray90",
#                                "absolute climate change winner" = "#d7191c",
#                                "relative climate change winner" = "#fdae61"), 
#                     labels = c("absolute change winner",
#                                "relative change winner",
#                                "no change",
#                                "relative change loser",
#                                "absolute change loser"),
#                     drop = FALSE, na.value = NA) +
#   new_scale_fill() +
#   scale_fill_manual(values = c("positive\ntrend" = "#EFCA08", "negative\ntrend" =  "#7DC4C9", "stable" =  "#EDE6F2"), na.value = NA) +
#   geom_stratum(width = 1/4, aes(fill = cfclim), show.legend = FALSE, colour = "grey20") + 
#   geom_stratum(width = 1/4, aes(fill = fact), show.legend = FALSE, colour = "grey20") +
#   guides(fill = "none") +
#   theme_bw() +
#   theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
#         plot.margin = margin(40,15,10,5), # trbl # top margin to later have label when arranging multiple plots
#         panel.border = element_blank(),
#         axis.ticks = element_blank(), axis.text = element_blank(),
#         text = element_text(size = 24),
#         plot.subtitle = element_text(margin = margin(0,0,30,0), size = 30),
#         legend.text = element_text(size = 24),
#         legend.position = "bottom", legend.direction = "vertical", 
#         legend.box.margin = margin(-50, 0, 0, 0),
#         legend.key.spacing.y = unit(0.2, "cm"), 
#         legend.title=element_blank())
# 
# # extract the legend:
# legend <- get_legend(leg_plot)
# 
# #  bird icon:
# uuid <- get_uuid(name = "Contopus virens", n = 1)
# img <- get_phylopic(uuid = uuid)
# icon <- ggplot() +
#   add_phylopic(x = 1, y = 1, img = img, remove_background = TRUE, vjust = 0.2) + # vjust for alignment with legend in final plot
#   theme_nothing()
# icon
# #get_attribution(uuid = uuid, text = TRUE)
# 
# # bottom row: icon + legend
# legend_plot <- plot_grid(NULL, icon, legend, NULL, 
#                          nrow = 1, ncol = 4,
#                          rel_widths = c(2.6, 0.4, 2.2, 3), # rel_widths = c(0.25, 0.05, 0.4, 0.3),
#                          scale = c(1, 1, 1, 1))
# legend_plot
# 
# # upper row: alluvial plots
# combined_plot_final <- plot_grid(combined_plot, legend_plot, nrow = 2, rel_heights = c(1, .2))
# 
# combined_plot_final
# 
# ggsave(filename = file.path(dir, "plots", "attribution", "alluvial_climate_lu_manuscript2.svg"), 
#        plot = combined_plot_final,
#        device = "svg",
#        width = 45,
#        height = 32,
#        units = "cm")

# session info:
writeLines(capture.output(sessionInfo()), file.path(dir, "results", "sessionInfo", "5_2_attribution_plots_trend_categories.txt"))
