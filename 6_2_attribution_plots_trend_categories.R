# alluvial plots for comparing factual and counterfactual occupancy trends
# categories: absolute and relative winners and losers of change

# packages: --------------------------------------------------------------------

library(dplyr)
library(ggplot2)
library(ggalluvial) # alluvial plots
library(ggnewscale) # to have multiple colour fill scales


# directories: -----------------------------------------------------------------

main_dir <- file.path("//NAS-2-P-SN-01.ibb.uni-potsdam.de", "daten$", "AG26", "Transfer",
                      "Schifferle_BBS_occupancy_models_2023")


# load data: -------------------------------------------------------------------

# species with fair model performance:

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

# attribution metrics:

load(file = file.path(main_dir, "results", "attribution", "attribution_metrics_final.RData")) # output of 6_1_attribution_metrics.R


# alluvial plots: --------------------------------------------------------------

# https://r-charts.com/flow/ggalluvial/
# change categories based on linear trend

# restructure data for plotting:


flow_df <- attr_metr_df %>% 
  # categorize relative impact (based on Langhammer et al. 2024):
  select(c(species, starts_with("slope"), starts_with("p_"))) %>% 
  # only species with good model performance:
  filter(species %in% spec_okay) %>% 
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

#save(flow_df, file =  file.path(main_dir, "results", "attribution", "trend_categories.RData"))

## climate change: ----

plot_df_clim <- flow_df %>%
  group_by(fact, cfclim, trend_change_clim) %>% 
  summarise(n = n()) %>% 
  arrange(cfclim, trend_change_clim) # change plotting order

# version for presentation:

jpeg(file = file.path("plots", "attribution", "summary_plots", "alluvial_climate_presentation.jpg"),
     width = 910, height = 1200, quality = 100)
svg(file = file.path("plots", "attribution", "summary_plots", "alluvial_climate_presentation.svg"),
    width = 14, height = 18)

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
                                              filter(species %in% spec_okay) %>% 
                                              group_by(cfclim) %>% 
                                              summarise(n = n()) %>% 
                                              arrange(cfclim) %>%  
                                              pull(n) %>% 
                                              rev, ")"),
                               no = paste0(stratum, " (",
                                           flow_df %>%
                                             filter(species %in% spec_okay) %>%
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
  annotate("text", x = 1, y = 85, label = "Counterfactual\nclimate", hjust = 0.5, size = 10) + # y = 85 for spec_okay
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

dev.off()

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
                                              filter(species %in% spec_okay) %>% 
                                              group_by(cfclim) %>% 
                                              summarise(n = n()) %>% 
                                              arrange(cfclim) %>%  
                                              pull(n) %>% 
                                              rev, ")"),
                               no = paste0(stratum, " (",
                                           flow_df %>%
                                             filter(species %in% spec_okay) %>%
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
  annotate("text", x = 1, y = 85, label = "Counterfactual\nclimate", hjust = 0.5, size = 9) + # y = 85 for spec_okay
  annotate("text", x = 2, y = 85, label = "Factual\nclimate", hjust = 0.5, size = 9) + # size = 9 for paper version
  coord_cartesian(clip = 'off') + # prevent annotations to get cut
  guides(fill = "none") +
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        plot.margin = margin(40,15,10,5), # trbl # top margin to later have label when arranging multiple plots
        panel.border = element_blank(),
        axis.ticks = element_blank(), axis.text = element_blank(),
        text = element_text(size = 24),
        plot.subtitle = element_text(margin = margin(0,0,30,0), size = 30),
        legend.text = element_text(size = 24),
        legend.position = "bottom", legend.direction = "vertical", 
        legend.box.margin = margin(-50, 0, 0, 0),
        legend.key.spacing.y = unit(0.2, "cm"), 
        legend.title=element_blank())


ggsave(filename = file.path("plots", "attribution", "summary_plots", "alluvial_climate_manuscript.svg"), 
       plot = clim_plot,
       device = "svg",
       width = 21,
       height = 29.7, # A4
       units = "cm")



# lodes form: ---

# rendering order unfortunate

plot_df_lodes <- to_lodes_form(plot_df_clim,
                               axes = c(1:2),
                               id = trend_change_clim) %>%
  mutate(x = factor(x, levels = c("cfclim", "fact"))) %>%
  # same order as in plot_df_clim: order in df determines rendering order
  mutate(trend_change_clim2 = rep(plot_df_clim$trend_change_clim, 2))


# jpeg(file = file.path("plots", "attribution", "summary_plots", "alluvial_climate_presentation2.jpg"),
#      width = 900, height = 1200, quality = 100)

ggplot(plot_df_lodes, aes(x = x, stratum = stratum,
                          alluvium = trend_change_clim,
                          y = n#,
                         # order = rep(c(1:9), 2)
                 )
                 ) + 
  geom_flow(stat = "alluvium", aes(fill = trend_change_clim2),
            color = "darkgray",width = 1/5, show.legend = TRUE) +
  scale_fill_manual(values = c("relative climate change loser" = "#abd9e9",
                               "absolute climate change loser" =  "#2c7bb6",
                               "no change" = "gray90",
                               "absolute climate change winner" = "#d7191c",
                               "relative climate change winner" = "#fdae61"),
                    drop = FALSE, na.value = NA, name = "impact of drivers") +
  new_scale_fill() +
  geom_stratum(aes(fill = stratum), width = 1/5, colour = "grey20") +
  scale_fill_manual(values = c("positive\ntrend" = "#EFCA08", "negative\ntrend" =  "#7DC4C9", "stable" =  "#EDE6F2"), na.value = NA) +
  # text in bars:
  geom_text(stat = "stratum",
            fontface = "italic",
            size = 7,
            aes(label = ifelse(test = after_stat(x) == "1",
                               yes = paste0(stratum, " (",
                                            plot_df_lodes %>%
                                              filter(x == "cfclim") %>%
                                              group_by(stratum) %>%
                                              summarise(n = sum(n)) %>%
                                              mutate(stratum = as.character(stratum)) %>%
                                              arrange(stratum) %>%  #xx seltsam
                                              pull(n), ")"),
                               no = paste0(stratum, " (",
                                           plot_df_lodes %>%
                                             filter(x == "fact") %>%
                                             group_by(stratum) %>%
                                             summarise(n = sum(n)) %>%
                                             mutate(stratum = as.character(stratum)) %>%
                                             arrange(desc(stratum)) %>% # xx seltsam
                                             pull(n), ")")))) +
  # species numbers in flow:
  geom_text(stat = "alluvium",
            fontface = "italic",
            aes(label = ifelse(test = after_stat(x) == "1",
                               yes = plot_df_lodes %>% filter(x == "cfclim") %>% pull(n),
                               no = "")),
            size = 6, nudge_x = 0.15) +
  scale_x_discrete(limits = c("cfclim", "fact"), expand = c(0.13, 0.0)) +
  labs(title = "Change in area of occupancy",
       subtitle = "numbers = number of species",
       x = "", y = "") +
  annotate("text", x = 1, y = 85, label = "Counterfactual\nclimate", hjust = 0.5, size = 8) + # y = 85 for spec_okay; 170 for all
  annotate("text", x = 2, y = 85, label = "Observed\nclimate change", hjust = 0.5, size = 8) + # xx
  guides(fill = "none") +
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        axis.ticks = element_blank(), axis.text = element_blank(),
        text = element_text(size = 26),
        plot.subtitle = element_text(margin = margin(0,0,30,0)),
        legend.text = element_text(size = 24),
        legend.position = "bottom", legend.direction = "vertical", legend.box.spacing = unit(-50, "pt"))

#dev.off()



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

jpeg(file = file.path("plots", "attribution", "summary_plots", "alluvial_lu_presentation.jpg"),
     width = 920, height = 1200, quality = 100)
svg(file = file.path("plots", "attribution", "summary_plots", "alluvial_lu_presentation.svg"),
    width = 14, height = 18)

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
  annotate("text", x = 1, y = 85, label = "Counterfactual\nland use", hjust = 0.5, size = 10) + # y = 85 for spec_okay; 170 for all
  annotate("text", x = 2, y = 85, label = "Observed\nland use change", hjust = 0.5, size = 10) + # xx
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

dev.off()


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
  annotate("text", x = 1, y = 85, label = "Counterfactual\nland use", hjust = 0.5, size = 9) + # y = 85 for spec_okay
  annotate("text", x = 2, y = 85, label = "Factual\nland use", hjust = 0.5, size = 9) + # size = 9 for paper version
  coord_cartesian(clip = 'off') + # prevent annotations to get cut
  guides(fill = "none") +
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        plot.margin = margin(40,15,10,5), # trbl
        panel.border = element_blank(),
        axis.ticks = element_blank(), axis.text = element_blank(),
        text = element_text(size = 24),
        plot.subtitle = element_text(margin = margin(0,0,30,0), size = 30),
        legend.text = element_text(size = 24),
        legend.position = "bottom", legend.direction = "vertical", 
        legend.box.margin = margin(-50, 0, 0, 0),
        legend.key.spacing.y = unit(0.2, "cm"), 
        legend.title=element_blank())

ggsave(filename = file.path("plots", "attribution", "summary_plots", "alluvial_lu_manuscript.svg"), 
       plot = lu_plot,
       device = "svg",
       width = 21,
       height = 29.7, # A4
       units = "cm")



## climate + land use change: ----

plot_df_climlu <- flow_df %>%
  group_by(fact, cfclimlu, trend_change_climlu) %>% 
  summarise(n = n()) %>% 
  arrange(cfclimlu, fact, trend_change_climlu) # change plotting order


# version for presentation:

jpeg(file = file.path("plots", "attribution", "summary_plots", "alluvial_climlu_presentation.jpg"),
     width = 900, height = 1200, quality = 100)

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
                                              filter(species %in% spec_okay) %>% 
                                              group_by(cfclimlu) %>% 
                                              summarise(n = n()) %>% 
                                              arrange(cfclimlu) %>%  
                                              pull(n) %>% 
                                              rev, ")"),
                               no = paste0(stratum, " (",
                                           flow_df %>%
                                             filter(species %in% spec_okay) %>%
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

dev.off()



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
                                              filter(species %in% spec_okay) %>% 
                                              group_by(cfclimlu) %>% 
                                              summarise(n = n()) %>% 
                                              arrange(cfclimlu) %>%  
                                              pull(n) %>% 
                                              rev, ")"),
                               no = paste0(stratum, " (",
                                           flow_df %>%
                                             filter(species %in% spec_okay) %>%
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
  annotate("text", x = 1, y = 85, label = "Counterfactual\nclimate & land use", hjust = 0.5, size = 9) + # y = 85 for spec_okay
  annotate("text", x = 2, y = 85, label = "Factual climate\n& land use", hjust = 0.5, size = 9) + # size = 9 for paper version
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

ggsave(filename = file.path("plots", "attribution", "summary_plots", "alluvial_climlu_manuscript.svg"), 
       plot = climlu_plot,
       device = "svg",
       width = 21,
       height = 29.7, # A4
       units = "cm")





# lodes form, seems more stable, but rendering order not optimal:
plot_df_lodes <- to_lodes_form(plot_df_climlu,
                               axes = c(1:2),
                               id = trend_change_climlu ) %>% 
  mutate(x = factor(x, levels = c("cfclimlu", "fact"))) %>% 
  # same order as in plot_df_climlu:
  mutate(trend_change_climlu2 = rep(plot_df_climlu$trend_change_climlu, 2))


# jpeg(file = file.path("plots", "attribution", "summary_plots", "alluvial_climlu.jpg"),
#      width = 900, height = 1200, quality = 100)

plot_df_lodes %>% # reihenfolge falsch xx
ggplot(aes(x = x, stratum = stratum, alluvium = trend_change_climlu,#rep(c(-1:-4, -7, -6, -5, -8:-10), 2),
           y = n)) + # change rendering order via order in alluvium argument
  geom_flow(stat = "alluvium", aes(fill = trend_change_climlu2),
            color = "darkgray", width = 1/5, show.legend = TRUE) +
  scale_fill_manual(values = c("relative global change loser" = "#abd9e9",
                               "absolute global change loser" =  "#2c7bb6",
                               "no change" = "gray90",
                               "absolute global change winner" = "#d7191c",
                               "relative global change winner" = "#fdae61"),
                    drop = FALSE, na.value = NA, name = "impact of drivers") +
  new_scale_fill() +
  geom_stratum(aes(fill = stratum), width = 1/5, colour = "grey20") +
  scale_fill_manual(values = c("positive\ntrend" = "#EFCA08", "negative\ntrend" =  "#7DC4C9", "stable" =  "#EDE6F2"), na.value = NA) +
  # text in bars:
  geom_text(stat = "stratum",
            fontface = "italic",
            size = 7,
            aes(label = ifelse(test = after_stat(x) == "1",
                               yes = paste0(stratum, " (",
                                            plot_df_lodes %>%
                                              filter(x == "cfclimlu") %>%
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
                               yes = plot_df_lodes %>% filter(x == "cfclimlu") %>% pull(n),
                               no = "")),
            size = 6, nudge_x = 0.15) +
  scale_x_discrete(limits = c("cfclimlu", "fact"), expand = c(0.15, 0.0)) +
  labs(title = "Change in area of occupancy",
       subtitle = "linear trend in 100 posterior draws of N. occ. routes per year\nnumbers = number of species",
       x = "", y = "") +
  annotate("text", x = 1, y = 85, label = "Counterfactual\nclimate & land use", hjust = 0.5, size = 8) +
  annotate("text", x = 2, y = 85, label = "Observed \nglobal change", hjust = 0.5, size = 8) + 
  guides(fill = "none") +
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        axis.ticks = element_blank(), axis.text = element_blank(),
        text = element_text(size = 26),
        plot.subtitle = element_text(margin = margin(0,0,30,0)),
        legend.text = element_text(size = 24),
        legend.position = "bottom", legend.direction = "vertical", legend.box.spacing = unit(-50, "pt"))

#dev.off()




# arrange plots for manuscript: ----

# library(ggpubr)
# 
# combined_plot <- ggarrange(clim_plot, lu_plot, ncol = 2, nrow = 1, 
#                            common.legend = TRUE, legend = "bottom",
#                            labels = "AUTO",
#                            font.label = list(size = 26),
#                            hjust = -0.3) +
#   theme(legend.box.spacing = unit(-100, "pt"))
# 



library(cowplot)
combined_plot2 <- plot_grid(clim_plot + theme(legend.position="none"),
                            lu_plot + theme(legend.position="none"),
                            align = 'vh',
                            labels = "AUTO",
                            label_size = 26)
combined_plot2

# plot to get legend from:
leg_plot <- plot_df_clim %>% 
  ggplot(aes(axis1 = cfclim, axis2 = fact, y = n)) +
  geom_alluvium(aes(fill = trend_change_clim), show.legend = TRUE, width = 1/5, colour = "grey50", alpha = 0.6) +
  scale_fill_manual(values = c("relative climate change loser" = "#abd9e9",
                               "absolute climate change loser" =  "#2c7bb6",
                               "no change" = "gray90",
                               "absolute climate change winner" = "#d7191c",
                               "relative climate change winner" = "#fdae61"), 
                    labels = c("absolute change winner",
                               "relative change winner",
                               "no change",
                               "relative change loser",
                               "absolute change loser"),
                    drop = FALSE, na.value = NA) +
  new_scale_fill() +
  scale_fill_manual(values = c("positive\ntrend" = "#EFCA08", "negative\ntrend" =  "#7DC4C9", "stable" =  "#EDE6F2"), na.value = NA) +
  geom_stratum(width = 1/4, aes(fill = cfclim), show.legend = FALSE, colour = "grey20") + 
  geom_stratum(width = 1/4, aes(fill = fact), show.legend = FALSE, colour = "grey20") +
  guides(fill = "none") +
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        plot.margin = margin(40,15,10,5), # trbl # top margin to later have label when arranging multiple plots
        panel.border = element_blank(),
        axis.ticks = element_blank(), axis.text = element_blank(),
        text = element_text(size = 24),
        plot.subtitle = element_text(margin = margin(0,0,30,0), size = 30),
        legend.text = element_text(size = 24),
        legend.position = "bottom", legend.direction = "vertical", 
        legend.box.margin = margin(-50, 0, 0, 0),
        legend.key.spacing.y = unit(0.2, "cm"), 
        legend.title=element_blank())

# extract the legend:
legend <- get_legend(
  # create some space to the left of the legend
  leg_plot
  #+ theme(legend.box.margin = margin(0, 0, 20, 0))
)


combined_plot_final <- plot_grid(combined_plot2, legend, 
                                 nrow = 2,
                                 rel_heights = c(1, .15))
combined_plot_final

# does make new plots with 
ggsave(filename = file.path("plots", "attribution", "summary_plots", "alluvial_climate_lu_manuscript.svg"), 
       plot = combined_plot_final,
       device = "svg",
       width = 45,
       height = 32,
       units = "cm")


# version with bird icon:

library(rphylopic)

uuid <- get_uuid(name = "Contopus virens", n = 1)
img <- get_phylopic(uuid = uuid)

icon <- ggplot() +
  add_phylopic(x = 1, y = 1, img = img, remove_background = TRUE, vjust = 0.2) + # vjust for alignment with legend in final plot
  theme_nothing()

icon
#get_attribution(uuid = uuid)

# bottom row:
legend_plot <- plot_grid(NULL, icon, legend, NULL, 
                         nrow = 1, ncol = 4,
                         rel_widths = c(2.6, 0.4, 2.2, 3), # rel_widths = c(0.25, 0.05, 0.4, 0.3),
                         scale = c(1, 1, 1, 1))
legend_plot

combined_plot_final2 <- plot_grid(combined_plot2, 
                                  legend_plot, 
                                  nrow = 2,
                                  rel_heights = c(1, .2))
combined_plot_final2
ggsave(filename = file.path("plots", "attribution", "summary_plots", "alluvial_climate_lu_manuscript2.svg"), 
       plot = combined_plot_final2,
       device = "svg",
       width = 45,
       height = 32,
       units = "cm")



# conceptual figures, replicated from Langhammer et al. 2024: ----


dummy_df <- expand.grid(impact = factor(c("absolute winners", "relative winners", 
                                          #"no change", 
                                          "relative losers", "absolute losers"), 
                                        levels = c("absolute winners", "relative winners", 
                                                   #"no change", 
                                                   "relative losers", "absolute losers")),
                        scenario = factor(c("observed change" ,"counterfactual scenario")),
                        x = c(0, 1)) %>% 
  as_tibble() %>% 
  arrange(impact) %>% 
  mutate(ymin = c(0, 0, 2, -2, # abs. winners
                  0, 0, -1, -6,  # rel. winners
                  #0, 0, 2, 1, # no change
                  0, 0, 1, 6, # rel. losers
                  0, 0, -6, -1), # abs. losers
         ymax = c(0, 0, 5, 1, # abs. winners
                  0, 0, -4, -9, # rel. winners
                  #0, 0, 5, 4, # no change
                  0, 0, 4, 9, # rel. losers
                  0, 0, -3, 3)) # abs. losers


# version for presentation:

conc_plot <- ggplot(dummy_df) +
  geom_ribbon(aes(x = x, ymin = ymin, ymax = ymax, fill = scenario)) +
  facet_wrap(~impact, ncol = 1, # change for horizontal version
             axes = "all", axis.labels = "all") + #, scales = "free"
  geom_hline(yintercept = 0, linetype = "longdash") +
  #geomtextpath::geom_labelpath(aes(x = x, y = ymin), size = 5, fill = "#F6F6FF") + # xx
  scale_x_continuous(expand = c(0, 0)) +
  scale_fill_manual(values = c(alpha("grey60", 0.6), alpha("grey10", 0.6))) +
  labs(x = "time", y = "occupied area", title = "Impact categories") +
  theme_bw() +
  theme(axis.ticks = element_blank(),
        axis.text = element_blank(),
        panel.grid = element_blank(),
        strip.background = element_rect(fill = NA),# colour added later
        legend.position = "bottom", legend.direction = "vertical",legend.title=element_blank(),
        text = element_text(size = 30), # paper version size = 25
        plot.title = element_text(margin=margin(0,0,20,0)),
        legend.key.spacing.y = unit(0.3, "cm"),
        axis.title = element_text(size = 25)
  ) 
conc_plot

# add colours for facet strips:
g <- ggplot_gtable(ggplot_build(conc_plot))
stripr <- which(grepl('strip', g$layout$name))
fills <- c(
  alpha("#2c7bb6", alpha = 0.6),
  alpha("#abd9e9", alpha = 0.6),
  #alpha("gray90", alpha = 0.6),
  alpha("#fdae61", alpha = 0.6),
  alpha("#d7191c", alpha = 0.6)
)

k <- 1
for (i in stripr) {
  j <- which(grepl('rect', g$grobs[[i]]$grobs[[1]]$childrenOrder))
  g$grobs[[i]]$grobs[[1]]$children[[j]]$gp$fill <- fills[k]
  k <- k+1
}

jpeg(file = file.path("plots", "attribution", "summary_plots",
                      "conceptual_trend_change_categories.jpg"),
     width = 350, height = 1000, quality = 100)

# svg(file = file.path("plots", "attribution", "summary_plots",
#                      "conceptual_trend_change_categories.svg"),
#     width = 9, height = 9)

grid::grid.draw(g)
dev.off()

# version for manuscript: