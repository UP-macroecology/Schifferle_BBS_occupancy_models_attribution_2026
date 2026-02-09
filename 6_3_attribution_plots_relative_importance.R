# plots regarding relative importance of climate and land use change for occupancy dynamics
# based on difference in mean absolute percentage error between factual and counterfcatual predictions

# packages: --------------------------------------------------------------------

library(dplyr)
library(ggplot2)
library(grid) # trait boxplots


# directories: -----------------------------------------------------------------

main_dir <- file.path("//NAS-2-P-SN-01.ibb.uni-potsdam.de", "daten$", "AG26", "Transfer",
                      "Schifferle_BBS_occupancy_models_2023")


# load data: -------------------------------------------------------------------

# species for which models worked fine:
load(file.path("data", "species_DOM_val_okay.RData")) # output of 6_2_attribution_plots_trend_categories.R
spec_okay

# attribution metrics:

load(file = file.path(main_dir, "results", "attribution", "attribution_metrics_final.RData")) # output of 6_1_attribution_metrics.R
attr_metr_df

# trend change categories: 
load(file.path(main_dir, "results", "attribution", "trend_categories.RData")) # output of 6_2_attribution_plots_trend_categories.R
flow_df

# traits:
load(file.path("data", "BBS_data_merged.RData")) # bbs_dt


# lollipop plot: ---------------------------------------------------------------
# relative importance of different drivers:

# reformat data:

rel_imp_df <- attr_metr_df %>% 
  filter(species %in% spec_okay) %>%
  select(-matches("(slope)|(p_)")) %>% 
  mutate(imp_clim = mape_cfclim - mape_fact,
         imp_lu = mape_cflu - mape_fact,
         imp_climlu = mape_cfclimlu - mape_fact) %>% 
  # add relative impact categories for ordering:
  left_join(flow_df) %>% 
  # add species orders etc.:
  left_join(bbs_dt %>% select(English_Common_Name, Scientific_Name, ORDER, Family) %>% distinct,
            by = c(species = "English_Common_Name")) %>% 
  mutate(ORDER = factor(ORDER, levels = names(sort(table(ORDER), decreasing = TRUE))))

plot_df <- rel_imp_df %>%
  select(species, fact, imp_clim, imp_lu, imp_climlu, trend_change_clim, trend_change_lu, trend_change_climlu, Scientific_Name, ORDER, Family) %>% 
  # trend categories: absolute / relative winner / loser without driver:
  mutate(trend_change_clim2 = gsub(pattern = " climate change", replacement = "", x = trend_change_clim),
         trend_change_lu2 = gsub(pattern = " land use change", replacement = "", x = trend_change_lu),
         trend_change_climlu2 = gsub(pattern = " global change", replacement = "", x = trend_change_climlu)) %>% 
  # convert to factor:
  mutate(trend_change_clim2 = factor(trend_change_clim2, levels = c("absolute winner", "relative winner", "no change", "relative loser", "absolute loser")),
         trend_change_lu2 = factor(trend_change_lu2, levels = c("absolute winner", "relative winner", "no change", "relative loser", "absolute loser")),
         trend_change_climlu2  = factor(trend_change_climlu2 , levels = c("absolute winner", "relative winner", "no change", "relative loser", "absolute loser"))) %>% 
  # convert to long format:
  tidyr::pivot_longer(cols = starts_with("imp"), names_to = "scenario", values_to = "value") %>% 
  mutate(scenario = factor(scenario, levels = c("imp_clim", "imp_lu", "imp_climlu")),
         scenario = recode(scenario, imp_clim = "climate", imp_lu = "land use", imp_climlu = "climate + land use")) %>% 
  rename("climate" = trend_change_clim2, "land use" = trend_change_lu2,  "climate + land use" = trend_change_climlu2) %>% 
  tidyr::pivot_longer(cols = c("climate", "land use", "climate + land use"), names_to = "scenario2", values_to = "impact") %>% 
  filter(scenario == scenario2) %>% 
  select(-c(scenario2))

plot_df2 <-   plot_df %>% 
  # order by global change impact:
  # and by value corresponding to global change impact:
  left_join(plot_df %>% filter(scenario == "climate + land use") %>% select(species, value_global = value)) %>% 
  #arrange(trend_change_climlu , desc(value_global), ORDER, Family) %>% 
  
  # alternative: order by factual linear occupancy trend:
  arrange(fact, trend_change_climlu, desc(value_global), ORDER, Family) %>% 
  
  mutate(plot_order = row_number(),
         species = factor(species)) %>% 
  mutate(value = ifelse(value < 0, 0, value)) %>% # one species with negative values, does not make sense for relative impact
  # linerange min:
  left_join(plot_df %>% group_by(species) %>% summarise(min_infl = min(value))) %>% 
  # linerange max:
  left_join(plot_df %>% group_by(species) %>% summarise(max_infl = max(value))) %>% 
  # change names for plotting:
  mutate(scenario = recode(scenario,  "climate" = "climate change", "land use" = "land use change", "climate + land use" = "climate & land use change"))


# version for presentation:

svg(file = file.path("plots", "attribution", "summary_plots", "lollipop_rel_importance_all_drivers.svg"),
    width = 12, height = 9)
# pdf(file = file.path("plots", "attribution", "summary_plots", "lollipop_rel_importance_all_drivers.pdf"),
#     width = 16, height = 9)

ggplot(plot_df2, aes(x = forcats::fct_reorder(Scientific_Name, desc(plot_order)), # Scientific_Name
                     y = value)) +
  geom_linerange(aes(x = forcats::fct_reorder(Scientific_Name, desc(plot_order)), 
                     ymin = min_infl, ymax = max_infl),
                 colour = "gray80",
                 linewidth = 0.3, show.legend = FALSE) +
  geom_point(aes(shape = scenario, fill = impact, size = scenario), color = "black", stroke = 0.3) +
  scale_size_manual(values = c("climate change" = 2.8, "land use change" = 2.4, "climate & land use change" = 2.3), guide = "none") +
  scale_shape_manual(values = c("climate change" = 21, "land use change" = 24, "climate & land use change" = 22)) +
  scale_fill_manual(values = c("relative loser" = "#abd9e9",
                               "absolute loser" =  "#2c7bb6",
                               "absolute winner" = "#d7191c",
                               "relative winner" = "#fdae61",
                               "no change" = "gray70"), drop = FALSE, na.value = NA) +
  labs(x = "Species", y = "Relative importance") +
  theme_light() +
  coord_flip() +
  theme(
    panel.grid.major.y = element_blank(),
    panel.border = element_blank(),
    axis.ticks.y = element_blank(),
    axis.text.y = element_text(family = "serif", vjust = 0.3, face = "italic", size = 8,
                               margin = margin(t = 0, r = -20, b = 0, l = 0)),
    axis.text.x = element_text(size = 14),
    text = element_text(size = 20),
    legend.justification = "top",
    legend.title = element_text(margin = margin(b = 15)),
    legend.key.spacing.y = unit(10, "pt"),
  ) +
  guides(shape = guide_legend("driver", 
                              override.aes = list(size = 5, stroke = 1.2), order = 1,
                              theme(legend.text = element_text(vjust = 1, size = 22),
                                    legend.title = element_text(size = 26, face = "bold"))),
         fill = guide_legend("relative impact", 
                             override.aes = list(shape = 21, size = 5), order = 2,
                             theme(legend.text = element_text(vjust = 0.6, size = 22),
                                   legend.title = element_text(size = 26, face = "bold"))))
dev.off()



# version for manuscript:

# # bird icon:
# library(rphylopic)
# uuid <- get_uuid(name = "Contopus virens", n = 1)
# img <- get_phylopic(uuid = uuid)


lollipop <- ggplot(plot_df2, aes(x = forcats::fct_reorder(Scientific_Name, desc(plot_order)), # Scientific_Name
                     y = value)) +

  
  geom_linerange(aes(x = forcats::fct_reorder(Scientific_Name, desc(plot_order)), 
                     ymin = pmax(0, min_infl), ymax = max_infl),
                 colour = "gray30",
                 linewidth = 0.3, show.legend = FALSE) +
  geom_point(aes(shape = scenario, fill = impact, size = scenario), color = "black", stroke = 0.3) +
  scale_size_manual(values = c("climate change" = 3.2, "land use change" = 2.8, "climate & land use change" = 2.7), guide = "none") +
  scale_shape_manual(values = c("climate change" = 21, "land use change" = 24, "climate & land use change" = 22)) +
  scale_fill_manual(values = c("relative loser" = "#abd9e9",
                               "absolute loser" =  "#2c7bb6",
                               "absolute winner" = "#d7191c",
                               "relative winner" = "#fdae61",
                               "no change" = "gray70"), drop = FALSE, na.value = NA) +
  scale_y_continuous(limits = c(-0.005, 0.45), expand = c(0,0)) +
  labs(x = "Species", y = "Relative importance") +
  theme_light() +
  coord_flip() +
  theme(
    panel.grid.major.y = element_blank(),
    #axis.ticks.y = element_blank(),
    axis.text.y = element_text(family = "serif", vjust = 0.3, face = "italic", size = 14,
                               margin = margin(t = 0, r = 0, b = 0, l = 0)),
    axis.title.y = element_text(margin = margin(0, 10, 0, 0)),
    axis.text.x = element_text(size = 16),
    text = element_text(size = 20),
    legend.box.background = element_rect(fill = "white", colour = NA),
    legend.position = "inside",
    legend.justification.inside = c(0.95, 0.99),
    legend.title = element_text(margin = margin(b = 10), size = 18, face = "bold"),
    legend.key.spacing.y = unit(3, "pt"),
  ) +
  guides(shape = guide_legend("driver", 
                              override.aes = list(size = 3), order = 1,
                              theme(legend.text = element_text(vjust = 0.8, size = 16))),
         fill = guide_legend("occupancy trend change", 
                             override.aes = list(shape = 21, size = 4), order = 2,
                             theme(legend.text = element_text(vjust = 0.6, size = 16)))) #+
  # add bird icon:
  #add_phylopic(x = 5, y = 0.42, img = img, height = 7)

  # geom_rect(data = data_bar,
  #           aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, fill = fill_var),
  #           alpha = 0.5, inherit.aes = FALSE)
  lollipop

ggsave(filename = file.path("plots", "attribution", "summary_plots", "lollipop_rel_importance_manuscript2.svg"), 
       plot = lollipop,
       device = "svg",
       width = 25,
       height = 35,
       units = "cm")


# add factual linear trend information:

lollipop2 <- ggplot(plot_df2, aes(x = forcats::fct_reorder(Scientific_Name, desc(plot_order)), # species
                                 y = value)) +
  
  # add rectangles depicting factual trend category:
  geom_rect(aes(ymin = -0.01, ymax = 0, xmin = 0, xmax = 45.5), fill = "#7DC4C9") +
  geom_rect(aes(ymin = -0.01, ymax = 0, xmin = 45.5, xmax = 49.5), fill = "#EDE6F2") +
  geom_rect(aes(ymin = -0.01, ymax = 0, xmin = 49.5, xmax = 80.5), fill = "#EFCA08") +
  
  # lines
  geom_linerange(aes(x = forcats::fct_reorder(Scientific_Name, desc(plot_order)),
                     ymin = pmax(0, min_infl), ymax = max_infl),
                 colour = "gray80",
                 linewidth = 0.3, show.legend = FALSE) +
  
  geom_linerange(aes(x = forcats::fct_reorder(Scientific_Name, desc(plot_order)), 
                     ymin = 0, ymax = min_infl),
                 colour = "gray80",
                 linewidth = 0.3, show.legend = FALSE, linetype = "dashed") +
  
  
  geom_point(aes(shape = scenario, fill = impact, size = scenario), color = "black", stroke = 0.3) +
  scale_size_manual(values = c("climate change" = 3.2, "land use change" = 2.8, "climate & land use change" = 2.7), guide = "none") +
  scale_shape_manual(values = c("climate change" = 21, "land use change" = 24, "climate & land use change" = 22)) +
  scale_fill_manual(values = c("relative loser" = "#abd9e9",
                               "absolute loser" =  "#2c7bb6",
                               "absolute winner" = "#d7191c",
                               "relative winner" = "#fdae61",
                               "no change" = "gray70"), drop = FALSE, na.value = NA) +
  
  scale_y_continuous(limits = c(-0.01, 0.45), expand = c(0,0)) +
  labs(x = "Species", y = "Relative importance") +
  theme_light() +
  coord_flip() +
  theme(
    panel.grid.major.y = element_blank(),
    #axis.ticks.y = element_blank(),
    axis.text.y = element_text(family = "serif", vjust = 0.3, face = "italic", size = 14,
                               margin = margin(t = 0, r = 0, b = 0, l = 0)),
    axis.title.y = element_text(margin = margin(0, 10, 0, 0)),
    axis.text.x = element_text(size = 16),
    text = element_text(size = 20),
    legend.box.background = element_rect(fill = "white", colour = NA),
    legend.position = "inside",
    legend.justification.inside = c(0.95, 0.01),
    legend.title = element_text(margin = margin(b = 10), size = 18, face = "bold"),
    legend.key.spacing.y = unit(3, "pt"),
  ) +
  guides(shape = guide_legend("driver",
                              override.aes = list(size = 3), order = 1,
                              theme(legend.text = element_text(vjust = 0.8, size = 16))),
         fill = guide_legend("occupancy trend change",
                             override.aes = list(shape = 21, size = 4), order = 2,
                             theme(legend.text = element_text(vjust = 0.6, size = 16))))


lollipop2
ggsave(filename = file.path("plots", "attribution", "summary_plots", "lollipop_rel_importance_manuscript2.svg"), 
       plot = lollipop2,
       device = "svg",
       width = 25,
       height = 35,
       units = "cm")


# add oberserved trend to legend:

lollipop3 <- lollipop2 +
  guides(a = guide_custom(title = "observed trend",
                          grid::rectGrob(gp = gpar(fill="#EFCA08", col=NA)),
                          width = unit(0.5, "cm"), height = unit(0.5, "cm"),
                          order = 3),
         b = guide_custom(title = NULL,
                          grid::rectGrob(gp = gpar(fill="#EDE6F2", col=NA)),
                          width = unit(0.5, "cm"), height = unit(0.5, "cm"),
                          order = 4),
         c = guide_custom(title = NULL,
                          grid::rectGrob(gp = gpar(fill="#7DC4C9", col=NA)),
                          width = unit(0.5, "cm"), height = unit(0.5, "cm"),
                          order = 5)) +
  labs(tag = "increase\n\nstable\n\ndecrease") +
  theme(plot.tag.position = c(0.825, 0.13),
        plot.tag = element_text(hjust = 0, size = 16))

# requires manual adjustments of legend item:
ggsave(filename = file.path("plots", "attribution", "summary_plots", "lollipop_rel_importance_manuscript3.svg"), 
       plot = lollipop3,
       device = "svg",
       width = 25,
       height = 35,
       units = "cm")






# boxplots summarising relative importance of climate vs. land use change: -----

# reformat data:

boxplot_overall_df <- rel_imp_df %>% 
  select(species, imp_clim, imp_lu, imp_climlu) %>% 
  tidyr::pivot_longer(cols = starts_with("imp"), names_to = "scenario", values_to = "influence", 
                      names_pattern = "_(.*)") %>% 
  mutate(scenario = recode(scenario, clim = "climate\nchange", lu = "land use\nchange", climlu = "climate & land use\nchange"),
         scenario = factor(scenario, levels = c("climate\nchange", "land use\nchange", "climate & land use\nchange")))


# raincloud plot:

# version for presentation:

p <- boxplot_overall_df %>%
  ggplot(aes(#x = forcats::fct_reorder(scenario, desc(scenario)), # horizontal
    x = scenario, # vertical
    y = influence, fill = scenario)) +
  # add half-violin from {ggdist} package
  ggdist::stat_halfeye(adjust = 1, # density: bandwidth multiplied with this
                       justification = -0.8,
                       .width = 0, point_colour = NA,
                       width = 0.4) + # does not work with ggplot2 version 4; downgrade via remotes::install_version("ggplot2", version = "3.5.2", repos = "https://cran.r-project.org")
  geom_boxplot(width = 0.4, alpha = 0.5, linewidth = 1) +
  #ggdist::stat_dots(side = "left", justification = 1.1) +
  geom_point(aes(fill = scenario), shape = 21, size = 4, alpha = .5, position = position_jitter(seed = 1, width = .15)) +
  viridis::scale_fill_viridis(discrete = TRUE) +
  theme_bw() +
  #coord_flip() + # horizontal
  theme_bw() +
  labs(y = "Relative importance", x = "") +
  theme(text = element_text(size = 40),
        #axis.title.x = element_text(size = 40, margin = margin(t = 20)),# horizontal
        #axis.text.y = element_text(size = 40, margin = margin(r = 10)) # horizontal
        axis.text.x = element_text(size = 40, margin = margin(t = 20)), # vertical
        axis.title.y = element_text(size = 40, margin = margin(r = 10)) # vertical
  ) +
  guides(fill = "none")
p

# jpeg(file = file.path("plots", "attribution", "summary_plots",
#                       "density_driver_raincloud_plot_vertical.jpeg"),
#      width = 1000, height = 1000, quality = 100)
# jpeg(file = file.path("plots", "attribution", "summary_plots",
#                       "density_driver_raincloud_plot_horizontal.jpeg"),
#      width = 1150, height = 750, quality = 100)
p
#dev.off()

# version for manuscript:

boxplot_rel_imp <- boxplot_overall_df %>%
  ggplot(aes(x = scenario, y = influence, fill = scenario)) +
  # add half-violin from {ggdist} package
  ggdist::stat_halfeye(adjust = 1, # density: bandwidth multiplied with this
                       justification = -0.8,
                       .width = 0, point_colour = NA,
                       width = 0.3) + # does not work with ggplot2 version 4; downgrade via remotes::install_version("ggplot2", version = "3.5.2", repos = "https://cran.r-project.org")
  geom_boxplot(width = 0.3, alpha = 0.5, linewidth = 1) +
  geom_point(aes(fill = scenario), shape = 21, size = 2, alpha = .5, position = position_jitter(seed = 1, width = .15)) +
  viridis::scale_fill_viridis(discrete = TRUE) +
  theme_bw() +
  labs(y = "Relative importance") +
  theme(axis.text.x = element_text(size = 14, margin = margin(t = 5)), 
        axis.title.y = element_text(size = 16), # vertical, margin = margin(r = 10)
        axis.text.y = element_text(size = 14),
        axis.title.x = element_blank(),
        plot.margin = margin(0, 0, 0, 1)
  ) +
  guides(fill = "none")

ggsave(filename = file.path("plots", "attribution", "summary_plots", "boxplot_rel_importance_manuscript.svg"), 
       plot = boxplot_rel_imp,
       device = "svg",
       width = 16,
       height = 16,
       units = "cm")


# plots traits and relative importance: ------------------------------------

# rel. influence ~ trait, grouped by winner/loser category:

# change order: same category next to each other, trait in legend

# reformat data:



boxplot_df %>% 
  select(species, Scientific_Name, ORDER, Family) %>% 
  distinct() %>%  View


rel_imp_df %>% 
  filter(species %in% c("Tufted Titmouse", "Carolina Chickadee", "Eastern Wood-Pewee", "Willow Flycatcher", "Acadian Flycatcher",
                        "Dickcissel", "Yellow-throated Vireo", "Fish Crow")) %>% # Currie & Venne 2017
  select(species, imp_clim, imp_lu, imp_climlu) %>% 
  arrange(desc(imp_clim))
