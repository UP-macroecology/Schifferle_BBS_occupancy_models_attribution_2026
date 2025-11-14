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
  select(species, imp_clim, imp_lu, imp_climlu, trend_change_clim, trend_change_lu, trend_change_climlu, Scientific_Name, ORDER, Family) %>% 
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
  arrange(trend_change_climlu , desc(value_global), ORDER, Family) %>% 
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
                       width = 0.4) +
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


# jpeg(file = file.path("plots", "attribution", "summary_plots",
#                       "density_driver_raincloud_plot_vertical.jpeg"),
#      width = 1000, height = 1000, quality = 100)
# jpeg(file = file.path("plots", "attribution", "summary_plots",
#                       "density_driver_raincloud_plot_horizontal.jpeg"),
#      width = 1150, height = 750, quality = 100)
p
#dev.off()

# version for manuscript:


# plots traits and relative importance: ------------------------------------

# rel. influence ~ trait, grouped by winner/loser category:

# change order: same category next to each other, trait in legend

# reformat data:

boxplot_df <- rel_imp_df %>%
  tidyr::pivot_longer(cols = starts_with("imp"), names_to = "scenario", values_to = "value") %>% 
  mutate(scenario = factor(scenario, levels = c("imp_clim", "imp_lu", "imp_climlu")),
         scenario = recode(scenario, imp_clim = "climate", imp_lu = "land use", imp_climlu = "climate + land use")) %>% 
  # add traits:
  left_join(bbs_dt %>% select(English_Common_Name, Habitat, Migration, Trophic.Level, Trophic.Niche, Primary.Lifestyle) %>% distinct,
            by = c(species = "English_Common_Name")) %>% 
  mutate(Migration = forcats::fct_recode(Migration, 
                                         "partially migratory" = "Part.migratory",
                                         "migratory" = "Migratory",
                                         "sedentary" = "Sedentary")) %>% 
  mutate(trend_change_clim = case_when(trend_change_clim == "absolute climate change winner" ~ "absolute climate\nchange winner",
                                       trend_change_clim == "relative climate change winner" ~ "relative climate\nchange winner",
                                       trend_change_clim == "no change" ~ "no change",
                                       trend_change_clim == "absolute climate change loser" ~ "absolute climate\nchange loser",
                                       trend_change_clim == "relative climate change loser" ~ "relative climate\nchange loser",
                                       .default = NA),
         trend_change_clim = factor(trend_change_clim, levels = c("absolute climate\nchange winner", "relative climate\nchange winner", "no change", "relative climate\nchange loser", "absolute climate\nchange loser")),
         trend_change_lu = case_when(trend_change_lu == "absolute land use change winner" ~ "absolute land use\nchange winner",
                                     trend_change_lu == "relative land use change winner" ~ "relative land use\nchange winner",
                                     trend_change_lu == "no change" ~ "no change",
                                     trend_change_lu == "absolute land use change loser" ~ "absolute land use\nchange loser",
                                     trend_change_lu == "relative land use change loser" ~ "relative land use\nchange loser",
                                .default = NA),
         trend_change_lu = factor(trend_change_lu, levels = c("absolute land use\nchange winner", "relative land use\nchange winner", "no change", "relative land use\nchange loser", "absolute land use\nchange loser")),
  ) 


## raincloud plot without differentiation between winners / losers: ----


# fct. for raincloud plots:
raincloud_plot_fct <- function(trait = "Migration", subset = NA, driver = "climate"){
  
  if(any(!is.na(subset))){
    plot_df <- plot_df %>% 
      filter(!!sym(trait) %in% subset)
  }
  
  labels <- plot_df %>% group_by(!!sym(trait)) %>% summarise(n = n()) %>% 
    mutate(label = paste0(!!sym(trait), "\nN = ", n)) %>% 
    pull(label)
  
  plot_df %>%
    ggplot(aes(x = !!sym(trait), y = value, fill = !!sym(trait))) +
    # add half-violin from {ggdist} package
    ggdist::stat_halfeye(adjust = .5, justification = -0.8, 
                         .width = 0, point_colour = NA,
                         width = 0.4) +
    geom_boxplot(width = 0.4, alpha = 0.5, linewidth = 1) +
    #ggdist::stat_dots(side = "left", justification = 1.1) +
    geom_point(aes(fill = !!sym(trait)), shape = 21, size = 4, alpha = .5, position = position_jitter(seed = 1, width = .2)) +
    viridis::scale_fill_viridis(discrete = TRUE) +
    theme_bw() +
    theme(
      legend.position = "none",
    ) +
    scale_x_discrete(labels = labels) +
    theme_bw() +
    labs(#title = "Climate change impact and migratory\nstrategy", 
      #y = expression("Relative impact" ~ (MAPE [counterfactual] - MAPE [factual])),
      y = paste("Relative", driver, "change impact"),
      x = "") +
    theme(text = element_text(size = 40),
          axis.title.y = element_text(margin = margin(r = 10))) +
    guides(fill = "none") +
    ggtitle(paste(trait))
}


### climate change: ----

plot_df <- boxplot_df %>% 
  filter(scenario == "climate")

# migration:

jpeg(file = file.path("plots", "attribution", "summary_plots", "trait_boxplots",
                      "rel_imp_climate_migration_raincloud_plot.jpg"),
     width = 850, height = 1000, quality = 100)
raincloud_plot_fct(trait = "Migration")
dev.off()

# habitat:

jpeg(file = file.path("plots", "attribution", "summary_plots", "trait_boxplots",
                      "rel_imp_climate_habitat_raincloud_plot.jpg"),
     width = 850, height = 1000, quality = 100)
raincloud_plot_fct(trait = "Habitat", subset = names(table(plot_df$Habitat)[which(table(plot_df$Habitat) >= 5)]))

dev.off()

# trophic level:

jpeg(file = file.path("plots", "attribution", "summary_plots", "trait_boxplots",
                      "rel_imp_climate_trophic_level_raincloud_plot.jpg"),
     width = 850, height = 1000, quality = 100)
raincloud_plot_fct(trait = "Trophic.Level", 
                   subset = names(table(plot_df$Trophic.Level)[which(table(plot_df$Trophic.Level) >= 5)]))
dev.off()

# trophic niche:

jpeg(file = file.path("plots", "attribution", "summary_plots", "trait_boxplots",
                      "rel_imp_climate_trophic_niche_raincloud_plot.jpg"),
     width = 850, height = 1000, quality = 100)
raincloud_plot_fct(trait = "Trophic.Niche", 
                   subset = names(table(plot_df$Trophic.Niche)[which(table(plot_df$Trophic.Niche) >= 5)]))
dev.off()

# primary lifestyle:
jpeg(file = file.path("plots", "attribution", "summary_plots", "trait_boxplots",
                      "rel_imp_climate_primary_lifestyle_raincloud_plot.jpg"),
     width = 850, height = 1000, quality = 100)
raincloud_plot_fct(trait = "Primary.Lifestyle", 
                   subset = names(table(plot_df$Primary.Lifestyle)[which(table(plot_df$Primary.Lifestyle) >= 5)]))
dev.off()


### land use change: ----

plot_df <- boxplot_df %>% 
  filter(scenario == "land use")


# migration:
jpeg(file = file.path("plots", "attribution", "summary_plots", "trait_boxplots",
                      "rel_imp_lu_migration_raincloud_plot.jpg"),
     width = 850, height = 1000, quality = 100)
raincloud_plot_fct(trait = "Migration", driver = "land use")
dev.off()

# habitat:

jpeg(file = file.path("plots", "attribution", "summary_plots", "trait_boxplots",
                      "rel_imp_lu_habitat_raincloud_plot.jpg"),
     width = 850, height = 1000, quality = 100)
raincloud_plot_fct(trait = "Habitat", driver = "land use", subset = names(table(plot_df$Habitat)[which(table(plot_df$Habitat) >= 5)]))

dev.off()

# trophic level:

jpeg(file = file.path("plots", "attribution", "summary_plots", "trait_boxplots",
                      "rel_imp_lu_trophic_level_raincloud_plot.jpg"),
     width = 850, height = 1000, quality = 100)
raincloud_plot_fct(trait = "Trophic.Level", driver = "land use",
                   subset = names(table(plot_df$Trophic.Level)[which(table(plot_df$Trophic.Level) >= 5)]))
dev.off()

# trophic niche:

jpeg(file = file.path("plots", "attribution", "summary_plots", "trait_boxplots",
                      "rel_imp_lu_trophic_niche_raincloud_plot.jpg"),
     width = 850, height = 1000, quality = 100)
raincloud_plot_fct(trait = "Trophic.Niche", driver = "land use",
                   subset = names(table(plot_df$Trophic.Niche)[which(table(plot_df$Trophic.Niche) >= 5)]))
dev.off()

# primary lifestyle:
jpeg(file = file.path("plots", "attribution", "summary_plots", "trait_boxplots",
                      "rel_imp_lu_primary_lifestyle_raincloud_plot.jpg"),
     width = 850, height = 1000, quality = 100)
raincloud_plot_fct(trait = "Primary.Lifestyle", driver = "land use",
                   subset = names(table(plot_df$Primary.Lifestyle)[which(table(plot_df$Primary.Lifestyle) >= 5)]))
dev.off()


## boxplots differentiating between winners / losers: ----


# fct. for boxplots differentiating between winners / losers:

# coloured stripe placement may need to be adjusted xx

boxplot_fct <- function(trait = "Migration", subset = NA, driver = "climate"){
  
  if(any(!is.na(subset))){
    plot_df <- plot_df %>% 
      filter(!!sym(trait) %in% subset)
  }
  
  # winner / loser category:
  if(driver == "climate"){x_axis <- "trend_change_clim"} 
  else if (driver == "land use") {x_axis <- "trend_change_lu"} 
  else {x_axis <- "trend_change_climlu"}
  
  plot_df %>% 
    ggplot(aes(x = !!sym(x_axis), 
               y = value, fill = !!sym(trait))) +
    geom_boxplot(width = 0.6, position = position_dodge(0.75)) +
    stat_summary(fun.data = get_box_stats, geom = "text", hjust = 0.5, vjust = 0.9, size = 7,
                 position = position_dodge(width = .75)) +
    viridis::scale_fill_viridis(discrete = TRUE, begin = 0.4, option = "D") +
    theme_bw() +
    labs(title = paste(driver, "change impact"), 
         y = expression(MAPE [counterfactual] - MAPE [factual]), 
         x = "",
         subtitle = "numbers = number of species") +
    theme(text = element_text(size = 26),
          axis.text.x = element_text(margin=margin(t=20), colour = "black"),
          legend.key.spacing.y = unit(0.3, 'cm'),
          legend.key.size = unit(1.2, 'cm'),
          legend.title = element_text(margin = margin(b = 10))) +
    guides(fill = guide_legend(trait)) +
    coord_cartesian(clip='off') +
    annotation_custom(
      grob=muh_grob, xmin = 0, xmax = 1, ymin = -0.03, ymax = -0.065)
}

### climate change: ----


# display boxes for winner / loser categories - climate change version:
muh_grob <- grid::rectGrob(
  x = 1:5, y= 0, gp=gpar(
    fill = c("#d7191c","#fdae61", "gray90", "#abd9e9","#2c7bb6"),
    col = "white",
    alpha=0.8))


# fct. to display number of species in each group:
get_box_stats <- function(y) {
  return(data.frame(y = 0.03 + max(y), # y position of label
                    label = length(y))
  )
}

plot_df <- boxplot_df %>% 
  filter(scenario == "climate")

# migration:

jpeg(file = file.path("plots", "attribution", "summary_plots", "trait_boxplots",
                      "rel_imp_climate_migration_boxplots.jpg"),
     width = 1500, height = 1000, quality = 100)
boxplot_fct(trait = "Migration", subset = NA, driver = "climate")
dev.off()

# habitat
jpeg(file = file.path("plots", "attribution", "summary_plots", "trait_boxplots",
                      "rel_imp_climate_habitat_boxplots.jpg"),
     width = 1500, height = 1000, quality = 100)
boxplot_fct(trait = "Habitat", 
            subset = names(table(plot_df$Habitat)[which(table(plot_df$Habitat) >= 5)]), 
            driver = "climate")
dev.off()

# trophic level:

jpeg(file = file.path("plots", "attribution", "summary_plots", "trait_boxplots",
                      "rel_imp_climate_trophic_level_boxplots.jpg"),
     width = 1500, height = 1000, quality = 100)
boxplot_fct(trait = "Trophic.Level", 
            subset = names(table(plot_df$Trophic.Level)[which(table(plot_df$Trophic.Level) >= 5)]), 
            driver = "climate")
dev.off()

# trophic niche:

jpeg(file = file.path("plots", "attribution", "summary_plots", "trait_boxplots",
                      "rel_imp_climate_trophic_niche_boxplots.jpg"),
     width = 1500, height = 1000, quality = 100)
boxplot_fct(trait = "Trophic.Niche", 
            subset = names(table(plot_df$Trophic.Niche)[which(table(plot_df$Trophic.Niche) >= 5)]), 
            driver = "climate")
dev.off()

# primary lifestyle:

jpeg(file = file.path("plots", "attribution", "summary_plots", "trait_boxplots",
                      "rel_imp_climate_primary_lifestyle_boxplots.jpg"),
     width = 1500, height = 1000, quality = 100)
boxplot_fct(trait = "Primary.Lifestyle", 
            subset = names(table(plot_df$Primary.Lifestyle)[which(table(plot_df$Primary.Lifestyle) >= 5)]), 
            driver = "climate")
dev.off()


### land use change: ----

# adjust:
get_box_stats <- function(y) {
  return(data.frame(y = 0.01 + max(y), # y position of label
                    label = length(y))
  )
}

# display boxes for winner / loser categories - land use change version:
muh_grob <- grid::rectGrob(
  x = 1:4, y= 1, gp=gpar(
    fill = c("#fdae61", "gray90", "#abd9e9","#2c7bb6"),
    col = "white",
    alpha=0.8)) # xx 

plot_df <- boxplot_df %>% 
  filter(scenario == "land use")

# migration:

jpeg(file = file.path("plots", "attribution", "summary_plots", "trait_boxplots",
                      "rel_imp_lu_migration_boxplots.jpg"),
     width = 1500, height = 1000, quality = 100)
boxplot_fct(trait = "Migration", subset = NA, driver = "land use")
dev.off()

# habitat
jpeg(file = file.path("plots", "attribution", "summary_plots", "trait_boxplots",
                      "rel_imp_lu_habitat_boxplots.jpg"),
     width = 1500, height = 1000, quality = 100)
boxplot_fct(trait = "Habitat", 
            subset = names(table(plot_df$Habitat)[which(table(plot_df$Habitat) >= 5)]), 
            driver = "land use")
dev.off()

# trophic level:

jpeg(file = file.path("plots", "attribution", "summary_plots", "trait_boxplots",
                      "rel_imp_lu_trophic_level_boxplots.jpg"),
     width = 1500, height = 1000, quality = 100)
boxplot_fct(trait = "Trophic.Level", 
            subset = names(table(plot_df$Trophic.Level)[which(table(plot_df$Trophic.Level) >= 5)]), 
            driver = "land use")
dev.off()

# trophic niche:

jpeg(file = file.path("plots", "attribution", "summary_plots", "trait_boxplots",
                      "rel_imp_lu_trophic_niche_boxplots.jpg"),
     width = 1500, height = 1000, quality = 100)
boxplot_fct(trait = "Trophic.Niche", 
            subset = names(table(plot_df$Trophic.Niche)[which(table(plot_df$Trophic.Niche) >= 5)]), 
            driver = "land use")
dev.off()

# primary lifestyle:

jpeg(file = file.path("plots", "attribution", "summary_plots", "trait_boxplots",
                      "rel_imp_lu_primary_lifestyle_boxplots.jpg"),
     width = 1500, height = 1000, quality = 100)
boxplot_fct(trait = "Primary.Lifestyle", 
            subset = names(table(plot_df$Primary.Lifestyle)[which(table(plot_df$Primary.Lifestyle) >= 5)]), 
            driver = "land use")
dev.off()
