# Script:   6_concept_figure.R
# Purpose:  Generate sub-plots for conceptual figure (Fig. 1)
# Inputs:   data/observed_time_series_1995_2019/<species>_obs_ts_sum_occ_routes.RData
#           results/fm_buffer750km/fact_pred_time_series_1995_2019/<species>_ts_sum_occ_routes_f_preds.RData (one file per species)
#           results/attribution/cfact_pred_time_series_1995_2019/<species>_ts_sum_occ_routes_cf_preds.RData (one file per species)
# Outputs:  plots/concept_fit_model.svg
#           plots/concept_apply_model.svg
#           plots/concept_winners_losers.svg
#           plots/concept_impact_attr.svg
# Runs on:  Local

source(file.path("scripts", "0_paths.R"))

# packages: --------------------------------------------------------------------

library(dplyr)
library(ggplot2)

# figure to illustrate overall steps: ------------------------------------------

example_species <- "Bell's Vireo"

# load time series: ------------------------------------------------------------

# observations time series:
load(file.path(dir, "data", "observed_time_series_1995_2019", 
               paste0(example_species, "_obs_ts_sum_occ_routes.RData"))) # output of 4_1_DOMs_predictions_time_series.R
ts_obs

# predicted time series for factual data:
load(file.path(dir, "results", "fm_buffer750km", "fact_pred_time_series_1995_2019",
               paste0(example_species, "_ts_sum_occ_routes_f_preds.RData"))) # output of 4_1_DOMs_predictions_time_series.R

ts_preds_fact

# predicted time series for counterfactual data:
load(file.path(dir, "results", "attribution", "cfact_pred_time_series_1995_2019", 
               paste0(example_species, "_ts_sum_occ_routes_cf_preds.RData"))) # output of 4_1_DOMs_predictions_time_series.R
ts_preds_cfact

ts_dt <- ts_obs %>% 
  left_join(ts_preds_fact %>% select(year, fact = median_Nocc)) %>% 
  left_join(ts_preds_cfact$cf_clim %>% select(year, cf_clim = median_Nocc)) %>% 
  left_join(ts_preds_cfact$cf_1995soc %>% select(year, cf_lu = median_Nocc)) %>% 
  left_join(ts_preds_cfact$cf_clim_1995soc %>% select(year, cf_climlu = median_Nocc))

colours <- c("Npres" = "black", "fact" = "#85CB33", 
             "cf_clim" = "#0D98BA", "cf_lu" = "#B7410E",
             "cf_climlu" = "#046865")

# fit models: ------------------------------------------------------------------

plot <- ggplot(ts_dt, aes(x = year)) +
  geom_line(aes(y = Npres, colour = "Npres"), linewidth = 0.5) +
  geom_point(aes(y = Npres, colour = "Npres"), size = 1) +
  geom_line(aes(y = fact, colour = "fact"), linewidth = 0.5) +
  #geom_line(aes(y = cf_clim, colour = "cf_clim"), linewidth = 0.5) +
  #geom_line(aes(y = cf_lu, colour = "cf_lu")) +
  #geom_line(aes(y = cf_climlu, colour = "cf_climlu")) +
  scale_x_continuous("Year", breaks = c(1995, 2019)) +
  scale_color_manual(values = colours, 
                     labels = c("factual simulation", "bird records")) +
  labs(x = "Year", y = "N occupied sites", color = "Legend") +
  theme_classic() +
  theme(axis.text.y = element_blank(),
        axis.ticks = element_blank(),
        legend.title = element_blank(),
        legend.key.spacing.y = unit(-0.2, "cm"),
        legend.direction = "vertical",
        legend.justification = "top",
        legend.text = element_text(size = 8),
        axis.title.y = element_text(size = 9),
        axis.title.x = element_text(size = 9, margin = margin(t = -5, r = 0, b = 0, l = 0)),
        legend.margin = margin(t = -5, r = 0, b = 0, l = 0),
        plot.margin = margin(t = 0, r = 0, b = 0, l = 0)) +
  guides(color = guide_legend(reverse = TRUE))

plot

ggsave(filename = file.path(dir, "plots", "concept_fit_model.svg"),
       plot = plot,
       device = "svg",
       width = 8.5,# 10
       height = 4, # A4 # 6
       units = "cm")


# apply models: ----------------------------------------------------------------

colours <- c("Npres" = "grey60", "fact" = "#85CB33", "cf" = "#0D98BA")

plot <- ts_dt %>% 
  tidyr::pivot_longer(cols = cf_clim:cf_climlu, names_to = "scenario", values_to = "cf_occ") %>% 
  mutate(scenario = case_when(scenario == "cf_clim" ~ "no climate change",
                              scenario == "cf_lu" ~ "no land use change",
                              scenario == "cf_climlu" ~ "no climate & land use change")) %>%
  mutate(scenario = factor(scenario, levels = c("no climate change", "no land use change", "no climate & land use change"))) %>% 
  ggplot(aes(x = year)) +
  geom_line(aes(y = Npres, colour = "Npres"), linewidth = 0.3) +
  geom_point(aes(y = Npres, colour = "Npres"), size = 0.8) +
  geom_line(aes(y = fact, colour = "fact"), linewidth = 0.3) +
  geom_line(aes(y = cf_occ, colour = "cf"), linewidth = 0.5) +
  facet_wrap(~scenario, nrow = 3) +
  scale_x_continuous("Year", breaks = c(1995, 2019)) +
  scale_color_manual(values = colours, 
                     labels = c("counterfactual simulation", "factual simulation", "bird records")) +
  labs(x = "Year", y = "N occupied sites", color = "Legend") +
  theme_classic() +
  theme(axis.text.y = element_blank(),
        axis.ticks = element_blank(),
        legend.title = element_blank(),
        legend.text = element_text(size = 8),
        legend.direction = "vertical",
        legend.position = "bottom",
        legend.justification = "left",
        legend.key.spacing.y = unit(-0.2, "cm"),
        legend.margin = margin(t = -5, r = 0, b = 0, l = 0),
        strip.text.x = element_text(size = 8), # 18
        strip.background = element_blank(),
        panel.border = element_rect(colour = "grey50", fill = NA, linewidth = 0.2),
        panel.widths = unit(4, "cm"),
        panel.spacing = unit(0.1, "cm"),
        axis.title.x = element_text(margin = margin(t = -5, r = 0, b = -5, l = 0), 
                                    size = 9),
        axis.title.y = element_text(size = 9),
        plot.margin = margin(t = 0, r = 5, b = 0, l = -90)) +
  guides(color = guide_legend(reverse = TRUE, ncol=2))

plot

ggsave(filename = file.path(dir, "plots", "concept_apply_model.svg"),
       plot = plot,
       device = "svg",
       width = 8,# 10
       height = 7.75, # A4 # 6
       units = "cm")

# ggsave(filename = file.path(dir, "plots", "concept_apply_model.png"), 
#        plot = plot,
#        device = "png",
#        width = 10,
#        height = 20, 
#        units = "cm",
#        dpi = 300)


# absolute and relative winner and loser categories: ---------------------------

dummy_df <- expand.grid(impact = factor(c("absolute winners", "relative winners", 
                                          "relative losers", "absolute losers"), 
                                        levels = c("absolute winners", "relative winners", 
                                                   "relative losers", "absolute losers")),
                        scenario = factor(c("                   factual" ,"              counterfactual")), # added spaces to get text placed in the wider part of the shaded area
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
  facet_wrap(~impact, ncol = 2, # change for horizontal version
             axes = "all", axis.labels = "all") +
  geom_hline(yintercept = 0, linetype = "dashed", linewidth = 0.5) +
  
  geomtextpath::geom_textpath(aes(x = x, y = ymean, label = scenario, colour = scenario), 
                              size = 2.5, text_only = TRUE) +
  
  scale_colour_manual(values = c("white", "black"), guide = "none") +
  scale_x_continuous(expand = c(0, 0)) +
  scale_fill_manual(values = c(alpha("grey10", 0.6), alpha("grey60", 0.6)), guide = "none" ) +
  labs(x = "Year", y = "N occupied sites") +
  theme_bw() +
  theme(axis.ticks = element_blank(),
        axis.text = element_blank(),
        axis.title.y = element_text(size = 9),
        axis.title.x = element_text(size = 9, margin = margin(t = 1, r = 0, b = 0, l = 0)),
        panel.grid = element_blank(),
        strip.background = element_rect(fill = NA),# colour added later
        strip.text.x = element_text(size = 8),
        plot.margin = margin(0,0,0,0)
  ) 
conc_plot_m

# add colours for facet strips:
g <- ggplot_gtable(ggplot_build(conc_plot_m))
stripr <- which(grepl('strip', g$layout$name))
fills <- c(
  alpha("#abd9e9", alpha = 0.6),
  alpha("#2c7bb6", alpha = 0.6),
  alpha("#d7191c", alpha = 0.6),
  alpha("#fdae61", alpha = 0.6)
)

k <- 1
for (i in stripr) {
  j <- which(grepl('rect', g$grobs[[i]]$grobs[[1]]$childrenOrder))
  g$grobs[[i]]$grobs[[1]]$children[[j]]$gp$fill <- fills[k]
  k <- k+1
}

grid::grid.draw(g)


svg(file = file.path(dir, "plots", "concept_winners_losers.svg"),
    width = 2.5, height = 2.1)

grid::grid.draw(g)
dev.off()

# png(file = file.path(dir, "plots", "concept_winners_losers.png"),
#     width = 10, height = 9, unit = "cm", res = 300, pointsize = 12)
# 
# grid::grid.draw(g)
# dev.off()


# quantification driver importance: --------------------------------------------

plot <- ggplot(ts_dt, aes(x = year)) +
  geom_ribbon(aes(ymin=Npres,ymax=fact), fill="blue") +
  geom_ribbon(aes(ymin=Npres,ymax=cf_clim), fill="yellow") +
  geom_line(aes(y = Npres), colour = "black", linewidth = 0.5) +
  geom_point(aes(y = Npres), colour = "black", size = 1) +
  geom_line(aes(y = fact, colour = "fact"), linewidth = 0.5) +
  geom_line(aes(y = cf_clim, colour = "cf"), linewidth = 0.5) +
  scale_x_continuous("Year", breaks = c(1995, 2019)) +
  scale_color_manual(values = colours, 
                     labels = c("counterfactual simulation", "factual simulation", "bird records")) +
  labs(x = "Year", y = "N occupied sites", color = "Legend") +
  theme_classic() +
  theme(axis.text.y = element_blank(),
        axis.ticks = element_blank(),
        legend.position = "bottom",
        legend.direction = "vertical",
        legend.justification = "left",
        legend.key.spacing.y = unit(-0.2, "cm"),
        legend.title = element_blank(),
        legend.text = element_text(size = 8),
        axis.title.y = element_text(size = 9),
        axis.title.x = element_text(size = 9, margin = margin(t = -5, r = 0, b = 0, l = 0)),
        legend.margin = margin(t = -5, r = 0, b = 0, l = 0),
        plot.margin = margin(t = 0, r = 45, b = 0, l = 0)) +
  guides(color = guide_legend(reverse = TRUE, ncol=2))

plot

ggsave(filename = file.path(dir, "plots", "concept_impact_attr.svg"), 
       plot = plot,
       device = "svg",
       width = 7.5,
       height = 5.5, # A4
       units = "cm")

# Inkscape:
# click on "Fill", change to "Pattern"
# Stripes 1:4, Polka dots medium
# use node tool to scale (square) and rotate
# scale: use ctrl for dots, but not for stripes


# session info:
writeLines(capture.output(sessionInfo()), file.path(dir, "results", "sessionInfo", "6_concept_figure.txt"))
