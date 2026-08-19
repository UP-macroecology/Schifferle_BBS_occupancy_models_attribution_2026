# Script:   5_4b_attribution_maps.R
# Purpose:  Maps on spatial patterns
# Inputs:   results/species_DOM_val_okay.RData
#           data/BBS_data_merged.RData
#           data/US_outline_ESRI102003.shp
#           data/Birdlife_range_maps/<species>.shp
#           results/attribution/attribution_metrics_final.RData
#           results/attribution/trend_categories.RData
# Outputs:  data/Birdlife_range_maps/US_clipped/<species>_US_clipped.shp
#           data/Birdlife_range_maps/tifs/<species>_US_clipped.tif
#           data/Birdlife_range_maps/tifs_attr/<species>_attr.tif
#           plots/attribution/maps/species_richness.svg (Fig. S4)
#           plots/attribution/maps/clim_lu_combined_all.svg (Fig. 5)
#           plots/attribution/maps/clim_change_impact_cats.svg (Fig. S5)
#           plots/attribution/maps/lu_change_impact_cats.svg (Fig. S6)
#           plots/attribution/maps/climlu_change_impact_cats.svg (Fig. S7)
# Runs on:  Local

# Steps:
# rasterize shapefiles
# add impact categories and driver importance as attributes
# generate maps on spatial patterns in:
# - community mean relative importance of climate change and land use change
# - occupancy trends (absolute / relative winners and losers)
# at the range level, based on overlaying BirdLife range 


source(file.path("scripts", "0_paths.R"))


# packages: --------------------------------------------------------------------

library(dplyr)
library(sf)
library(terra)
library(ggplot2)
library(tidyterra)
library(cowplot)


# directories: -----------------------------------------------------------------

# directory with extracted BirdLife ranges:
output_dir <- file.path(dir, "data", "Birdlife_range_maps")

# directory to store maps:
plot_dir <- file.path(dir, "plots", "attribution", "maps")
if(!dir.exists(plot_dir)){dir.create(plot_dir, recursive = TRUE)}

# prepare data: ----------------------------------------------------------------

# species for attribution:
load(file = file.path(dir, "results", "species_DOM_val_okay.RData")) # output of 4_0_DOMs_predictions_y_routes_scenarios.R
spec_attr

# add scientific names:
load(file.path(dir, "data", "BBS_data_merged.RData")) # bbs_dt; output of 1_0_dataprep_BBS_bird_data.R

# outline conterminous US
US_albers_sf <- read_sf(file.path(dir, "data", "US_outline_ESRI102003.shp")) # output of 1_0_dataprep_climate.R

# attribution metrics:
load(file = file.path(dir, "results", "attribution", "attribution_metrics_final.RData")) # output of 5_1_attribution_metrics.R
attr_metr_df

# impact categories:
load(file =  file.path(dir, "results", "attribution", "trend_categories.RData")) # output of 5_2_attribution_plots_trend_categories.R
flow_df

## rasterize shapefile for easier summary calculations: ------------------------

# clip to conterminous USA, specify extent to ensure same grid for all, rasterize
# extent USA:
extent_US <- st_bbox(US_albers_sf)

# iterate over species:
spec_names <- bbs_dt %>% 
  select(English_Common_Name, Scientific_Name) %>% 
  distinct %>% 
  filter(English_Common_Name %in% spec_attr)

for(s in 1:nrow(spec_names)){
  
  spec <- spec_names$English_Common_Name[s]
  
  print(paste(s, spec))
  
  range_sf <- read_sf(file.path(output_dir, paste0(spec, ".shp"))) %>% 
    # add species common name:
    mutate(species = spec_names$English_Common_Name[s])
  
  # clip to USA:
  print("clip")
  
  range_sf_US <- range_sf %>%
    st_transform(crs = "ESRI:102003") %>% 
    # simplify geometry to avoid processing issues:
    st_make_valid() %>% 
    # union breeding-only and year-round presence:
    st_union() %>% 
    st_intersection(US_albers_sf)
  
  #plot(st_geometry(range_sf_US), main = paste(s, spec))
  
  st_write(range_sf_US, file.path(output_dir, "US_clipped",  paste0(spec, "_US_clipped.shp")), append = FALSE)
  
  # rasterize:
  print("rasterize")
  
  gdalUtilities::gdal_rasterize(src_datasource = file.path(output_dir, "US_clipped",  paste0(spec, "_US_clipped.shp")),
                                dst_filename = file.path(output_dir, "tifs",  paste0(spec, "_US_clipped.tif")),
                                burn = 1,
                                te = c(extent_US[1],extent_US[2],extent_US[3],extent_US[4]), # target extent of the resulting raster
                                tr = c(50000, 50000),
                                at = TRUE,
                                a_nodata = -99999) # value for cells with missing data
  }


## add impact categories and driver importance as attributes: ------------------

rel_imp_df <- attr_metr_df %>%
  select(-matches("(slope)|(p_)")) %>%
  mutate(imp_clim = mape_cfclim - mape_fact,
         imp_lu = mape_cflu - mape_fact,
         imp_climlu = mape_cfclimlu - mape_fact)

imp_cat_df <- flow_df %>% 
  select(species, starts_with("trend")) %>% 
  mutate(trend_change_clim = case_match(trend_change_clim, "absolute climate change winner" ~ 1,
                                        "relative climate change winner" ~ 2,
                                        "no change" ~ 3,
                                        "relative climate change loser" ~ 4,
                                        "absolute climate change loser" ~ 5),
         trend_change_lu = case_match(trend_change_lu, "absolute land use change winner" ~ 1,
                                      "relative land use change winner" ~ 2,
                                      "no change" ~ 3,
                                      "relative land use change loser" ~ 4,
                                      "absolute land use change loser" ~ 5),
         trend_change_climlu = case_match(trend_change_climlu, "absolute global change winner" ~ 1,
                                          "relative global change winner" ~ 2,
                                          "no change" ~ 3,
                                          "relative global change loser" ~ 4,
                                          "absolute global change loser" ~ 5))


# species ranges as multilayer tifs with attribution data:

for(s in 1:nrow(spec_names)){
  
  spec <- spec_names$English_Common_Name[s]
  
  print(paste(s, spec))
  
  range_tif <- rast(file.path(output_dir, "tifs",  paste0(spec, "_US_clipped.tif")))
  
  range_tif$imp_clim <- ifel(range_tif == 1, rel_imp_df %>% filter(species == spec) %>% pull(imp_clim), NA)
  range_tif$imp_lu <- ifel(range_tif == 1, rel_imp_df %>% filter(species == spec) %>% pull(imp_lu), NA) 
  range_tif$imp_climlu <- ifel(range_tif == 1, rel_imp_df %>% filter(species == spec) %>% pull(imp_climlu), NA)
  
  range_tif$cat_clim <- ifel(range_tif == 1, imp_cat_df %>% filter(species == spec) %>% pull(trend_change_clim), NA) 
  range_tif$cat_lu <- ifel(range_tif == 1, imp_cat_df %>% filter(species == spec) %>% pull(trend_change_lu), NA) 
  range_tif$cat_climlu <- ifel(range_tif == 1, imp_cat_df %>% filter(species == spec) %>% pull(trend_change_climlu), NA) 
  
  writeRaster(range_tif, file.path(output_dir, "tifs_attr",  paste0(spec, "_attr.tif")), overwrite = TRUE)
}


# maps: ------------------------------------------------------------------------


## 1) general species richness map (Fig. S 3): ----

raster_list <- list.files(file.path(output_dir, "tifs"), pattern = "tif", full.names = TRUE)
ranges <- rast(raster_list)
sr <- sum(ranges, na.rm = TRUE)

sr_map <- ggplot() +
  geom_spatraster(data = sr) +
  scale_fill_whitebox_c(palette = "viridi") +
  labs(
    fill = "N species",
    title = "Species richness") +
  theme_bw() +
  theme(text = element_text(size = 24),
        legend.key.height = unit(1, "cm"),
        plot.margin = margin(0,0,0,0))
sr_map

ggsave(filename = file.path(plot_dir, "species_richness.svg"),
       plot = sr_map,
       device = "svg",
       width = 29.7,
       height = 18, # A4
       units = "cm")


## 2) climate change mean importance: ----

raster_list <- list.files(file.path(output_dir, "tifs_attr"), pattern = ".tif$", full.names = TRUE)

ranges_clim <- rast(raster_list, lyrs = seq(2, 560, by = 7)) # imp_clim
mean_clim_imp <- mean(ranges_clim, na.rm = TRUE)

# use color scale matching boxplots:
# viridis purple for climate change importance
clim_col <- viridis::viridis(n = 3)[1]
clim_imp_map <- ggplot() +
  geom_spatraster(data = mean_clim_imp * 100) + # as percentage
  scale_fill_gradient(low = "#FFFFFF", high = clim_col, na.value = "transparent") +
  labs(fill = "community mean importance [%]") +
  theme_bw() +
  guides(colour = "none") +
  theme(text = element_text(size = 24),
        legend.key.height = unit(1, "cm"),
        legend.key.width = unit(.2, "npc"),
        legend.position = "bottom",
        legend.title.position = "top",
        legend.box.margin = margin(0,0,0,0),
        legend.margin = margin(0,0,0,0),
        plot.margin = unit(c(2, 0, 0, 0), "cm")) # trbl # top margin to later have label when arranging multiple plots)
clim_imp_map

ggsave(filename = file.path(plot_dir, "mean_climate_change_importance.svg"),
       plot = clim_imp_map,
       device = "svg",
       width = 29.7,
       height = 21, # A4
       units = "cm")


## 3) land use change mean importance: ----

ranges_lu <- rast(raster_list, lyrs = seq(3, 560, by = 7))
mean_lu_imp <- mean(ranges_lu, na.rm = TRUE)

# use color scale matching boxplots:
# viridis petrol for land use change importance
lu_col <- viridis::viridis(n = 3)[2]

lu_imp_map <- ggplot() +
  geom_spatraster(data = mean_lu_imp * 100) + # as percentage
  scale_fill_gradient(low = "#FFFFFF", high = lu_col, na.value = "transparent") +
  labs(fill = "community mean importance [%]") +
  guides(colour = "none") +
  theme_bw() +
  theme(text = element_text(size = 24),
        legend.key.height = unit(1, "cm"),
        legend.key.width = unit(.2, "npc"),
        legend.position = "bottom",
        legend.title.position = "top",
        legend.box.margin = margin(0,0,0,0),
        legend.margin = margin(0,0,0,0),
        plot.margin = unit(c(2, 0, 0, 0), "cm"))
lu_imp_map

ggsave(filename = file.path(plot_dir, "mean_lu_change_importance.svg"),
       plot = lu_imp_map,
       device = "svg",
       width = 29.7,
       height = 21, # A4
       units = "cm")




## 4) dominant driver (driver combination): ----

clim_col <- viridis::viridis(n = 3)[1]
lu_col <- viridis::viridis(n = 3)[2]
climlu_col <- viridis::viridis(n = 3)[3]

# dominant driver:

## climate + land use change mean importance: ----
ranges_climlu <- rast(raster_list, lyrs = seq(4, 560, by = 7))
mean_climlu_imp <- mean(ranges_climlu, na.rm = TRUE)

rast_stack <- c(mean_clim_imp, mean_lu_imp, mean_climlu_imp)
dom_driver_rast <- app(rast_stack, which.max)
names(dom_driver_rast) <- "dom_driver"

# convert to factor:
cls <- data.frame(id = 1:3, driver = c("climate change", "land use change", "climate & land use change"))
levels(dom_driver_rast) <- cls

# plot as binary map:
dom_map <- ggplot() +
  geom_spatraster(data = dom_driver_rast, aes(fill = driver)) + # colour
  scale_fill_viridis_d(na.translate = FALSE, breaks = c("climate change", "climate & land use change", "land use change")) +
  guides(fill = guide_legend(order = 1, ncol = 2), 
         colour = guide_legend(order = 2)) +
  theme_bw() +
  theme(text = element_text(size = 24),
        legend.position = "bottom",
        legend.title = element_blank(),
        legend.key.spacing.y = unit(0.3, "cm"),
        legend.box.margin = margin(0,0,0,-50),
        legend.margin = margin(0,0,0,0),
        legend.text = element_text(size = 24, hjust = 0),
        plot.margin = margin(0, 0, 0, 0, unit = "cm")
        ) #trbl
dom_map


## -> arrange climate, land use importance and dominant driver (Fig. 5): ----

combined_plot <- plot_grid(clim_imp_map, lu_imp_map, dom_map, 
                           align = 'vh',
                           labels = c("A: climate change", "B: land use change", "C: dominating driver at community level"),
                           label_size = 26, 
                           nrow = 2,
                           hjust = c(-0.14,-0.14,-0.05),
                           vjust = 2
                           )

combined_plot

ggsave(filename = file.path(plot_dir, "clim_lu_combined_all.svg"),
       plot = combined_plot,
       device = "svg",
       width = 40,
       height = 35,
       units = "cm")


## 5) climate change impact categories (Fig. S 4): ----

raster_list <- list.files(file.path(output_dir, "tifs_attr"), pattern = ".tif$", full.names = TRUE)

abs_winners <- imp_cat_df %>% filter(trend_change_clim == 1) %>% pull(species)
rel_winners <- imp_cat_df %>% filter(trend_change_clim == 2) %>% pull(species)
rel_losers <- imp_cat_df %>% filter(trend_change_clim == 4) %>% pull(species)
abs_losers <- imp_cat_df %>% filter(trend_change_clim == 5) %>% pull(species)

aw_raster_subset <- grep(paste(abs_winners, collapse="|"),  raster_list, value = TRUE)
rw_raster_subset <- grep(paste(rel_winners, collapse="|"),  raster_list, value = TRUE)
rl_raster_subset <- grep(paste(rel_losers, collapse="|"),  raster_list, value = TRUE)
al_raster_subset <- grep(paste(abs_losers, collapse="|"),  raster_list, value = TRUE)

# load rasters:
aw_ranges <- rast(aw_raster_subset, lyrs = seq(1, 560, by = 7))
rw_ranges <- rast(rw_raster_subset, lyrs = seq(1, 560, by = 7))
rl_ranges <- rast(rl_raster_subset, lyrs = seq(1, 560, by = 7))
al_ranges <- rast(al_raster_subset, lyrs = seq(1, 560, by = 7))

# sum across species:
aw_sum_cat <- sum(aw_ranges, na.rm = TRUE)
rw_sum_cat <- sum(rw_ranges, na.rm = TRUE)
rl_sum_cat <- sum(rl_ranges, na.rm = TRUE)
al_sum_cat <- sum(al_ranges, na.rm = TRUE)

clim_imp_cat <- c(aw_sum_cat, rw_sum_cat, rl_sum_cat, al_sum_cat)
names(clim_imp_cat) <- c("absolute winners", "relative winners", "relative losers", "absolute losers")

clim_cat_map <- ggplot() +
  geom_sf(data = US_albers_sf, colour = "black", fill = NA) +
  geom_spatraster(data = clim_imp_cat) +
  facet_wrap(~lyr, ncol = 2) +
  scale_fill_whitebox_c(palette = "viridi") +
  labs(fill = "N species",
       title = "Climate change impact") +
  theme_bw() +
  theme(text = element_text(size = 21),
        legend.key.height = unit(1, "cm"))
clim_cat_map

# add colours for facet strips:
g <- ggplot_gtable(ggplot_build(clim_cat_map))
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

ggsave(filename = file.path(plot_dir, "clim_change_impact_cats.svg"),
       plot = g,
       device = "svg",
       width = 29.7,
       height = 21, # A4
       units = "cm")


## 6) land use change impact categories (Fig. S 5): ----

#abs_winners <- imp_cat_df %>% filter(trend_change_lu == 1) %>% pull(species) # none
rel_winners <- imp_cat_df %>% filter(trend_change_lu == 2) %>% pull(species)
rel_losers <- imp_cat_df %>% filter(trend_change_lu == 4) %>% pull(species)
abs_losers <- imp_cat_df %>% filter(trend_change_lu == 5) %>% pull(species)

rw_raster_subset <- grep(paste(rel_winners, collapse="|"),  raster_list, value = TRUE)
rl_raster_subset <- grep(paste(rel_losers, collapse="|"),  raster_list, value = TRUE)
al_raster_subset <- grep(paste(abs_losers, collapse="|"),  raster_list, value = TRUE)

# load rasters:
rw_ranges <- rast(rw_raster_subset, lyrs = seq(1, 560, by = 7))
rl_ranges <- rast(rl_raster_subset, lyrs = seq(1, 560, by = 7))
al_ranges <- rast(al_raster_subset, lyrs = seq(1, 560, by = 7))

# sum across species:
rw_sum_cat <- sum(rw_ranges, na.rm = TRUE)
rl_sum_cat <- sum(rl_ranges, na.rm = TRUE)
al_sum_cat <- sum(al_ranges, na.rm = TRUE)

# add empty raster for absolute winners:
aw_sum_cat <- subst(rw_sum_cat, 1, NA)

lu_imp_cat <- c(aw_sum_cat, rw_sum_cat, rl_sum_cat, al_sum_cat)
names(lu_imp_cat) <- c("absolute winners", "relative winners", "relative losers", "absolute losers")

lu_cat_map <- ggplot() +
  geom_sf(data = US_albers_sf, colour = "black", fill = NA) +
  geom_spatraster(data = lu_imp_cat) +
  facet_wrap(~lyr, ncol = 2) +
  scale_fill_whitebox_c(palette = "viridi") +
  labs(fill = "N species",
       title = "Land use change impact") +
  theme_bw() +
  theme(text = element_text(size = 21),
        legend.key.height = unit(1, "cm"))
lu_cat_map

# add colours for facet strips:
g <- ggplot_gtable(ggplot_build(lu_cat_map))
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

ggsave(filename = file.path(plot_dir, "lu_change_impact_cats.svg"),
       plot = g,
       device = "svg",
       width = 29.7,
       height = 21, # A4
       units = "cm")


## 7) climate + land use change impact categories (Fig. S 6): ----

abs_winners <- imp_cat_df %>% filter(trend_change_climlu == 1) %>% pull(species)
rel_winners <- imp_cat_df %>% filter(trend_change_climlu == 2) %>% pull(species)
rel_losers <- imp_cat_df %>% filter(trend_change_climlu == 4) %>% pull(species)
abs_losers <- imp_cat_df %>% filter(trend_change_climlu == 5) %>% pull(species)

aw_raster_subset <- grep(paste(abs_winners, collapse="|"),  raster_list, value = TRUE)
rw_raster_subset <- grep(paste(rel_winners, collapse="|"),  raster_list, value = TRUE)
rl_raster_subset <- grep(paste(rel_losers, collapse="|"),  raster_list, value = TRUE)
al_raster_subset <- grep(paste(abs_losers, collapse="|"),  raster_list, value = TRUE)

# load rasters:
aw_ranges <- rast(aw_raster_subset, lyrs = seq(1, 560, by = 7))
rw_ranges <- rast(rw_raster_subset, lyrs = seq(1, 560, by = 7))
rl_ranges <- rast(rl_raster_subset, lyrs = seq(1, 560, by = 7))
al_ranges <- rast(al_raster_subset, lyrs = seq(1, 560, by = 7))

# sum across species:
aw_sum_cat <- sum(aw_ranges, na.rm = TRUE)
rw_sum_cat <- sum(rw_ranges, na.rm = TRUE)
rl_sum_cat <- sum(rl_ranges, na.rm = TRUE)
al_sum_cat <- sum(al_ranges, na.rm = TRUE)

climlu_imp_cat <- c(aw_sum_cat, rw_sum_cat, rl_sum_cat, al_sum_cat)
names(climlu_imp_cat) <- c("absolute winners", "relative winners", "relative losers", "absolute losers")

climlu_cat_map <- ggplot() +
  geom_sf(data = US_albers_sf, colour = "black", fill = NA) +
  geom_spatraster(data = climlu_imp_cat) +
  facet_wrap(~lyr, ncol = 2) +
  scale_fill_whitebox_c(palette = "viridi") +
  labs(fill = "N species",
       title = "Climate & land use change impact") +
  theme_bw() +
  theme(text = element_text(size = 21),
        legend.key.height = unit(1, "cm"))
climlu_cat_map

# add colours for facet strips:
g <- ggplot_gtable(ggplot_build(climlu_cat_map))
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

ggsave(filename = file.path(plot_dir, "climlu_change_impact_cats.svg"),
       plot = g,
       device = "svg",
       width = 29.7,
       height = 21, # A4
       units = "cm")

# session info:
writeLines(capture.output(sessionInfo()), file.path(dir, "results", "sessionInfo", "5_4b_attribution_maps.txt"))
