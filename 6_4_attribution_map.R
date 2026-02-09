# make maps based on overlaying BirdLife range maps of species used in attribution, 
# add relative importance values of drivers to see hotspots of climate change / land use change impact:

# packages: --------------------------------------------------------------------

library(sf)
library(dplyr)
library(doParallel)
library(gdalUtilities)
library(terra)
library(ggplot2)
library(tidyterra)

# directories: -----------------------------------------------------------------

main_dir <- file.path("//NAS-2-P-SN-01.ibb.uni-potsdam.de", "daten$", "AG26", "Transfer",
                      "Schifferle_BBS_occupancy_models_2023")

# Birdlife range maps:
datashare_Birdlife <- file.path("//NAS-2-P-SN-01.ibb.uni-potsdam.de/daten$", "AG26", "Arbeit", "datashare", "data", "biodat", "distribution", "Birdlife", "BOTW_2022")
#datashare_Birdlife <- file.path("/mnt", "ibb_share", "zurell", "biodat", "distribution", "Birdlife", "BOTW_2022")

# output directory:
# output_dir <- file.path("/mnt", "ibb_share", "zurell_transfer", "Schifferle_BBS_occupancy_models_2023", "data", "Birdlife_range_maps")
#output_dir <- file.path("data", "Birdlife_range_maps")
output_dir <- file.path(main_dir, "data", "Birdlife_range_maps")


# prepare data: ----------------------------------------------------------------

# selected species:
load(file.path("data", "species_DOM_val_okay.RData")) # output of 6_2_attribution_plots_trend_categories.R
spec_okay

# add scientific names:
load(file.path("data", "BBS_data_merged.RData")) # bbs_dt; output of: 
spec_names <- bbs_dt %>% 
  select(English_Common_Name, Scientific_Name) %>% 
  distinct %>% 
  filter(English_Common_Name %in% spec_okay)


## shapefile extraction from Birdlife gdb: ----

# register cores for parallel computation:
registerDoParallel(cores = 20)

# only parts of range where species is considered extant (presence = 1) and which
# is used either throughout the whole year (seasonal = 1) or during the breeding season (seasonal = 2)

# # account for taxonomic changes between BBS and BL range maps:
# # (suitable resources:
# # https://www.iucnredlist.org/
# # https://explorer.natureserve.org/Search)
spec_name_change_df <- data.frame("BBS_name" = "Dryobates villosus", "BL_name" = "Leuconotopicus villosus")
spec_name_change_df[2,] <- c("Dryocopus pileatus", "Hylatomus pileatus")

foreach(s = 2:nrow(spec_names),
        .packages = c("gdalUtilities"),
        .verbose = TRUE,
        .errorhandling = "remove",
        .inorder = FALSE) %dopar% {
          
          #spec <- sub("_", " ", spec_names$Scientific_Name[s])
          spec <- spec_names$Scientific_Name[s]
          
          # check if species was already processed:
          exists <- file.exists(file.path("data", "Birdlife_range_maps", paste0(spec_names$English_Common_Name[s], ".shp")))
          if(exists) {
            print(paste(spec, "ran already."))
            next
          }
          
          # adjust species name if BL uses other taxonomy than BBS:
          if(spec %in% spec_name_change_df$BBS_name){
            spec <- spec_name_change_df$BL_name[which(spec_name_change_df$BBS_name == spec)]
          }
          
          # extract shapefile:
          gdalUtilities::ogr2ogr(src_datasource_name = file.path(datashare_Birdlife, "BOTW.gdb"),
                                 layer = "All_Species",
                                 where = paste0("sci_name = '", spec, "' AND (SEASONAL = '1' OR SEASONAL = '2') AND PRESENCE = '1'"), # SEASONAL = 1: resident throughout the year, SEASONAL = 2: breeding season
                                 dst_datasource_name = file.path(output_dir, paste0(spec_names$English_Common_Name[s], ".shp")),
                                 overwrite = TRUE)
        }

## rasterize for easier summary calculations: ----

# clip to US, specify extent to ensure same grid for all, rasterize

# conterminous US:
US_albers_sf <- read_sf(file.path("data", "US_outline_ESRI102003.shp"))
# extent:
extent_US <- st_bbox(US_albers_sf)


for(s in 1:nrow(spec_names)){
  
  spec <- spec_names$English_Common_Name[s]
  
  print(paste(s, spec))
  
  range_sf <- read_sf(file.path(output_dir, paste0(spec, ".shp"))) %>% 
    # add species common name:
    mutate(species = spec_names$English_Common_Name[s])
  
  # clip to US:
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
                                at = TRUE, # xx to decide
                                a_nodata = -99999) # value for cells with missing data
  }


## add impact categories and driver importance as attributes: ----


# attribution metrics:

load(file = file.path(main_dir, "results", "attribution", "attribution_metrics_final.RData")) # output of 6_1_attribution_metrics.R
attr_metr_df

rel_imp_df <- attr_metr_df %>%
  select(-matches("(slope)|(p_)")) %>%
  mutate(imp_clim = mape_cfclim - mape_fact,
         imp_lu = mape_cflu - mape_fact,
         imp_climlu = mape_cfclimlu - mape_fact)

# impact categories:
load(file =  file.path(main_dir, "results", "attribution", "trend_categories.RData"))
flow_df

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


# species ranges:
# add info as multilayer tifs

for(s in 1:nrow(spec_names)){
  
  spec <- spec_names$English_Common_Name[s]
  
  print(paste(s, spec))
  
  range_tif <- rast(file.path(output_dir, "tifs",  paste0(spec, "_US_clipped.tif")))
  #plot(range_tif)
  
  range_tif$imp_clim <- ifel(range_tif == 1, rel_imp_df %>% filter(species == spec) %>% pull(imp_clim), NA)
  range_tif$imp_lu <- ifel(range_tif == 1, rel_imp_df %>% filter(species == spec) %>% pull(imp_lu), NA) 
  range_tif$imp_climlu <- ifel(range_tif == 1, rel_imp_df %>% filter(species == spec) %>% pull(imp_climlu), NA)
  
  range_tif$cat_clim <- ifel(range_tif == 1, imp_cat_df %>% filter(species == spec) %>% pull(trend_change_clim), NA) 
  range_tif$cat_lu <- ifel(range_tif == 1, imp_cat_df %>% filter(species == spec) %>% pull(trend_change_lu), NA) 
  range_tif$cat_climlu <- ifel(range_tif == 1, imp_cat_df %>% filter(species == spec) %>% pull(trend_change_climlu), NA) 
  
  writeRaster(range_tif, file.path(output_dir, "tifs_attr",  paste0(spec, "_attr.tif")), overwrite = TRUE)
  
}


# maps: ------------------------------------------------------------------------

# elevation:
US_contours_subset_simpl <- st_read(file.path("data", "US_elev_contours_simplified.shp")) # output of 5_0_DOM_predictions_USA.R

## 1) general species richness map (locations of ranges of the species in our selection): ----

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

ggsave(filename = file.path("plots", "attribution", "maps", "species_richness.svg"), 
       plot = sr_map,
       device = "svg",
       width = 29.7,
       height = 18, # A4
       units = "cm")

# load(file = file.path("data", "species_ecoregions.RData"))
# spec_eco_df %>% 
#   filter(species %in% spec_okay) %>% 
#   group_by(eco) %>% 
#   summarise(n = n())


## 2) climate change mean importance: ----

# use color scale matching boxplots:
# viridis purple for climate change importance
clim_col <- viridis::viridis(n = 3)[1]


raster_list <- list.files(file.path(output_dir, "tifs_attr"), pattern = ".tif$", full.names = TRUE)
ranges_clim <- rast(raster_list, lyrs = seq(2, 560, by = 7)) # imp_clim

mean_clim_imp <- mean(ranges_clim, na.rm = TRUE)

clim_imp_map <- ggplot() +
  geom_spatraster(data = mean_clim_imp) +
  geom_sf(data = US_contours_subset_simpl, aes(colour = Contour), linewidth = 0.2) +
  scale_colour_gradientn(colours = c(terrain.colors(8)[-9], "grey80"),
                         transform = "sqrt", 
                         breaks = c(100, 200, 500, 700, 1000, 1500, 2100, 2500, 2900)) +
  scale_fill_gradient(low = "#FFFFFF", high = clim_col, na.value = "transparent") +
  labs(fill = "mean range-level importance") +
  theme_bw() +
  guides(colour = "none") +
  theme(text = element_text(size = 24),
        legend.key.height = unit(1, "cm"),
        legend.key.width = unit(.2, "npc"),
        plot.margin = margin(unit(c(3,3,0,3), "lines")),
        legend.position = "bottom",
        legend.title.position = "top") # trbl # top margin to later have label when arranging multiple plots)
clim_imp_map



ggsave(filename = file.path("plots", "attribution", "maps", "mean_climate_change_importance.svg"), 
       plot = clim_imp_map,
       device = "svg",
       width = 29.7,
       height = 21, # A4
       units = "cm")

## 3) land use change mean importance: ----

# use color scale matching boxplots:
# viridis petrol for land use change importance
lu_col <- viridis::viridis(n = 3)[2]

ranges_lu <- rast(raster_list, lyrs = seq(3, 560, by = 7))

mean_lu_imp <- mean(ranges_lu, na.rm = TRUE)

lu_imp_map <- ggplot() +
  geom_spatraster(data = mean_lu_imp) +
  geom_sf(data = US_contours_subset_simpl, aes(colour = Contour), linewidth = 0.2) +
  scale_colour_gradientn(colours = c(terrain.colors(8)[-9], "grey80"),
                         transform = "sqrt", 
                         breaks = c(100, 200, 500, 700, 1000, 1500, 2100, 2500, 2900),
                         name = "elevation [m]") +
  scale_fill_gradient(low = "#FFFFFF", high = lu_col, na.value = "transparent") +
  labs(fill = "mean range-level importance") +
  guides(colour = "none") +
  theme_bw() +
  theme(text = element_text(size = 24),
        legend.key.height = unit(1, "cm"),
        legend.key.width = unit(.2, "npc"),
        plot.margin = margin(unit(c(3,3,0,10), "lines")),
        legend.position = "bottom",
        legend.title.position = "top")
lu_imp_map

ggsave(filename = file.path("plots", "attribution", "maps", "mean_lu_change_importance.svg"), 
       plot = lu_imp_map,
       device = "svg",
       width = 29.7,
       height = 21, # A4
       units = "cm")

## 3.1) is mean importance of land use change anywhere larger than mean importance of climate change?: ----

# difference:
diff_clim_lu_rast <- mean_clim_imp - mean_lu_imp
diff_clim_lu_rast$binary <- mean_lu_imp > mean_clim_imp
# mark raster cells where lu > clim:
v <- terra::as.polygons(diff_clim_lu_rast, aggregate = FALSE) %>% 
  filter(binary == 1)

diff_clim_lu_rast_map <- ggplot() +
  geom_spatraster(data = diff_clim_lu_rast, aes(fill = mean), ) +
  geom_spatvector(data = v, fill = NA) +
  scale_fill_gradient2(name = "difference mean range-level importance (climate - land use change)", 
                       low = lu_col, mid = "white",high = clim_col, na.value = "transparent") +
  theme_bw() +
  theme(text = element_text(size = 24),
        legend.key.height = unit(1, "cm"),
        legend.key.width = unit(.1, "npc"),
        plot.margin = margin(unit(c(3,3,0,10), "lines")),
        legend.position = "bottom",
        legend.title.position = "top",
        legend.box = "vertical",
        legend.box.just = "left",
        plot.tag.position = c(unit(0.08, "npc"), unit(0.03, "npc")),
        plot.tag = element_text(hjust = 0, size = 22)
        ) +
  guides(test = guide_custom(grob = grid::rectGrob(width = unit(1, "cm"), height = unit(1, "cm"),
                                                   x = unit(0.5, "npc"), y = unit(0.5, "npc")))) +
  labs(tag = "land use change more important than climate change") # xx wording not ideal
diff_clim_lu_rast_map

ggsave(filename = file.path("plots", "attribution", "maps", "clim_vs_lu_importance.svg"), 
       plot = diff_clim_lu_rast_map,
       device = "svg",
       width = 29.7,
       height = 21, # A4
       units = "cm")


## 3.2) for manuscript: combine climate and land use in one plot: ----

library(cowplot)
combined_plot <- plot_grid(clim_imp_map, lu_imp_map,
                           align = 'vh',
                           labels = c("A: climate change", "B: land use change"),
                           label_size = 26, 
                           nrow = 1,
                           hjust = -0.01)

combined_plot

# plot to get legend from:
leg_plot <- ggplot() +
  geom_spatraster(data = mean_clim_imp) +
  geom_sf(data = US_contours_subset_simpl, aes(colour = as.factor(Contour)), linewidth = 0.05) +
  scale_colour_manual(values = c(terrain.colors(8)[-9], "grey80"),
                      name = "elevation [m]") +
  theme_bw() +
  theme(text = element_text(size = 24),
        legend.position = "bottom",
        legend.box.margin = margin(0, 0, 0, 0),
        legend.margin = margin(0, 0, 0, 0)) +
  guides(colour = guide_legend(override.aes = list(linewidth = 4), nrow = 1, title.position = "top"),
         fill = "none")
leg_plot

# extract the legend:
legend <- get_legend(leg_plot)

combined_plot_final <- plot_grid(combined_plot, legend, 
                                 nrow = 2,
                                 rel_heights = c(1, .1))
combined_plot_final

# does make new plots with 
ggsave(filename = file.path("plots", "attribution", "maps", "clim_lu_combined.svg"), 
       plot = combined_plot_final,
       device = "svg",
       width = 45,
       height = 23,
       units = "cm")

## 4) climate + land use change mean importance: ----

ranges_climlu <- rast(raster_list, lyrs = seq(4, 560, by = 7))

mean_climlu_imp <- mean(ranges_climlu, na.rm = TRUE)

climlu_imp_map <- ggplot() +
  geom_spatraster(data = mean_climlu_imp) +
  scale_fill_whitebox_c(palette = "viridi") +
  labs(fill = "",
       title = "Mean climate & land use change importance",
       subtitle = "based on BirdLife ranges of 80 species used in attribution") +
  theme_bw() +
  theme(text = element_text(size = 21),
        legend.key.height = unit(1, "cm"))
climlu_imp_map

ggsave(filename = file.path("plots", "attribution", "maps", "mean_climlu_change_importance.svg"), 
       plot = climlu_imp_map,
       device = "svg",
       width = 29.7,
       height = 21, # A4
       units = "cm")


## 5) climate change impact categories: ----

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

ggsave(filename = file.path("plots", "attribution", "maps", "clim_change_impact_cats.svg"), 
       plot = g,
       device = "svg",
       width = 29.7,
       height = 21, # A4
       units = "cm")

## 6) land use change impact: ----

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

ggsave(filename = file.path("plots", "attribution", "maps", "lu_change_impact_cats.svg"), 
       plot = g,
       device = "svg",
       width = 29.7,
       height = 21, # A4
       units = "cm")

## 7) climate + land use change impact categories: ----

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


ggsave(filename = file.path("plots", "attribution", "maps", "climlu_change_impact_cats.svg"), 
       plot = g,
       device = "svg",
       width = 29.7,
       height = 21, # A4
       units = "cm")
