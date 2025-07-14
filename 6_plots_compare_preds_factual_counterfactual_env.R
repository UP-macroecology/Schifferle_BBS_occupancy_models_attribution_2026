# time series (and maps) to compare DOM predictions for factual as well as for
# counterfactual environmental data:

# packages: --------------------------------------------------------------------

library(dplyr)
library(sf)
library(terra)
library(foreach)
library(doParallel)
library(ggplot2)
library(tidyterra)
library(ggrepel)
library(ggtext)


# functions: -------------------------------------------------------------------

source("0_functions.R")


# directories: -----------------------------------------------------------------

res_dir <- file.path("T:", "Schifferle_BBS_occupancy_models_2023", "results")
#res_dir <- file.path("//NAS-2-P-SN-01.ibb.uni-potsdam.de", "daten$", "AG26", "Transfer", "Schifferle_BBS_occupancy_models_2023", "results")


# load data: -------------------------------------------------------------------

# selected species:
load(file = file.path("data", "species_set_analysis.RData")) # output of 3_1_DOM_CV_evaluation_metrics.R
final_species

# selected variables:
load(file = file.path("data", "selected_variables.RData")) # selvar_final; output of 1_2_variable_selection.R
selvar_final

# routes-years:
load(file = file.path("data", "BBS_for_occ_selection.RData")) # route_sel_dt; output of 1_3_match_BBS_to_env_data.R 

# selected routes spatial data (to buffer presences):
routes_sel_sf <- st_read(file.path("data", "route_selection_1995_2019_surv_beg_end_max_5y_miss_v2_spat_thin_100km_max_30_r_per_BCR_centroids.shp")) # output of 1_1_route_selection.R

# route-year-species information (only surveyed)
load(file = file.path("data", "BBS_for_occ_spec_records.RData")) # bbs_dt_occ; output of 1_0_reformat_BBS_data.R

# environmental data for each US grid cell:
load(file = file.path("data", "US_grid_env_data.RData")) # output of 1_3_env_data_df_contUS.R; clim_lu_cells_sf

# US cells:
US_cells_sf <- read_sf(file.path("data", "cell_centroids_US_ESRI102003.shp")) # output of 1_3_env_data_df_contUS.R 
US_cells_sf # 4149 cells


# compare time series of sum of occupied routes: -------------------------------
# observations vs. factual vs. counterfactual:

#res_dir <- file.path("T:", "Schifferle_BBS_occupancy_models_2023", "results")
res_dir <- file.path("//NAS-2-P-SN-01.ibb.uni-potsdam.de", "daten$", "AG26", "Transfer", 
                     "Schifferle_BBS_occupancy_models_2023", "results")

#lu_scen <- "1995soc"
#lu_scen <- "1901soc"

# directory counterfactual predictions:
#cf_dir <- file.path(res_dir, "attribution", paste0("fm_y_preds_routes_cf_", lu_scen))
cf_dir <- file.path(res_dir, "attribution", "fm_y_preds_routes_cf_1995_all")

# directory to save plots:
#plot_dir <- file.path("plots", "attribution", "time_series_y_preds_cf")
#plot_dir <- file.path("plots", "attribution", paste0("fm_y_preds_routes_cf_", lu_scen))
plot_dir <- file.path("plots", "attribution", "fm_y_preds_routes_cf_1995_all")
if(!dir.exists(plot_dir)){dir.create(plot_dir, recursive = TRUE)} # xx


for(s in 108:length(final_species)){ # 108!
  
  spec <- final_species[s]
  
  print(paste(s, spec))
  
  # observations:
  rel_routes <- training_routes(species = spec, buffer_km = 750, output = "RTENOs")
  occ_dt_spec <- BBS_pres_abs_spec(species = spec) %>%
    filter(RTENO %in% rel_routes)
  # route-level presence:
  # sum all routes for each year (temporal trend)
  obs_temp_trend <- occ_dt_spec %>% #xx
    group_by(Year) %>% #xx
    summarise(pres_sum = sum(presence, na.rm = TRUE)) #xx

  
  # factual predictions:
  
  # check where to look for model output (did MCMC fitting work with less or only with more iterations?)
  
  if(file.exists(file.path(res_dir, "fm_buffer750km", "refit_2000_2000", paste0("out_", spec, "_fm_buffer_750.RData")))){
    output_dir <- file.path(res_dir, "fm_buffer750km", "refit_2000_2000")
  } else {
    output_dir <- file.path(res_dir, "fm_buffer750km")
  }
  
  #load(file.path(output_dir, paste0("postproc_", spec, "_fm_buffer750.RData")))
  # sum across route sections:
  #preds_routes <- apply(res_list$y_preds, MAR = c(1,3,4), FUN = max) # as soon as one predicted presence -> 1
  # save this:
  #save(preds_routes, file = file.path(res_dir, "fm_buffer750km", "y_preds_route_level_section_sum", paste0(spec, "_y_preds_route_level_section_sum.RData")))

  load(file.path(res_dir, "fm_buffer750km", "y_preds_route_level_section_sum", paste0(spec, "_y_preds_route_level_section_sum.RData")))
  
  # sum across routes for each year:
  preds_years <- apply(preds_routes, MAR = c(2,3), FUN = sum, na.rm = TRUE) # xx
  # mean:
  preds_years_mean_f <- apply(preds_years, MAR = 1, FUN = mean)
  
  
  # calculate this time series for all beforehand and save xx
  
  # counterfactual predictions:
  
  time_series_df <- data.frame(spec = spec, year = 1995:2019)
  time_series_df$obs <- obs_temp_trend$pres_sum
  time_series_df$fact <- preds_years_mean_f
  
  # iterate over counterfactual scenarios:

  cf_files <- list.files(cf_dir, pattern = spec)
  
  for(v in 1:length(cf_files)){
    
    print(v)
    
    # extract scenario from file name:
    scen <- gsub(cf_files[v], pattern = paste0("(", paste0(spec, "_y_preds_cf_"), ")|(.RData)"), replacement = "")
    #scen <- gsub(cf_files[v], pattern = paste0("(", paste0(spec, "_y_preds_cf_"), ")|(_occ1fact.RData)"), replacement = "")
    
    load(file.path(cf_dir, cf_files[v]))
    
    # sum across routes for each year:
    preds_years2 <- apply(y_preds_route_cf, MAR = c(2,3), FUN = sum, na.rm = TRUE) # xx
    # mean:
    preds_years_mean2 <- apply(preds_years2, MAR = 1, FUN = mean)
    
    time_series_df[[scen]] <- preds_years_mean2
  }
  
  time_series_df
  
  # compare with using predictions for initial occupancy based on factual data
  # ggplot, mark lines:
  # https://r-graph-gallery.com/web-line-chart-with-labels-at-end-of-line.html
  
  
  time_series_df_lf <- time_series_df %>% 
    tidyr::pivot_longer(cols = -c(spec, year), names_to = "scenario", values_to = "N_routes") %>% 
    mutate(scenario = factor(scenario),
           name_lab = if_else(year == 2019, scenario, NA_character_),
           name_lab = case_match(name_lab,
                                 "obs" ~ "observations",
                                 "fact" ~ "all factual",
                                 "counterclim" ~ "all climate counterfact.",
                                 "1995soc" ~ "all land use counterfact.",
                                 "counterclim_1995soc" ~ "all counterfact."))
  

  ## linetype:
  line_type1 <- rep(1, 2)
  names(line_type1) <- c("obs", "fact")
  line_type2 <- rep(5, 3)
  names(line_type2) <- c("counterclim", "1995soc", "counterclim_1995soc")
  line_type <- c(line_type1, line_type2)

  ## colours:
  cols <- c("counterclim" = "#0D98BA", 
            "1995soc" = "#B7410E", 
            "counterclim_1995soc" = "#046865", 
            "obs" = "black", 
            "fact" = "#85CB33")  
  
  
  p <- ggplot(data = time_series_df_lf, 
              aes(x = year, y = N_routes, group = scenario)) +
    # geometric annotations that play the role of grid lines (to avoud grid lines on the right where labels are)
    geom_vline(xintercept = seq(1995, 2020, by = 5), color = "grey90", linewidth = .6) +
    geom_segment(data = tibble(y = seq(round(min(time_series_df_lf$N_routes), digits = -1), 
                                       round(max(time_series_df_lf$N_routes), digits = -1), length = 5), # xx
                               x1 = 1995, x2 = 2020),
                 aes(x = x1, xend = x2, y = y, yend = y), inherit.aes = FALSE, color = "grey90", linewidth = .6) +
    # data lines:
    geom_line(aes(colour = scenario, linetype = scenario), linewidth = 1) +
    geom_point(data = time_series_df_lf %>%  filter(scenario == "obs"), 
               aes(colour = scenario), size = 3) +
    ylab("N routes with presences") +
    theme_bw() +
    ggtitle(spec) +
    geom_text_repel(
      aes(color = scenario, label = name_lab),
      fontface = "bold", 
      size = 6, #6
      direction = "y", xlim = c(2020, NA), hjust = 0,
      segment.size = .7, segment.alpha = .5, segment.linetype = 1,
      #box.padding = .4,
      segment.curvature = -0.1, segment.ncp = 3, segment.angle = 20,
      max.overlaps = 30
    ) +
    scale_x_continuous(expand = c(0, 0), limits = c(1995, 2027)) +
    theme(legend.position = "none", panel.grid = element_blank(),
          text = element_text(size = 25)) +
    scale_linetype_manual(values = line_type) +
    scale_colour_manual(values = cols)
   p 
   
   jpeg(file = file.path(plot_dir, paste0(spec, "_ts_y_preds_cf_1995_all.jpg")), 
        width = 1000, height = 700, quality = 100)
   print(p)
   dev.off()
  
    # other versions:
    # 
    # time_series_df_lf <- time_series_df %>% 
    #   tidyr::pivot_longer(cols = -c(spec, year), names_to = "scenario", values_to = "N_routes") %>% 
    #   mutate(scenario = factor(scenario),
    #          name_lab = if_else(year == 2019, scenario, NA_character_),
    #          name_lab = case_match(name_lab,
    #                            "obs" ~ "observations",
    #                            "fact" ~ "all factual",
    #                            "bio1" ~ "annual mean temp.",
    #                            "bio2" ~ "mean diurnal range",
    #                            "bio3" ~ "isothermality",
    #                            "bio7" ~ "temp. annual range",
    #                            "bio14" ~ "prec. driest month",
    #                            "bio15" ~ "prec. seasonality",
    #                            "pr_mean_spring" ~ "prec. spring",
    #                            "pr_mean_summer" ~ "prec. summer",
    #                            "pr_mean_autumn" ~ "prec. autumn",
    #                            "pr_mean_winter" ~ "prec. winter",
    #                            "counterclim" ~ "all climate counterfact.",
    #                            "urbanareas" ~ "urban area",
    #                            "sum_annual_crops" ~ "annual crops",
    #                            "primary_nonforests" ~ "primary non-forest",
    #                            "secondary_nonforests" ~ "secondary non-forest",
    #                            "managed_pastures" ~ "pasture",
    #                            lu_scen ~ "all land use counterfact.",
    #                            paste0("counterclim_", lu_scen) ~ "all counterfact."))
    # 
    # # specify values for aesthetics:
    # 
    # ## colours:
    # cols_clim <- c("bio1" = "#E76254FF", "bio2" = "#F7AA58FF", "#FFD06FFF", 
    #                "bio7" = "#EF8A47FF", "bio14" = "#1E88E5FF", "bio15" = "#1976D2FF", 
    #                "pr_mean_spring" = "#AADCE0FF", "pr_mean_summer" = "#72BCD5FF",
    #                "pr_mean_autumn" = "#528FADFF", "pr_mean_winter" = "#376795FF",  
    #                "counterclim" = "#1565C0FF")
    # cols_lu <- c("managed_pastures" = "#88AB38FF", "sum_annual_crops" = "#B4BF3AFF", 
    #              "urbanareas" = "#7F793CFF", "secondary_nonforests" = "#3B7D31FF", 
    #              "primary_nonforests" = "#225F2FFF",
    #              "lu_scen" = "#775B24FF",
    #              "counterclim_lu" = "#B86092FF")
    # names(cols_lu)[which(names(cols_lu) == "lu_scen")] <- lu_scen
    # names(cols_lu)[which(names(cols_lu) == "counterclim_lu")] <- paste0("counterclim_", lu_scen)
    # 
    # cols <- c(cols_clim, cols_lu, "obs" = "black", "fact" = "#FFCD12FF")
    # 
    # ## linewidth:
    # line_width1 <- rep(0.8, 5)
    # names(line_width1) <- c("obs", "fact", "counterclim", lu_scen, paste0("counterclim_", lu_scen))
    # line_width2 <- rep(0.2, 15)
    # names(line_width2) <- selvar_final
    # line_width <- c(line_width1, line_width2)
    # 
    # ## linetype:
    # line_type1 <- rep(1, 2)
    # names(line_type1) <- c("obs", "fact")
    # line_type2 <- rep(5, 18)
    # names(line_type2) <- c(selvar_final, "counterclim", lu_scen, paste0("counterclim_", lu_scen))
    # line_type <- c(line_type1, line_type2)
    # 
    # # plot:
    # p <- ggplot(data = time_series_df_lf, 
    #             aes(x = year, y = N_routes, group = scenario)) +
    #   # geometric annotations that play the role of grid lines (to avoud grid lines on the right where labels are)
    #   geom_vline(xintercept = seq(1995, 2020, by = 5), color = "grey90", linewidth = .6) +
    #   geom_segment(data = tibble(y = seq(round(min(time_series_df_lf$N_routes), digits = -1), 
    #                                      round(max(time_series_df_lf$N_routes), digits = -1), length = 5), # xx
    #                              x1 = 1995, x2 = 2020),
    #                aes(x = x1, xend = x2, y = y, yend = y), inherit.aes = FALSE, color = "grey90", linewidth = .6) +
    #   # data lines:
    #   geom_line(aes(colour = scenario, linewidth = scenario, linetype = scenario)) +
    #   ylab("N routes with presences") +
    #   theme_bw() +
    #   ggtitle(paste(spec, "land use", lu_scen)) +
    #   geom_text_repel(
    #     aes(color = scenario, label = name_lab),
    #     fontface = "bold", 
    #     size = 6,
    #     direction = "y", xlim = c(2020, NA), hjust = 0,
    #     segment.size = .7, segment.alpha = .5, segment.linetype = 1,
    #     #box.padding = .4,
    #     segment.curvature = -0.1, segment.ncp = 3, segment.angle = 20,
    #     max.overlaps = 30
    #   ) +
    #   scale_x_continuous(expand = c(0, 0), limits = c(1995, 2027)) +
    #   theme(legend.position = "none", panel.grid = element_blank(),
    #         text = element_text(size = 20))
  
}


# play around with style:

# specify values for aesthetics:

## colours:
cols_clim <- c("bio1" = "#E76254FF", "bio2" = "#F7AA58FF", "#FFD06FFF", 
               "bio7" = "#EF8A47FF", "bio14" = "#1E88E5FF", "bio15" = "#1976D2FF", 
               "pr_mean_spring" = "#AADCE0FF", "pr_mean_summer" = "#72BCD5FF",
               "pr_mean_autumn" = "#528FADFF", "pr_mean_winter" = "#376795FF",  
               "counterclim" = "#1565C0FF")
cols_lu <- c("managed_pastures" = "#88AB38FF", "sum_annual_crops" = "#B4BF3AFF", 
             "urbanareas" = "#7F793CFF", "secondary_nonforests" = "#3B7D31FF", 
             "primary_nonforests" = "#225F2FFF",
             "1901soc" = "#775B24FF",
             "counterclim_1901soc" = "#B86092FF")
cols <- c(cols_clim, cols_lu, "obs" = "black", "fact" = "#FFCD12FF")

#c("#AADCE0FF", "#72BCD5FF", "#528FADFF", "#376795FF", "#1E466EFF")
#c("#E76254FF", "#EF8A47FF", "#F7AA58FF", "#FFD06FFF") # hiroshige
#c("#7D9D33FF", "#CED38CFF", "#DCC949FF", "#BCA888FF", "#CD8862FF", "#775B24FF") # kakapo
#c("#EBCF2EFF", "#B4BF3AFF", "#88AB38FF", "#5E9432FF", "#3B7D31FF", "#225F2FFF", "#244422FF", "#252916FF") # alkalay2


# versions for presentation:

# specify values for aesthetics:

## colours:
cols_clim <- c("bio1" = "#E76254FF", "bio2" = "#F7AA58FF", "#FFD06FFF",
               "bio7" = "#EF8A47FF", "bio14" = "#1E88E5FF", "bio15" = "#1976D2FF",
               "pr_mean_spring" = "#AADCE0FF", "pr_mean_summer" = "#72BCD5FF",
               "pr_mean_autumn" = "#528FADFF", "pr_mean_winter" = "#376795FF",
               "counterclim" = "#0D98BA")

cols_lu <- c("managed_pastures" = "#88AB38FF", "sum_annual_crops" = "#B4BF3AFF", 
             "urbanareas" = "#7F793CFF", "secondary_nonforests" = "#3B7D31FF", 
             "primary_nonforests" = "#225F2FFF",
             "lu_scen" = "#B7410E",
             "counterclim_lu" = "#046865")

# reduced version:
cols_clim <- c("bio1" = "#87CEFA", "bio2" = "#87CEFA", "#87CEFA",
               "bio7" = "#87CEFA", "bio14" = "#87CEFA", "bio15" = "#87CEFA",
               "pr_mean_spring" = "#87CEFA", "pr_mean_summer" = "#87CEFA",
               "pr_mean_autumn" = "#87CEFA", "pr_mean_winter" = "#87CEFA",
               "counterclim" = "#0D98BA")  #6EB4D1
cols_lu <- c("managed_pastures" = "#F3885A", "sum_annual_crops" = "#F3885A",
             "urbanareas" = "#F3885A", "secondary_nonforests" = "#F3885A",
             "primary_nonforests" = "#F3885A",
             "lu_scen" = "#B7410E",
             "counterclim_lu" = "#046865") #C47335
names(cols_lu)[which(names(cols_lu) == "lu_scen")] <- lu_scen
names(cols_lu)[which(names(cols_lu) == "counterclim_lu")] <- paste0("counterclim_", lu_scen)

cols <- c(cols_clim, cols_lu, "obs" = "black", "fact" = "#85CB33")

## linewidth:
line_width1 <- rep(1.2, 5)
names(line_width1) <- c("obs", "fact", "counterclim", lu_scen, paste0("counterclim_", lu_scen))
line_width2 <- rep(0.2, 15)
names(line_width2) <- selvar_final
line_width <- c(line_width1, line_width2)

## linetype:
line_type1 <- rep(1, 5)
names(line_type1) <- c("obs", "fact", lu_scen, "counterclim", paste0("counterclim_", lu_scen))
line_type2 <- rep(5, 15)
names(line_type2) <- c(selvar_final)
line_type <- c(line_type1, line_type2)

# plot:
p <- time_series_df_lf %>% 
  #filter(scenario %in% c("counterclim", "1901soc", "obs", "fact", "counterclim_1901soc")) %>% 
ggplot(aes(x = year, y = N_routes, group = scenario)) +
  # geometric annotations that play the role of grid lines (to avoud grid lines on the right where labels are)
  geom_vline(xintercept = seq(1995, 2020, by = 5), color = "grey90", linewidth = .6) +
  geom_segment(data = tibble(y = seq(round(min(time_series_df_lf$N_routes), digits = -1), 
                                     round(max(time_series_df_lf$N_routes), digits = -1), length = 5), # xx
                             x1 = 1995, x2 = 2020),
               aes(x = x1, xend = x2, y = y, yend = y), inherit.aes = FALSE, color = "grey90", linewidth = .6) +
  # data lines:
  geom_line(aes(colour = scenario, linewidth = scenario, linetype = scenario)) +
  geom_point(data = time_series_df_lf %>%  filter(scenario == "obs"), 
             aes(colour = scenario), size = 4) +
  ylab("N routes with presences") +
  theme_bw() +
  ggtitle(paste0(spec, ", counterfactual data ", substr(lu_scen, start = 1, stop = 4))) +
  geom_text_repel(
    # reduced version:
    data = time_series_df_lf %>% filter(scenario %in% c("counterclim", lu_scen, "obs", "fact", paste0("counterclim_", lu_scen))),
    aes(color = scenario, label = name_lab),
    fontface = "plain", 
    size = 9,
    direction = "y", xlim = c(2020, NA), hjust = 0,
    segment.size = .7, segment.alpha = .5, segment.linetype = 1,
    #box.padding = .4,
    segment.curvature = -0.1, segment.ncp = 3, segment.angle = 20,
    max.overlaps = 30
  ) +
  scale_x_continuous(expand = c(0, 0), limits = c(1995, 2027)) +
  theme(legend.position = "none", panel.grid = element_blank(),
        text = element_text(size = 40)) +
  scale_colour_manual(values = cols) +
  scale_linewidth_manual(values = line_width) +
  scale_linetype_manual(values = line_type)
p

jpeg(file = file.path(plot_dir, paste0(spec, "_ts_y_preds_cf_", lu_scen, "2.jpg")), 
     width = 1400, height = 1000, quality = 100)
print(p)
dev.off()



# compare maps of predictions for factual data vs. counterfactual data: --------

# occupancy map 1995 factual, counterfactual, difference

spec <- final_species[123]

# find cells within 750 km buffer around presences:
spec_pres_buffer_sf <- training_routes(species = spec, buffer_km = 750, output = "buffer")

cells_in_buffer_sf <- clim_lu_cells_sf %>%
  st_filter(., spec_pres_buffer_sf)

# predictions to factual data:
load(file.path(res_dir, "fm_preds_US", paste0(spec, "_US_occ_preds.RData")))

# mean occupancy:
occ_start_stop_sf <- cells_in_buffer_sf %>% 
  filter(year == 1995) %>% # same cellIDs in each year, doesn't matter which
  select(cellID) %>% 
  mutate(occ_mean_1995 = apply(Z[,1,], MARGIN = 1, FUN = mean),
         occ_mean_2019 = apply(Z[,25,], MARGIN = 1, FUN = mean))

# predictions to counterfactual data:

for(i in selvar_final){
  
  print(i)
  
  load(file.path(res_dir, "fm_preds_US", "fm_preds_US_counterfactual", 
                 paste0(spec, "_occ_preds_cf_", i, ".RData")))
  
  occ_start_stop_sf <- occ_start_stop_sf %>% 
    mutate("occ_mean_cf_{i}_1995" := apply(Z_cf[,1,], MARGIN = 1, FUN = mean) / 100,
           "occ_mean_cf_{i}_2019" := apply(Z_cf[,25,], MARGIN = 1, FUN = mean) / 100,
           # difference between factal and counterfactual:
           "diff_cf_{i}_1995" := occ_mean_1995 - !!as.name(paste0("occ_mean_cf_",i,"_1995")),
           "diff_cf_{i}_2019" := occ_mean_2019 - !!as.name(paste0("occ_mean_cf_",i,"_2019")))
}

occ_start_stop_sf


# convert to raster:
occ_start_stop_rast <- occ_start_stop_sf %>% 
  mutate(x = st_coordinates(.)[,1],
         y = st_coordinates(.)[,2]) %>% 
  select(x, y, everything(), -cellID) %>% 
  st_drop_geometry() %>% 
  rast(., type='xyz', crs=crs(clim_lu_cells_sf))

# occ. and difference in two plots:

# 1) plot mean occupancy factual vs. counterfactual scenarios:

ggplot() +
  geom_spatraster(data = occ_start_stop_rast %>% select(starts_with("occ"))) +
  #geom_sf(data = elev_dt, aes(colour = Contour), linewidth = 0.2) +
  facet_wrap(~lyr, ncol = 7) +
  scale_fill_viridis_c(na.value = "transparent", option = "plasma", limits = c(0, 1)) +
  scale_colour_gradientn(colours = c(terrain.colors(8)[-8], "grey80"),
                         transform = "sqrt", guide = "none") +
  labs(fill = "occ. prob.") +
  theme_bw() +
  theme(plot.margin = unit(c(0.25, 0.25, 0.25, 0.25), "cm"))

# 2) plot differences:

limit_df <- occ_start_stop_sf %>% 
  select(starts_with("diff")) %>% 
  st_drop_geometry()
limit <-abs(max(range(limit_df)))  * c(-1, 1) # to center colour scale

ggplot() +
  geom_spatraster(data = occ_start_stop_rast %>% select(starts_with("diff"))) +
  #geom_sf(data = elev_dt, aes(colour = Contour), linewidth = 0.2) +
  facet_wrap(~lyr, ncol = 7) + 
  scale_fill_distiller("diff.", type = "div", 
                       palette = "RdBu", limit = limit, na.value = "transparent",
                       transform = LightLogR::symlog_trans(base = 2, thr = 0.01),
                       breaks = c(-0.5, -0.1, 0, 0.1, 0.5)) +
  #scale_colour_gradientn(colours = c(terrain.colors(8)[-8], "grey80"),
  #                      transform = "sqrt", guide = "none") +
  labs(fill = "difference") +
  theme_bw() +
  ggtitle(spec) +
  theme(plot.margin = unit(c(0.25, 0.25, 0.25, 0.25), "cm"))

