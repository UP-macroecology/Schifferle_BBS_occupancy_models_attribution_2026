# use fitted DOMs to predict occupancy (colonisation and extinction) probability
# across conterminous US (maps):


# packages: --------------------------------------------------------------------

library(doParallel)
library(dplyr)
library(flocker)
library(cmdstanr)
set_cmdstan_path(path = NULL)
#set_cmdstan_path("C:/Users/schifferle1/Documents/cmdstan-2.34.1")
library(sf)
library(terra)
library(tidyterra)
library(ggplot2)
library(cowplot)

# functions: -----

source("0_functions.R")

# register cores for parallel computation:
ncores <- 2
cl <- makeCluster(ncores, setup_timeout = 0.5)
registerDoParallel(cl)

# directories: ----

print(tempdir())

#dir <- file.path("/import", "ecoc9z", "data-zurell", "schifferle", "BBS_occupancy_models_2023")
dir <- file.path("M:", "Documents", "DEBTs", "analysis", "Schifferle_BBS_occupancy_models_2023")
#dir <- getwd()

# logfiles:
log_dir <- file.path("logfiles", "fm_preds_US")

# results_dir <- file.path("//NAS-2-P-SN-01.ibb.uni-potsdam.de", "users$", "schifferle1", "Documents", "DEBTs", "analysis", "Schifferle_BBS_occupancy_models_2023",
#                      "results")
results_dir <- file.path(dir, "results", "fm_buffer750km")


# load data: -------------------------------------------------------------------

# routes-years:
load(file = file.path("data", "BBS_for_occ_selection.RData")) # route_sel_dt; output of 1_3_match_BBS_to_env_data.R 

# route-year-species information (only surveyed)
load(file = file.path("data", "BBS_for_occ_spec_records.RData")) # bbs_dt_occ; output of 1_0_reformat_BBS_data.R

# selected routes spatial data (to buffer presences):
routes_sel_sf <- st_read(file.path("data", "route_selection_1995_2019_surv_beg_end_max_5y_miss_v2_spat_thin_100km_max_30_r_per_BCR_centroids.shp")) # output of 1_1_route_selection.R

# selected variables:
load(file = file.path("data", "selected_variables.RData")) # selvar_final; output of 1_2_variable_selection.R

# selected species:
load(file = file.path("data", "species_set_analysis.RData"))
final_species

# environmental data for each US grid cell:

load(file = file.path("data", "US_grid_env_data.RData")) # output of 1_3_env_data_df_contUS.R; clim_lu_cells_sf
clim_lu_cells_df <- sf::st_drop_geometry(clim_lu_cells_sf)


# scale env. data: ----

# load mean and sd with which training data were scaled:
load(file.path("data", "route_env_dt_scale_pars.RData")) # output of 2_1_DOM_flocker_fit_fm.R; env_scale_pars

clim_lu_cells_df_scaled <- clim_lu_cells_df

# iterate over columns:
for(c in colnames(clim_lu_cells_df[, 3:ncol(clim_lu_cells_df)])){
  print(c)
  clim_lu_cells_df_scaled[, c] <- as.numeric(scale(clim_lu_cells_df[, c], 
                                                   center = as.numeric(env_scale_pars$center[c]),
                                                  scale = as.numeric(env_scale_pars$scale[c])))
  }

# # check:
# round(apply(clim_lu_cells_df_scaled, 2, mean), 1)
# round(apply(clim_lu_cells_df_scaled, 2, sd), 1)


# # prepare US elevation contours to add to prediction maps: ----
# 
# # to clip to conterminous US:
# US_albers_sf <- read_sf(file.path("data", "US_outline_ESRI102003.shp"))
# 
# # elevation data, downloaded from https://www.sciencebase.gov/catalog/item/581d051de4b08da350d523c3
# US_contours <- read_sf(file.path("data", "USGS_contours_small_scale", "cont49l010a.shp")) %>% 
#   st_transform(crs = "ESRI:102003") %>% 
#   st_intersection(y = US_albers_sf) # clip to conterminous US
# 
# sort(unique(US_contours$Contour))
# table(US_contours$Contour)
# 
# contour_subset <- c(100, 200, 500, 700, 1000, 1500, 2100, 2500, 2900)
# 
# US_contours_subset <- US_contours %>% 
#   filter(Contour %in% contour_subset)
# US_contours_subset
# 
# rm(US_contours)
# #plot(st_geometry(US_contours_subset))
# 
# library(rmapshaper) # works better for simplifying these than sf::st_simplify
# 
# US_contours_subset_simpl_all <- vector(mode = "list", length = length(contour_subset))
# 
# for(i in 1:length(contour_subset)){
#   
#   print(contour_subset[i])
#   
#   US_contours_subset_simpl <- ms_simplify(US_contours_subset %>%  filter(Contour == contour_subset[i]), # fatal error if trying all at once
#                                           keep = 0.0005,
#                                           keep_shapes = FALSE) %>% 
#     filter(!st_is_empty(.))
#   
#   US_contours_subset_simpl_all[[i]] <- US_contours_subset_simpl
#   
#   print(US_contours_subset_simpl_all)
#   
# }
# 
# US_contours_subset_simpl_all_sf <- bind_rows(US_contours_subset_simpl_all)
# US_contours_subset_simpl <- US_contours_subset_simpl_all_sf
# 
# plot(US_contours_subset_simpl["Contour"])
# rm(US_contours_subset_simpl_all_sf)
# 
# st_write(US_contours_subset_simpl, 
#          file.path("data", "US_elev_contours_simplified.shp"), append = FALSE)

US_contours_subset_simpl <- st_read(file.path("data", "US_elev_contours_simplified.shp"))


# load model for species and predict: ----

# log overall progress:
prog_log_file <- file(file.path(log_dir, "fm_preds_US_750_progress.txt"), open = "wt") # write console output here
sink(prog_log_file, type = "message")
sink(prog_log_file, type = "output")

# species:
foreach(spec = final_species[43:length(final_species)],
        .packages = c("dplyr", "collapse", "flocker", "cmdstanr", "brms", "sf", "terra", "tidyterra", "ggplot2"),
        .errorhandling = "pass", #"remove",
        .verbose = TRUE) %dopar% {

#for(spec in final_species[1:length(final_species)]){ # xx}
  
#  print(spec)

  # # check where to look for model output (did MCMC fitting work with less or only with more iterations?)
  # if(file.exists(file.path(results_dir, "refit_2000_2000", paste0("out_", spec, "_CV_fold1.RData")))){
  #   output_dir <- file.path(results_dir, "refit_2000_2000")
  # } else {
  #   output_dir <- results_dir
  # }
  # 
  # print(output_dir)
  # 
  # # load fitted model:
  # 
  # skip_to_next <- FALSE
  # tryCatch(print(load(file.path(output_dir, paste0("out_", spec, "_fm_buffer750.RData")))), # output of 2_1_DOM_flocker_single_model.R
  #          error = function(e) { skip_to_next <<- TRUE})
  # if(skip_to_next) { next }
  # 
  # 
  # # buffers?!
  # 
  # # cells in buffer:
  spec_pres_buffer_sf <- training_routes(species = spec, buffer_km = 750, output = "buffer")

  cells_in_buffer_sf <- clim_lu_cells_sf %>%
    st_filter(., spec_pres_buffer_sf)
  # 
  # #plot(st_geometry(spec_pres_buffer_sf))
  # #plot(st_geometry(clim_lu_cells_sf), add = TRUE)
  # #plot(st_geometry(cells_in_buffer_sf), add = TRUE, col = "red")
  # 
  # # env. data within buffer:
  # 
  # clim_lu_cells_df_scaled_buff <- clim_lu_cells_df_scaled %>%
  #   filter(cellID %in% cells_in_buffer_sf$cellID)
  # 
  # # predict:
  # 
  # starttime <- Sys.time()
  # Z <- get_Z_mod(flocker_fit = out,
  #                draw_ids = seq(1, 4000, 4),
  #                new_data = clim_lu_cells_df_scaled_buff) # from 0_functions.R
  # endtime <- Sys.time()
  # endtime - starttime # 1.2 min
  # dim(Z)
  # 
  # # save predictions:
  # # round to 3 decimal numbers to save storage space:
  # Z <- round(Z, 3)
  # save(Z, file = file.path("results", "fm_preds_US", paste0(spec, "_US_occ_preds.RData")))
  # #save(Z, file = file.path(dir, "results", "fm_preds_US", paste0(spec, "_US_occ_preds.RData")))
  # 
  load(file.path(dir, "results", "fm_preds_US", paste0(spec, "_US_occ_preds.RData"))) # xx
  
  
  # compare Z with output of flocker::get_Z(): ----
  
  ## make predictions with flocker::get_Z():
  
  # # dummy data:
  # 
  # years <- unique(clim_lu_cells_df$year)
  # nyears <- length(years)
  # 
  # # dummy obs. as array sites x surveys x years:
  # nsites <- length(unique(clim_lu_cells_df_scaled_buff$cellID))
  # print(nsites)
  # 
  # nsurveys <- 5
  # 
  # y_array <- array(NA, dim = c(nsites, nsurveys, nyears))
  # for (t in 1:nyears){
  #   y_array[1:nsites, 1:nsurveys, t] <- as.matrix(0)
  # }
  # 
  # # environmental data:
  # env_cov <- vector("list", length = nyears)
  # for (t in 1:nyears){
  #   env_cov[[t]] <- clim_lu_cells_df_scaled_buff[which(clim_lu_cells_df_scaled_buff$year == years[t]), 
  #                                           c(selvar_final, paste0(selvar_final, "_3yrs"))]
  # }
  # 
  # # covariate for detection probability:
  # det_cov <- vector("list", length = 1)
  # names(det_cov) <- "route_section"
  # det_cov$route_section <- array(NA, dim = c(nsites, nsurveys, nyears))
  # det_cov$route_section[ , , 1:nyears] <- matrix(rep(c("Sect1", "Sect2", "Sect3", "Sect4", "Sect5"), nsites), nsites, byrow = TRUE)
  # 
  # starttime <- Sys.time()
  # 
  # # make flocker data:
  # fd_new <- make_flocker_data_dynamic(
  #   obs = y_array,
  #   unit_covs = env_cov, 
  #   event_covs = det_cov, 
  #   quiet = FALSE
  # )
  # save(fd_new, file = file.path("results", paste0(spec, "_US_flockerdt.RData")))
  # #load(file.path("results", paste0(spec, "_US_flockerdt.RData")))
  # rm(y_array, env_cov, det_cov) 
  # 
  # endtime <- Sys.time()
  # print(endtime-starttime) # 32 min
  # 
  # # predict occupancy probability:
  # starttime_pred <- Sys.time()
  # occ_posterior_US <- get_Z(out, 
  #                        new_data = fd_new, # requires flocker data object
  #                        history_condition = FALSE,
  #                        draw_ids = seq(1, 4000, 4)) # use 1000 draws
  # 
  # endtime_pred <- Sys.time()
  # print(endtime_pred-starttime_pred) # > 5 hours
  # 
  # save(occ_posterior_US, file = file.path(res_dir, "preds_US_fm", paste0(spec, "_US_occ_preds.RData")))
  
  
  ## compare predictions of flocker::get_Z() and get_Z_mod():
  
  # load(file.path("results", "fm_preds_US", paste0(spec, "_US_occ_preds.RData")))
  # 
  # # check whether flocker::get_Z() and get_Z_mod() yield same output: yes
  # dim(occ_posterior_US) # 2872, 25, 1000: US cells within buffer, years, draws
  # identical(Z, occ_posterior_US) # FALSE
  # all.equal(Z, occ_posterior_US) # class and attributes differ
  # occ_posterior_US_arr <- array(occ_posterior_US, dim = dim(Z))
  # all.equal(Z, occ_posterior_US_arr) # TRUE
  # identical(Z, occ_posterior_US_arr) # FALSE
  # range(Z - occ_posterior_US_arr) # very small numerical differences, otherwise the same
  # 
  
  
  
  
  
  
  # plots: ----
  
  
  # plot predicted occupancy probability across US in first and last year + uncertainty + difference:
  
  occ_start_stop_sf <- cells_in_buffer_sf %>% 
    filter(year == 1995) %>% # same cellIDs in each year, doesn't matter which
    select(cellID) %>% 
    mutate(occ_mean_1995 = apply(Z[,1,], MARGIN = 1, FUN = mean),
           occ_sd_1995 = apply(Z[,1,], MARGIN = 1, FUN = sd),
           occ_mean_2019 = apply(Z[,25,], MARGIN = 1, FUN = mean),
           occ_sd_2019 = apply(Z[,25,], MARGIN = 1, FUN = sd),
           occ_diff = occ_mean_2019 - occ_mean_1995)
  
  
  # plot(occ_start_stop_sf["occ_mean_1995"])
  # plot(occ_start_stop_sf["occ_sd_1995"])
  
  # convert to raster:
  occ_start_stop_rast <- occ_start_stop_sf %>% 
    mutate(x = st_coordinates(.)[,1],
           y = st_coordinates(.)[,2]) %>% 
    select(x, y, occ_mean_1995, occ_sd_1995, occ_mean_2019, occ_sd_2019, occ_diff) %>% 
    st_drop_geometry() %>% 
    rast(., type='xyz', crs=crs(clim_lu_cells_sf))
  
  # adjust elevation data:
  elev_dt <- st_intersection(US_contours_subset_simpl, spec_pres_buffer_sf) %>% 
    st_crop(y = occ_start_stop_rast) %>% 
    filter(Contour != 100)
  
  # observations: xx
  # species presences-absences:
  occ_dt_spec <- BBS_pres_abs_spec(species = spec)

  # match to spatial data:
  routes_obs_sf <- routes_sel_sf %>%
    left_join(occ_dt_spec, by = c(RTENO_BBS = "RTENO")) %>%
    select(-c(BCR, ObsN, doy, paste0("Count", seq(10, 50, 10)))) %>%
    filter(Year == 1995 | Year == 2019) %>% 
    mutate(lyr = paste0("occ_mean_", Year)) # needed to plot spatrasters and detections in same facets, same as name of spatraster layers

  
  # plot occ. prob. first and last year:
  occ_plot <- ggplot() +
    geom_spatraster(data = occ_start_stop_rast %>% 
                      select(-c(occ_sd_1995, occ_sd_2019, occ_diff))) +
    geom_sf(data = elev_dt, aes(colour = Contour), linewidth = 0.2) +
    geom_sf(data = routes_obs_sf %>% filter(presence == 1) %>% st_filter(., spec_pres_buffer_sf), color = "black", size = 0.5) +
    geom_sf(data = routes_obs_sf %>% filter(presence == 0) %>% st_filter(., spec_pres_buffer_sf), color = "grey90", size = 0.5) +
    facet_wrap(~lyr, ncol = 1) + # wrap by layer
    scale_fill_viridis_c(na.value = "transparent", option = "plasma", limits = c(0, 1)) +
    scale_colour_gradientn(colours = c(terrain.colors(8)[-8], "grey80"),
                           transform = "sqrt", guide = "none") +
    labs(fill = "occ. prob.") +
    theme_bw() +
    theme(plot.margin = unit(c(0.25, 0.25, 0.25, 0.25), "cm"))
  
  
  # plot standard dev. occ. prob. first and last year:
  occ_sd_plot <- ggplot() +
    geom_spatraster(data = occ_start_stop_rast %>% 
                      select(-c(occ_mean_1995, occ_mean_2019, occ_diff))) +
    geom_sf(data = elev_dt, aes(colour = Contour),linewidth = 0.2) +
    facet_wrap(~lyr, ncol = 1) +
    scale_fill_viridis_c(na.value = "transparent", option = "mako") +
    scale_colour_gradientn(colours = c(terrain.colors(8)[-8], "grey80"),
                           transform = "sqrt", guide = "none") +
    labs(fill = "sd. \nocc. prob.") +
    theme_bw() +
    theme(plot.margin = unit(c(0.25, 0.25, 0.25, 0.25), "cm"))
  

  # plot difference in occ. prob. between first and last year:
  limit <- max(abs(minmax(occ_start_stop_rast$occ_diff))) * c(-1, 1) # to center colour scale
  
  occ_diff <- ggplot() +
    geom_spatraster(data = occ_start_stop_rast %>% 
                      select(occ_diff)) +
    geom_sf(data = elev_dt, aes(colour = Contour), linewidth = 0.2) +
    scale_fill_distiller("diff. occ. prob. \n(2019 - 1995)", type = "div", 
                         palette = "RdBu", limit = limit, na.value = "transparent") +
    scale_colour_gradientn(colours = c(terrain.colors(8)[-9], "grey80"),
                           transform = "sqrt", 
                           breaks = c(100, 200, 500, 700, 1000, 1500, 2100, 2500, 2900)) +
    labs(fill = "diff. occ. prob. \n(2019 - 1995)", col = "elevation \n[m]") +
    guides(fill = guide_colourbar(order = 1),
           colour = guide_colourbar(order = 2)) +
    theme_dark() +
    theme(panel.grid.major = element_line(colour = "grey80"),
          panel.background = element_rect(colour = NA, fill = '#282C33'),
          legend.box = 'horizontal',
          plot.margin = unit(c(0.7, 0.25, 0.5, 0.25), "cm"))
  
  #plot_dir <- file.path(dir, "plots", "fm_preds_US") 
  plot_dir <- file.path("plots", "fm_preds_US") 
  
  if(!dir.exists(plot_dir)){dir.create(plot_dir)}
  
  jpeg(file = file.path(plot_dir, paste0("occ_preds_", spec, ".jpg")), 
        width = 1600, height = 1400, quality = 100, res = 170)

  # arrange plots:
  title <- ggdraw() + draw_label(spec, x = 0, hjust = 0) + theme(plot.margin = margin(7,0,7,3))
  
  second_row <- plot_grid(occ_plot, occ_sd_plot, labels = c("A", "B"), 
                         label_size = 12, nrow = 1, rel_widths = c(1, 1), axis = "l")
  
  print(plot_grid(title, second_row, occ_diff, 
                  label_size = 12, nrow = 3, labels = c("", "", "C"),
                  rel_heights = c(0.1, 2, 1)))

  dev.off()

}

sink(type="message")
sink(type="output")

stopCluster(cl)

rm(list=ls())
gc()

# # species observations (to add to maps): ----
# 
# # selected routes spatial data (to buffer presences):
# routes_sel_sf <- st_read(file.path("data", "route_selection_1995_2019_surv_beg_end_max_5y_miss_v2_spat_thin_100km_max_30_r_per_BCR_centroids.shp")) # output of 1_1_route_selection.R
# 
# # routes and years:
# load(file = file.path("data", "route_year_env_data.RData")) # merged route, year, environmental data
# 
# # species presences-absences:
# occ_dt_spec <- BBS_pres_abs_spec(species = spec)
# 
# # match to spatial data:
# routes_obs_sf <- routes_sel_sf %>% 
#   left_join(occ_dt_spec, by = c(RTENO_BBS = "RTENO")) %>% 
#   select(-c(BCR, ObsN, doy, paste0("Count", seq(10, 50, 10)))) %>% 
#   mutate(lyr = factor(Year)) # needed to plot spatrasters and detections in same facets, same as name of spatraster layers
# 



# ----

# # add observations:
# detection_colors <- c("1" = "black", "0" = "grey80")
# ggplot() +
#   geom_spatraster(data = test_rast, aes(fill = pred_mean_1995)) +
#   # observations:
#   geom_sf(data = routes_obs_sf[which(!is.na(routes_obs_sf$presence) & routes_obs_sf$Year == 1995),],
#           aes(colour = as.character(presence)),
#           size = 1, shape = 21) +
#   #facet_wrap(~lyr)+
#   scale_fill_viridis_c(na.value = "transparent") +
#   labs(fill = "mean occ. prob.") +
#   scale_color_manual(values = detection_colors, name = "detection") +
#   theme_bw() +
#   theme(panel.grid.major = element_line(colour = "transparent"),
#         text = element_text(size = 25)) +
#   ggtitle(spec)


# # fitted initial occupancy, colonisation and extinction probs. for new env. data:
# fitted_values_US <- fitted_flocker(out,
#                                    new_data = as.data.frame(clim_lu_cells_df_scaled), # tibble doesn't work
#                                    components = c("occ", "col", "ex")) # no det. because we don't have new data for this; quite fast
# 
# str(fitted_values_US)
# dim(fitted_values_US$linpred_occ) # 4000 draws for each cell-year combination
# 
# 
# # derive occupancy prob. for each year based on initial occ. prob. and colonisation
# # and extinction probabilities of each year:
# 
# # df with each US grid cell - year combination:
# post_df <- clim_lu_cells_df %>% 
#   select(cellID, year)
# 
# # occupancy prob. for each year and all cells:
# # save as array [cells, years, draws]
# pred_occ <- array(dim = c(length(unique(post_df$cellID)), nyears, ncol(fitted_values_US$linpred_occ)))
# 
# for(t in 1:nyears){
#   
#   print(years[t])
#   
#   # initial occupancy for first year:
#   if(t == 1){
#     
#     pred_occ[, 1, ] <- fitted_values_US$linpred_occ[which(post_df$year == years[t]), ] # posterior draws initial occupancy
#     
#   } else {
#     
#     # following years:
#     
#     # use colonisation and extinction prob. of previous year to calculate this years occupancy prob.:
#     year_ind <- which(post_df$year == years[t-1])
#     eps <- fitted_values_US$linpred_ex[year_ind, ]
#     col <- fitted_values_US$linpred_col[year_ind, ]
#     
#     pred_occ[, t, ] <- pred_occ[, t-1, ] * (1 - eps) + (1 - pred_occ[, t-1, ]) * col # 4000 draws for occupancy at all sites in year t
#   }
# }
# 
# 
# # maps of predicted values: ----------------------------------------------------
# 
# # add coordinates as columns, needed for conversion to raster:
# clim_lu_cells_sf <- clim_lu_cells_sf %>% 
#   mutate(x = st_coordinates(.)[,1]) %>% 
#   mutate(y = st_coordinates(.)[,2])
# 
# 
# ## mean colonisation prob.: ----
# 
# # mean col. prob.:
# clim_lu_cells_sf$mean_col <- apply(fitted_values_US$linpred_col, FUN = mean, MAR = 1) # mean for each cell
# 
# # convert to raster:
# col_preds_rast <- clim_lu_cells_sf %>% 
#   select(x, y, year, mean_col) %>% 
#   st_drop_geometry() %>% 
#   tidyr::pivot_wider(names_from = year, values_from = mean_col) %>% 
#   rast(., type='xyz', crs=crs(clim_lu_cells_sf))
# 
# # map predicted mean col. prob. per year:
# p_col <- ggplot() +
#   geom_spatraster(data = col_preds_rast) +
#   facet_wrap(~lyr) + # wrap by layer
#   #scale_fill_viridis_c(na.value = "transparent", transform = "log") +
#   scale_fill_viridis_c(na.value = "transparent", transform = "sqrt") +
#   #scale_fill_viridis_c(na.value = "transparent") +
#   # scale_fill_whitebox_c(palette = "muted",na.value = "white") +
#   labs(fill = "mean col. prob.") +
#   theme_bw() +
#   theme(panel.grid.major = element_line(colour = "transparent")) +
#   ggtitle(spec)
# 
# dir.create(file.path("plots", "DOM_preds_US"))
# jpeg(file = file.path("plots", "DOM_preds_US", paste0("col_", spec, ".jpg")), 
#      width = 1200, height = 900, quality = 100)
# p_col
# dev.off()
# 
# 
# ## mean extinction prob.: ----
# 
# # mean ex. prob.:
# clim_lu_cells_sf$mean_ex <- apply(fitted_values_US$linpred_ex, FUN = mean, MAR = 1) # mean for each cell
# 
# # convert to raster:
# ex_preds_rast <- clim_lu_cells_sf %>% 
#   select(x, y, year, mean_ex) %>% 
#   st_drop_geometry() %>% 
#   tidyr::pivot_wider(names_from = year, values_from = mean_ex) %>% 
#   rast(., type='xyz', crs=crs(clim_lu_cells_sf))
# 
# # map predicted mean col. prob. per year:
# p_ex <- ggplot() +
#   geom_spatraster(data = ex_preds_rast) +
#   facet_wrap(~lyr) + 
#   scale_fill_viridis_c(na.value = "transparent") +
#   labs(fill = "mean ex. prob.") +
#   theme_bw() +
#   theme(panel.grid.major = element_line(colour = "transparent")) +
#   ggtitle(spec)
# 
# jpeg(file = file.path("plots", "DOM_preds_US", paste0("ex_", spec, ".jpg")), 
#      width = 1200, height = 900, quality = 100)
# p_ex
# dev.off()
# 
# 
# ## occupancy probability: ----
# 
# # mean occ. prob.:
# mean_occ <- numeric(length = nrow(clim_lu_cells_sf))
# 
# for(t in 1:nyears){
#   print(years[t])
#   year_ind <- which(post_df$year == years[t])
#   mean_occ[year_ind] <- apply(pred_occ[, t, ], FUN = mean, MAR = 1)
# }
# 
# clim_lu_cells_sf$mean_occ <- mean_occ
# rm(mean_occ)
# 
# # convert to raster:
# occ_preds_rast <- clim_lu_cells_sf %>% 
#   select(x, y, year, mean_occ) %>% 
#   st_drop_geometry() %>% 
#   tidyr::pivot_wider(names_from = year, values_from = mean_occ) %>% 
#   rast(., type='xyz', crs=crs(clim_lu_cells_sf))
# 
# # map predicted mean occupancy per year:
# 
# p_occ_mean <- ggplot() +
#   geom_spatraster(data = occ_preds_rast) +
#   facet_wrap(~lyr) +
#   scale_fill_viridis_c(na.value = "transparent") +
#   labs(fill = "mean occ. prob.") +
#   theme_bw() +
#   theme(panel.grid.major = element_line(colour = "transparent")) +
#   ggtitle(spec)
# 
# jpeg(file = file.path("plots", "DOM_preds_US", paste0("occ_", spec, ".jpg")), 
#      width = 1200, height = 900, quality = 100)
# p_occ_mean
# dev.off()
# 
# # add observations:
# detection_colors <- c("1" = "black", "0" = "grey80")
# p_occ_mean_obs <- ggplot() +
#   geom_spatraster(data = occ_preds_rast) +
#   # observations:
#   geom_sf(data = routes_obs_sf[which(!is.na(routes_obs_sf$occ_route)),],
#           aes(colour = as.character(occ_route)),
#           size = 1, shape = 21) + 
#   facet_wrap(~lyr)+
#   scale_fill_viridis_c(na.value = "transparent") +
#   labs(fill = "mean occ. prob.") +
#   scale_color_manual(values = detection_colors, name = "detection") +
#   theme_bw() +
#   theme(panel.grid.major = element_line(colour = "transparent"),
#         text = element_text(size = 25)) +
#   ggtitle(spec)
# 
# jpeg(file = file.path("plots", "DOM_preds_US", paste0("obs_occ_", spec, ".jpg")), 
#      width = 1800, height = 1200, quality = 100)
# p_occ_mean_obs
# dev.off()
# 
# 
# # uncertainty: sd of estimated occ. prob.:
# sd_occ <- numeric(length = nrow(clim_lu_cells_sf))
# 
# for(t in 1:nyears){
#   print(years[t])
#   year_ind <- which(post_df$year == years[t])
#   sd_occ[year_ind] <- apply(pred_occ[, t, ], FUN = sd, MAR = 1)
# }
# 
# clim_lu_cells_sf$sd_occ <- sd_occ
# rm(sd_occ)
# 
# # convert to raster:
# occ_sd_preds_rast <- clim_lu_cells_sf %>% 
#   select(x, y, year, sd_occ) %>% 
#   st_drop_geometry() %>% 
#   tidyr::pivot_wider(names_from = year, values_from = sd_occ) %>% 
#   rast(., type='xyz', crs=crs(clim_lu_cells_sf))
# 
# # map sd predicted occupancy per year:
# 
# p_occ_sd_obs <- ggplot() +
#   geom_spatraster(data = occ_sd_preds_rast) +
#   # observations:
#   geom_sf(data = routes_obs_sf[which(!is.na(routes_obs_sf$occ_route)),],
#           aes(colour = as.character(occ_route)),
#           size = 1, shape = 21) + 
#   facet_wrap(~lyr)+
#   scale_fill_whitebox_c(palette = "muted", na.value = "white") +  labs(fill = "mean occ. prob.") +
#   labs(fill = "sd occ. prob.") +
#   scale_color_manual(values = detection_colors, name = "detection") +
#   theme_bw() +
#   theme(panel.grid.major = element_line(colour = "transparent"),
#         text = element_text(size = 25)) +
#   ggtitle(spec)
# 
# jpeg(file = file.path("plots", "DOM_preds_US", paste0("obs_sd_occ_", spec, ".jpg")), 
#      width = 1800, height = 1200, quality = 100)
# p_occ_sd_obs
# dev.off()

