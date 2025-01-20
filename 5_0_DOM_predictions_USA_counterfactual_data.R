# use fitted DOMs to predict occupancy probability across conterminous US 
# attribution: 
# use counterfactual data for one environmental variable at a time
# predict to these data with the fitted DOMs
# compare predictions to predictions using factual data only


# packages: --------------------------------------------------------------------

library(dplyr)
library(flocker)
library(cmdstanr)
#set_cmdstan_path(path = NULL)
set_cmdstan_path("C:/Users/schifferle1/Documents/cmdstan-2.34.1")
library(sf)
library(terra)
library(tidyterra)
library(ggplot2)


# functions: -----

source("0_functions.R")


# directories: ----

# res_dir <- file.path("M:", "Documents", "DEBTs", "analysis", "Schifferle_BBS_occupancy_models_2023",
#                            "results")
res_dir <- file.path("results")


# load data: -------------------------------------------------------------------

# routes-years:
load(file = file.path("data", "BBS_for_occ_selection.RData")) # route_sel_dt; output of 1_3_match_BBS_to_env_data.R 

# route-year-species information (only surveyed)
load(file = file.path("data", "BBS_for_occ_spec_records.RData")) # bbs_dt_occ; output of 1_0_reformat_BBS_data.R

# selected variables:
load(file = file.path("data", "selected_variables.RData")) # selvar_final; output of 1_2_variable_selection.R
selvar_final

# selected routes spatial data (to buffer presences):
routes_sel_sf <- st_read(file.path("data", "route_selection_1995_2019_surv_beg_end_max_5y_miss_v2_spat_thin_100km_max_30_r_per_BCR_centroids.shp")) # output of 1_1_route_selection.R

# selected species:
# sorted by ecoregion:
load(file = file.path("data", "final_species_selection_eco_sorted.RData")) # final_species_eco_sorted; output of 1_2_species_selection.R
final_species_eco_sorted

# factual environmental data for each US grid cell:

load(file = file.path("data", "US_grid_env_data.RData")) # output of 1_3_env_data_df_contUS.R; clim_lu_cells_sf
fact_env_df <- sf::st_drop_geometry(clim_lu_cells_sf)

# counterfactual environmental data for each US grid cell:
load(file = file.path("data", "US_grid_env_data_cf.RData")) # output of 1_3_env_data_df_contUS.R; clim_lu_cells_sf
counterfact_env_df <- sf::st_drop_geometry(clim_lu_cells_sf)


# scale env. data: ----

# load mean and sd with which training data were scaled:
load(file.path("data", "route_env_dt_scale_pars.RData")) # output of 2_1_DOM_flocker_fit_fm.R; env_scale_pars

# factual data:
fact_env_df_scaled <- fact_env_df
# iterate over columns:
for(c in colnames(fact_env_df[, 3:ncol(fact_env_df)])){
  print(c)
  fact_env_df_scaled[, c] <- as.numeric(scale(fact_env_df[, c], 
                                                   center = as.numeric(env_scale_pars$center[c]),
                                                   scale = as.numeric(env_scale_pars$scale[c])))
}

# counterfactual data:
counterfact_env_df_scaled <- counterfact_env_df
# iterate over columns:
for(c in colnames(counterfact_env_df[, 3:ncol(counterfact_env_df)])){
  print(c)
  counterfact_env_df_scaled[, c] <- as.numeric(scale(counterfact_env_df[, c], 
                                              center = as.numeric(env_scale_pars$center[c]),
                                              scale = as.numeric(env_scale_pars$scale[c])))
}

# # check:
# round(apply(fact_env_df_scaled, 2, mean), 1)
# round(apply(fact_env_df_scaled, 2, sd), 1)
# round(apply(counterfact_env_df_scaled, 2, mean), 1)
# round(apply(counterfact_env_df_scaled, 2, sd), 1)


# load model for species and predict: ----

# directory to store predictions:
preds_dir <- file.path("results", "fm_preds_US", "fm_preds_US_counterfactual")
if(!dir.exists(preds_dir)){
  dir.create(preds_dir)
}

# species:
spec <- "Carolina Wren"

for(spec in final_species_eco_sorted){ # continue with Ovenbird, 43
  
  print(spec)
  
  # load fitted model:
  load(file.path(res_dir, "fm_buffer750km", paste0("out_", spec, "_fm_buffer750.RData"))) # output of 2_1_DOM_flocker_single_model.R
  
  # buffers?!
  
  # cells in buffer:
  spec_pres_buffer_sf <- training_routes(species = spec, buffer_km = 750, output = "buffer")
  
  cells_in_buffer_sf <- clim_lu_cells_sf %>% 
    st_filter(., spec_pres_buffer_sf)
  
  #plot(st_geometry(spec_pres_buffer_sf))
  #plot(st_geometry(clim_lu_cells_sf), add = TRUE)
  #plot(st_geometry(cells_in_buffer_sf), add = TRUE, col = "red")
  
  # env. data within buffer: 
  
  fact_env_df_scaled_buff <- fact_env_df_scaled %>% 
    filter(cellID %in% cells_in_buffer_sf$cellID)
  counterfact_env_df_scaled_buff <- counterfact_env_df_scaled %>% 
    filter(cellID %in% cells_in_buffer_sf$cellID)
  
  # iterate over variables, predict to counterfactual data for focal variable and factual data for all others:
  for(v in selvar_final){
    
    print(v)
    
    # gather new data:
    new_data <- fact_env_df_scaled_buff %>% 
      select(-all_of(c(v, paste0(v, "_3yrs")))) %>% # xx counterfcat bio1!
      left_join(counterfact_env_df_scaled_buff[, c("cellID", "year", v, paste0(v, "_3yrs"))])
      
    # predict:
    Z_cf <- get_Z_mod(flocker_fit = out, 
                   draw_ids = seq(1, 4000, 4),
                   new_data = new_data) # from 0_functions.R
    
    # save predictions:
    # round to 3 decimal numbers to save storage space:
    Z_cf <- round(Z_cf, 3)
    #Z_cf <- as.integer(round(Z_cf * 100))
    save(Z_cf, file = file.path(preds_dir,  paste0(spec, "_occ_preds_cf_", v, ".RData")))

  }

  
  # sum difference in occupancy probability in 2019 between factual and counterfactual data
  # or difference in difference between 1995 and 2019?
  # load predictions to counterfactual data:
  var_imp_test1 <- rep(NA, length(selvar_final))
  for(i in 1:length(selvar_final)){
    print(i)
    load(file.path(preds_dir, paste0(spec, "_occ_preds_cf_", selvar_final[i], ".RData")))
    var_imp_test1[i] <- sum(apply(Zf[,25,], MARGIN = 1, FUN = mean) - apply(Z_cf[,25,], MARGIN = 1, FUN = mean))
  }

  var_imp_test1
  
  var_imp_test1_scal <- round(var_imp_test1 / sum(abs(var_imp_test1)) * 100)
  names(var_imp_test1_scal) <- selvar_final
  sort(abs(var_imp_test1_scal), decreasing = TRUE)
  # reflects that there is much more change in land use than in climate?
  # account for uncertainty! xx
  
  # or difference in difference between 1995 and 2019?:
  var_imp_test2 <- rep(NA, length(selvar_final))
  for(i in 1:length(selvar_final)){
    print(i)
    load(file.path(preds_dir, paste0(spec, "_occ_preds_cf_", selvar_final[i], ".RData")))
    
    diff_f <- apply(Zf[,25,], MARGIN = 1, FUN = mean) - apply(Zf[,1,], MARGIN = 1, FUN = mean)
    diff_cf <- apply(Z_cf[,25,], MARGIN = 1, FUN = mean) - apply(Z_cf[,1,], MARGIN = 1, FUN = mean)
    
    var_imp_test2[i] <- sum(diff_f - diff_cf)
  }
  
  var_imp_test2_scal <- round(var_imp_test2 / sum(abs(var_imp_test2)) * 100)
  names(var_imp_test2_scal) <- selvar_final
  sort(abs(var_imp_test2_scal), decreasing = TRUE)
  # looks quite different
  
  # plots: ----
  
  # load predictions to factual data:
  load(file.path("results", "fm_preds_US", paste0(spec, "_US_occ_preds.RData")))
  Zf <- round(Z, 3)
  
  # load predictions to counterfactual data:
  v <- "urbanareas"
  load(file.path(preds_dir, paste0(spec, "_occ_preds_cf_", v, ".RData")))
  


  # plot predicted occupancy probability across US in first and last year + uncertainty + difference:
  
  occ_start_stop_sf <- cells_in_buffer_sf %>% 
    filter(year == 1995) %>% # same cellIDs in each year, doesn't matter which
    select(cellID) %>% 
    mutate(occ_mean_1995_f = apply(Zf[,1,], MARGIN = 1, FUN = mean),
           occ_sd_1995_f = apply(Zf[,1,], MARGIN = 1, FUN = sd),
           occ_mean_2019_f = apply(Zf[,25,], MARGIN = 1, FUN = mean),
           occ_sd_2019_f = apply(Zf[,25,], MARGIN = 1, FUN = sd),
           occ_diff_f = occ_mean_1995_f - occ_mean_2019_f,
           occ_mean_1995_cf = apply(Z_cf[,1,], MARGIN = 1, FUN = mean),
           occ_sd_1995_cf = apply(Z_cf[,1,], MARGIN = 1, FUN = sd),
           occ_mean_2019_cf = apply(Z_cf[,25,], MARGIN = 1, FUN = mean),
           occ_sd_2019_cf = apply(Z_cf[,25,], MARGIN = 1, FUN = sd),
           occ_diff_cf = occ_mean_1995_cf - occ_mean_2019_cf,
           occ_diff_f_cf_1995 = occ_mean_1995_f - occ_mean_1995_cf,
           occ_diff_f_cf_2019 = occ_mean_2019_f - occ_mean_2019_cf)
  

  
  # plot(occ_start_stop_sf["occ_mean_1995"])
  # plot(occ_start_stop_sf["occ_sd_1995"])
  # plot(occ_start_stop_sf["occ_diff_f_cf_2019"])
  
  occ_start_stop_rast <- occ_start_stop_sf %>% 
    mutate(x = st_coordinates(.)[,1],
           y = st_coordinates(.)[,2]) %>% 
    select(x, y, 
           occ_mean_1995_f, occ_sd_1995_f, occ_mean_2019_f, occ_sd_2019_f, occ_diff_f, 
           occ_mean_1995_cf, occ_sd_1995_cf, occ_mean_2019_cf, occ_sd_2019_cf, occ_diff_cf, 
           occ_diff_f_cf_1995, occ_diff_f_cf_2019) %>%
    st_drop_geometry() %>% 
    rast(., type = 'xyz', crs = crs(clim_lu_cells_sf))
  
  # ggplot() +
  #   geom_spatraster(data = occ_start_stop_rast, aes(fill = occ_diff_f_cf_2019)) +
  #   #facet_wrap(~lyr, ncol = 2) + # wrap by layer
  #   #scale_fill_viridis_c(na.value = "transparent", transform = "log") +
  #   #scale_fill_viridis_c(na.value = "transparent", transform = "sqrt") +
  #   scale_fill_viridis_c(na.value = "transparent") +
  #   #scale_fill_whitebox_c(palette = "muted",na.value = "white") +
  #   labs(fill = "occ. prob.") +
  #   theme_bw() +
  #   theme(#panel.grid.major = element_line(colour = "transparent"),
  #     plot.margin = unit(c(0.1, 0.1, 0.1, 0.2), "null")) +
  #   ggtitle(spec)

  
  # difference in occ. prob. between first and last year:
  limit <- max(abs(minmax(occ_start_stop_rast$occ_diff_f_cf_2019))) * c(-1, 1) # to center colour scale
  
  ggplot() +
    geom_spatraster(data = occ_start_stop_rast, aes(fill = occ_diff_f_cf_2019)) +
    #scale_fill_viridis_c(na.value = "transparent", transform = "sqrt") +
    #scale_fill_whitebox_c(palette = "muted", na.value = "transparent") +
    scale_fill_distiller(paste0("diff. occ. prob. 2019 \n(fact. - cf ", v, ")"), type = "div", 
                         palette = "RdBu", limit = limit, na.value = "transparent") +
    labs(fill = "occ. prob.") +
    #theme_bw() +
    theme_dark() +
    # ggthemes::theme_solarized(light=FALSE) +
    theme(panel.grid.major = element_line(colour = "grey80"),
          panel.background = element_rect(colour = NA, fill = '#282C33'),
          #       axis.text = element_text(colour = 'grey70'),
          #       legend.text = element_text(color = 'white'),
          #       legend.title = element_text(colour = 'white'),
          #       plot.title = element_text(colour = 'white'),
          plot.margin = unit(c(0.1, 0.1, 0.1, 0.2), "null")
    ) +
  ggtitle(spec)
  
  
  
  plot_dir <- file.path("plots", "fm_preds_US")
  if(!dir.exists(plot_dir)){
    dir.create(plot_dir)
  }
  jpeg(file = file.path(plot_dir, paste0("occ_preds_", spec, ".jpg")), 
       width = 1200, height = 1200, quality = 100, res = 200)
  print(cowplot::plot_grid(occ_plot, occ_diff, labels = c('A', 'B'), label_size = 12,
                           nrow = 2, rel_heights = c(1, 0.5), axis = "l"))
  dev.off()
  
}
