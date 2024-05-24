# dynamic occupancy model predictions of occupancy, colonisation and extinction probability
# for whole conterminous US (maps):


# packages: --------------------------------------------------------------------

library(dplyr)
library(flocker)
#install.packages("cmdstanr", repos = c("https://mc-stan.org/r-packages/", getOption("repos")))
library(cmdstanr)
#set_cmdstan_path(path = NULL)
set_cmdstan_path("C:/Users/schifferle1/Documents/cmdstan-2.34.1")
library(sf)
library(terra)
library(tidyterra)
library(ggplot2)


# load data: -------------------------------------------------------------------

# species:

spec <- "Eurasian Collared-Dove"

# load fitted model:

load(file.path("results", paste0("out_flocker_", spec, ".RData"))) # output of 2_1_DOM_flocker_single_model.R
out

# environmental data for each US grid cell:

load(file = file.path("data", "US_grid_env_data.RData")) # output of 1_3_env_data_df_contUS.R
clim_lu_cells_df <- sf::st_drop_geometry(clim_lu_cells_sf)
## scale:
clim_lu_cells_df_scaled <- clim_lu_cells_df %>% # use same mean and sd for scaling as for training data? xx
  mutate(across(bio2:pr_winter_3yrs, ~ (scale(.)) %>% as.vector()))

years <- unique(clim_lu_cells_df$year)
nyears <- length(years)

# species observations (to add to maps):

# routes:
routes_sel_sf <- sf::st_read(file.path("data", "route_selection_1991_2015_surv_beg_end_max_5y_miss_v2_spat_thin_100km_max_30_r_per_BCR_centroids.shp")) # output of 1_1_route_selection.R
# routes and years:
load(file = file.path("data", "route_year_env_data.RData")) # merged route, year, environmental data
# route-year-species information (only surveyed)
load(file = file.path("data", "BBS_for_occ_spec_records.RData")) # output of 1_0_reformat_BBS_data.R

presences_spec <- bbs_dt_occ %>% 
  select(c(English_Common_Name, RTENO, Year, paste0("Count", seq(10, 50, 10)))) %>% 
  filter(English_Common_Name == spec)

# match observations to routes-year-env:
occ_dt_spec <- route_sel_env_dt_final %>% 
  # add observations:
  collapse::join(presences_spec, on = c("RTENO", "Year"), how = "left") %>% 
  # if route was surveyed but species not observed, replace NA with 0:
  mutate(across(Count10:Count50, ~ 
                  case_when(Surveyed == 1 & is.na(.) ~ 0,
                            .default = .))) %>%
  # convert bird counts to presence / absence:
  mutate(across(Count10:Count50, ~ 
                  case_when(. > 1 ~ 1,
                            .default = .))) %>% 
  # sum presence/absence on route:
  mutate(occ_route = rowSums(across(Count10:Count50))) %>%
  # convert presence absence:
  mutate(occ_route = ifelse(occ_route > 0, 1, 0)) %>% 
  select(RTENO, Year, occ_route)

# match to spatial data:
routes_obs_sf <- routes_sel_sf %>% 
  left_join(occ_dt_spec, by = c(RTENO_BBS = "RTENO")) %>% 
  mutate(lyr = factor(Year)) # needed to plot spatrasters and detections in same facets, same as name of spatraster layers


# model predictions: -----------------------------------------------------------

# fitted initial occupancy, colonisation and extinction probs. for new env. data:
fitted_values_US <- fitted_flocker(out,
                                   new_data = as.data.frame(clim_lu_cells_df_scaled), # tibble doesn't work
                                   components = c("occ", "col", "ex")) # no det. because we don't have new data for this; quite fast

str(fitted_values_US)
dim(fitted_values_US$linpred_occ) # 4000 draws for each cell-year combination


# derive occupancy prob. for each year based on initial occ. prob. and colonisation
# and extinction probabilities of each year:

# df with each US grid cell - year combination:
post_df <- clim_lu_cells_df %>% 
  select(cellID, year)

# occupancy prob. for each year and all cells:
# save as array [cells, years, draws]
pred_occ <- array(dim = c(length(unique(post_df$cellID)), nyears, ncol(fitted_values_US$linpred_occ)))

for(t in 1:nyears){
  
  print(years[t])
  
  # initial occupancy for first year:
  if(t == 1){
    
    pred_occ[, 1, ] <- fitted_values_US$linpred_occ[which(post_df$year == years[t]), ] # posterior draws initial occupancy
    
  } else {
    
    # following years:
    
    # use colonisation and extinction prob. of previous year to calculate this years occupancy prob.:
    year_ind <- which(post_df$year == years[t-1])
    eps <- fitted_values_US$linpred_ex[year_ind, ]
    col <- fitted_values_US$linpred_col[year_ind, ]
    
    pred_occ[, t, ] <- pred_occ[, t-1, ] * (1 - eps) + (1 - pred_occ[, t-1, ]) * col # 4000 draws for occupancy at all sites in year t
  }
}


# maps of predicted values: ----------------------------------------------------

# add coordinates as columns, needed for conversion to raster:
clim_lu_cells_sf <- clim_lu_cells_sf %>% 
  mutate(x = st_coordinates(.)[,1]) %>% 
  mutate(y = st_coordinates(.)[,2])


## mean colonisation prob.: ----

# mean col. prob.:
clim_lu_cells_sf$mean_col <- apply(fitted_values_US$linpred_col, FUN = mean, MAR = 1) # mean for each cell

# convert to raster:
col_preds_rast <- clim_lu_cells_sf %>% 
  select(x, y, year, mean_col) %>% 
  st_drop_geometry() %>% 
  tidyr::pivot_wider(names_from = year, values_from = mean_col) %>% 
  rast(., type='xyz', crs=crs(clim_lu_cells_sf))

# map predicted mean col. prob. per year:
p_col <- ggplot() +
  geom_spatraster(data = col_preds_rast) +
  facet_wrap(~lyr) + # wrap by layer
  #scale_fill_viridis_c(na.value = "transparent", transform = "log") +
  scale_fill_viridis_c(na.value = "transparent", transform = "sqrt") +
  #scale_fill_viridis_c(na.value = "transparent") +
  # scale_fill_whitebox_c(palette = "muted",na.value = "white") +
  labs(fill = "mean col. prob.") +
  theme_bw() +
  theme(panel.grid.major = element_line(colour = "transparent")) +
  ggtitle(spec)

dir.create(file.path("plots", "DOM_preds_US"))
jpeg(file = file.path("plots", "DOM_preds_US", paste0("col_", spec, ".jpg")), 
     width = 1200, height = 900, quality = 100)
p_col
dev.off()


## mean extinction prob.: ----

# mean ex. prob.:
clim_lu_cells_sf$mean_ex <- apply(fitted_values_US$linpred_ex, FUN = mean, MAR = 1) # mean for each cell

# convert to raster:
ex_preds_rast <- clim_lu_cells_sf %>% 
  select(x, y, year, mean_ex) %>% 
  st_drop_geometry() %>% 
  tidyr::pivot_wider(names_from = year, values_from = mean_ex) %>% 
  rast(., type='xyz', crs=crs(clim_lu_cells_sf))

# map predicted mean col. prob. per year:
p_ex <- ggplot() +
  geom_spatraster(data = ex_preds_rast) +
  facet_wrap(~lyr) + 
  scale_fill_viridis_c(na.value = "transparent") +
  labs(fill = "mean ex. prob.") +
  theme_bw() +
  theme(panel.grid.major = element_line(colour = "transparent")) +
  ggtitle(spec)

jpeg(file = file.path("plots", "DOM_preds_US", paste0("ex_", spec, ".jpg")), 
     width = 1200, height = 900, quality = 100)
p_ex
dev.off()


## occupancy probability: ----

# mean occ. prob.:
mean_occ <- numeric(length = nrow(clim_lu_cells_sf))

for(t in 1:nyears){
  print(years[t])
  year_ind <- which(post_df$year == years[t])
  mean_occ[year_ind] <- apply(pred_occ[, t, ], FUN = mean, MAR = 1)
}

clim_lu_cells_sf$mean_occ <- mean_occ
rm(mean_occ)

# convert to raster:
occ_preds_rast <- clim_lu_cells_sf %>% 
  select(x, y, year, mean_occ) %>% 
  st_drop_geometry() %>% 
  tidyr::pivot_wider(names_from = year, values_from = mean_occ) %>% 
  rast(., type='xyz', crs=crs(clim_lu_cells_sf))

# map predicted mean occupancy per year:

p_occ_mean <- ggplot() +
  geom_spatraster(data = occ_preds_rast) +
  facet_wrap(~lyr) +
  scale_fill_viridis_c(na.value = "transparent") +
  labs(fill = "mean occ. prob.") +
  theme_bw() +
  theme(panel.grid.major = element_line(colour = "transparent")) +
  ggtitle(spec)

jpeg(file = file.path("plots", "DOM_preds_US", paste0("occ_", spec, ".jpg")), 
     width = 1200, height = 900, quality = 100)
p_occ_mean
dev.off()

# add observations:
detection_colors <- c("1" = "black", "0" = "grey80")
p_occ_mean_obs <- ggplot() +
  geom_spatraster(data = occ_preds_rast) +
  # observations:
  geom_sf(data = routes_obs_sf[which(!is.na(routes_obs_sf$occ_route)),],
          aes(colour = as.character(occ_route)),
          size = 1, shape = 21) + 
  facet_wrap(~lyr)+
  scale_fill_viridis_c(na.value = "transparent") +
  labs(fill = "mean occ. prob.") +
  scale_color_manual(values = detection_colors, name = "detection") +
  theme_bw() +
  theme(panel.grid.major = element_line(colour = "transparent"),
        text = element_text(size = 25)) +
  ggtitle(spec)

jpeg(file = file.path("plots", "DOM_preds_US", paste0("obs_occ_", spec, ".jpg")), 
     width = 1800, height = 1200, quality = 100)
p_occ_mean_obs
dev.off()


# uncertainty: sd of estimated occ. prob.:
sd_occ <- numeric(length = nrow(clim_lu_cells_sf))

for(t in 1:nyears){
  print(years[t])
  year_ind <- which(post_df$year == years[t])
  sd_occ[year_ind] <- apply(pred_occ[, t, ], FUN = sd, MAR = 1)
}

clim_lu_cells_sf$sd_occ <- sd_occ
rm(sd_occ)

# convert to raster:
occ_sd_preds_rast <- clim_lu_cells_sf %>% 
  select(x, y, year, sd_occ) %>% 
  st_drop_geometry() %>% 
  tidyr::pivot_wider(names_from = year, values_from = sd_occ) %>% 
  rast(., type='xyz', crs=crs(clim_lu_cells_sf))

# map sd predicted occupancy per year:

p_occ_sd_obs <- ggplot() +
  geom_spatraster(data = occ_sd_preds_rast) +
  # observations:
  geom_sf(data = routes_obs_sf[which(!is.na(routes_obs_sf$occ_route)),],
          aes(colour = as.character(occ_route)),
          size = 1, shape = 21) + 
  facet_wrap(~lyr)+
  scale_fill_whitebox_c(palette = "muted", na.value = "white") +  labs(fill = "mean occ. prob.") +
  labs(fill = "sd occ. prob.") +
  scale_color_manual(values = detection_colors, name = "detection") +
  theme_bw() +
  theme(panel.grid.major = element_line(colour = "transparent"),
        text = element_text(size = 25)) +
  ggtitle(spec)

jpeg(file = file.path("plots", "DOM_preds_US", paste0("obs_sd_occ_", spec, ".jpg")), 
     width = 1800, height = 1200, quality = 100)
p_occ_sd_obs
dev.off()