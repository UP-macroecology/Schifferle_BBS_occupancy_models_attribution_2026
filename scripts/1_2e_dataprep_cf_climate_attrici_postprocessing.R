# Script:   1_2d_dataprep_cf_climate_attrici_postprocessing.R
# Purpose:  Generate counterfactual variables from ATTRICI output to simulate counterfactual occupancy dynamics with dynamic occupancy models
# Inputs:   data/Counterfactual_env_data/ISIMIP_GSWP3_W5E5/attrici_detrending/output/<var>_detrended/US_<var>_detrended_1994_2019.nc
#           data/US_outline_ESRI102003.shp
#           data/selected_variables.RData
#           data/route_selection_1995_2019_surv_beg_end_max_5y_miss_v2_spat_thin_100km_max_30_r_per_BCR_centroids.shp
# Outputs:  data/Counterfactual_env_data/ISIMIP_GSWP3_W5E5/attrici_detrending/output/ATTRICI_CLIM_ESRI102003_tifs/<var>_<yyyymm>_ESRI102003.tif
#           data/Counterfactual_env_data/ISIMIP_GSWP3_W5E5/attrici_detrending/output/ATTRICI_CLIM_ESRI102003_tifs/bioclim/<var>_<year>.tif
#           data/Counterfactual_env_data/ISIMIP_GSWP3_W5E5/attrici_detrending/output/ATTRICI_CLIM_ESRI102003_tifs/seasonal/<var>_<year>.tif
#           plots/attrici/gmt_raw_smoothed.svg
#           plots/attrici/logp.svg
#           plots/attrici/<var>_sel_routes_monthly.svg
#           plots/attrici/<var>_sel_routes_daily.svg
# Runs on:  Local

# Steps:
# 1) postprocess temperature: detrended tasskew and tasrange to tasmin and tasmax
# 2) convert netCDF output to tifs
# 3) calculate counterfactual version of the climate variables used to fit dynamic occupancy models
# 4) plots to check detrending

source(file.path("scripts", "0_paths.R"))


# packages: --------------------------------------------------------------------

library(ncdf4)
library(dplyr)
library(sf)
library(terra)
library(ggplot2)
library(tidyterra)
library(doParallel)
library(patchwork)
library(cowplot)
library(lubridate)


# directories: -----------------------------------------------------------------

# directory where ATTRICI stores output:
attrici_out <- file.path(dir, "data", "Counterfactual_env_data", 
                         "ISIMIP_GSWP3_W5E5", "attrici_detrending", "output")

# ATTRICI input files:
attrici_in <- file.path(dir, "data", "Counterfactual_env_data", 
                         "ISIMIP_GSWP3_W5E5", "attrici_detrending", "input_files")

# directory to store postprocessed ATTRICI output:
res_dir_proj <- file.path(attrici_out, "ATTRICI_CLIM_ESRI102003_tifs")
if(!dir.exists(res_dir_proj)){dir.create(res_dir_proj)}

# directory to store detrended bioclimatic rasters:
bioclim_folder <- file.path(attrici_out, "ATTRICI_CLIM_ESRI102003_tifs", "bioclim")
if(!dir.exists(bioclim_folder)){dir.create(bioclim_folder, recursive = TRUE)}

# directory to store detrended seasonal variables:
seasonal_folder <- file.path(attrici_out, "ATTRICI_CLIM_ESRI102003_tifs", "seasonal")
if(!dir.exists(seasonal_folder)){dir.create(seasonal_folder, recursive = TRUE)}

# plots to check detrending:
plot_dir <- file.path(dir, "plots", "attrici")
if(!dir.exists(plot_dir)){dir.create(plot_dir)}


# load data: -------------------------------------------------------------------

# US outline:
US_albers_sf <- read_sf(file.path(dir, "data", "US_outline_ESRI102003.shp")) # output of 1_0_dataprep_climate.R

# selected variables:
load(file = file.path(dir, "data", "selected_variables.RData")) # selvar_final; output of 1_2a_dataprep_env_variable_selection.R
selvar_final

# BBS route selection (route centroids):
routes_sel_sf <- st_read(file.path(dir, "data", "route_selection_1995_2019_surv_beg_end_max_5y_miss_v2_spat_thin_100km_max_30_r_per_BCR_centroids.shp")) # output of 1_1_dataprep_BBS_route_selection.R
nrow(routes_sel_sf) # 539


# 1) temperature postprocessing: -----------------------------------------------

## get from tas, tasskew, tasrange to tasmin and tasmax

# calculate tasmin and tasmax from detrended tasrange and tasskew:
# using attrici postprocess-tas function throws errors coming from cdo,
# I do equivalent postprocessing here:

# see https://github.com/ISI-MIP/attrici/tree/main:
# tasskew = (tas - tasmin) / tasrange
# tasrange = tasmax - tasmin
# ->
# tasmin = tas - (tasskew * tasrange)
# tasmax = tasmin + tasrange


# counterfactual data:
nc_cfact_tas <- ncdf4::nc_open(file.path(attrici_out, "tas_detrended", "US_tas_detrended_1994_2019.nc"))
nc_cfact_tasskew <- ncdf4::nc_open(file.path(attrici_out, "tasskew_detrended", "US_tasskew_detrended_1994_2019.nc"))
nc_cfact_tasrange <- ncdf4::nc_open(file.path(attrici_out, "tasrange_detrended", "US_tasrange_detrended_1994_2019.nc"))

tas_data <- ncvar_get(nc_cfact_tas, "cfact")
tasskew_data <- ncvar_get(nc_cfact_tasskew, "cfact")
tasrange_data <- ncvar_get(nc_cfact_tasrange, "cfact")

# tasmin = tas - (tasskew * tasrange)
tasmin_data <- tas_data - (tasskew_data * tasrange_data)

dims <- nc_cfact_tas$dim

# define new variable:
tasmin_var <- ncvar_def(name = "cfact", units = "K", dim = list(dims$lon, dims$lat, dims$time), missval = -9999)

# create NetCDF file:
nc_out <- nc_create(file.path(attrici_out, "US_tasmin_detrended_1994_2019.nc"), 
                    list(tasmin_var))

# write data into the NetCDF file:
ncvar_put(nc_out, tasmin_var, tasmin_data)

# tasmax = tasmin + tasrange:

nc_cfact_tasmin <- ncdf4::nc_open(file.path(attrici_out, "US_tasmin_detrended_1994_2019.nc"))
tasmin_data <- ncvar_get(nc_cfact_tasmin, "cfact")

tasmax_data <- tasmin_data + tasrange_data 

# define new variable:
tasmax_var <- ncvar_def(name = "cfact", units = "K", dim = list(dims$lon, dims$lat, dims$time), missval = -9999)

# create NetCDF file:
nc_out <- nc_create(file.path(attrici_out, "US_tasmax_detrended_1994_2019.nc"), 
                    list(tasmax_var))

# write data into the NetCDF file:
ncvar_put(nc_out, tasmax_var, tasmax_data)

# close all NetCDF files when done
nc_close(nc_cfact_tas)
nc_close(nc_cfact_tasskew)
nc_close(nc_cfact_tasrange)
nc_close(nc_cfact_tasmin)


# 2) extract monthly tifs from nc files: ---------------------------------------

start <- 1994
end <- 2019

vars <- c("tas", "tasrange", "tasskew", "pr", "tasmin", "tasmax")

for(i in 1:length(vars)){
  
  var <- vars[i]
  print(paste(i, var))
  
  # output attrici:
  if(var %in% c("tas", "tasrange", "tasskew", "pr")){
    nc_cfact <- ncdf4::nc_open(file.path(attrici_out, paste0(var, "_detrended"),
                                         paste0("US_", var, "_detrended_1994_2019.nc"))) 
  } else {
    nc_cfact <- ncdf4::nc_open(file.path(attrici_out,
                                         paste0("US_", var, "_detrended_1994_2019.nc"))) 
  }


  # extract time and convert it to date:
  time <- ncvar_get(nc_cfact, "time")
  time_units <- ncatt_get(nc_cfact, "time", "units")$value
  time_origin <- lubridate::as_date(time_units)
  time_date <- time_origin + time
  range(time_date)
  
  # extract latitude and longitude
  lat <- ncvar_get(nc_cfact, "lat")
  lon <- ncvar_get(nc_cfact, "lon")
  
  #nc_cfact$var$cfact # detrended, counterfactual data
  #nc_cfact$var$y # factual data
  #nc_cfact$var$logp # how well do estimated model parameters match prior distribution
  
  var_arr <- ncvar_get(nc_cfact, varid = "cfact")
  #dim(var_arr) # lon, lat, time in days
  
  # to match coordinates, time and variable data in a data frame:
  lonlattime_df <- as_tibble(expand.grid(lon, lat, time_date)) %>% 
    rename("lon" = Var1, "lat" = Var2, "date" = Var3) %>% 
    mutate(year = lubridate::year(date),
           month = lubridate::month(date))
  
  # add variable values to lon-lat-time:
  var_vec_long <- as.vector(var_arr) # reshape variable data

  # add to data frame:
  var_df <- lonlattime_df %>% 
    mutate(cfact = var_vec_long) %>% 
    filter(date >= lubridate::as_date("1994-01-01"))
  
  # clean up:
  rm(var_arr)
  rm(lonlattime_df)
  
  # extract one tif per month:

  # example raster to align to:
  ex_rast <- rast(file.path(clim_path, "ISIMIP_CLIM_ESRI102003", "pr_199501_ESRI102003.tif"))
  ex_rast
  
  for(y in start:end){  
    
    print(y)
    
    for(m in 1:12){
      
      print(m)
      
      dt_export <- var_df %>% 
        filter(year == y & month == m) %>% 
        # aggregate to monthly mean:
        group_by(lon, lat) %>% 
        summarise({{var}} := mean(cfact, na.rm = TRUE)) %>%
        select(c(lon, lat, {{var}})) %>% 
        rast(crs = "EPSG:4326") %>% 
        project(y = "ESRI:102003", method = "average") %>%
        resample(y = ex_rast, method = "average") %>%
        mask(US_albers_sf) 

      writeRaster(dt_export,
                  filename = file.path(res_dir_proj, paste0(var, "_", y, stringr::str_pad(m, 2, pad = 0),  "_ESRI102003.tif")),
                  overwrite = TRUE)
    }
  }
  rm(var_vec_long)
  
  nc_close(nc_cfact)
}


# 3) calculate variables selected as covariates in DOMs: -----------------------

# register cores for parallel computation:
ncores <- 2 
cl <- makeCluster(ncores, setup_timeout = 0.5)
registerDoParallel(cl)

## bioclimatic variables: ----

start <- 1995
end <- 2019

# iterate over years:
foreach(year = start:end, 
        .packages = c("raster", "terra", "dismo") , 
        .verbose = TRUE) %dopar% {
          
          for(var in c("tasmin", "tasmax", "pr")){
            
            # corresponding filenames, June previous year to May current year:
            files <- c(paste0(var, "_", year-1, stringr::str_pad(c(6:12), width = 2, pad = "0"), "_ESRI102003.tif"),
                       paste0(var, "_", year, stringr::str_pad(c(1:5), width = 2, pad = "0"), "_ESRI102003.tif"))
            
            # bricks to store the month-wise values:
            if(var == "tasmin"){
              tasmin_June_May <- raster::stack(x = file.path(res_dir_proj, files))
            } else if(var == "tasmax"){
              tasmax_June_May <- raster::stack(x = file.path(res_dir_proj, files))
            } else {
              pr_June_May <- raster::stack(x = file.path(res_dir_proj, files))
            }
          }
          
          # calculate bioclimatic variables:
          biovars <- dismo::biovars(prec = pr_June_May,  # monthly = 12 variables for each variable; Raster brick
                                    tmin = tasmin_June_May, 
                                    tmax = tasmax_June_May)
          
          biovars_rast <- terra::rast(biovars) # convert to terra object
          
          # save tifs:
          terra::writeRaster(biovars_rast,
                             filename = file.path(bioclim_folder, 
                                                  paste0(names(biovars), "_", year, ".tif")), 
                             overwrite = TRUE)
        }


## seasonal temperature and precipitation summaries: ----

selvar_final # pr_mean_spring, pr_mean_summer, pr_mean_autumn, pr_mean_winter

months_seasons_ls <- list(spring = stringr::str_pad(c(3:5), width = 2, pad = "0"),
                          summer = stringr::str_pad(c(6:8), width = 2, pad = "0"),
                          autumn = stringr::str_pad(c(9:11), width = 2, pad = "0"),
                          winter = stringr::str_pad(c(12, 1:2), width = 2, pad = "0"))

# iterate over years:
foreach(year = 1995:2019, 
        .packages = c("raster", "terra") , 
        .verbose = TRUE) %dopar% {
          
          # mean precipitation:
          var <- "pr"
          
          # iterate over seasons:
          for(s in 1:4){
            
            # corresponding filenames:
            if(s < 4){
              # spring, summer, autumn:
              files <- paste0(var, "_", year, months_seasons_ls[[s]], "_ESRI102003.tif")
            } else {
              # winter: December previous year and Jan. and Febr. current year:
              files <- c(paste0(var, "_", year-1, "12", "_ESRI102003.tif"),
                         paste0(var, "_", year, c("01", "02"), "_ESRI102003.tif"))
            }
            
            # mean values for season:
            rast_dt <- terra::rast(x = file.path(res_dir_proj, files))
            mean_rast <- terra::mean(rast_dt)
            names(mean_rast) <- paste0(var, "_mean_", names(months_seasons_ls)[[s]])
            
            # save tifs:
            terra::writeRaster(mean_rast,
                               filename = file.path(seasonal_folder, 
                                                    paste0(var, "_mean_", names(months_seasons_ls)[s], "_", year, ".tif")), 
                               overwrite = TRUE)
          }
        }

stopCluster(cl)


# 4) explorative plots: --------------------------------------------------------

## check smoothed global mean temperature file: ----

# load raw global mean temperature time series:
gmt_raw <- nc_open(file.path(dir, "data", "Counterfactual_env_data", "ISIMIP_GSWP3_W5E5", "attrici_detrending",
                             "input_files", "gswp3-w5e5_gmt_raw_1994_2019.nc"), verbose = TRUE)

# extract the time variable and convert it to date:
time <- ncvar_get(gmt_raw, "time")
time_units <- ncatt_get(gmt_raw, "time", "units")$value
time_origin <- lubridate::as_date(time_units)
time_date <- time_origin + time

# extract mean temperature:
gmt <- ncvar_get(gmt_raw, "tas", verbose = TRUE)
gmt_C <- gmt - 273.15 # Kelvin to Celsius

# compile data frame:
df <- data.frame("date" = time_date, "gmt_raw" = gmt_C)

# plot:
ggplot(df) +
  geom_line(aes(x = date, y = gmt_raw)) +
  labs(y = "mean global temperature [°C]") +
  theme_bw()

nc_close(gmt_raw)

# load smoothed global mean temperature time series:

gmt_smoothed <- nc_open(file.path(dir, "data", "Counterfactual_env_data", "ISIMIP_GSWP3_W5E5", "attrici_detrending",
                                  "input_files", "gswp3-w5e5_ssa_gmt_1994_2019.nc"))

print(gmt_smoothed)

# extract the time variable and convert it to date:
time <- ncvar_get(gmt_smoothed, "time")
time_units <- ncatt_get(gmt_smoothed, "time", "units")$value
time_origin <- lubridate::as_date(time_units)
time_date2 <- time_origin + time

# extract mean temperature:
gmt_ssa <- ncvar_get(gmt_smoothed, "tas")
gmt_ssa_C <- gmt_ssa - 273.15 # Kelvin to Celsius

# compile data frame:
df_smoothed <- data.frame("date_smoothed" = time_date2, "gmt_smoothed" = gmt_ssa_C)

# plot:
gmt_plot <- ggplot() +
  geom_line(aes(x = date, y = gmt_raw), data = df) +
  geom_line(aes(x = date_smoothed, y = gmt_smoothed), data = df_smoothed, colour = "pink", alpha = 0.7, linewidth = 2) +
  labs(y = "mean global temperature [°C]") +
  theme_bw() +
  ggtitle("Smoothed global mean temperature, input for ATTRICI")
gmt_plot

ggsave(filename = file.path(plot_dir, "gmt_raw_smoothed.svg"),
       plot = gmt_plot,
       width = 8, height = 5)

nc_close(gmt_smoothed)


## compare time series for single locations: ----

# functions:

# compile data to plot factual and counterfactual time series:
compile_plot_data <- function(var, timestep = "month", BBS_routes = routes_sel_sf){
  
  if(timestep == "month"){
    
    # months:
    d <- paste0(rep(1995:2019, each = length(stringr::str_pad(1:12, 2, pad = 0))),
                stringr::str_pad(1:12, 2, pad = 0))
    d <- lubridate::as_date(d, format = "%Y%m")
    
  } else if (timestep == "year"){
    # years:
    d <- lubridate::as_date(paste(1995:2019, "01-01", sep="-"), format = "%Y-%m-%d")
  }
  
  
  if(grepl("bio", var)){
    
    # factual time series:
    files_f <- list.files(file.path(clim_path, "bioclim"), pattern = paste0(var, "_(", paste0(1995:2019, collapse = "|"), ").tif"), full.names = TRUE)
    # counterfactual time series:
    files_cf <- list.files(file.path(res_dir_proj, "bioclim"), pattern = paste0(var, "_(", paste0(1995:2019, collapse = "|"), ").tif"), full.names = TRUE)
    
  } else if(grepl("pr_mean", var)){
    
    # factual time series:
    files_f <- list.files(file.path(clim_path, "seasonal"), pattern = paste0(var, "_(", paste0(1995:2019, collapse = "|"), ").tif"), full.names = TRUE)
    # counterfactual time series:
    files_cf <- list.files(file.path(res_dir_proj, "seasonal"), pattern = paste0(var, "_(", paste0(1995:2019, collapse = "|"), ").tif"), full.names = TRUE)
    
  } else if(var %in% c("pr", "tasmin", "tasmax", "tas")){
    
    # factual:
    files_f <- list.files(file.path(clim_path, "ISIMIP_CLIM_ESRI102003"),
                          pattern = paste0(var, "_(", paste0(1995:2019, collapse = "|"), ")"), full.names = TRUE)
    
    # counterfactual:
    files_cf <- list.files(res_dir_proj,
                           pattern = paste0(var, "_(", paste0(1995:2019, collapse = "|"), ")"), full.names = TRUE)
  
    } else if (var %in% c("tasrange", "tasskew")){
    
    # factual:
    files_f <- list.files(file.path(attrici_in, "ATTRICI_CLIM_ESRI102003_tifs"),
                          pattern = paste0(var, "_(", paste0(1995:2019, collapse = "|"), ")"), full.names = TRUE)
    
    # counterfactual:
    files_cf <- list.files(res_dir_proj,
                           pattern = paste0(var, "_(", paste0(1995:2019, collapse = "|"), ")"), full.names = TRUE)
    }
  
  # load climatic data:
  # factual:
  factual_rast <- rast(files_f)
  
  # counterfactual time series:
  cfactual_rast <- rast(files_cf)
  
  # plot time series for BBS routes:
  
  # extract data at route locations:
  routes_f <- extract(x = factual_rast, y = BBS_routes, cells = TRUE, ID = FALSE)
  routes_f <- cbind(routes_f, "RTENO" = BBS_routes$RTENO_BBS)
  routes_cf <- extract(x = cfactual_rast, y = BBS_routes, cells = TRUE, ID = FALSE)
  routes_cf <- cbind(routes_cf, "RTENO" = BBS_routes$RTENO_BBS)
  
  # assemble data frame for plotting:
  plot_df <- routes_f %>%
    tidyr::pivot_longer(cols = !c(cell, RTENO), values_to = "factual") %>%
    cbind("date" = rep(d, times = nrow(routes_sel_sf))) %>%
    cbind("counterfactual" = tidyr::pivot_longer(data = routes_cf, cols = !c(cell, RTENO), values_to = "counterfactual")$counterfactual) %>%
    select(-c(name)) %>%
    tidyr::pivot_longer(cols = c(factual, counterfactual), values_to = "value", names_to = "scenario") %>%
    mutate(cell = factor(cell))
  
  return(plot_df)
  
}


# plot time series:
plot_ts_fun <- function(var, route_subset, plot_data){
  
  plot_list <- vector(mode = "list", length = length(route_subset))
  
  for(i in 1:length(route_subset)){
    
    print(i)
    
    route_id <- route_subset[i]
    
    p_ts <- plot_data %>%
      filter(RTENO == route_id) %>%
      ggplot(aes(x = date, y = value, colour = scenario, fill = scenario, linewidth = scenario)) +
      geom_smooth(method = "lm", alpha = 0.3, linewidth = 0.2) +
      geom_line() +
      scale_colour_manual(values = c("counterfactual" = "#0D98BA",
                                     "factual" = "#85CB33")) +
      scale_fill_manual(values = c("counterfactual" = "#0D98BA",
                                   "factual" = "#85CB33")) +
      scale_linewidth_manual(values = c("counterfactual" = 0.2,
                                        "factual" = 0.1), guide = "none") +
      guides(fill = guide_legend(title = "Scenario"),
             colour = guide_legend(title = "Scenario")) +
      ylab(var) +
      ggtitle(paste("Route:", route_id)) +
      theme_bw() +
      theme(text = element_text(size = 11),
            plot.margin = unit(c(0.5, 2, 0.5, 0.5), "cm"),
            axis.title.x = element_blank())
    
    p_map <- ggplot() +
      geom_sf(data = routes_sel_sf, size = 0.1) +
      geom_sf(data = routes_sel_sf %>% filter(RTENO_BBS == route_id), colour = "red2", size = 1.5) +
      theme_void() +
      theme(panel.border = element_rect(fill = NA, color = "black"))
    
    gg_inset_map <- ggdraw() +
      draw_plot(p_ts) +
      draw_plot(p_map, x = 0.78, y = 0.7, width = 0.2, height = 0.2)
    
    plot_list[[i]] <- gg_inset_map
  }
  
  return(plot_list)
}


## logp maps: ----

vars <- c("pr", "tas", "tasrange", "tasskew")

logp_list <- vector(mode = "list", length = length(vars))
logp_plot_list <- vector(mode = "list", length = length(vars))

# example raster to align to:
ex_rast <- rast(file.path(clim_path, "ISIMIP_CLIM_ESRI102003", "pr_199501_ESRI102003.tif"))

for(i in 1:length(vars)){
  
  var <- vars[i]
  print(paste(i, var))
  
  # output attrici:
  nc_cfact <- ncdf4::nc_open(file.path(attrici_out, paste0(var, "_detrended"), 
                                       paste0("US_", var, "_detrended_1994_2019.nc")))
  
  # extract latitude and longitude
  lat <- ncvar_get(nc_cfact, "lat")
  lon <- ncvar_get(nc_cfact, "lon")
  
  # plot log-probability of the observed data given the estimated model parameters:
  
  var_arr_logp <- ncvar_get(nc_cfact, varid = "logp")
  dim(var_arr_logp) # lon, lat
  
  logp_list[[i]] <- as_tibble(expand.grid(lon, lat)) %>% 
    rename("lon" = Var1, "lat" = Var2) %>% 
    mutate(logp = as.vector(var_arr_logp)) %>% 
    rast(crs = "EPSG:4326") %>% 
    project(y = "ESRI:102003", method = "average") %>%
    resample(y = ex_rast, method = "average") %>%
    mask(US_albers_sf) 
  
  logp_plot_list[[i]] <- ggplot() +
    geom_spatraster(data = logp_list[[i]]) +
    # add selected routes:
    geom_sf(data = routes_sel_sf, colour = "white", shape = 4, size = 0.5) +
    scale_fill_viridis_c(na.value = NA) + 
    theme_bw() + 
    labs(title = paste(var, "logp"), fill = "logp") 

  nc_close(nc_cfact)
}

# plot:
svg(file.path(plot_dir, "logp.svg"), onefile = TRUE, width = 10, height = 5)
wrap_plots(logp_plot_list, ncol = 2, nrow = 2)
dev.off()


## plot time series for subset of routes: ----

# manually selected routes:
routes_US_sel <- c(84014006, 84089004, 84053011, 84044036, 84090001, 84025081, 84052122, 84083315, 84080002, 84061052,
                   84085022, 84014106, 84055011, 84033226, 84006013, 84091055, 84081029, 84060080, 84051113, 84038038)

# plot selected routes:
ggplot() +
  geom_sf(data = routes_sel_sf, size = 0.1) +
  geom_sf(data = routes_sel_sf %>% filter(RTENO_BBS %in% routes_US_sel), colour = "red2", size = 1.5) +
  theme_bw()

# # iterate over climate variables:
# names_long <- c("annual mean temperature", "diurnal temperature range", "isothermality", "annual temperature range",
#                 "precipitation driest month", "precipitation seasonality",
#                 "precipitation spring", "precipitation summer", "precipitation autumn", "precipitation winter")
# 
# plot_data_list <- vector(mode = "list", length = 10)
# 
# for(v in 1:length(names_long)){
#   
#   var <- selvar_final[v]
#   
#   print(var)
#   
#   plot_data_list[[v]] <- compile_plot_data(var = var, timestep = "year")
#   
#   plot_list <- plot_ts_fun(var = var, route_subset = routes_US_sel, plot_data = plot_data_list[[v]])
#   
#   svg(file.path(dir, "plots", "attrici", paste0(names_long[v], ".svg")), onefile = TRUE, width = 10, height = 20)
#   print(wrap_plots(plot_list, ncol = 2, nrow = 10))
#   dev.off()
#   
# } 


#### monthly time series extracted from netCDFs (not from tifs): -----

# route selection:
coords_sf <- routes_sel_sf %>% 
  filter(RTENO_BBS %in% routes_US_sel) %>% 
  st_transform(crs = "EPSG:4326")

coords_sf

# US mask:
US_mask <- rast(file.path(attrici_in, "US_mask.nc"))
coords_mask <- extract(x = US_mask, y = coords_sf, xy = TRUE)

vars <- c("pr", "tas", "tasskew", "tasrange", "tasmin", "tasmax")

for(v in 1:length(vars)){
  
  var <- vars[v]
  
  print(var)
  
  # factual data:
  
  if(var %in% c("tasskew", "tasrange")){
    nc_fact <- nc_open(file.path(attrici_in, paste0("US_", var, "_1994_2019.nc")))
  } else {
    nc_fact <- nc_open(file.path(clim_path, paste0("US_", var, "_1994_2019.nc")))
    
  }
  if(var %in% c("tasmin", "tasmax")){
    # counterfactual data:
    nc_cfact <- nc_open(file.path(attrici_out,  paste0("US_", var, "_detrended_1994_2019.nc")))
  } else {
    # counterfactual data:
    nc_cfact <- nc_open(file.path(attrici_out, paste0(var, "_detrended"), paste0("US_", var, "_detrended_1994_2019.nc")))
  }

  # get dates:
  # extract time and convert it to date:
  time <- ncvar_get(nc_fact, "time")
  #time2 <- ncvar_get(nc_cfact, "time") # identical
  time_units <- ncatt_get(nc_fact, "time", "units")$value
  #time_units2 <- ncatt_get(nc_cfact, "time", "units")$value # same
  time_origin <- lubridate::as_date(time_units)
  #time_origin2 <- lubridate::as_date(time_units2) # same
  
  time_date <- time_origin + time
  dates_subset <- subset(time_date, time_date >= lubridate::as_date("1994-01-01"))
  
  # iterate over routes:
  
  data_routes <- vector(mode = "list", length = nrow(coords_mask))
  
  for(i in 1:nrow(coords_mask)){
    
    print(i)
    
    lon <- coords_mask$x[i]  # longitude of location
    lat <- coords_mask$y[i] # latitude  of location
    
    # get values at location lonlat
    fact_out <- ncvar_get(nc_fact, varid = var,
                          start = c(which.min(abs(nc_fact$dim$lon$vals - lon)), # look for closest long
                                    which.min(abs(nc_fact$dim$lat$vals - lat)),  # look for closest lat
                                    1), 
                          count = c(1,1, length(dates_subset)))
    
    fact_out
    
    cfact_out <- ncvar_get(nc_cfact, varid = "cfact",
                           start = c(which.min(abs(nc_cfact$dim$lon$vals - lon)), # look for closest long
                                     which.min(abs(nc_cfact$dim$lat$vals - lat)),  # look for closest lat
                                     1), 
                           count = c(1,1, length(dates_subset))) 
    
    cfact_out
    
    
    # create dataframe
    data_routes[[i]] <- data.frame(date = dates_subset, factual_d = fact_out, counterfactual_d = cfact_out) %>% 
      mutate(month = month(dates_subset),
             year = year(dates_subset)) %>% 
      group_by(month, year) %>% 
      summarise(factual = mean(factual_d),
                counterfactual = mean(counterfactual_d)) %>% 
      mutate(date = as_date(paste0(year, "-", month, "-01"))) %>% 
      tidyr::pivot_longer(cols = c(factual, counterfactual), names_to = "scenario", values_to = "value") %>% 
      mutate(RTENO = coords_sf$RTENO_BBS[i])
    data_routes
    
  }
  data_routes_all <- do.call(rbind, data_routes)
  
  plot_list <- plot_ts_fun(var = var, route_subset = routes_US_sel, plot_data = data_routes_all)
  
  svg(file.path(dir, "plots", "attrici", paste0(var, "_sel_routes_monthly.svg")), onefile = TRUE, width = 13, height = 20)
  print(wrap_plots(plot_list, ncol = 2, nrow = 10))
  dev.off()
  
  nc_close(nc_fact)
  nc_close(nc_cfact)
}


### daily time series: ---- 

# plot first few years or last few years, adjust code

vars <- c("pr", "tas", "tasskew", "tasrange","tasmin", "tasmax") 

for(v in 1:length(vars)){
  
  var <- vars[v]
  print(var)
  
  # factual data:
  
  if(var %in% c("tasskew", "tasrange")){
    nc_fact <- nc_open(file.path(attrici_in, paste0("US_", var, "_1994_2019.nc")))
  } else {
    nc_fact <- nc_open(file.path(clim_path, paste0("US_", var, "_1994_2019.nc")))
    
  }
  
  if(var %in% c("tasmin", "tasmax")){
    # counterfactual data:
    nc_cfact <- nc_open(file.path(attrici_out,  paste0("US_", var, "_detrended_1994_2019.nc")))
  } else {
    # counterfactual data:
    nc_cfact <- nc_open(file.path(attrici_out, paste0(var, "_detrended"), paste0("US_", var, "_detrended_1994_2019.nc")))
  }
  
  # get dates:
  # extract time and convert it to date:
  time <- ncvar_get(nc_fact, "time")
  #time2 <- ncvar_get(nc_cfact, "time") # identical
  time_units <- ncatt_get(nc_fact, "time", "units")$value
  #time_units2 <- ncatt_get(nc_cfact, "time", "units")$value # same
  time_origin <- lubridate::as_date(time_units)
  #time_origin2 <- lubridate::as_date(time_units2) # same
  
  time_date <- time_origin + time
  
  #dates_subset <- subset(time_date, time_date <= lubridate::as_date("1997-12-31")) # first few years
  dates_subset <- subset(time_date, time_date >= lubridate::as_date("2016-01-01")) # first few years
  
  # iterate over routes:
  
  data_routes <- vector(mode = "list", length = nrow(coords_mask))
  
  for(i in 1:nrow(coords_mask)){
    
    print(i)
    
    lon <- coords_mask$x[i]  # longitude of location
    lat <- coords_mask$y[i] # latitude  of location
    
    # get values at location lonlat
    # fact_out <- ncvar_get(nc_fact, varid = var,
    #                       start = c(which.min(abs(nc_fact$dim$lon$vals - lon)), # look for closest long
    #                                which.min(abs(nc_fact$dim$lat$vals - lat)),  # look for closest lat
    #                                1),
    #                       count = c(1,1, length(dates_subset))) # first few years
    
    fact_out <- ncvar_get(nc_fact, varid = var,
                          start = c(which.min(abs(nc_fact$dim$lon$vals - lon)), # look for closest long
                                    which.min(abs(nc_fact$dim$lat$vals - lat)),  # look for closest lat
                                    length(time_date) - length(dates_subset)),
                          count = c(1,1, length(dates_subset))) # last few years
    
    fact_out
    
    # cfact_out <- ncvar_get(nc_cfact, varid = "cfact",
    #                        start = c(which.min(abs(nc_cfact$dim$lon$vals - lon)), # look for closest long
    #                                 which.min(abs(nc_cfact$dim$lat$vals - lat)),  # look for closest lat
    #                                 1),
    #                        count = c(1,1, length(dates_subset))) # first few years
    
    cfact_out <- ncvar_get(nc_cfact, varid = "cfact",
                           start = c(which.min(abs(nc_cfact$dim$lon$vals - lon)), # look for closest long
                                     which.min(abs(nc_cfact$dim$lat$vals - lat)),  # look for closest lat
                                     length(time_date) - length(dates_subset)),
                           count = c(1,1, length(dates_subset))) # last few years
    cfact_out
    
    
    # create dataframe
    data_routes[[i]] <- data.frame(date = dates_subset, factual = fact_out, counterfactual = cfact_out) %>% 
      tidyr::pivot_longer(cols = c(factual, counterfactual), names_to = "scenario", values_to = "value") %>% 
      mutate(RTENO = coords_sf$RTENO_BBS[i])
    data_routes
    
  }
  data_routes_all <- do.call(rbind, data_routes)
  
  plot_list <- plot_ts_fun(var = var, route_subset = routes_US_sel, plot_data = data_routes_all)
  
  # svg(file.path(plot_dir, paste0(var, "_sel_routes_daily_1994_1997.svg")), onefile = TRUE, # first few years
  #     width = 18, height = 40)
  svg(file.path(plot_dir, paste0(var, "_sel_routes_daily_2016_2019.svg")), onefile = TRUE, # last few years
      width = 18, height = 40)
  print(wrap_plots(plot_list, ncol = 1, nrow = 20))
  dev.off()

  nc_close(nc_fact)
  nc_close(nc_cfact)
}


# session info:
writeLines(capture.output(sessionInfo()), file.path(dir, "results", "sessionInfo", "1_2e_dataprep_cf_climate_attrici_postprocessing.txt"))
