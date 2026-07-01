# Script:   1_2d_dataprep_cf_climate_attrici_postprocessing.R
# Purpose:  Generate counterfactual variables from ATTRICI output to simulate counterfactual occupancy dynamics with dynamic occupancy models
# Inputs:   data/Counterfactual_env_data/ISIMIP_GSWP3_W5E5/attrici_detrending/output/<var>_detrended/US_<var>_detrended_1901_2019.nc
#           data/US_outline_ESRI102003.shp
#           data/selected_variables.RData
# Outputs:  data/Counterfactual_env_data/ISIMIP_GSWP3_W5E5/attrici_detrending/output/ATTRICI_CLIM_ESRI102003_tifs/<var>_<yyyymm>_ESRI102003.tif
#           data/Counterfactual_env_data/ISIMIP_GSWP3_W5E5/attrici_detrending/output/ATTRICI_CLIM_ESRI102003_tifs/bioclim/<var>_<year>.tif
#           data/Counterfactual_env_data/ISIMIP_GSWP3_W5E5/attrici_detrending/output/ATTRICI_CLIM_ESRI102003_tifs/seasonal/<var>_<year>.tif
# Runs on:  Local

# Steps:
# 1) convert ATTRICI output to tifs, postprocess temperature
# 2) calculate counterfactual version of the climate variables used to fit dynamic occupancy models

source(file.path("scripts", "0_paths.R"))


# packages: --------------------------------------------------------------------

library(ncdf4)
library(dplyr)
library(sf)
library(terra)
library(ggplot2)


# directories: -----------------------------------------------------------------

# directory where ATTRICI stores output:
attrici_out <- file.path(dir, "data", "Counterfactual_env_data", 
                         "ISIMIP_GSWP3_W5E5", "attrici_detrending", "output")

# directory to store postprocessed ATTRICI output:
res_dir_proj <- file.path(attrici_out, "ATTRICI_CLIM_ESRI102003_tifs")
if(!dir.exists(res_dir_proj)){dir.create(res_dir_proj)}

# directory to store detrended bioclimatic rasters:
bioclim_folder <- file.path(attrici_out, "ATTRICI_CLIM_ESRI102003_tifs", "bioclim")
if(!dir.exists(bioclim_folder)){dir.create(bioclim_folder, recursive = TRUE)}

# directory to store detrended seasonal variables:
seasonal_folder <- file.path(attrici_out, "ATTRICI_CLIM_ESRI102003_tifs", "seasonal")
if(!dir.exists(seasonal_folder)){dir.create(seasonal_folder, recursive = TRUE)}


# load data: -------------------------------------------------------------------

# US outline:
US_albers_sf <- read_sf(file.path(dir, "data", "US_outline_ESRI102003.shp")) # output of 1_0_dataprep_climate.R

# selected variables:
load(file = file.path(dir, "data", "selected_variables.RData")) # selvar_final; output of 1_2a_dataprep_env_variable_selection.R
selvar_final


# 1) postprocessing ATTRICI output: --------------------------------------------

## extract one tif per month from nc files: ----

start <- 1994
end <- 2019

vars <- c("tas", "tasrange", "tasskew", "pr")

for(i in 1:length(vars)){
  
  var <- vars[i]
  print(paste(i, var))
  
  # output attrici:
  nc_cfact <- ncdf4::nc_open(file.path(attrici_out, paste0(var, "_detrended"), 
                                       paste0("US_", var, "_detrended_1901_2019.nc")))
  
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
  #nc_cfact$var$logp # log posterior predictive density? lppd, to see predictive accuracy of model relating to GMT
  
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
  nc_close(nc_cfact)
  rm(lonlattime_df)
  
  # extract one tif per month:
  for(y in intersect((start-3):end, unique(lonlattime_df$year))){
  
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
}


## postprocess temperature: ----------------------------------------------------
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


# load files:
tasskew_files <- list.files(file.path(res_dir_proj), pattern = "tasskew", full.names = TRUE)
tasskew_rast <- rast(tasskew_files)

tasrange_files <- list.files(file.path(res_dir_proj), pattern = "tasrange", full.names = TRUE)
tasrange_rast <- rast(tasrange_files)

tas_files <- list.files(file.path(res_dir_proj), pattern = "tas_", full.names = TRUE)
tas_rast <- rast(tas_files)

# convert to tasmin, tasmax:
tasmin_rast <- tas_rast - (tasskew_rast * tasrange_rast)
tasmin_rast
tasmax_rast <- tasmin_rast + tasrange_rast
tasmax_rast

# save tasmin as separate tifs:
tasmin_names <- paste0("tasmin", "_",
                       paste0(rep(1994:2019, each = length(stringr::str_pad(1:12, 2, pad = 0))),
                              stringr::str_pad(1:12, 2, pad = 0)), "_ESRI102003.tif")

writeRaster(tasmin_rast,
            filename = file.path(res_dir_proj, tasmin_names),
            overwrite = TRUE)

# save tasmax as separate tifs:
tasmax_names <- paste0("tasmax", "_",
                       paste0(rep(1994:2019, each = length(stringr::str_pad(1:12, 2, pad = 0))),
                              stringr::str_pad(1:12, 2, pad = 0)), "_ESRI102003.tif")

writeRaster(tasmax_rast,
            filename = file.path(res_dir_proj, tasmax_names),
            overwrite = TRUE)


# 2) calculate variables selected as covariates in DOMs: -----------------------


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
              tasmin_June_May <- raster::stack(x = file.path(clim_path, files))
            } else if(var == "tasmax"){
              tasmax_June_May <- raster::stack(x = file.path(clim_path, files))
            } else {
              pr_June_May <- raster::stack(x = file.path(clim_path, files))
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
            rast_dt <- terra::rast(x = file.path(clim_path, files))
            mean_rast <- terra::mean(rast_dt)
            names(mean_rast) <- paste0(var, "_mean_", names(months_seasons_ls)[[s]])
            
            # save tifs:
            terra::writeRaster(mean_rast,
                               filename = file.path(seasonal_folder, 
                                                    paste0(var, "_mean_", names(months_seasons_ls)[s], "_", year, ".tif")), 
                               overwrite = TRUE)
          }
        }

# session info:
writeLines(capture.output(sessionInfo()), file.path(dir, "results", "sessionInfo", "1_2d_dataprep_cf_climate_attrici_postprocessing.txt"))

# # explorative plots: -----------------------------------------------------------
# 
# # compare time series for single locations:
# 
# library(rts) # to generate raster time series
# 
# # variable
# var <- "pr" #"tasmax" #"tasmin"
# 
# # months:
# d <- paste0(rep(1994:2019, each = length(stringr::str_pad(1:12, 2, pad = 0))),
#             stringr::str_pad(1:12, 2, pad = 0))
# d <- lubridate::as_date(d, format = "%Y%m")
# 
# # load climatic data:
# # factual time series:
# files_factual <- list.files(file.path(fclim_path, "ISIMIP_CLIM_ESRI102003"), 
#                             pattern = paste0(var, "_"), full.names = TRUE)
# # for focal time period:
# files_factual_s <- grep(x = files_factual, pattern = paste0(paste0("(", 1994:2019, "[0-9])", collapse = "|")), value = TRUE)
# factual_rast <- rast(files_factual_s)
# factual_ts <- rts(factual_rast, d)
# 
# # counterfactual time series:
# files_cfactual <- list.files(res_dir_proj, pattern = paste0(var, "_"), full.names = TRUE)
# cfactual_rast <- rast(files_cfactual)
# cfactual_ts <- rts(cfactual_rast, d)
# 
# # 1) choose location on map and plot factual and counterfactual time series:
# plot(factual_rast[[1]])
# c <- click(factual_rast[[1]], n=1, cell=TRUE)$cell # click on map
# 
# # restructure data:
# plot_df <- data.frame("date" = d, "counterfactual" = cfactual_ts[c][,1])
# plot_df$factual <- factual_ts[c][,1]
# 
# # plot time series:
# plot(x = plot_df$date, y = plot_df$factual, type = "l", main = var)
# points(x = plot_df$date, y = plot_df$counterfactual, type = "l", col = "red")
# legend("bottomleft", col = c("black", "red"), legend = c("factual", "counterfactual"), lty = 1)
# 
# 
# # 2) plot time series for BBS routes:
# 
# # selected routes:
# routes_sel_sf <- st_read(file.path("data", "route_selection_1995_2019_surv_beg_end_max_5y_miss_v2_spat_thin_100km_max_30_r_per_BCR_centroids.shp")) # output of 1_1_route_selection.R
# 
# # extract data:
# factual_routes <- extract(x = factual_rast, y = routes_sel_sf, cells = TRUE, ID = FALSE)
# factual_routes <- cbind(factual_routes, "RTENO" = routes_sel_sf$RTENO_BBS)
# cfactual_routes <- extract(x = cfactual_rast, y = routes_sel_sf, cells = TRUE, ID = FALSE)
# cfactual_routes <- cbind(cfactual_routes, "RTENO" = routes_sel_sf$RTENO_BBS)
# 
# # assemble data frame for plotting:
# plot_df_2 <- factual_routes %>% 
#   tidyr::pivot_longer(cols = !c(cell, RTENO), values_to = "factual") %>% 
#   cbind("date" = rep(d, times = nrow(routes_sel_sf))) %>%
#   cbind("counterfactual" = tidyr::pivot_longer(data = cfactual_routes, cols = !c(cell, RTENO), values_to = "counterfactual")$counterfactual) %>% 
#   select(-c(name)) %>% 
#   tidyr::pivot_longer(cols = c(factual, counterfactual), values_to = "value", names_to = "scenario") %>% 
#   mutate(cell = factor(cell))
# 
# # plot:
# nrow(routes_sel_sf)/12
# page <- 8
# 
# plot_df_2 %>% 
#   #filter(cell %in% unique(plot_df_2$cell)[528:539]) %>% 
#   #filter(cell == 5790) %>% 
#   ggplot(aes(x = date, y = value, colour = scenario, fill = scenario)) +
#   ggforce::facet_wrap_paginate(~RTENO, nrow = 3, ncol = 4, page = page, scales = "free_y") +
#   #facet_wrap(~cell, nrow = 3, ncol = 4) +
#   geom_line() +
#   geom_smooth(method = "lm") +
#   ggtitle(var) +
#   theme_bw() +
#   theme(text = element_text(size = 18))
# 
# # map corresponding routes:
# RTENOS_plot <- unique(plot_df_2$RTENO)[(page*12-11):(page*12)]
# plot(st_geometry(routes_sel_sf))
# points(routes_sel_sf %>% filter(RTENO_BBS %in% RTENOS_plot), col = "red")