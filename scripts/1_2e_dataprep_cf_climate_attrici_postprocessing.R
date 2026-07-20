# Script:   1_2d_dataprep_cf_climate_attrici_postprocessing.R
# Purpose:  Generate counterfactual variables from ATTRICI output to simulate counterfactual occupancy dynamics with dynamic occupancy models
# Inputs:   data/Counterfactual_env_data/ISIMIP_GSWP3_W5E5/attrici_detrending/output/<var>_detrended/US_<var>_detrended_1901_2019.nc
#           data/US_outline_ESRI102003.shp
#           data/selected_variables.RData
#           data/route_selection_1995_2019_surv_beg_end_max_5y_miss_v2_spat_thin_100km_max_30_r_per_BCR_centroids.shp
# Outputs:  data/Counterfactual_env_data/ISIMIP_GSWP3_W5E5/attrici_detrending/output/ATTRICI_CLIM_ESRI102003_tifs/<var>_<yyyymm>_ESRI102003.tif
#           data/Counterfactual_env_data/ISIMIP_GSWP3_W5E5/attrici_detrending/output/ATTRICI_CLIM_ESRI102003_tifs/bioclim/<var>_<year>.tif
#           data/Counterfactual_env_data/ISIMIP_GSWP3_W5E5/attrici_detrending/output/ATTRICI_CLIM_ESRI102003_tifs/seasonal/<var>_<year>.tif
#           plots/attrici/<var>_USA_logp.svg
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
library(tidyterra)
library(doParallel)
library(patchwork)


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

# session info:
writeLines(capture.output(sessionInfo()), file.path(dir, "results", "sessionInfo", "1_2d_dataprep_cf_climate_attrici_postprocessing.txt"))


# explorative plots: -----------------------------------------------------------

# compare time series for single locations:

# functions:

# compile data to plot factual and counterfcatual time series:
compile_plot_data <- function(var, timestep = "month", BBS_routes = routes_sel_sf){
  
  if(timestep == "month"){
    
    # months:
    d <- paste0(rep(1994:2019, each = length(stringr::str_pad(1:12, 2, pad = 0))),
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
    
  } else if(var %in% c("pr", "tasmin", "tasmax")){
    
    # factual:
    files_f <- list.files(file.path(clim_path, "ISIMIP_CLIM_ESRI102003"),
                          pattern = paste0(var, "_(", paste0(1994:2019, collapse = "|"), ")"), full.names = TRUE)
    
    # counterfactual:
    files_cf <- list.files(res_dir_proj,
                           pattern = paste0(var, "_(", paste0(1994:2019, collapse = "|"), ")"), full.names = TRUE)
  }
  
  # load climatic data:
  # factual:
  factual_rast <- rast(files_f)
  
  # counterfactual time series:
  cfactual_rast <- rast(files_cf)
  
  # plot time series for BBS routes:
  
  # extract data at route locations:
  routes_f <- extract(x = factual_rast, y = routes_sel_sf, cells = TRUE, ID = FALSE)
  routes_f <- cbind(routes_f, "RTENO" = routes_sel_sf$RTENO_BBS)
  routes_cf <- extract(x = cfactual_rast, y = routes_sel_sf, cells = TRUE, ID = FALSE)
  routes_cf <- cbind(routes_cf, "RTENO" = routes_sel_sf$RTENO_BBS)
  
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
      scale_linewidth_manual(values = c("counterfactual" = 0.3,
                                        "factual" = 0.2), guide = "none") +
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
                                       paste0("US_", var, "_detrended_1901_2019.nc")))
  
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

#save(logp_list, file = file.path(attrici_out, "logp_tifs_list.RData"))


# plot:

svg(file.path(dir, "plots", "attrici", "logp.svg"), onefile = TRUE, width = 10, height = 5)
wrap_plots(logp_plot_list, ncol = 2, nrow = 2)
dev.off()


## selvars representative routes: ----

# subset of representative routes for different logp values:
logp_rast <- rast(logp_list)
names(logp_rast) <- c("pr", "tas", "tasrange", "tasskew")

# logp at route locations:
logp_routes <- extract(x = logp_rast, y = routes_sel_sf, cells = TRUE, ID = FALSE)
logp_routes <- cbind(logp_routes, "RTENO" = routes_sel_sf$RTENO_BBS)

summary(logp_routes$pr)
summary(logp_routes$tas)
summary(logp_routes$tasrange)
summary(logp_routes$tasskew)

quantiles <- logp_routes %>% 
  tidyr::pivot_longer(cols = vars, names_to = "var", values_to = "logp") %>% 
  group_by(var) %>% 
  summarise(max = max(logp),
            q75 = quantile(logp, 0.75),
            q50 = quantile(logp, 0.5),
            q25 = quantile(logp, 0.25),
            min = min(logp))

pr_routes <- c(logp_routes %>% filter(pr <= quantiles$min[1]) %>% arrange(desc(pr)) %>% slice(1) %>% pull(RTENO),
               logp_routes %>% filter(pr <= quantiles$q25[1]) %>% arrange(desc(pr)) %>% slice(1) %>% pull(RTENO),
               logp_routes %>% filter(pr <= quantiles$q50[1]) %>% arrange(desc(pr)) %>% slice(1) %>% pull(RTENO),
               logp_routes %>% filter(pr <= quantiles$q75[1]) %>% arrange(desc(pr)) %>% slice(1) %>% pull(RTENO),
               logp_routes %>% filter(pr <= quantiles$max[1]) %>% arrange(desc(pr)) %>% slice(1) %>% pull(RTENO))

tas_routes <- c(logp_routes %>% filter(tas <= quantiles$min[2]) %>% arrange(desc(tas)) %>% slice(1) %>% pull(RTENO),
               logp_routes %>% filter(tas <= quantiles$q25[2]) %>% arrange(desc(tas)) %>% slice(1) %>% pull(RTENO),
               logp_routes %>% filter(tas <= quantiles$q50[2]) %>% arrange(desc(tas)) %>% slice(1) %>% pull(RTENO),
               logp_routes %>% filter(tas <= quantiles$q75[2]) %>% arrange(desc(tas)) %>% slice(1) %>% pull(RTENO),
               logp_routes %>% filter(tas <= quantiles$max[2]) %>% arrange(desc(tas)) %>% slice(1) %>% pull(RTENO))

tasrange_routes <- c(logp_routes %>% filter(tasrange <= quantiles$min[3]) %>% arrange(desc(tasrange)) %>% slice(1) %>% pull(RTENO),
                logp_routes %>% filter(tasrange <= quantiles$q25[3]) %>% arrange(desc(tasrange)) %>% slice(1) %>% pull(RTENO),
                logp_routes %>% filter(tasrange <= quantiles$q50[3]) %>% arrange(desc(tasrange)) %>% slice(1) %>% pull(RTENO),
                logp_routes %>% filter(tasrange <= quantiles$q75[3]) %>% arrange(desc(tasrange)) %>% slice(1) %>% pull(RTENO),
                logp_routes %>% filter(tasrange <= quantiles$max[3]) %>% arrange(desc(tasrange)) %>% slice(1) %>% pull(RTENO))

tasskew_routes <- c(logp_routes %>% filter(tasskew <= quantiles$min[4]) %>% arrange(desc(tasskew)) %>% slice(1) %>% pull(RTENO),
                     logp_routes %>% filter(tasskew <= quantiles$q25[4]) %>% arrange(desc(tasskew)) %>% slice(1) %>% pull(RTENO),
                     logp_routes %>% filter(tasskew <= quantiles$q50[4]) %>% arrange(desc(tasskew)) %>% slice(1) %>% pull(RTENO),
                     logp_routes %>% filter(tasskew <= quantiles$q75[4]) %>% arrange(desc(tasskew)) %>% slice(1) %>% pull(RTENO),
                     logp_routes %>% filter(tasskew <= quantiles$max[4]) %>% arrange(desc(tasskew)) %>% slice(1) %>% pull(RTENO))

route_subset <- c(pr_routes, tas_routes, tasrange_routes, tasskew_routes)

# plot selected routes:
ggplot() +
  geom_sf(data = routes_sel_sf, size = 0.1) +
  geom_sf(data = routes_sel_sf %>% filter(RTENO_BBS %in% route_subset), colour = "red2", size = 1.5) +
  theme_void() +
  theme(panel.border = element_rect(fill = NA, color = "black"))


# iterate over climate variables:
names_long <- c("annual mean temperature", "diurnal temperature range", "isothermality", "annual temperature range",
                "precipitation driest month", "precipitation seasonality",
                "precipitation spring", "precipitation summer", "precipitation autumn", "precipitation winter")

plot_data_list <- vector(mode = "list", length = 10)

for(v in 1:length(names_long)){
  
  var <- selvar_final[v]
  
  print(var)
  
  plot_data_list[[v]] <- compile_plot_data(var = var, timestep = "year")
  
  plot_list <- plot_ts_fun(var = var, route_subset = route_subset, plot_data = plot_data_list[[v]])
  
  svg(file.path(dir, "plots", "attrici", paste0(names_long[v], ".svg")), onefile = TRUE, width = 10, height = 20)
  print(wrap_plots(plot_list, ncol = 2, nrow = 10))
  dev.off()
  
}


## monthly pr, tas, tasmin, tasmax: ----

# at representative routes:

vars <- c("pr", "tasmin", "tasmax")

plot_data_list <- vector(mode = "list", length = length(vars))

for(v in 1:length(vars)){
  
  var <- vars[v]
  
  print(var)
  
  plot_data_list[[v]] <- compile_plot_data(var = var, timestep = "month")
  
  plot_list <- plot_ts_fun(var = var, route_subset = route_subset, plot_data = plot_data_list[[v]])
  
  svg(file.path(dir, "plots", "attrici", paste0(var, ".svg")), onefile = TRUE, width = 13, height = 20)
  print(wrap_plots(plot_list, ncol = 2, nrow = 10))
  dev.off()
  
}
save(plot_data_list, file = file.path(attrici_out, "attrici_check_monthly_data.RData"))


### low logp routes: ----

routes_lowlogp <- unlist(lapply(X = c("pr", "tas", "tasrange", "tasskew"), 
               FUN = function(.x){logp_routes %>% arrange(!!sym(.x)) %>% slice(1:5) %>% pull(RTENO)}))

# plot selected routes:
ggplot() +
  geom_sf(data = routes_sel_sf, size = 0.1) +
  geom_sf(data = routes_sel_sf %>% filter(RTENO_BBS %in% routes_lowlogp), colour = "red2", size = 1.5) +
  theme_bw()

for(v in 1:length(vars)){
  
  var <- vars[v]
  
  print(var)
  
  plot_list <- plot_ts_fun(var = var, route_subset = routes_lowlogp, plot_data = plot_data_list[[v]])
  
  save(plot_data_list, file = file.path(attrici_out, paste0(var, ".RData")))
  
  svg(file.path(dir, "plots", "attrici", paste0(var, "_low_logp_routes.svg")), onefile = TRUE, width = 13, height = 20)
  print(wrap_plots(plot_list, ncol = 2, nrow = 10))
  dev.off()
  
}

### high logp routes: ----

routes_highlogp <- unlist(lapply(X = c("pr", "tas", "tasrange", "tasskew"), 
                                FUN = function(.x){logp_routes %>% arrange(desc(!!sym(.x))) %>% slice(1:5) %>% pull(RTENO)}))

# plot selected routes:
ggplot() +
  geom_sf(data = routes_sel_sf, size = 0.1) +
  geom_sf(data = routes_sel_sf %>% filter(RTENO_BBS %in% routes_highlogp), colour = "red2", size = 1.5) +
  theme_bw()

for(v in 1:length(vars)){
  
  var <- vars[v]
  
  print(var)
  
  plot_list <- plot_ts_fun(var = var, route_subset = routes_highlogp, plot_data = plot_data_list[[v]])
  
  save(plot_data_list, file = file.path(attrici_out, paste0(var, ".RData")))
  
  svg(file.path(dir, "plots", "attrici", paste0(var, "_high_logp_routes.svg")), onefile = TRUE, width = 13, height = 20)
  print(wrap_plots(plot_list, ncol = 2, nrow = 10))
  dev.off()
  
}

# manual selection across USA:

routes_US_sel <- c(84014006, 84089004, 84053011, 84044036, 84090001, 84025081, 84052122, 84083315, 84080002, 84061052,
                   84085022, 84014106, 84055011, 84033226, 84006013, 84091055, 84081029, 84060080, 84051113, 84038038)

# plot selected routes:
ggplot() +
  geom_sf(data = routes_sel_sf, size = 0.1) +
  geom_sf(data = routes_sel_sf %>% filter(RTENO_BBS %in% routes_US_sel), colour = "red2", size = 1.5) +
  theme_bw()

for(v in 1:length(vars)){
  
  var <- vars[v]
  
  print(var)
  
  plot_list <- plot_ts_fun(var = var, route_subset = routes_US_sel, plot_data = plot_data_list[[v]])
  
  save(plot_data_list, file = file.path(attrici_out, paste0(var, ".RData")))
  
  svg(file.path(dir, "plots", "attrici", paste0(var, "_sel_routes.svg")), onefile = TRUE, width = 13, height = 20)
  print(wrap_plots(plot_list, ncol = 2, nrow = 10))
  dev.off()
  
}
