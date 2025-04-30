# Prepare environmental data later used as covariates in occupancy models:

# climate data:

# steps: 
# 1) get daily climate data (tasmax, tasmin, pr) for 1990-2019 from ISIMIP: 
## https://data.isimip.org/search/tree/ISIMIP3a/InputData/climate/atmosphere/gswp3-w5e5/obsclim/
## bounding box: South: 24 North: 50 West: -126 East: -66
# 2) aggregate to monthly data
# 3) calculate bioclim variables
# 4) as alternative to bioclim quarter/month variables I calculate seasonal variables

# 5) for the attribution part I repeat the same steps using counterfactual climate data from ISIMIP:
# https://data.isimip.org/search/tree/ISIMIP3a/InputData/climate/atmosphere/gswp3-w5e5/counterclim/
# (adjust directory)


# packages:
library(ncdf4)
library(doParallel)
library(dplyr)
library(sf)
library(terra)

# directories:
# data stored as "pr.zip", "tasmin.zip", "tasmax.zip" here:
#clim_path <- file.path("data", "Env_data", "ISIMIP_GSWP3_W5E5") # factual data
#clim_path <- file.path("data", "Counterfactual_env_data", "ISIMIP_GSWP3_W5E5") # counterfactual data


# register cores for parallel computation:
ncores <- 3 
cl <- makeCluster(ncores, setup_timeout = 0.5)
registerDoParallel(cl)

start <- 1995
end <- 2019

# load data:

# outline conterminous US, to later mask SpatRasters
# library(spData)
# if (requireNamespace("sf", quietly = TRUE)) {
#   data(us_states)
# }
# US_albers_sf <- us_states %>%
#   st_union() %>%
#   st_transform(crs = "ESRI:102003")
# # save as shp:
# write_sf(US_albers_sf, file.path("data", "US_outline_ESRI102003.shp"))
US_albers_sf <- read_sf(file.path("data", "US_outline_ESRI102003.shp"))


# daily data to monthly means: -----

# iterate over variables (tasmin, tasmax, precipitation):
foreach(var = c("tasmin", "tasmax", "pr"), 
        .packages = c("terra", "ncdf4", "dplyr") , 
        .verbose = TRUE) %dopar% {
          
          print(var)
  
          # list files:
          zipfiles <- utils::unzip(file.path(clim_path, paste0(var, ".zip")), list = TRUE)
          nc_files <- grep(x = zipfiles$Name, pattern = ".nc", value = TRUE) # every nc file contains daily values of 10 years
          
          # unzip files:
          utils::unzip(file.path(clim_path, paste0(var, ".zip")), exdir = clim_path) # unzip the top directory
          
          # iterate over files:
          
          for(f in 1:length(nc_files)){
            
            print(paste(f, nc_files[f]))
            
            nc_dt <- ncdf4::nc_open(file.path(clim_path, nc_files[f])) # open nc file for reading      
                    
            # extract the time variable and convert it to date:
            time <- ncvar_get(nc_dt, "time")
            time_units <- ncatt_get(nc_dt, "time", "units")$value
            time_origin <- lubridate::as_date(time_units)
            time_date <- time_origin + time
          
            # extract latitude and longitude variables
            lat <- ncvar_get(nc_dt, "lat")
            lon <- ncvar_get(nc_dt, "lon")
            
            # to match coordinates, time and variable data in a data frame:
            lonlattime_df <- as_tibble(expand.grid(lon, lat, time_date)) %>% 
              rename("lon" = Var1, "lat" = Var2, "date" = Var3) %>% 
              mutate(year = lubridate::year(date),
                     month = lubridate::month(date))
          
            # add variable values to lon-lat-time:
            var_arr <- ncvar_get(nc_dt, var) # dim: lon, lat, time
            var_vec_long <- as.vector(var_arr) # reshape variable data
            
            # add to data frame:
            var_df <- lonlattime_df %>% 
              mutate({{var}} := var_vec_long)
            
            # extract one tif per month:
            for(y in intersect((start-3):end, unique(lonlattime_df$year))){
              
              print(y)
              
              for(m in 1:12){
                
                print(m)
                
                dt_export <- var_df %>% 
                  filter(year == y & month == m) %>% 
                  # aggregate to monthly mean:
                  group_by(lon, lat) %>% 
                    summarise({{var}} := mean(.data[[var]], na.rm = TRUE)) %>%
                  select(c(lon, lat, {{var}})) %>% 
                  rast(crs = "EPSG:4326") %>% 
                  crop(ext(c(-126, -66, 24, 50))) %>% # cut extent
                  project(y = "ESRI:102003", method = "average") %>% # project
                  mask(US_albers_sf) 

                # directory to store results:
                res_dir_proj <- file.path(clim_path, "ISIMIP_CLIM_ESRI102003")
                if(!dir.exists(res_dir_proj)){dir.create(res_dir_proj)}
                
                writeRaster(dt_export,
                            filename = file.path(res_dir_proj, paste0(var, "_", y, stringr::str_pad(m, 2, pad = 0),  "_ESRI102003.tif")),
                            overwrite = TRUE)
              }
              }
            nc_close(nc_dt)
            }
        }


# bioclim vars: ----

# for each year and
# for three first years (later used as covariate for initial occupancy)

# use 12 months before survey started:

#---
# in which month were most routes surveyed?
# load(file = file.path("data", "BBS_data_merged.RData"))
# # route selection:
# load(file = file.path("data", "route_selection_1995_2019_surv_beg_end_max_5y_miss_v2_spat_thin_100km_max_30_r_per_BCR.RData")) # output of 1_1_route_selection.R
# bbs_dt %>%
#   filter(RTENO %in% sel_routes_final) %>% 
#   group_by(Month) %>% 
#   summarise(n = n()) %>% 
#   mutate(percent = round(n/sum(n) * 100, 1))# 85 % of selected routes surveyed in June
#---


# -> use months from June previous year to May current year to calculate bioclim vars

# folder to store bioclim rasters:
bioclim_folder <- file.path(clim_path, "bioclim")
if(!dir.exists(bioclim_folder)){dir.create(bioclim_folder, recursive = TRUE)}

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
              tasmin_June_May <- raster::stack(x = file.path(clim_path, "ISIMIP_CLIM_ESRI102003", files))
            } else if(var == "tasmax"){
              tasmax_June_May <- raster::stack(x = file.path(clim_path, "ISIMIP_CLIM_ESRI102003", files))
            } else {
              pr_June_May <- raster::stack(x = file.path(clim_path, "ISIMIP_CLIM_ESRI102003", files))
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


# bioclims summarising three years before focal period as predictor for initial occupancy: ----


# template to store output (raster brick required by dismo):
out <- raster::brick(raster::raster(list.files(file.path(clim_path, "ISIMIP_CLIM_ESRI102003"), full.names = TRUE)[1]),
                     values = FALSE)
pr_mean <- tasmin_mean <- tasmax_mean <- out

# 3 year-mean for each month:

for(var in c("tasmin", "tasmax", "pr")){
  
  print(var)
  
  # corresponding filenames, May start year to June three years earlier:
  files <- c(paste0(var, "_", start-3, stringr::str_pad(c(6:12), width = 2, pad = "0"), "_ESRI102003.tif"),
             paste0(var, "_", start-2, stringr::str_pad(c(1:12), width = 2, pad = "0"), "_ESRI102003.tif"),
             paste0(var, "_", start-1, stringr::str_pad(c(1:12), width = 2, pad = "0"), "_ESRI102003.tif"),
             paste0(var, "_", start, stringr::str_pad(c(1:5), width = 2, pad = "0"), "_ESRI102003.tif"))
  
  # iterate over months:
  for(m in 1:12){
    
    print(m)
    
    # monthly means:
    out[[m]] <- raster::stack(x = file.path(clim_path, "ISIMIP_CLIM_ESRI102003",
                                            files[which(grepl(paste0(stringr::str_pad(m, width = 2, pad = "0"), "_"), files))])) %>% 
      raster::calc(fun = mean)
    
    if(var == "tasmin"){
      tasmin_mean <- out
    } else if(var == "tasmax"){
      tasmax_mean <- out
    } else {
      pr_mean <- out
    }
  }
}

# reorder June to May:
tasmin_mean_June_May <- tasmin_mean[[c(6:12, 1:5)]]
tasmax_mean_June_May <- tasmax_mean[[c(6:12, 1:5)]]
pr_mean_June_May <- pr_mean[[c(6:12, 1:5)]]

# raster::plot(tasmin_mean_June_May)
# raster::plot(tasmax_mean_June_May)
# raster::plot(pr_mean_June_May)

# calculate bioclimatic variables:
biovars_init_occ <- dismo::biovars(prec = pr_mean_June_May,
                                   tmin = tasmin_mean_June_May, 
                                   tmax = tasmax_mean_June_May)

biovars_rast <- terra::rast(biovars_init_occ) # convert to terra object

# save tifs:
terra::writeRaster(biovars_rast,
                   filename = file.path(bioclim_folder, 
                                        paste0(names(biovars_rast), "_",start-3, "_", start, ".tif")), 
                   overwrite = TRUE)

terra::plot(biovars_rast)


# explore/plot bioclimatic variables:
files <- list.files(bioclim_folder, full.names = TRUE)
bioclims_rast <- terra::rast(files[which(grepl("2000", files))]) # 2000
bioclims_rast_scaled <- terra::scale(bioclims_rast)
dir.create("plots/bioclim_vars_scaled_2000")
for(i in 1:19){
  jpeg(file = file.path("plots", "bioclim_vars_scaled_2000", paste0(names(bioclims_rast_scaled)[i], ".jpg")), 
      width = 800, height = 500, quality = 100)
  terra::plot(bioclims_rast_scaled[[i]], main = names(bioclims_rast_scaled)[i])
  dev.off()
}



# seasonal temperature and precipitation summaries (fixed months): -----

# spring = March, April, May
# summer = June, July, August
# autumn = September, October, November
# winter = December, January, February

# folder to store rasters:
seasonal_folder <- file.path(clim_path, "seasonal")
if(!dir.exists(seasonal_folder)){dir.create(seasonal_folder, recursive = TRUE)}

# months considered for each season:

months_seasons_ls <- list(spring = stringr::str_pad(c(3:5), width = 2, pad = "0"),
     summer = stringr::str_pad(c(6:8), width = 2, pad = "0"),
     autumn = stringr::str_pad(c(9:11), width = 2, pad = "0"),
     winter = stringr::str_pad(c(12, 1:2), width = 2, pad = "0"))

# iterate over years:
foreach(year = 1995:2019, 
        .packages = c("raster", "terra") , 
        .verbose = TRUE) %dopar% {
          
          # mean min. and max. temperature and precipitation:
          
          for(var in c("tasmin", "tasmax", "pr")){
            
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
              rast_dt <- terra::rast(x = file.path(clim_path, "ISIMIP_CLIM_ESRI102003", files))
              #terra::plot(rast_dt)
              mean_rast <- terra::mean(rast_dt)
              names(mean_rast) <- paste0(var, "_mean_", names(months_seasons_ls)[[s]])
              #terra::plot(mean_rast)
              
              # save tifs:
              terra::writeRaster(mean_rast,
                                 filename = file.path(seasonal_folder, 
                                                      paste0(var, "_mean_", names(months_seasons_ls)[s], "_", year, ".tif")), 
                                 overwrite = TRUE)
            }
          }
          
          # mean temperature:
          
          # iterate over seasons:
          for(s in 1:4){
            
            # corresponding filenames:
            if(s < 4){
              # spring, summer, autumn:
              files <- c(paste0("tasmin", "_", year, months_seasons_ls[[s]], "_ESRI102003.tif"),
                         paste0("tasmax", "_", year, months_seasons_ls[[s]], "_ESRI102003.tif"))
            } else {
              # winter: December previous year and Jan. and Febr. current year:
              files <- c(paste0("tasmin", "_", year-1, "12", "_ESRI102003.tif"),
                         paste0("tasmin", "_", year, c("01", "02"), "_ESRI102003.tif"),
                         paste0("tasmax", "_", year-1, "12", "_ESRI102003.tif"),
                         paste0("tasmax", "_", year, c("01", "02"), "_ESRI102003.tif"))
            }
            
            # mean values for season:
            rast_dt <- terra::rast(x = file.path(clim_path, "ISIMIP_CLIM_ESRI102003", files))
            #terra::plot(rast_dt)
            mean_rast <- terra::mean(rast_dt)
            names(mean_rast) <- paste0("tas_mean_", names(months_seasons_ls)[[s]])
            #terra::plot(mean_rast)
            
            # save tifs:
            terra::writeRaster(mean_rast,
                               filename = file.path(seasonal_folder, 
                                                    paste0("tas_mean_", names(months_seasons_ls)[s], "_", year, ".tif")), 
                               overwrite = TRUE)
          }
        }



# seasonal climate summarising three years before focal period as predictor for initial occupancy: ----

for(var in c("tasmin", "tasmax", "pr")){
  
  # iterate over seasons:
  for(s in 1:4){
    
    # corresponding filenames:
    
    if(s ==1){ # spring (before survey = 1995, 1994, 1993)
        
      files <- c(paste0(var, "_", start, months_seasons_ls[[s]], "_ESRI102003.tif"),
                 paste0(var, "_", start-1, months_seasons_ls[[s]], "_ESRI102003.tif"),
                 paste0(var, "_", start-2, months_seasons_ls[[s]], "_ESRI102003.tif"))
      
    } else if(s >= 2 & s <= 3){ # summer, autumn(before survey = 1994, 1993, 1992):
        
        files <- c(paste0(var, "_", start-1, months_seasons_ls[[s]], "_ESRI102003.tif"),
                   paste0(var, "_", start-2, months_seasons_ls[[s]], "_ESRI102003.tif"),
                   paste0(var, "_", start-3, months_seasons_ls[[s]], "_ESRI102003.tif"))
       
         } else {   # winter: December previous year and Jan. and Febr. current year (1994/1995, 1993/1994, 1992/1993):
          
          files <- c(paste0(var, "_", start-1, "12", "_ESRI102003.tif"),
                     paste0(var, "_", start, c("01", "02"), "_ESRI102003.tif"),
                     paste0(var, "_", start-2, "12", "_ESRI102003.tif"),
                     paste0(var, "_", start-1, c("01", "02"), "_ESRI102003.tif"),
                     paste0(var, "_", start-3, "12", "_ESRI102003.tif"),
                     paste0(var, "_", start-2, c("01", "02"), "_ESRI102003.tif"))
    }
    
    # mean values for season:
    rast_dt <- terra::rast(x = file.path(clim_path, "ISIMIP_CLIM_ESRI102003", files))
    #terra::plot(rast_dt)
    mean_rast <- terra::mean(rast_dt)
    names(mean_rast) <- paste0(var, "_mean_", names(months_seasons_ls)[[s]])
    #terra::plot(mean_rast)
    
    # save tifs:
    terra::writeRaster(mean_rast,
                       filename = file.path(seasonal_folder, 
                                            paste0(var, "_mean_", names(months_seasons_ls)[s], "_1992_1995.tif")), 
                       overwrite = TRUE)
  }
}

# mean temperature:

# iterate over seasons:
for(s in 1:4){
  
  # corresponding filenames:
  if(s == 1){
    # spring, summer, autumn:
    files <- c(paste0("tasmin", "_", start, months_seasons_ls[[s]], "_ESRI102003.tif"),
               paste0("tasmin", "_", start-1, months_seasons_ls[[s]], "_ESRI102003.tif"),
               paste0("tasmin", "_", start-2, months_seasons_ls[[s]], "_ESRI102003.tif"),
               paste0("tasmax", "_", start, months_seasons_ls[[s]], "_ESRI102003.tif"),
               paste0("tasmax", "_", start-1, months_seasons_ls[[s]], "_ESRI102003.tif"),
               paste0("tasmax", "_", start-2, months_seasons_ls[[s]], "_ESRI102003.tif"))
  } else if(s >= 2 & s <= 3) {
    
    files <- c(paste0("tasmin", "_", start-1, months_seasons_ls[[s]], "_ESRI102003.tif"),
               paste0("tasmin", "_", start-2, months_seasons_ls[[s]], "_ESRI102003.tif"),
               paste0("tasmin", "_", start-3, months_seasons_ls[[s]], "_ESRI102003.tif"),
               paste0("tasmax", "_", start-1, months_seasons_ls[[s]], "_ESRI102003.tif"),
               paste0("tasmax", "_", start-2, months_seasons_ls[[s]], "_ESRI102003.tif"),
               paste0("tasmax", "_", start-3, months_seasons_ls[[s]], "_ESRI102003.tif")) 
    
  } else {
    
    # winter: December previous year and Jan. and Febr. current year:
    files <- c(paste0("tasmin", "_", start-1, "12", "_ESRI102003.tif"),
               paste0("tasmin", "_", start, c("01", "02"), "_ESRI102003.tif"),
               paste0("tasmin", "_", start-2, "12", "_ESRI102003.tif"),
               paste0("tasmin", "_", start-1, c("01", "02"), "_ESRI102003.tif"),
               paste0("tasmin", "_", start-3, "12", "_ESRI102003.tif"),
               paste0("tasmin", "_", start-2, c("01", "02"), "_ESRI102003.tif"),
               paste0("tasmax", "_", start-1, "12", "_ESRI102003.tif"),
               paste0("tasmax", "_", start, c("01", "02"), "_ESRI102003.tif"),
               paste0("tasmax", "_", start-2, "12", "_ESRI102003.tif"),
               paste0("tasmax", "_", start-1, c("01", "02"), "_ESRI102003.tif"),
               paste0("tasmax", "_", start-3, "12", "_ESRI102003.tif"),
               paste0("tasmax", "_", start-2, c("01", "02"), "_ESRI102003.tif"))
  }
  
  # mean values for season:
  rast_dt <- terra::rast(x = file.path(clim_path, "ISIMIP_CLIM_ESRI102003", files))
  #terra::plot(rast_dt)
  mean_rast <- terra::mean(rast_dt)
  names(mean_rast) <- paste0("tas_mean_", names(months_seasons_ls)[[s]])
  #terra::plot(mean_rast)
  
  # save tifs:
  terra::writeRaster(mean_rast,
                     filename = file.path(seasonal_folder, 
                                          paste0("tas_mean_", names(months_seasons_ls)[s], "_1992_1995.tif")), 
                     overwrite = TRUE)
}



# explorations: ----
# plot seasonal variables: 

files <- list.files(seasonal_folder, full.names = TRUE)
sclims_rast <- terra::rast(files[which(grepl("2000", files))]) # 2000
sclims_rast_scaled <- terra::scale(sclims_rast)
dir.create("plots/seasonal_vars_scaled_2000")
for(i in 1:16){
  jpeg(file = file.path("plots", "seasonal_vars_scaled_2000", paste0(names(sclims_rast_scaled)[i], ".jpg")), 
       width = 800, height = 500, quality = 100)
  terra::plot(sclims_rast_scaled[[i]], main = names(sclims_rast_scaled)[i])
  dev.off()
}

# plot overall trend:
var <- "bio15"
files <- list.files(bioclim_folder, full.names = TRUE)
cl_rast <- terra::rast(files[which(grepl(paste0(var, "_[0-9]{4}.tif"), files))])
values_df <- values(cl_rast, dataframe = TRUE) 
#dim(values_df) # each column = one year

plot(x = 1995:2019, y = colSums(values_df, na.rm = TRUE)/(colSums(values_df, na.rm = TRUE)[1]), type = "o", 
     main = var#, 
     #ylim = c(-1, 1)
     )


# explore/plot selected climatic variables of each year:

files1 <- list.files(file.path("data", "Env_data", "ISIMIP_GSWP3_W5E5", "bioclim"), full.names = TRUE)
files2 <- list.files(file.path("data", "Env_data", "ISIMIP_GSWP3_W5E5", "seasonal"), full.names = TRUE)

dir.create("plots/clim_1995_2019")

# variable selection:
load(file = file.path("data", "selected_variables.RData")) # output of 1_2_variable_selection.R
selvar_final
sel_clim_var <- grep(pattern = "(bio)|(pr_mean)", x = selvar_final, value = TRUE) 

sel_clim_files1 <- grep(pattern = paste0(paste0(sel_clim_var, "_"), collapse = "|"), x = files1, value = TRUE)
sel_clim_files1_years <- grep(pattern = paste0(1995:2019, collapse = "|"), x = sel_clim_files1, value = TRUE)
sel_clim_files1_years <- grep(pattern = "_1992_1995", x = sel_clim_files1_years, value = TRUE, invert = TRUE)

sel_clim_files2 <- grep(pattern = paste0(sel_clim_var, collapse = "|"), x = files2, value = TRUE)
sel_clim_files2_years <- grep(pattern = paste0(1995:2019, collapse = "|"), x = sel_clim_files2, value = TRUE)
sel_clim_files2_years <- grep(pattern = "_1992_1995", x = sel_clim_files2_years, value = TRUE, invert = TRUE)

sel_clim_files <- c(sel_clim_files1_years, sel_clim_files2_years)


# scale variables:
load(file = file.path("data", "route_env_dt_scale_pars.RData")) # output of 2_1_DOM_flocker_fit_fm.R
env_scale_pars

for(v in sel_clim_var){
  
  print(v)
  
  cl_rast <- terra::rast(grep(pattern = paste0(v, "_"), x = sel_clim_files, value = TRUE))
  
  # scale:
  cl_rast_scaled <- scale(cl_rast, center = as.numeric(env_scale_pars$center[v]), scale = as.numeric(env_scale_pars$scale[v]))
  
  range <- range(values(cl_rast_scaled), na.rm = TRUE)
  
  for(y in 1:length(1995:2019)){
    
    print(y)
    
    jpeg(file = file.path("plots", "clim_1995_2019", paste0(names(cl_rast_scaled)[1], "_", (1995:2019)[y], ".jpg")), 
         width = 800, height = 500, quality = 100)
    terra::plot(cl_rast_scaled[[y]], main = paste(names(cl_rast_scaled)[1], (1995:2019)[y]), range = range)
    dev.off()
    
  }
}
