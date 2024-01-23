# Prepare environmental data later used as covariates in occupancy models:

# climate data:

# steps: 
# 1) get climate data for 1995-2019 from Chelsa or ISIMIP, mask for conterminous US, transform to equal area (use square or hexagone grid?)
# 2) from ISIMIP I get daily data -> aggregate to monthly data
# 3) calculate bioclim variables

# packages:
library(stars)
library(doParallel)
library(gdalUtilities)
library(dplyr)

# register cores for parallel computation:
ncores <- 3 
cl <- makeCluster(ncores, setup_timeout = 0.5)
registerDoParallel(cl)

# climate data: daily ISIMIP Chelsa data to monthly means and project: -----

# xx so far only historical, add future!

# download ISIMIP Chelsa observed climate:

# CHELSA-W5E5v1.0
# e.g. https://data.isimip.org/datasets/4014abe4-32fb-46c9-b9b7-53b13f12166f/ -> configure download
# chelsa-w5e5v1.0_obsclim_tasmax_1800arcsec_global_daily
# bounding box: South: 24 North: 50 West: -126 East: -66
# 1990-2016

chelsa_path <- file.path("data", "Env_data", "ISIMIP_CHELSA-W5E5v1.0")

# iterate over variables (tasmin, tasmax, precipitation):

for(var in c("tasmin", "tasmax", "pr")) {
  
  print(var)
  
  # list files:
  zipfiles <- utils::unzip(file.path(chelsa_path, paste0(var, ".zip")), 
                           list = TRUE)
  nc_files <- grep(x = zipfiles$Name, pattern = ".nc", value = TRUE) # every nc file contains daily values of one month
  
  # unzip files:
  utils::unzip(file.path(chelsa_path, paste0(var, ".zip")), exdir = chelsa_path) # unzip the top directory
  
  
  # iterate over files / months:
  
  foreach(i = 1:length(nc_files), 
          .packages = c("gdalUtilities", "stars") , 
          .verbose = TRUE) %dopar% {
    
    nc_dt <- read_stars(file.path(chelsa_path, nc_files[i])) # uses GDAL driver for netCDF files
    # delta = cell size = 0.5°; offset: the start coordinate (or time) value of the first pixel
    names(nc_dt) <- var
    
    # set coordinate system:
    # "All global CHELSA products are in a geographic coordinate system referenced to the WGS 84 horizontal datum"
    st_crs(nc_dt) <- "EPSG:4326"
    
    # monthly mean:
    monthly_mean <- st_apply(nc_dt, c("x", "y"), mean)
    #plot(monthly_mean)
    
    # save as tif:
    filename <- stringr::str_sub(nc_files[i],-9,-4) # yearmonth
    write_stars(monthly_mean, file.path(chelsa_path, "monthly_EPSG4326", paste0(var, "_", filename, ".tif")))
    
    # reproject with gdal: (reprojecting with stars and the saving did somehow not work)
    gdalUtilities::gdalwarp(srcfile = file.path(chelsa_path, "monthly_EPSG4326", paste0(var, "_", filename, ".tif")),
                            dstfile = file.path(chelsa_path, "monthly_albers_proj", paste0(var, "_", filename, "_ESRI102003.tif")),
                            overwrite = TRUE,
                            #tr = ,
                            r = "bilinear", # resampling method
                            t_srs = "ESRI:102003",
                            dstnodata = -9999 # no data value in destination file ###
                            )
  } # xx edges of US are not covered any more after reprojection -> xx (hexagonal grids?)
}


# bioclim vars: ----

# for each year and
# for three first years (later used as covariate for initial occupancy)

# use 12 months before survey started:

#---
# when were most routes surveyed?
load(file = file.path("data", "BBS_data_merged.RData"))
# route selection:
load(file = file.path("data", "route_selection_25ys_surv_beg_end_max_5y_miss_max_30_r_per_BCR_v2_spat_thin_100km.RData"))
bbs_dt %>%
  filter(RTENO %in% sel_routes_final) %>% 
  group_by(Month) %>% 
  summarise(n = n()) %>% 
  mutate(percent = round(n/sum(n) * 100, 1))# 85 % of selected routes surveyed in June
#---


# use months from June previous year to May current year to calculate bioclim vars

# folder to store bioclim rasters:
bioclim_folder <- file.path(chelsa_path, "bioclim")
if(!dir.exists(bioclim_folder)){dir.create(bioclim_folder, recursive = TRUE)}

# iterate over years:
foreach(year = 1995:2016, 
        .packages = c("raster", "terra", "dismo") , 
        .verbose = TRUE) %dopar% {
          
          for(var in c("tasmin", "tasmax", "pr")){

            # corresponding filenames, June previous year to May current year:
            files <- c(paste0(var, "_", year-1, stringr::str_pad(c(6:12), width = 2, pad = "0"), "_ESRI102003.tif"),
                       paste0(var, "_", year, stringr::str_pad(c(1:5), width = 2, pad = "0"), "_ESRI102003.tif"))
            
            # bricks to store the month-wise values:
            if(var == "tasmin"){
              tasmin_June_May <- raster::stack(x = file.path(chelsa_path, "monthly_albers_proj", files))
            } else if(var == "tasmax"){
              tasmax_June_May <- raster::stack(x = file.path(chelsa_path, "monthly_albers_proj", files))
            } else {
              pr_June_May <- raster::stack(x = file.path(chelsa_path, "monthly_albers_proj", files))
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



# predictor for initial occupancy: bioclims summarizing three years before start:

year <- 1995 # start year

# template to store output:
out <- raster::brick(raster::raster(list.files(file.path(chelsa_path, "monthly_albers_proj"), full.names = TRUE)[1]),
                     values = FALSE)
pr_mean <- tasmin_mean <- tasmax_mean <- out

# mean for each month within 3 years prior to survey:

for(var in c("tasmin", "tasmax", "pr")){
  
  print(var)
  
  # corresponding filenames, June previous year to May three years earlier:
  files <- c(paste0(var, "_", year-3, stringr::str_pad(c(6:12), width = 2, pad = "0"), "_ESRI102003.tif"),
             paste0(var, "_", year-2, stringr::str_pad(c(1:12), width = 2, pad = "0"), "_ESRI102003.tif"),
             paste0(var, "_", year-1, stringr::str_pad(c(1:12), width = 2, pad = "0"), "_ESRI102003.tif"),
             paste0(var, "_", year, stringr::str_pad(c(1:5), width = 2, pad = "0"), "_ESRI102003.tif"))
  
  # iterate over months:
  for(m in 1:12){
    
    print(m)
    
    # monthly means:
    out[[m]] <- raster::stack(x = file.path(chelsa_path, "monthly_albers_proj",
                                            files[which(grepl(paste0(stringr::str_pad(m, width = 2, pad = "0"), "_"), files))])) %>% 
      raster::calc(fun = mean)
    
    names(out[[m]]) <- paste0("month", stringr::str_pad(m, width = 2, pad = "0"))
    
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

raster::plot(tasmin_mean_June_May)
raster::plot(tasmax_mean_June_May)
raster::plot(pr_mean_June_May)

# calculate bioclimatic variables:
biovars_init_occ <- dismo::biovars(prec = pr_mean_June_May,
                                   tmin = tasmin_mean_June_May, 
                                   tmax = tasmax_mean_June_May)

biovars_rast <- terra::rast(biovars_init_occ) # convert to terra object

# save tifs:
terra::writeRaster(biovars_rast,
                   filename = file.path(bioclim_folder, 
                                        paste0(names(biovars_init_occ), "_",year-3, "_", year, ".tif")), 
                   overwrite = TRUE)

terra::plot(biovars_rast)


# explore/plot bioclimatic variables:
files <- list.files(bioclim_folder, full.names = TRUE)
bioclims_rast <- terra::rast(files[which(grepl("2000", files))]) # 2000
bioclims_rast_scaled <- scale(bioclims_rast)
dir.create("plots/bioclim_vars_scaled_2000")
for(i in 1:19){
  jpeg(file = file.path("plots", "bioclim_vars_scaled_2000", paste0(names(bioclims_rast_scaled)[i], ".jpg")), 
      width = 800, height = 500, quality = 100)
  terra::plot(bioclims_rast_scaled[[i]], main = names(bioclims_rast_scaled)[i])
  dev.off()
}



# ---
library(ggplot2)
library(viridis)
ggplot() + 
  geom_stars(data = nc_dt2) +
  facet_wrap(~time)+
  theme_void() +
  scale_fill_viridis() #+
#scale_x_discrete(expand = c(0, 0)) +
#scale_y_discrete(expand = c(0, 0))

ggplot() + 
  geom_stars(data = monthly_mean2) +
  theme_void() +
  scale_fill_viridis()