# Prepare environmental data later used as covariates in occupancy models:

# land use data:

# steps: 
# 1) get annual land use data for 1985-2019 from ISIMIP: 
## https://data.isimip.org/search/tree/ISIMIP3a/InputData/socioeconomic/landuse/histsoc/ (based on LUH2 land use states)
# 2) extract data from nc files, convert to spatial data format
# 3) transform to equal area (use square or hexagone grid?)
# 4) mask for conterminous US

# ISIMIP files:
## landuse-5crops_histsoc_annual_1901_2021.nc 
## landuse-forests-and-natural-vegetation_histsoc_annual_1901_2021.nc 
## landuse-pastures_histsoc_annual_1901_2021.nc 
## landuse-urbanareas_histsoc_annual_1901_2021.nc 

# variables:
# we discard C3 and C4 perennial crops as they occur only in few regions in the US (C3: trees, vine, C4: sugar cane)
# also I keep managed pasture and rangeland seperate for now xx
# we summarise annual crops follwoing Naimi et al. 2022:
# "The land classes considered were primary forested land (primf), primary non-forested land (primn), 
# potentially forested secondary land (secdf), potentially non-forested secondary land (secdn), 
# pasture resulting from the combination of classes of managed pasture (pastr) and rangelands (range), urban land (urban), 
# annual crops resulting from the combination of the classes of C3 and C4 annual crops (c3ann and c4 ann) with C3 nitrogen-fixing crops (c3nfx), 
# perennial crops resulting from the combination of C3 and C4 perennial crops (c3per and c4per)."


# packages:

library(ncdf4)
library(terra)
library(doParallel)
library(stars)
library(dplyr)
library(sf)

# load data: ----

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

# ISIMIP land use data:
lu_path <- file.path("data", "Env_data", "ISIMIP_land_use_and_irrigation")
# files:
lu_files <- list.files(lu_path, full.names = TRUE, pattern = ".nc")


# extract one tif for every land use variable and year: ----

# iterate over files:
for(f in 1:length(lu_files)){
  
  print(paste("file", f, lu_files[f]))
  
  nc_dt <- ncdf4::nc_open(lu_files[f]) # open nc file for reading
  
  # extract the time variable and convert it to years
  time <- ncvar_get(nc_dt, "time")
  time_units <- ncatt_get(nc_dt, "time", "units")$value
  time_origin <- lubridate::year(sub(".* since ", "", time_units))
  time_years <- time_origin + time
  
  # extract latitude and longitude variables
  lat <- ncvar_get(nc_dt, "lat")
  lon <- ncvar_get(nc_dt, "lon")
  
  # to match coordinates, time and variable data in a data frame:
  lonlattime <- as.matrix(expand.grid(lon, lat, time_years))
  var_df <- data.frame(lonlattime)
  colnames(var_df) <- c("lon", "lat", "year")
     
  # variables in file:
  vars <- names(nc_dt$var)
                          
  # for crops, summarise annual crops (following Naimi et al. 2022):
  if(grepl(pattern = "5crops", x = lu_files[f])){
    
    ann_crops_vars <- c("c3ann_irrigated", "c3ann_rainfed", "c3nfx_irrigated", "c3nfx_rainfed",
                        "c4ann_irrigated", "c4ann_rainfed")
    
    # get data for one variable:
    for(v in ann_crops_vars){
      
      print(v)
      
      var_arr <- ncvar_get(nc_dt, varid = v) # dim: lon, lat, time
      var_vec_long <- as.vector(var_arr) # reshape variable data
      
      # add to data frame:
      
      var_df <- var_df %>% 
        mutate({{v}} := var_vec_long)
    }
    
    rm(var_arr, var_vec_long)
    
    var_df2 <- var_df %>% 
      mutate(sum_annual_crops = rowSums(across(!c("lon", "lat", "year")), na.rm = TRUE),
             sum_annual_crops = if_else(if_all(!c("lon", "lat", "year", "sum_annual_crops"), is.na), NA, sum_annual_crops)) %>% # otherwise NAs become zero
      select(-c("c3ann_irrigated", "c3ann_rainfed", "c3nfx_irrigated", "c3nfx_rainfed",
                  "c4ann_irrigated", "c4ann_rainfed"))
    
    v <- "sum_annual_crops"
    
    # extract one .tif-file per focal year:
    
    for(y in 1990:2019){
      
      print(y)
      
      dt_export <- var_df2 %>% 
        filter(year == y) %>% 
        select(-year) %>% 
        rast(crs = "EPSG:4326") %>% 
        crop(ext(c(-126, -66, 24, 50))) %>% # cut extent
        terra::project(y = "ESRI:102003", method = "average") %>% # project
        terra::mask(US_albers_sf) # cut out US
      
      #plot(dt_export)
      
      # directory to store results:
      res_dir_proj <- file.path(lu_path, "ISIMIP_LU_ESRI102003")
      if(!dir.exists(res_dir_proj)){dir.create(res_dir_proj)}
      
      writeRaster(dt_export,
                  filename = file.path(res_dir_proj, paste0(v, "_", y, "_ESRI102003.tif")),
                  overwrite = TRUE)
    }
    
    rm(var_df2)
    
  } else { # all non-crop land use data:
    

      # get data for one variable:
      for(v in vars){
        
        print(v)
        
        var_arr <- ncvar_get(nc_dt, varid = v) # dim: lon, lat, time
        var_vec_long <- as.vector(var_arr) # reshape variable data
        
        # add to data frame with lon, lat and year:
        var_df2 <- var_df %>% 
          mutate({{v}} := var_vec_long)
        
        rm(var_arr, var_vec_long)
        
        # extract one .tif-file per focal year:
        
        for(y in 1990:2019){
          
          print(y)
          
          dt_export <- var_df2 %>% 
            filter(year == y) %>% 
            select(-year) %>% 
            rast(crs = "EPSG:4326") %>% 
            crop(ext(c(-126, -66, 24, 50))) %>% # cut extent
            terra::project(y = "ESRI:102003", method = "average") %>% # project
            terra::mask(US_albers_sf) # cut out US
          
          #plot(dt_export)
          
          # directory to store results:
          res_dir_proj <- file.path(lu_path, "ISIMIP_LU_ESRI102003")
          if(!dir.exists(res_dir_proj)){dir.create(res_dir_proj)}
          
          writeRaster(dt_export,
                      filename = file.path(res_dir_proj, paste0(v, "_", y, "_ESRI102003.tif")),
                      overwrite = TRUE)
          
        }
        rm(var_df2)
      }
  }
  rm(var_df)
}


# summarise land use variables for three years before focal period as predictor for initial occupancy: ----

start <- 1995 # start year

lu_classes <- unique(gsub(x = list.files(res_dir_proj),
                          pattern = "_[0-9]{4}_ESRI102003.tif", 
                          replacement = ""))

# iterate over variables:
for(luc in lu_classes){
  
  print(luc)
  
  luc_3yrs <- rast(file.path(res_dir_proj, paste0(luc, "_", c((start-1):(start-3)), "_ESRI102003.tif")))
  names(luc_3yrs) <- c((start-1):(start-3))
  luc_3yrs_mean <- mean(luc_3yrs)
  names(luc_3yrs_mean) <- luc
  
  writeRaster(luc_3yrs_mean, file.path(res_dir_proj, paste0(luc, "_mean_", start-3, "_", start-1, "_ESRI102003.tif")),
              overwrite = TRUE)
  
}


# explorations: ----

# explore/plot land use variables:
files <- list.files(res_dir_proj, full.names = TRUE)
lu_rast <- terra::rast(files[which(grepl("2010", files))]) # 2010
dir.create("plots/land_use_classes_2010_2")
for(i in 1:nlyr(lu_rast)){
  jpeg(file = file.path("plots", "land_use_classes_2010_2", paste0(names(lu_rast)[i], ".jpg")), 
       width = 800, height = 500, quality = 100)
  terra::plot(lu_rast[[i]], main = names(lu_rast)[i])
  dev.off()
}

# plot overall trend:
var <- "secondary_nonforests"
lu_rast <- terra::rast(files[which(grepl(paste0(var, "_[0-9]{4}_ESRI102003"), files))])
values_df <- values(lu_rast, dataframe = TRUE) 
#dim(values_df) # each column = one year

plot(x = 1990:2019, y = colSums(values_df, na.rm = TRUE)/(colSums(values_df, na.rm = TRUE)[1]), type = "o", 
     main = var, ylim = c(-2, 2))


# future land use