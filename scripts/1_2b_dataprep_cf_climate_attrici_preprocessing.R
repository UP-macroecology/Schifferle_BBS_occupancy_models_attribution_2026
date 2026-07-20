# Script:   1_2b_dataprep_cf_climate_attrici_preprocessing.R
# Purpose:  Create US mask as netCDF file as input for detrending climate data with ATTRICI
# Inputs:   data/US_outline_ESRI102003.shp
#           <clim_path>/gswp3-w5e5_obsclim_<var>_lat24.0to50.0lon-126.0to-66.0_daily_<year1>_<year2>.nc (example file)
# Outputs:  data/Counterfactual_env_data/ISIMIP_GSWP3_W5E5/attrici_detrending/input_files/US_mask.nc
# Runs on:  Local
# Notes:    Detrending is done with the external command-line tool ATTRICI (Mengel et al. 2021), 
#           data preparation for this is happening here.


source(file.path("scripts", "0_paths.R"))


# packages: --------------------------------------------------------------------

library(ncdf4)
library(dplyr)
library(sf)
library(terra)
library(ggplot2)


# directories: -----------------------------------------------------------------

# directory to store ATTRICI input:
output_path <- file.path(dir, "data", "Counterfactual_env_data",
                         "ISIMIP_GSWP3_W5E5", "attrici_detrending", "input_files")
if(!dir.exists(output_path)){dir.create(output_path, recursive = TRUE)}


# 1) prepare input for ATTRICI: ------------------------------------------------

## create mask for the conterminous USA: -----

# US outline:
US_albers_sf <- read_sf(file.path(dir, "data", "US_outline_ESRI102003.shp")) # output of 1_0_dataprep_climate.R
# transform to WGS84 (to match ISIMIP data):
US_albers_sf_wgs84 <- st_transform(US_albers_sf, crs = 4326)

# example ISIMIP file:
pr_files <- list.files(clim_path, pattern = "_pr_.*\\.nc$", full.names = TRUE)
pr1 <- ncdf4::nc_open(pr_files[1])

# create mask:
nc_rast <- rast(pr_files[1])
mask_test <- terra::mask(nc_rast$pr_1, US_albers_sf_wgs84)
mask_test <-  ifel(is.na(mask_test), 0, 1) # 1s for cells within US, 0s otherwise
names(mask_test) <- "mask"

# write as netCDF-file:

# extract values:
mask_values <- as.integer(values(mask_test))

# define dimensions:
lon <- xFromCol(mask_test, seq_len(ncol(mask_test)))  
lat <- yFromRow(mask_test, seq_len(nrow(mask_test)))  
londim <- ncdf4::ncdim_def(name = "lon", units = "degrees_east", vals = lon, longname = "Longitude")
latdim <- ncdf4::ncdim_def(name = "lat", units = "degrees_north", vals = lat, longname = "Latitude")

# define variable:
vardef <- ncdf4::ncvar_def(name = "mask", units = "", dim = list(londim, latdim), prec = "integer")

# create NetCDF file:
ncpath <- file.path(output_path, "US_mask.nc")
ncout <- ncdf4::nc_create(ncpath, vardef)
# fill in values:
ncdf4::ncvar_put(nc = ncout, varid = vardef, vals = mask_values)
# close file:
ncdf4::nc_close(ncout)

# # check:
# plot(rast(file.path(output_path, "US_mask.nc")))
# length(which(ncvar_get(mask_test, "mask") == 1))


# ## merge climate data from 1901 to 2019 into single nc file:  --------
# # -> now instead done with CDO in 1_2c_attrici_input_preps.sh

# vars <- c("tas", "tasmin", "tasmax", "pr") # tas needed to calculate tasrange and tasskew for ATTRICI, which are then detrended and converted back to tasmin, tasmax
# 
# for(i in 1:length(vars)){
#   
#   var_name <- vars[i]
#   
#   print(var_name)
#   
#   var_files <- list.files(clim_path, pattern = paste0("obsclim_", var_name, "_.*\\.nc$"), full.names = TRUE)
#   
#   # extract data of all files:
#   data_list <- lapply(var_files, function(f) {
#     print(f)
#     nc <- nc_open(f)
#     data <- ncvar_get(nc, var_name)
#     nc_close(nc)
#     return(data)
#   })
#   
#   # merge data:
#   data_merged <- abind::abind(
#     data_list, along = 3)
#   dim(data_merged) # 120, 52, 43464
#   
#   # get variable definition and dimension from one file:
#   example_file <- nc_open(var_files[1])
#   
#   # # note: date origin not the same across nc files!
#   # example_file$dim$time # days since 1900-01-01
#   # # 1-8 (1901-1980): days since 1860-1-1
#   # # 9-12 (1981-2019): days since 1900-1-1
#   
#   var <- example_file$var[[var_name]]
#   var$dim
#   nc_close(example_file)
#   
#   var$dim[[3]]$vals <- seq(from = var$dim[[3]]$vals[1], length = dim(data_merged)[3]) # 58438
#   var$dim[[3]]$len <- dim(data_merged)[3]
#   
#   # define variable:
#   var_def <- ncvar_def(name = var$name, units = var$units, dim = var$dim,
#                        longname = var$longname)
#   
#   # create new file:
#   file_merged <- nc_create(
#     filename = file.path(output_path, paste0("US_", var_name, "_1901_2019.nc")), 
#     vars = var_def)
#   
#   # fill in values:
#   ncvar_put(nc = file_merged, varid = var_name, vals = data_merged)
#   
#   # close file:
#   nc_close(file_merged)
#   
# }

# session info:
writeLines(capture.output(sessionInfo()), file.path(dir, "results", "sessionInfo", "1_2b_dataprep_cf_climate_attrici_preprocessing.txt"))
