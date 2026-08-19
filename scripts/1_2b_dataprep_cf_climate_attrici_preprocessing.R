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


# directories: -----------------------------------------------------------------

# directory to store ATTRICI input:
output_path <- file.path(dir, "data", "Counterfactual_env_data",
                         "ISIMIP_GSWP3_W5E5", "attrici_detrending", "input_files")
if(!dir.exists(output_path)){dir.create(output_path, recursive = TRUE)}


# 1) prepare input for ATTRICI: ------------------------------------------------

## create mask for the conterminous USA: -----

# example ISIMIP file:
pr_files <- list.files(clim_path, pattern = "_pr_.*\\.nc$", full.names = TRUE)
pr1 <- ncdf4::nc_open(pr_files[1])

# US outline:
US_albers_sf <- read_sf(file.path(dir, "data", "US_outline_ESRI102003.shp")) # output of 1_0_dataprep_climate.R
US_albers_sf_wgs84 <- st_transform(US_albers_sf, crs = 4326)

# create mask:
nc_rast <- rast(pr_files[1])
mask_ATTR <- terra::mask(nc_rast$pr_1, US_albers_sf_wgs84)
mask_ATTR <-  ifel(is.na(mask_ATTR), 0, 1) # 1s for cells within US, 0s otherwise
names(mask_ATTR) <- "mask"

# write as netCDF-file:

# extract values:
mask_values <- as.integer(values(mask_ATTR))

# define dimensions:
lon <- xFromCol(mask_ATTR, seq_len(ncol(mask_ATTR)))  
lat <- yFromRow(mask_ATTR, seq_len(nrow(mask_ATTR)))  
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
new_mask <- rast(file.path(output_path, "US_mask.nc"))
plot(new_mask)
ext(new_mask)
length(which(mask_values == 1)) # 3514

# session info:
writeLines(capture.output(sessionInfo()), file.path(dir, "results", "sessionInfo", "1_2b_dataprep_cf_climate_attrici_preprocessing.txt"))
