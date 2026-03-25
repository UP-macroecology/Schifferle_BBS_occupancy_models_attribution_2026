# Generate counterfactual climate data by removing climate change from 1995,
# the beginning of our time series, onwards with the ATTRICI approach (Mengel et al. 2021)

# 1) prepare climate .nc-files from ISIMIP data that will then be used as input for
# ATTRICI (Mengel et al. 2021) to calculate counterfactual climate data
# 2) ATTRICI used via Windows command line tool, code pasted here
# 3) convert ATTRICI output to tifs, postprocess temperature
# 4) calculate variables that will then be used to simulate occupancy dynamics with fitted dynamic occupancy models


# packages: --------------------------------------------------------------------

library(ncdf4)
library(dplyr)
library(sf)
library(terra)
library(ggplot2)


# directories: -----------------------------------------------------------------

# factual climate data downloaded from ISIMIP (see 1_0_dataprep_climate.R)
fclim_path <- file.path("data", "Env_data", "ISIMIP_GSWP3_W5E5")

# directory to store ATTRICI input:
output_path <- file.path("T:", "Schifferle_BBS_occupancy_models_2023", "data", "Counterfactual_env_data",
                         "ISIMIP_GSWP3_W5E5", "attrici_detrending", "input_files")


# 1) prepare input for ATTRICI: ------------------------------------------------

## create mask for the conterminous USA: -----

# US outline:
US_albers_sf <- read_sf(file.path("data", "US_outline_ESRI102003.shp")) # output of 1_0_dataprep_climate.R
# transform to WGS84 (to match ISIMIP data):
US_albers_sf_wgs84 <- st_transform(US_albers_sf, crs = 4326)

# example ISIMIP file:
pr_files <- list.files(fclim_path, pattern = "_pr_.*\\.nc$", full.names = TRUE)
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


## merge climate data from 1901 to 2019 into single nc file: --------

vars <- c("tas", "tasmin", "tasmax", "pr") # tas needed to calculate tasrange and tasskew for ATTRICI, which are then detrended and converted back to tasmin, tasmax

for(i in 1:length(vars)){
  
  var_name <- vars[i]
  
  print(var_name)
  
  var_files <- list.files(fclim_path, pattern = paste0("obsclim_", var_name, "_.*\\.nc$"), full.names = TRUE)
  
  # extract data of all files:
  data_list <- lapply(var_files, function(f) {
    print(f)
    nc <- nc_open(f)
    data <- ncvar_get(nc, var_name)
    nc_close(nc)
    return(data)
  })
  
  # merge data:
  data_merged <- abind::abind(
    data_list, along = 3)
  dim(data_merged) # 120, 52, 43464
  
  # get variable definition and dimension from one file:
  example_file <- nc_open(var_files[1])
  
  # # note: date origin not the same across nc files!
  # example_file$dim$time # days since 1900-01-01
  # # 1-8 (1901-1980): days since 1860-1-1
  # # 9-12 (1981-2019): days since 1900-1-1
  
  var <- example_file$var[[var_name]]
  var$dim
  nc_close(example_file)
  
  var$dim[[3]]$vals <- seq(from = var$dim[[3]]$vals[1], length = dim(data_merged)[3]) # 58438
  var$dim[[3]]$len <- dim(data_merged)[3]
  
  # define variable:
  var_def <- ncvar_def(name = var$name, units = var$units, dim = var$dim,
                       longname = var$longname)
  
  # create new file:
  file_merged <- nc_create(
    filename = file.path(output_path, paste0("US_", var_name, "_1901_2019.nc")), 
    vars = var_def)
  
  # fill in values:
  ncvar_put(nc = file_merged, varid = var_name, vals = data_merged)
  
  # close file:
  nc_close(file_merged)
  
}


# 2) detrending using ATTRICI tool (Mengel et al. 2021): -----------------------

# used command line:

# (requires CDO installation on the HPC)
# cloned bootstrapping branch (recommended by Matthias Mengel): https://github.com/ISI-MIP/attrici/tree/bootstrapping
# create and activate virtual environment:
# schifferle@ecoc9:~/DEBTs/attrici$ python3 -m venv env
# schifferle@ecoc9:~/DEBTs/attrici$ source env/bin/activate
# install ATTRICI:
# (env) schifferle@ecoc9:~/DEBTs/attrici$ pip install -e .[dev]

# 1) preprocessing:

# attrici preprocess-tas /mnt/ibb_share/zurell_transfer/Schifferle_BBS_occupancy_models_2023/data/Counterfactual_env_data/ISIMIP_GSWP3_W5E5/attrici_detrending/input_files/US_tas_1901_2019.nc
# /mnt/ibb_share/zurell_transfer/Schifferle_BBS_occupancy_models_2023/data/Counterfactual_env_data/ISIMIP_GSWP3_W5E5/attrici_detrending/input_files/US_tasmin_1901_2019.nc
# /mnt/ibb_share/zurell_transfer/Schifferle_BBS_occupancy_models_2023/data/Counterfactual_env_data/ISIMIP_GSWP3_W5E5/attrici_detrending/input_files/US_tasmax_1901_2019.nc 
#/mnt/ibb_share/zurell_transfer/Schifferle_BBS_occupancy_models_2023/data/Counterfactual_env_data/ISIMIP_GSWP3_W5E5/attrici_detrending/input_files/US_tasrange_1901_2019.nc 
#/mnt/ibb_share/zurell_transfer/Schifferle_BBS_occupancy_models_2023/data/Counterfactual_env_data/ISIMIP_GSWP3_W5E5/attrici_detrending/input_files/US_tasskew_1901_2019.nc 

# 2) detrending, see scripts:
# 1_2c_attrici_US_pr.sh
# 1_2c_attrici_US_tas.sh
# 1_2c_attrici_US_tasrange.sh
# 1_2c_attrici_US_tasskew.sh

# 3) merge outputs:

# (env) schifferle@ecoc9z:~/DEBTs/attrici$ attrici merge-output /mnt/ibb_share/zurell_transfer/Schifferle_BBS_occupancy_models_2023/data/Counterfactual_env_data/ISIMIP_GSWP3_W5E5/attrici_detrending/output/pr_detrended/timeseries/pr /mnt/ibb_share/zurell_transfer/Schifferle_BBS_occupancy_models_2023/data/Counterfactual_env_data/ISIMIP_GSWP3_W5E5/attrici_detrending/output/pr_detrended/US_pr_detrended_1901_2019.nc  
# (env) schifferle@ecoc9z:~/DEBTs/attrici$ attrici merge-output /mnt/ibb_share/zurell_transfer/Schifferle_BBS_occupancy_models_2023/data/Counterfactual_env_data/ISIMIP_GSWP3_W5E5/attrici_detrending/output/tas_detrended/timeseries/tas /mnt/ibb_share/zurell_transfer/Schifferle_BBS_occupancy_models_2023/data/Counterfactual_env_data/ISIMIP_GSWP3_W5E5/attrici_detrending/output/tas_detrended/US_tas_detrended_1901_2019.nc 
# (env) schifferle@ecoc9z:~/DEBTs/attrici$ attrici merge-output /mnt/ibb_share/zurell_transfer/Schifferle_BBS_occupancy_models_2023/data/Counterfactual_env_data/ISIMIP_GSWP3_W5E5/attrici_detrending/output/tasskew_detrended/timeseries/tasskew /mnt/ibb_share/zurell_transfer/Schifferle_BBS_occupancy_models_2023/data/Counterfactual_env_data/ISIMIP_GSWP3_W5E5/attrici_detrending/output/tasskew_detrended/US_tasskew_detrended_1901_2019.nc
# (env) schifferle@ecoc9z:~/DEBTs/attrici$ attrici merge-output /mnt/ibb_share/zurell_transfer/Schifferle_BBS_occupancy_models_2023/data/Counterfactual_env_data/ISIMIP_GSWP3_W5E5/attrici_detrending/output/tasrange_detrended/timeseries/tasrange /mnt/ibb_share/zurell_transfer/Schifferle_BBS_occupancy_models_2023/data/Counterfactual_env_data/ISIMIP_GSWP3_W5E5/attrici_detrending/output/tasrange_detrended/US_tasrange_detrended_1901_2019.nc
