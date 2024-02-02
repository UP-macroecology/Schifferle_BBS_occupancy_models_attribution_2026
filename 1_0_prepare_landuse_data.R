# LUH2 land use data
# states (include transitions? xx)
# downloaded to datashare
# reproject/resample to same grid as climate data
# summarise annual crops 

# packages:

library(ncdf4)
library(terra)
library(doParallel)

# register cores for parallel computation:
ncores <- 3 
cl <- makeCluster(ncores, setup_timeout = 0.5)
registerDoParallel(cl)

# load data: ----
luh2_path <- file.path("//ibb-fs01.ibb.uni-potsdam.de", "daten$", "AG26", "Arbeit", "datashare", "data", "envidat", 
                       "socioeconomic", "LUH2_v2h", "global")

nc_dt <- ncdf4::nc_open(file.path(luh2_path, "states.nc"))

# names of the land-use classes:
lu_classes <- names(nc_dt$var)[1:14]

# years for which data should be extracted:
years <- c(1987:2015) # historical period for LUH2 until 2015

# reproject and resample: ----
# as climate data, crop to same extent
# load example climate file:
clim_ex <- rast(file.path("data", "Env_data", "ISIMIP_CHELSA-W5E5v1.0", "monthly_albers_proj", "pr_199001_ESRI102003.tif"))
extent_clim <- as.vector(ext(clim_ex)) # terra objects are not exportable for parallel processing
res_clim <- res(clim_ex) # terra objects are not exportable for parallel processing

# path to store data:
luh2_output_path <- file.path("data", "Env_data", "LUH2")
# create destination folder if it doesn't exist yet:
if(!dir.exists(luh2_output_path)){
  dir.create(luh2_output_path, recursive = TRUE)
}
# path to store projected and cropped output:
luh2_output_path2 <- file.path("data", "Env_data", "LUH2", "albers_proj")
# create destination folder if it doesn't exist yet:
if(!dir.exists(luh2_output_path2)){
  dir.create(luh2_output_path2, recursive = TRUE)
}

# iterate over land use classes:
foreach(luc = lu_classes, 
        .packages = c("gdalUtilities", "terra") , 
        .verbose = TRUE) %dopar% {

  # load all LUH2 data from current class (= each year from 850-2015, for whole globe)
  luh2_luc <- rast(file.path(luh2_path, "states.nc"), subds = luc)
  names(luh2_luc) <- as.character(850:2015) # years

  # focal years:
  luh2_luc_out <- luh2_luc[[as.character(years)]]
  # change name to land use class:
  names(luh2_luc_out) <- rep(luc, length(years))

  # write the layers to file to use them with gdalUtilities:
  writeRaster(luh2_luc_out, file.path(luh2_output_path , paste0(luc, "_", years, "_global.tif")),
              overwrite = TRUE,
              NAflag = -9999)
  
  # iterate over years:
  for(year in years){
    
    gdalUtilities::gdalwarp(srcfile = file.path(luh2_output_path , paste0(luc, "_", year, "_global.tif")), 
                            dstfile = file.path(luh2_output_path2, paste0(luc, "_", year, "_ESRI102003_ave.tif")),
                            overwrite = TRUE,
                            tr = res_clim,# target resolution in meters (same as unit of target srs)
                            r = "average", # continuous cover data, resampling method
                            t_srs = "ESRI:102003", # target spatial reference
                            dstnodata = -9999,
                            srcnodata = -9999,
                            te = extent_clim[c(1,3,2,4)] # extent of output file to be created
    )
  }
}

# Naimi et al 2022: 
# The land classes considered were primary forested land (primf), primary non-forested land (primn), 
# potentially forested secondary land (secdf), potentially non-forested secondary land (secdn), 
# pasture resulting from the combination of classes of managed pasture (pastr) and rangelands (range), urban land (urban), 
# annual crops resulting from the combination of the classes of C3 and C4 annual crops (c3ann and c4 ann) with C3 nitrogen-fixing crops (c3nfx), 
# perennial crops resulting from the combination of C3 and C4 perennial crops (c3per and c4per).

# we keep C3 and C4 perennial crops distinct as they may be relevant for different species groups (C3: trees, vine, C4: sugar cane)
# also I keep managed pasture and rangeland seperate for now xx

# summarise annual crops: ----
# c3ann, c4ann, c3nfx (following Naimi et al. 2022):
for(year in years){
  
  print(year)
  
  annual_crops <- rast(file.path(luh2_output_path2, paste0(c("c3ann", "c4ann", "c3nfx"), "_", year, "_ESRI102003_ave.tif")))
  annual_crops_sum <- sum(annual_crops)
  names(annual_crops_sum) <- "sum_annual_crops"
  writeRaster(annual_crops_sum, file.path(luh2_output_path2, paste0("sum_annual_crops_", year, "_ESRI102003_ave.tif")),
              overwrite = TRUE)
}

# explore/plot land use variables:
files <- list.files(luh2_output_path2, full.names = TRUE)
lu_rast <- terra::rast(files[which(grepl("2010", files))]) # 2010
dir.create("plots/land_use_classes_2010")
for(i in 1:nlyr(lu_rast)){
  jpeg(file = file.path("plots", "land_use_classes_2010", paste0(names(lu_rast)[i], ".jpg")), 
       width = 800, height = 500, quality = 100)
  terra::plot(lu_rast[[i]], main = names(lu_rast)[i])
  dev.off()
}


# predictor for initial occupancy: land use variables summarizing three years before start: ----

year <- 1991 # start year

for(luc in c(lu_classes, "sum_annual_crops")){
  
  print(luc)
  
  luc_3yrs <- rast(file.path(luh2_output_path2, paste0(luc, "_", c((year-1):(year-3)), "_ESRI102003_ave.tif")))
  names(luc_3yrs) <- c((year-1):(year-3))
  luc_3yrs_mean <- mean(luc_3yrs)
  names(luc_3yrs_mean) <- luc
  
  writeRaster(luc_3yrs_mean, file.path(luh2_output_path2, paste0(luc, "_mean_", year-3, "_", year-1, "_ESRI102003_ave.tif")),
              overwrite = TRUE)

}



# future land use?
