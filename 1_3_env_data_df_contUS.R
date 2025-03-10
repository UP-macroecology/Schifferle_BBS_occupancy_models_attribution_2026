# Extract values of selected environmental variables at each grid cell to predict with DOMs 
# to whole conterminous US

# 1) extract values of env. variables at each grid cell for each year
# 2) extract values of env. variables summarising 3 years before the focal period (as predictors of initial occupancy)

# -> df: grid cell id | x coord | y coord | year | bio1 | ...

# packages: --------------------------------------------------------------------

library(sf)
library(dplyr)
library(terra)

# directories:

# environmental data:
#env_dir <- file.path("data", "Env_data") # factual data used to fit DOMs
env_dir <- file.path("data", "Counterfactual_env_data") # counterfactual data for prediction / attribution

# output (env. data for each US grid cell):

#out_RData <- file.path("data", "US_grid_env_data.RData")
#out_csv <- file.path("data", "US_grid_env_data.csv")
out_RData <- file.path("data", "US_grid_env_data_cf.RData")
out_csv <- file.path("data", "US_grid_env_data_cf.csv")

# load data: -------------------------------------------------------------------

# selected variables:
load(file = file.path("data", "selected_variables.RData")) # selvar_final; output of 1_2_variable_selection.R
selvar_final

# outline conterminous US:
US_albers_sf <- read_sf(file.path("data", "US_outline_ESRI102003.shp")) # output of 1_2_variable_selection.R

# env. data files:
bioclim_files <- list.files(file.path(env_dir, "ISIMIP_GSWP3_W5E5", "bioclim"), full.names = TRUE)
lu_files <- list.files(file.path(env_dir, "ISIMIP_land_use_and_irrigation", "ISIMIP_LU_ESRI102003"), full.names = TRUE)
sclim_files <- list.files(file.path(env_dir, "ISIMIP_GSWP3_W5E5", "seasonal"), full.names = TRUE) #xx

years <- 1995:2019


# grid cell coordinates: -------------------------------------------------------

# extract grid cell centroid from climatic and land use rasters,
# then use grid cells for which data for climatic as well as land use variables are available
# (no land use over water, e.g. Great Lakes, but land use data for outside of cont. US)

# bioclim. and land use rasters:
i = 1
bioclim_year <- rast(bioclim_files[which(grepl(paste0("bio.{1,2}_", years[i], ".tif"), bioclim_files))])
lu_year <- rast(lu_files[which(grepl(paste0(years[i], "_ESRI102003.tif$"), lu_files))])



# extract coordinates of cell centroids within conterminous US:

plot(lu_year[[1]])

lu_coords <- lu_year %>% 
  crds(df = TRUE) %>% 
  round(digits = 2)
  
bioclim_coords <- bioclim_year %>% 
  crds(df = TRUE) %>% 
  round(digits = 2)

nrow(bioclim_coords) # 4151

# subset of cell centroids for which climatic and land use data are available:
# for climatic variables two cells more along Florida Keys, drop them
clim_lu_cells <- bind_rows(bioclim_coords, lu_coords) %>% 
  group_by(x, y) %>% 
  summarise(n = n()) %>% 
  filter(n == 2) %>% 
  select(-n) %>% 
  st_as_sf(coords = c("x", "y"), crs = "ESRI:102003") %>% 
  mutate(cellID = row_number())

ggplot(clim_lu_cells) +
  geom_sf(aes(fill = cellID, colour = cellID)) +
  theme_bw()
# cellID goes up to down from left to right
#write_sf(clim_lu_cells, file.path("data", "cell_centroids_US_ESRI102003.shp"))



# add environmental data: ------------------------------------------------------

# yearly data::

# iterate over years:
for(i in 1:length(years)){

  print(years[i])
  
  clim_lu_cells_year <- clim_lu_cells
  clim_lu_cells_year$year <- years[i]
  
  # bioclimatic variables of year i:
  bioclim_year <- rast(bioclim_files[which(grepl(paste0("bio.{1,2}_", years[i], ".tif"), bioclim_files))])
  # reduce to selected variables:
  bioclim_year_sel <- bioclim_year[[selvar_final[grepl(pattern = "bio", x = selvar_final)]]]
  # extract values of each bioclimatic variable at each grid cell:
  for(biovar in names(bioclim_year_sel)){
    clim_lu_cells_year[, biovar] <- bioclim_year_sel[[biovar]] %>% 
      terra::extract(y = clim_lu_cells_year) %>% 
      pull(biovar)
  }
  
  # land use variables of year i:
  lu_year <- rast(lu_files[which(grepl(paste0(years[i], "_ESRI102003.tif$"), lu_files))])
  # reduce to selected variables:
  lu_year_sel <- lu_year[[selvar_final[!grepl(pattern = "bio|pr_", x = selvar_final)]]] # xx
  # extract values of each land use class at each grid cell:
  for(luvar in names(lu_year_sel)){
    clim_lu_cells_year[, luvar] <- lu_year_sel[[luvar]] %>% 
      terra::extract(y = clim_lu_cells_year) %>% 
      pull(luvar)
  }
  
  # seasonal climate variables of year i:
  sclim_year <- rast(sclim_files[which(grepl(paste0("(spring|summer|autumn|winter)", "_", years[i], ".tif"), sclim_files))])
  # reduce to selected variables:
  sclim_year_sel <- sclim_year[[selvar_final[grepl(pattern = "(spring|summer|autumn|winter)", x = selvar_final)]]]
  # extract values of each variable at each grid cell:
  for(sclimvar in names(sclim_year_sel)){
    clim_lu_cells_year[, sclimvar] <- sclim_year_sel[[sclimvar]] %>% 
      terra::extract(y = clim_lu_cells_year) %>% 
      pull(sclimvar)
  }
  
  # merge data for all years:
  if(i == 1){
    clim_lu_cells_all <- clim_lu_cells_year
  } else{
    clim_lu_cells_all <- rbind(clim_lu_cells_all, clim_lu_cells_year)
  }
}


# 3 year summaries before focal period of bioclims, seasonal clims and land use:

# bioclimatic variables:
bioclim_3yrs_sp <- rast(bioclim_files[which(grepl("1992_1995", bioclim_files))])
# reduce to selected variables:
bioclim_3yrs_sp_sel <- bioclim_3yrs_sp[[selvar_final[grepl(pattern = "bio", x = selvar_final)]]]
# extract value of each bioclimatic variable at each grid cell:
for(biovar in names(bioclim_3yrs_sp_sel)){
  clim_lu_cells[, paste0(biovar, "_3yrs")] <- bioclim_3yrs_sp_sel[[biovar]] %>% 
    terra::extract(y = clim_lu_cells) %>% 
    pull(biovar)
}

# land use variables:
lu_3yrs_sp <- rast(lu_files[which(grepl("1992_1994", lu_files))])
# reduce to selected variables:
lu_3yrs_sp_sel <- lu_3yrs_sp[[selvar_final[!grepl(pattern = "bio|pr_", x = selvar_final)]]] # xx
# extract value of each lu variable at each route centroid:
for(luvar in names(lu_3yrs_sp_sel)){
  clim_lu_cells[, paste0(luvar, "_3yrs")] <- lu_3yrs_sp_sel[[luvar]] %>% 
    terra::extract(y = clim_lu_cells) %>% 
    pull(luvar)
}

# seasonal climate variables:
sclim_3yrs_sp <- rast(sclim_files[which(grepl("1992_1995", sclim_files))])
# reduce to selected variables:
sclim_3yrs_sp_sel <- sclim_3yrs_sp[[selvar_final[grepl(pattern = "(spring|summer|autumn|winter)", x = selvar_final)]]]
# extract value of each bioclimatic variable at each route centroid:
for(sclimvar in names(sclim_3yrs_sp_sel)){
  clim_lu_cells[, paste0(sclimvar, "_3yrs")] <- sclim_3yrs_sp_sel[[sclimvar]] %>% 
    terra::extract(y = clim_lu_cells) %>% 
    pull(sclimvar)
}

# match to df with yearly data:
clim_lu_cells_sf <- clim_lu_cells_all %>% 
  left_join(st_drop_geometry(clim_lu_cells), by = "cellID")

# write to file:
save(clim_lu_cells_sf, file = out_RData)
write.csv(st_drop_geometry(clim_lu_cells_sf), file = out_csv,
          row.names = FALSE)
