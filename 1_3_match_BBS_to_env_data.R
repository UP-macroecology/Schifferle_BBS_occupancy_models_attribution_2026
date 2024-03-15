# Extract values of selected bioclimatic and land-use variables at route centroids for occupancy models:

# 1) extract values of bioclim. and land use variables at each route for each year
# 2) extract bioclim. and land use values summarising 3 years before the focal period (as predictors of initial occupancy)

# packages: --------------------------------------------------------------------

library(sf)
library(dplyr)
library(terra)

# load BBS data: ---------------------------------------------------------------

# BBS data (which routes surveyed in which years) formatted for occupancy modelling:
load(file = file.path("data", "BBS_for_occ.RData")) # output of 1_0_reformat_BBS_data.R
route_dt
nrow(route_dt) # 224124

# BBS route selection (route centroids):
routes_sel_sf <- st_read(file.path("data", "route_selection_1991_2015_surv_beg_end_max_5y_miss_v2_spat_thin_100km_max_30_r_per_BCR_centroids.shp")) # output of 1_1_route_selection.R
nrow(routes_sel_sf) # 476
plot(st_geometry(routes_sel_sf))

# BBS data, only selected routes and focal time period:
route_sel_dt <- route_dt %>%
  filter(RTENO %in% routes_sel_sf$RTENO_BBS) %>%
  filter(Year >= 1991 & Year <= 2015)
nrow(route_sel_dt) # 11900

# selected variables (output of 1_2_variable_selection.R):
load(file = file.path("data", "selected_variables.RData"))
selvar
selvar <- c(selvar, "bio1", "bio4", "bio7", "bio10", "bio11", "bio12", "bio15", "bio16", "bio17", "secdf", "primf") # xx

load(file = file.path("data", "selected_variables_seasonal.RData"))
selvar_seasonal
selvar_seasonal2 <- unique(c(selvar_seasonal, selvar)) # xx
selvar <- selvar_seasonal2
# xx

# 1) bioclim and land use variables for each route-year combination: -----------

# files:
bioclim_files <- list.files(file.path("data", "Env_data", "ISIMIP_CHELSA-W5E5v1.0", "bioclim"), full.names = TRUE)
lu_files <- list.files(file.path("data", "Env_data", "LUH2", "albers_proj"), full.names = TRUE)
sclim_files <- list.files(file.path("data", "Env_data", "ISIMIP_CHELSA-W5E5v1.0", "seasonal"), full.names = TRUE) #xx


# iterate over years:
## extract BBS data matching year i, load bioclim and land use data of year i, extract values at route centroids:
years <- seq(min(route_sel_dt$Year), max(route_sel_dt$Year))

for(i in 1:length(years)){
  
  print(i)
  
  # routes surveyed in year i:
  route_IDs_year <- route_sel_dt %>% 
    filter(Year == years[i]) %>% 
    pull(RTENO)
  routes_sel_year_sf <- routes_sel_sf %>% 
    filter(RTENO_BBS %in% route_IDs_year) %>% 
    mutate(Year = years[i])
  
  # bioclimatic variables of year i:
  bioclim_year <- rast(bioclim_files[which(grepl(paste0("bio.{1,2}_", years[i], ".tif"), bioclim_files))])
  
  # reduce to selected variables:
  bioclim_year_sel <- bioclim_year[[selvar[grepl(pattern = "bio", x = selvar)]]]
  
  # extract values of each bioclimatic variable at each relevant route location:
  for(biovar in names(bioclim_year_sel)){
    routes_sel_year_sf[, biovar] <- bioclim_year_sel[[biovar]] %>% 
      terra::extract(y = routes_sel_year_sf) %>% 
      pull(biovar)
  }
  
  # land use variables of year i:
  lu_year <- rast(lu_files[which(grepl(paste0(years[i], "_ESRI102003_ave.tif$"), lu_files))])
  
  # reduce to selected variables:
  lu_year_sel <- lu_year[[selvar[!grepl(pattern = "bio|pr_", x = selvar)]]] # xx
  
  # extract values of each land use class at each relevant route location:
  for(luvar in names(lu_year_sel)){
    routes_sel_year_sf[, luvar] <- lu_year_sel[[luvar]] %>% 
      terra::extract(y = routes_sel_year_sf) %>% 
      pull(luvar)
  }
  
  # seasonal climate variables of year i: xx
  sclim_year <- rast(sclim_files[which(grepl(paste0("(spring|summer|autumn|winter)", "_", years[i], ".tif"), sclim_files))])
  # reduce to selected variables:
  sclim_year_sel <- sclim_year[[selvar[grepl(pattern = "pr_", x = selvar)]]]
  # extract values of each variable at each relevant route location:
  for(sclimvar in names(sclim_year_sel)){
    routes_sel_year_sf[, sclimvar] <- sclim_year_sel[[sclimvar]] %>% 
      terra::extract(y = routes_sel_year_sf) %>% 
      pull(sclimvar)
  }
  
  # merge data for all years:
  if(i == 1){
    routes_sel_all_sf <- routes_sel_year_sf
  } else{
    routes_sel_all_sf <- rbind(routes_sel_all_sf, routes_sel_year_sf)
  }
}
routes_sel_all_sf

# merge with BBS data:
route_sel_env_dt1 <- route_sel_dt %>% 
  left_join(routes_sel_all_sf, by = c(RTENO = "RTENO_BBS", Year = "Year")) %>% 
  select(-geometry)


# 2) variables summarising 3 years before focal period: ------------------------

# bioclimatic variables:

bioclim_3yrs_sp <- rast(bioclim_files[which(grepl("1988_1991", bioclim_files))])

# reduce to selected variables:
bioclim_3yrs_sp_sel <- bioclim_3yrs_sp[[selvar[grepl(pattern = "bio", x = selvar)]]]

# extract value of each bioclimatic variable at each route centroid:
for(biovar in names(bioclim_3yrs_sp_sel)){
  
  routes_sel_sf[, paste0(biovar, "_3yrs")] <- bioclim_3yrs_sp_sel[[biovar]] %>% 
    terra::extract(y = routes_sel_sf) %>% 
    pull(biovar)
}

# land use variables:

lu_3yrs_sp <- rast(lu_files[which(grepl("1988_1990", lu_files))])
lu_3yrs_sp

# reduce to selected variables:
lu_3yrs_sp_sel <- lu_3yrs_sp[[selvar[!grepl(pattern = "bio|pr_", x = selvar)]]] # xx

# extract value of each lu variable at each route centroid:
for(luvar in names(lu_3yrs_sp_sel)){
  
  routes_sel_sf[, paste0(luvar, "_3yrs")] <- lu_3yrs_sp_sel[[luvar]] %>% 
    terra::extract(y = routes_sel_sf) %>% 
    pull(luvar)
}

# seasonal climate variables: xx

sclim_3yrs_sp <- rast(sclim_files[which(grepl("1988_1991", sclim_files))])

# reduce to selected variables:
sclim_3yrs_sp_sel <- sclim_3yrs_sp[[selvar[grepl(pattern = "(spring|summer|autumn|winter)", x = selvar)]]]

# extract value of each bioclimatic variable at each route centroid:
for(sclimvar in names(sclim_3yrs_sp_sel)){
  
  routes_sel_sf[, paste0(sclimvar, "_3yrs")] <- sclim_3yrs_sp_sel[[sclimvar]] %>% 
    terra::extract(y = routes_sel_sf) %>% 
    pull(sclimvar)
}


# match to BBS data:
route_sel_env_dt_final <- route_sel_env_dt1 %>% 
  left_join(routes_sel_sf, by = c(RTENO = "RTENO_BBS")) %>% 
  select(-geometry)

# write to file:
save(route_sel_env_dt_final, file = file.path("data", "route_year_env_data.RData"))
write.csv(route_sel_env_dt_final, file = file.path("data", "route_year_env_data.csv"),
          row.names = FALSE)
