# Prepare covariates for occupancy models:
# extract climate and land-use data at routes centroid

## 1) extract values of bioclimatic variables and land use classes at each route for each year
## 2) extract 3 years bioclim summary for each route

# packages:

library(sf)
library(dplyr)
library(terra)


# load BBS data: ---------------------------------------------------------------

## BBS data formatted for occupancy modelling:
load(file = file.path("data", "BBS_for_occ.RData")) # output of 1_0_reformat_BBS_data.R
route_dt
nrow(route_dt) # 224124

## BBS route selection (route centroids):
routes_sel_sf <- st_read(file.path("data", "route_selection_25ys_surv_beg_end_max_5y_miss_max_30_r_per_BCR_v2_spat_thin_100km_centroids.shp")) # output of 1_1_route_selection.R
nrow(routes_sel_sf) # 539
plot(st_geometry(routes_sel_sf))

## BBS data, only selected routes and focal time period:
route_sel_dt <- route_dt %>% 
  filter(RTENO %in% routes_sel_sf$RTENO_BBS) %>% 
  filter(Year >= 1995 & Year <= 2019)
nrow(route_sel_dt) # 13475


# 1) bioclim variables and land use classes for each route-year combination: ----

# files:
bioclim_files <- list.files(file.path("data", "Env_data", "ISIMIP_CHELSA-W5E5v1.0", "bioclim"), full.names = TRUE)
lu_files <- list.files(file.path("data", "Env_data", "LUH2", "albers_proj"), full.names = TRUE)

# iterate over years:
## extract BBS data matching year i, load bioclim and land use data of year i, extract values at route centroids:
years <- seq(min(route_sel_dt$Year), 2016) # Chelsa data for historic period only until 2016, LUH2 until 2015 xx max(route_sel_dt$Year))

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
  
  # extract values of each bioclimatic variable at each relevant route location:
  for(biovar in names(bioclim_year)){
    routes_sel_year_sf[, biovar] <- bioclim_year[[biovar]] %>% 
      terra::extract(y = routes_sel_year_sf) %>% 
      pull(biovar)
  }
  
  if(years[i] <= 2015){ # LUH2 historic period only until 2015
    
    # land use variables of year i:
    lu_year <- rast(lu_files[which(grepl(paste0(years[i], "_ESRI102003_ave.tif$"), lu_files))])
    
    # extract values of each land use class at each relevant route location:
    for(luvar in names(lu_year)){
      routes_sel_year_sf[, luvar] <- lu_year[[luvar]] %>% 
        terra::extract(y = routes_sel_year_sf) %>% 
        pull(luvar)
    }
  } else {
    for(luvar in names(lu_year)){
      routes_sel_year_sf[, luvar] <- NA
    }
  }

  # merge data for all years:
  if(i == 1){
    routes_sel_all_sf <- routes_sel_year_sf
  } else{
    routes_sel_all_sf <- rbind(routes_sel_all_sf, routes_sel_year_sf)
  }
}

# xx 4 routes on the coasts not covered by bioclim raster -> recalculate raster with buffer?

# merge with BBS data:
route_sel_dt[which(route_sel_dt$Year == years[i]),]

route_sel_dt2 <- route_sel_dt %>% 
  left_join(routes_sel_all_sf, by = c(RTENO = "RTENO_BBS", Year = "Year")) %>% 
  select(-geometry)


# 2) bioclimatic variables summarising 3 years before focal period at each route: ----

bioclim_3yrs_sp <- rast(bioclim_files[which(grepl("1992_1995", bioclim_files))])
bioclim_3yrs_sp

# extract value of each bioclimatic variable at each route centroid:

for(biovar in names(bioclim_3yrs_sp)){
  
  routes_sel_sf[, paste0(biovar, "_3yrs")] <- bioclim_3yrs_sp[[biovar]] %>% 
    terra::extract(y = routes_sel_sf) %>% 
    pull(biovar)
}

# match to BBS data:
route_sel_dt3 <- route_sel_dt2 %>% 
  left_join(routes_sel_sf, by = c(RTENO = "RTENO_BBS")) %>% 
  select(-geometry)

View(route_sel_dt3)

# write to file:
save(route_sel_dt3, file = file.path("data", "route_year_env_data.RData"))
write.csv(route_sel_dt3, file = file.path("data", "route_year_env_data.csv"),
          row.names = FALSE)
