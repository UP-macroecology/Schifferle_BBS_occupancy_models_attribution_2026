# assemble environmental data at BBS route locations for predictions for ISIMIP, 1901 - 2019;
# predictors for initial occupancy: 1901-1903, for colonisation and extinction: 1904-2019;
# factual data: obsclim, counterfactual data: counterclim;
# used to predict with DOMs fitted to 1995-2019

# note: since we'll do experiments obsclim + histsoc & counterclim + histsoc
# I use factual land use data within the counterfactual data


# packages: --------------------------------------------------------------------

library(dplyr)
library(sf)
library(terra)

# obsclim or counterclim:
#env_path <- file.path("data", "Env_data")
env_path <- file.path("data", "Counterfactual_env_data")


# load data: -------------------------------------------------------------------


# BBS route selection (route centroids):
routes_sel_sf <- st_read(file.path("data", "route_selection_1995_2019_surv_beg_end_max_5y_miss_v2_spat_thin_100km_max_30_r_per_BCR_centroids.shp")) # output of 1_1_route_selection.R
nrow(routes_sel_sf) # 539
plot(st_geometry(routes_sel_sf))

# selected variables (output of 1_2_variable_selection.R):
load(file = file.path("data", "selected_variables.RData"))
selvar_final

# for predictions for ISIMIP from 1904 - 2019 (1901 - 1903 to predict initial occupancy):
years <- 1904:2019


# env. data for each route-year combination: -----------------------------------


## data for each year: ------------------------

# files:
bioclim_files <- list.files(file.path(env_path, "ISIMIP_GSWP3_W5E5", "bioclim"), full.names = TRUE)
sclim_files <- list.files(file.path(env_path, "ISIMIP_GSWP3_W5E5", "seasonal"), full.names = TRUE)
lu_files <- list.files(file.path("data", "Env_data", "ISIMIP_land_use_and_irrigation", "ISIMIP_LU_ESRI102003"), full.names = TRUE) # ! xx

# iterate over years:
## load bioclim and land use data of year i, extract values at route centroids:

for(i in 1:length(years)){
  
  print(years[i])
  
  # store values of year i in sf:
  routes_sel_year_sf <- routes_sel_sf %>% 
    mutate(Year = years[i])
  
  # bioclimatic variables of year i:
  bioclim_year <- rast(bioclim_files[which(grepl(paste0("bio.{1,2}_", years[i], ".tif"), bioclim_files))])
  # reduce to selected variables:
  bioclim_year_sel <- bioclim_year[[selvar_final[grepl(pattern = "bio", x = selvar_final)]]]
  
  # extract values of each bioclimatic variable at each relevant route location:

  for(biovar in names(bioclim_year_sel)){
    routes_sel_year_sf[, biovar] <- bioclim_year_sel[[biovar]] %>% 
      terra::extract(y = routes_sel_sf) %>% 
      pull(biovar)
  }
  
  # seasonal climate variables of year i:
  sclim_year <- rast(sclim_files[which(grepl(paste0("(spring|summer|autumn|winter)", "_", years[i], ".tif"), sclim_files))])
  # reduce to selected variables:
  sclim_year_sel <- sclim_year[[selvar_final[grepl(pattern = "pr_", x = selvar_final)]]]
  # extract values of each variable at each relevant route location:
  for(sclimvar in names(sclim_year_sel)){
    routes_sel_year_sf[, sclimvar] <- sclim_year_sel[[sclimvar]] %>% 
      terra::extract(y = routes_sel_year_sf) %>% 
      pull(sclimvar)
  }
  
  # land use variables of year i:
  lu_year <- rast(lu_files[which(grepl(paste0("[a-z]_", years[i], "_ESRI102003.tif$"), lu_files))])
  # reduce to selected variables:
  lu_year_sel <- lu_year[[selvar_final[!grepl(pattern = "bio|pr_|mean", x = selvar_final)]]] # xx
  # extract values of each land use class at each relevant route location:
  for(luvar in names(lu_year_sel)){
    routes_sel_year_sf[, luvar] <- lu_year_sel[[luvar]] %>% 
      terra::extract(y = routes_sel_year_sf) %>% 
      pull(luvar)
  }
  
  # merge data for all years:
  if(i == 1){
    routes_sel_all_sf <- routes_sel_year_sf
  } else{
    routes_sel_all_sf <- rbind(routes_sel_all_sf, routes_sel_year_sf)
  }
}
routes_sel_all_sf

# as data frame:
route_sel_env_dt1 <- routes_sel_all_sf %>% 
  st_drop_geometry()


## add variables summarising 3 years before focal period: ----------------------

# bioclimatic variables:
bioclim_3yrs_sp <- rast(bioclim_files[which(grepl(paste0(years[1]-3, "_", years[1]), bioclim_files))])

# reduce to selected variables:
bioclim_3yrs_sp_sel <- bioclim_3yrs_sp[[selvar_final[grepl(pattern = "bio", x = selvar_final)]]]
# extract value of each bioclimatic variable at each route centroid:
for(biovar in names(bioclim_3yrs_sp_sel)){
  routes_sel_sf[, paste0(biovar, "_3yrs")] <- bioclim_3yrs_sp_sel[[biovar]] %>% 
    terra::extract(y = routes_sel_sf) %>% 
    pull(biovar)
}

# seasonal climate variables:
sclim_3yrs_sp <- rast(sclim_files[which(grepl(paste0(years[1]-3, "_", years[1]), sclim_files))])
# reduce to selected variables:
sclim_3yrs_sp_sel <- sclim_3yrs_sp[[selvar_final[grepl(pattern = "(spring|summer|autumn|winter)", x = selvar_final)]]]
# extract value of each bioclimatic variable at each route centroid:
for(sclimvar in names(sclim_3yrs_sp_sel)){
  routes_sel_sf[, paste0(sclimvar, "_3yrs")] <- sclim_3yrs_sp_sel[[sclimvar]] %>% 
    terra::extract(y = routes_sel_sf) %>% 
    pull(sclimvar)
}


# land use variables:
lu_3yrs_sp <- rast(lu_files[which(grepl(paste0(years[1]-3, "_", years[1]-1), lu_files))])
# reduce to selected variables:
lu_3yrs_sp_sel <- lu_3yrs_sp[[selvar_final[!grepl(pattern = "bio|pr_|mean", x = selvar_final)]]]
# extract value of each lu variable at each route centroid:
for(luvar in names(lu_3yrs_sp_sel)){
  routes_sel_sf[, paste0(luvar, "_3yrs")] <- lu_3yrs_sp_sel[[luvar]] %>% 
    terra::extract(y = routes_sel_sf) %>% 
    pull(luvar)
}

# match to yearly env. data:
route_sel_env_dt_ISIMIP <- route_sel_env_dt1 %>%
  left_join(routes_sel_sf, by = "RTENO_BBS") %>%
  dplyr::select(-geometry) %>% 
  arrange(RTENO_BBS)


# write to file:
if(env_path  == file.path("data", "Env_data")){
  save(route_sel_env_dt_ISIMIP, file = file.path("data", "route_sel_env_dt_ISIMIP_obsclim.RData"))
  write.csv(route_sel_env_dt_ISIMIP, file = file.path("data", "route_sel_env_dt_ISIMIP_obsclim.csv"),
            row.names = FALSE)
} else {
  save(route_sel_env_dt_ISIMIP, file = file.path("data", "route_sel_env_dt_ISIMIP_counterclim.RData"))
  write.csv(route_sel_env_dt_ISIMIP, file = file.path("data", "route_sel_env_dt_ISIMIP_counterclim.csv"),
            row.names = FALSE)
}



# 1901soc - extract land use data: ----


# store values of year i in sf:
routes_sel_year_sf <- routes_sel_sf %>%
  mutate(Year = 1901)

# land use variables of year i:
lu_year <- rast(lu_files[which(grepl(paste0("[a-z]_1901_ESRI102003.tif$"), lu_files))])
# reduce to selected variables:
lu_year_sel <- lu_year[[selvar_final[!grepl(pattern = "bio|pr_|mean", x = selvar_final)]]] # xx
# extract values of each land use class at each relevant route location:
for(luvar in names(lu_year_sel)){
  routes_sel_year_sf[, luvar] <- lu_year_sel[[luvar]] %>% 
    terra::extract(y = routes_sel_year_sf) %>% 
    pull(luvar)
}
routes_sel_year_sf 
# as data frame:
route_sel_env_dt1901 <- routes_sel_year_sf %>% 
  st_drop_geometry()

# match to obsclim data:
load(file = file.path("data", "route_sel_env_dt_ISIMIP_obsclim.RData"))
route_sel_env_dt_ISIMIP
# land use variables:
lu_vars <- selvar_final[!grepl(pattern = "bio|pr_|mean", x = selvar_final)]

route_sel_env_dt_ISIMIP_obsclim_1901soc <- route_sel_env_dt_ISIMIP %>% 
  # replace histsoc with 1901soc:
  select(-matches(lu_vars)) %>% 
  left_join(route_sel_env_dt1901, by = "RTENO_BBS") %>% 
  mutate(urbanareas_3yrs = urbanareas,
         managed_pastures_3yrs = managed_pastures,
         primary_nonforests_3yrs = primary_nonforests,
         secondary_nonforests_3yrs = secondary_nonforests,
         sum_annual_crops_3yrs = sum_annual_crops) %>% 
  select(-Year.y) %>% 
  rename(Year = Year.x)
colnames(route_sel_env_dt_ISIMIP_obsclim_1901soc)

# save:
save(route_sel_env_dt_ISIMIP_obsclim_1901soc, file = file.path("data", "route_sel_env_dt_ISIMIP_obsclim_1901soc.RData"))
write.csv(route_sel_env_dt_ISIMIP, file = file.path("data", "route_sel_env_dt_ISIMIP_obsclim_1901soc.csv"),
          row.names = FALSE)

# some explorations: ----

library(ggplot2)

route_sel_env_dt_ISIMIP %>% 
  filter(RTENO_BBS == "84092077") %>% 
  ggplot(aes(x = Year, y = urbanareas)) +
  geom_line() +
  geom_vline(xintercept = 1995, linetype = "dashed") +
  theme_bw()
colnames(route_sel_env_dt_ISIMIP)
unique(route_sel_env_dt_ISIMIP$RTENO_BBS)
