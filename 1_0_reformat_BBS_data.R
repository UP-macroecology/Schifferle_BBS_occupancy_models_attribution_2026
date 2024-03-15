# reformat BBS data for occupancy modelling
# one row per route-year combination, routes either surveyed or not
# centroids of routes as shapefile

# packages: ----

library(dplyr)
library(sf)

# reformat BBS data: -----------------------------------------------------------

load(file = file.path("data", "BBS_data_merged.RData")) # output of DEBTs\analysis\Schifferle_BBS_explorations_2023\BBS_data_prep.R

# reduce dataset to necessary columns:
bbs_dt_occ <- bbs_dt %>% 
  select(English_Common_Name, AOU, RTENO, Latitude, Longitude, BCR, Year, paste0("Count", seq(10, 50, 10)), Month, Day, ObsN) %>% 
  # convert month and date to day of year:
  mutate(date = lubridate::ymd(paste(Year, Month, Day, sep = "/"))) %>% 
  mutate(doy = lubridate::yday(date)) %>% 
  select(-c(Month, Day, date)) %>% 
  # add column on whether route was surveyed (needed later):
  mutate(Surveyed = 1) %>% 
  # site needs to be numeric:
  mutate(RTENO = as.numeric(RTENO))

save(bbs_dt_occ, file = file.path("data", "BBS_for_occ_spec_records.RData"))

# expand data to have one row per route and year:
route_dt <- tidyr::expand_grid(RTENO = unique(bbs_dt_occ$RTENO),
                               Year = min(bbs_dt_occ$Year):max(bbs_dt_occ$Year)) %>% # 224'124
  # join route data:
  collapse::join(bbs_dt_occ[, c("RTENO", "Latitude", "Longitude", "BCR")], on = c("RTENO"), how = "left") %>% 
  
  # add observer and date when route was surveyed:
  collapse::join(bbs_dt_occ[, c("RTENO", "Year", "ObsN", "doy")], on = c("RTENO", "Year"), how = "left") %>%
  # all route-year combinations without date / observer haven't been surveyed:
  mutate(Surveyed = if_else(is.na(doy), 0, 1))

save(route_dt, file = file.path("data", "BBS_for_occ.RData"))


# save shapefile of route centroids: -------------------------------------------

datashare_BBS <- file.path("//ibb-fs01.ibb.uni-potsdam.de", "daten$", "AG26", "Arbeit", "datashare", "data", "biodat", "distribution", "BBS")

routes_sf <- read_sf(file.path(datashare_BBS, "bbs_routes", "bbsrtsl020.shp")) %>% 
  st_transform(crs = "ESRI:102003") %>% # Albers Equal Area projection
  mutate(RTENO_BBS = as.integer(paste0("840", stringr::str_pad(RTENO, width = 5, side = "left", pad = "0")))) # reformat RTENO to match RTENO from BBS data imported with the bbsAssistant package

# route centroids:

routes_sf2 <- routes_sf %>% 
  group_by(RTENO_BBS) %>% summarise %>% # merge lines if routes in shapefile consist of multiple adjacent lines
  mutate(centroid_X = NA) %>% 
  mutate(centroid_Y = NA)

# coordinates of route centroids, if route geometry is available:
for(j in 1:nrow(routes_sf2)){
  
  print(j)
  
  centroid_coords <- routes_sf2[j,] %>% 
    st_bbox %>% 
    st_as_sfc %>% 
    st_centroid %>% 
    st_coordinates
  
  routes_sf2$centroid_X[j] <- centroid_coords[1]
  routes_sf2$centroid_Y[j] <- centroid_coords[2]
}
# change geometry from route line to centroid:
routes_sf_centr <- routes_sf2 %>% 
  st_drop_geometry %>% 
  st_as_sf(coords = c("centroid_X", "centroid_Y"), crs = "ESRI:102003")

# save routes as shapefile - centroids:
routes_sf_centr %>% 
  st_write(file.path("data", "BBS_routes_all_centroids.shp"),
           append = FALSE)
