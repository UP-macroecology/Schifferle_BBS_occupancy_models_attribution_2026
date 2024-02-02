# select BBS routes to use in occupancy models:

# packages: ----

library(dplyr)
library(tidyr)
library(collapse)
library(sf)

# load data: ----

## BBS data:

# BBS data formatted for occupancy modelling:
load(file = file.path("data", "BBS_for_occ.RData")) # output of 1_0_reformat_BBS_data.R

# spatial data on routes:
datashare_BBS <- file.path("//ibb-fs01.ibb.uni-potsdam.de", "daten$", "AG26", "Arbeit", "datashare", "data", "biodat", "distribution", "BBS")
routes_sf <- read_sf(file.path(datashare_BBS, "bbs_routes", "bbsrtsl020.shp")) %>% 
  st_transform(crs = "ESRI:102003") %>% # Albers Equal Area projection
  mutate(RTENO_BBS = as.integer(paste0("840", stringr::str_pad(RTENO, width = 5, side = "left", pad = "0")))) # reformat RTENO to match RTENO from BBS data imported with the bbsAssistant package


# route selection: ----

# first remove routes too rarely sampled, then remove routes to get approximately equal coverage of Bird Conservation Regions (BCRs)

## 1) temporal coverage: ----

# routes surveyed in considered time period:
route_dt_time <- route_dt %>% 
  # time period to consider:
  filter(Year >= 1991 & Year <= 2015) # no data in 2020 due to covid

nyears <- length(seq(min(route_dt_time$Year), max(route_dt_time$Year), 1))

# test different versions:

# keep routes that are sampled:
subset_versions <- c("each year", 
                   "min. every 2nd y.", 
                   "min. every 3rd y.", 
                   "1 y. missing", 
                   "2 y. missing", 
                   "5 y. missing", 
                   "10 y. missing",
                   "surveyed beginning & end, max. 5 years missing",
                   "surveyed beginning & end, max. 10 years missing")

sel_routes_temp_coverage <- function(route_dt, nyears, version) {
  
  # route_dt_time: df with columns RTENO (route ID), Year, Surveyed (1: route surveyed in this year, 0: route not surveyed)
  # nyears: length of considered time period (years)
  # version: how often must route be surveyed to be used
  # output: IDs of routes to use
  
  # number of years routes are surveyed within the time period and max. number of years missed between consecutive surveys:
  routes_survey_inf <- route_dt %>% 
    select(RTENO, Year, Surveyed) %>%
    arrange(RTENO, Year) %>%
    # number of missing years between consecutive surveys:
    group_by(RTENO, grp = with(rle(Surveyed), rep(seq_along(lengths), times = lengths))) %>% 
    mutate(counter = seq_along(grp)) %>% 
    mutate(gap = if_else(Surveyed == 1, 0, counter)) %>%
    # maximum number of missing years between consecutive surveys and number of years route was surveyed in total:
    ungroup() %>% 
    group_by(RTENO) %>% 
    summarise(max_gap = max(gap), years_surveyed = sum(Surveyed))
  
  if(version == "each year"){
    routes_subset <- routes_survey_inf %>% 
      filter(years_surveyed == nyears) %>% 
      pull(RTENO)
  }
  
  if(version == "min. every 2nd y."){
    routes_subset <- routes_survey_inf %>% 
      filter(max_gap <= 1) %>% 
      pull(RTENO)
  }
  
  if(version == "min. every 3rd y."){
    routes_subset <- routes_survey_inf %>% 
      filter(max_gap <= 2) %>% 
      pull(RTENO)
  }
  
  if(version == "1 y. missing"){
    routes_subset <- routes_survey_inf %>% 
      filter(years_surveyed >= nyears-1) %>% 
      pull(RTENO)
  }
  
  if(version == "2 y. missing"){
    routes_subset <- routes_survey_inf %>% 
      filter(years_surveyed >= nyears-2) %>% 
      pull(RTENO)
  }
  
  if(version == "5 y. missing"){
    routes_subset <- routes_survey_inf %>% 
      filter(years_surveyed >= nyears-5) %>% 
      pull(RTENO)
  }
  
  if(version == "10 y. missing"){
    routes_subset <- routes_survey_inf %>% 
      filter(years_surveyed >= nyears-10) %>% 
      pull(RTENO)
  }
  
  if(version == "surveyed beginning & end, max. 5 years missing"){
    routes_subset <- route_dt %>% 
      select(RTENO, Year, Surveyed) %>% 
      group_by(RTENO) %>% 
      # only keep routes / groups that was surveyed in start year and end year:
      filter(any(Year == min(route_dt$Year) & Surveyed == 1)) %>% 
      filter(any(Year == max(route_dt$Year) & Surveyed == 1)) %>%
      # how many years are missing in between:
      summarise(years_surveyed = sum(Surveyed), years_total = n()) %>% 
      mutate(missing_years = years_total - years_surveyed) %>%
      filter(missing_years <= 5) %>% 
      pull(RTENO)
  }
  
  if(version == "surveyed beginning & end, max. 10 years missing"){
    routes_subset <- route_dt %>% 
      select(RTENO, Year, Surveyed) %>% 
      group_by(RTENO) %>% 
      # only keep routes / groups that was surveyed in start year and end year:
      filter(any(Year == min(route_dt$Year) & Surveyed == 1)) %>% 
      filter(any(Year == max(route_dt$Year) & Surveyed == 1)) %>%
      # how many years are missing in between:
      summarise(years_surveyed = sum(Surveyed), years_total = n()) %>% 
      mutate(missing_years = years_total - years_surveyed) %>%
      filter(missing_years <= 10) %>% 
      pull(RTENO)
  }
  return(routes_subset)
}

routes_sel <- sel_routes_temp_coverage(route_dt = route_dt_time, 
                         nyears = nyears, 
                         version = "surveyed beginning & end, max. 5 years missing") # 1990-2015: 781

# save selected routes as shapefile - lines:
routes_sf %>% 
  filter(RTENO_BBS %in% routes_sel) %>% 
  st_write(file.path("data", "route_selection_1991_2015_surv_beg_end_max_5y_miss.shp"),
           append = FALSE)


## 2) spatial thinning: ----

# require minimum distance of route centroids of 100 km
# equal area projection of route centroids

# route centroids:

routes_sf_sel <- routes_sf %>% 
  filter(RTENO_BBS %in% routes_sel) %>% 
  group_by(RTENO_BBS) %>% summarise %>% # merge lines if routes in shapefile consist of multiple adjacent lines
  mutate(centroid_X = NA) %>% 
  mutate(centroid_Y = NA)

# coordinates of route centroids, if route geometry is available:
for(j in 1:nrow(routes_sf_sel)){
  
  print(j)
  
  centroid_coords <- routes_sf_sel[j,] %>% 
    st_bbox %>% 
    st_as_sfc %>% 
    st_centroid %>% 
    st_coordinates
  
  routes_sf_sel$centroid_X[j] <- centroid_coords[1]
  routes_sf_sel$centroid_Y[j] <- centroid_coords[2]
}
# change geometry from route line to centroid:
routes_sf_centr <- routes_sf_sel %>% 
  st_drop_geometry %>% 
  st_as_sf(coords = c("centroid_X", "centroid_Y"), crs = "ESRI:102003")

# load thinning function:
source("0_functions.R") # to get thin()

routes_sf_centr_thinned <- thin(sf = routes_sf_centr, thin_dist = 50000, runs = 1, ncores = 1) # 100 km apart -> 50 km thinning distance (radius)
nrow(routes_sf_centr)
nrow(routes_sf_centr_thinned) # 1991-2015: 533, 1995-2019: 632

# save selected routes as shapefile - lines:
routes_sf %>% 
  filter(RTENO_BBS %in% routes_sf_centr_thinned$RTENO_BBS) %>% 
  st_write(file.path("data", "route_selection_1991_2015_surv_beg_end_max_5y_miss_v2_spat_thin_100km.shp"),
           append = FALSE)

# save selected routes as shapefile:
routes_sf_centr_thinned %>% 
  st_write(file.path("data", "route_selection_1991_2015_surv_beg_end_max_5y_miss_v2_spat_thin_100km_centroids.shp"),
           append = FALSE)


## 3) spatial coverage: ----

# spatial subsampling as Jarzyna et al. 2017 & 2018:
# "We removed routes from BCRs with more than 30 routes (in order of
# proximity to remaining routes) until all BCRs had only 30 or fewer routes."

# I test two versions:

# 1) first remove one of the two closest points, 
# then calculate distances again and again remove one of the two closest points...

subsample1 <- function(data, npoints_keep = 30){
  
  # data = sf object
  
  # number of points to remove:
  remove_n <- nrow(data) - npoints_keep
  
  output_data <- data
  
  for(r in 1:remove_n){
    
    #print(r)
    
    # distance matrix:
    dist_mat <- output_data %>% 
      st_distance()
    
    # which pair of points are closest:
    dist_mat[dist_mat==units::as_units(0, "m")] <- NA  # exclude 0 (distance to itself)
    ind <- arrayInd(which.min(dist_mat), dim(dist_mat))
    
    # corresponding RTENO:
    remove_this <- output_data %>% 
      slice(ind[sample(c(1,2), 1)]) %>% # arbitrarily remove one of two closest points
      pull(RTENO)
    
    # remove that route:
    output_data <- output_data %>% 
      filter(RTENO != remove_this)
  }
  return(output_data)
}

# 2) from: https://gist.github.com/danlwarren/271288d5bab45d2da549:

# randomly choose one point, calculate distance to all others, choose the most distant one,
# then find distance of all not chosen points to the chosen points and choose again the most distant point to all chosen points, repeat
# input:
# x, a data frame containing the columns to be used to calculate distances along with whatever other data you need
# cols, a vector of column names or indices to use for calculating distances
# npoints, the number of rarefied points to spit out
# e.g. subsample2(my.data, c("latitude", "longitude"), 200)

# advantage: random drawing -> repeat to get different sets -> what to use? xx

subsample2 <- function(x, cols, npoints){
  
  # create empty vector for output:
  inds <- vector(mode = "numeric")
  
  # create distance matrix:
  this.dist <- as.matrix(dist(x[, cols], upper=TRUE))
  
  # draw first index at random:
  inds <- c(inds, as.integer(runif(1, 1, length(this.dist[,1]))))
  
  # get second index from maximally distant point from first one (necessary for apply ftc)
  inds <- c(inds, which.max(this.dist[,inds]))
  
  while(length(inds) < npoints){
    # for each point, find its distance to the closest point that's already been selected:
    min.dists <- apply(this.dist[,inds], 1, min)
    
    # select the point that is furthest from everything already selected:
    this.ind <- which.max(min.dists)
    
    # get rid of ties, if they exist:
    if(length(this.ind) > 1){
      print("Breaking tie...")
      this.ind <- this.ind[1]
    }
    inds <- c(inds, this.ind)
  }
  return(x[inds,])
}

# apply both versions:

route_dt_temp_ss <- route_dt %>% 
  # keep routes left after checking temporal coverage and after spatial thinning:
  filter(RTENO %in% routes_sf_centr_thinned$RTENO_BBS)

# number of routes per BCR:
routes_per_BCR_dt <- route_dt_temp_ss %>% 
  select(RTENO, BCR, Latitude, Longitude) %>% 
  distinct() %>% 
  # convert to sf:
  left_join(routes_sf_centr[, c("RTENO_BBS")], by = c("RTENO" = "RTENO_BBS")) %>% 
  st_as_sf() %>% 
  # number of routes per BCR:
  group_by(BCR) %>% 
  mutate(n = n()) %>% 
  mutate(thin = if_else(n > 30, 1, 0))

# keep all routes with less than 30 routes per BCR:
routes_keep1 <- routes_per_BCR_dt %>% 
  filter(thin == 0) %>% 
  pull(RTENO)

# subsample routes in each BCR with > 30 routes:
subsample_dt <- routes_per_BCR_dt %>% 
  filter(thin == 1)

BCRs_to_subsample <- subsample_dt %>% 
  pull(BCR) %>% 
  unique

## version 1:

# iterate over BCRs:
for (b in 1:length(BCRs_to_subsample)){

  subsample_dt_BCR <- subsample_dt %>% 
    filter(BCR == BCRs_to_subsample[b])
  
  subsample_BCR <- subsample1(subsample_dt_BCR, 30)
  routes_keep_BCR <- subsample_BCR %>% pull(RTENO)

  # store:
  if(b == 1){
    routes_keep2 <- routes_keep_BCR
  } else {
    routes_keep2 <- c(routes_keep2, routes_keep_BCR)
  }
  plot(st_geometry(subsample_dt_BCR), main = BCRs_to_subsample[b])
  plot(st_geometry(subsample_BCR), add = TRUE, col = "red", pch = 20)
}

# subset data with selected routes:
route_dt_spat_temp_ss1 <- route_dt_temp_ss %>% 
  filter(RTENO %in% c(routes_keep1, routes_keep2))

length(unique(route_dt_spat_temp_ss1$RTENO)) # 1991-2015: 476 routes, 1995-2019: 539 routes


## version 2:

# iterate over BCRs:
for (b in 1:length(BCRs_to_subsample)){
  
  subsample_dt_BCR <- subsample_dt %>% 
    filter(BCR == BCRs_to_subsample[b]) %>% 
    st_drop_geometry()

  subsample_BCR <- subsample2(x = subsample_dt_BCR, cols = c("Latitude", "Longitude"), npoints = 30)
  routes_keep_BCR <- subsample_BCR %>% pull(RTENO)
  
  # store:
  if(b == 1){
    routes_keep2 <- routes_keep_BCR
  } else {
    routes_keep2 <- c(routes_keep2, routes_keep_BCR)
  }
  
  # plots:
  subsample_dt_BCR %>%
    left_join(routes_sf_centr_thinned[, "RTENO_BBS"], by = c("RTENO" = "RTENO_BBS")) %>%
    st_as_sf() %>%
    st_geometry() %>%
    plot(main = paste("version 2:", BCRs_to_subsample[b]))

  subsample_dt_BCR %>%
    left_join(routes_sf_centr_thinned[, "RTENO_BBS"], by = c("RTENO" = "RTENO_BBS")) %>%
    filter(RTENO %in% routes_keep_BCR) %>%
    st_as_sf() %>%
    st_geometry() %>%
    plot(add = TRUE, col = "red", pch = 20)
}

# subset data with selected routes:
route_dt_spat_temp_ss2 <- route_dt_temp_ss %>% 
  filter(RTENO %in% c(routes_keep1, routes_keep2))
length(unique(route_dt_spat_temp_ss2$RTENO)) # 1995-2019: 539; 1991-2015: 476

# save output:
sel_routes_final <- route_dt_spat_temp_ss2 %>% 
  pull(RTENO) %>% 
  unique
save(sel_routes_final, 
     file = file.path("data", "route_selection_1991_2015_surv_beg_end_max_5y_miss_v2_spat_thin_100km_max_30_r_per_BCR.RData"))

# save selected routes as shapefile - lines:
routes_sf %>% 
  filter(RTENO_BBS %in% sel_routes_final) %>% 
  st_write(file.path("data", "route_selection_1991_2015_surv_beg_end_max_5y_miss_v2_spat_thin_100km_max_30_r_per_BCR.shp"),
           append = FALSE)
# save selected routes as shapefile - centroids:
routes_sf_centr %>% 
  filter(RTENO_BBS %in% sel_routes_final) %>% 
  st_write(file.path("data", "route_selection_1991_2015_surv_beg_end_max_5y_miss_v2_spat_thin_100km_max_30_r_per_BCR_centroids.shp"),
           append = FALSE)
