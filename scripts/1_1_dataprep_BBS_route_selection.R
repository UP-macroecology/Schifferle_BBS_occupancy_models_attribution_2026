# select BBS routes used to fit dynamic occupancy models:

# routes selected based on:
# 1) time period, number of years a route was surveyed (test different versions)
# 2) spatial thinning 1: minimum distance 100 km
# 3) spatial thinning 2: maximum 30 routes per Bird Conservation Region (BCR)


# packages: --------------------------------------------------------------------

library(dplyr)
library(tidyr)
library(collapse)
library(sf)
library(ggplot2)


# functions: -------------------------------------------------------------------

source("0_functions.R") # to get thin() and thin.max()


# load data: -------------------------------------------------------------------

# BBS data formatted for occupancy modelling:
load(file = file.path("data", "BBS_for_occ.RData")) # output of 1_0_dataprep_BBS_bird_data.R
route_dt

# spatial data on routes (centroids):
routes_sf <- read_sf(file.path("data", "BBS_routes_all_centroids.shp")) # output of 1_0_dataprep_BBS_bird_data.R


# route selection: -------------------------------------------------------------

# first remove routes too rarely sampled

start_year <- 1995 # tested different start and end years
end_year <- 2019


## 1) temporal coverage: ----

# routes surveyed in considered time period:
route_dt_time <- route_dt %>% 
  # time period to consider:
  filter(Year >= start_year & Year <= end_year) # no data in 2020 due to Covid-19

nyears <- length(seq(min(route_dt_time$Year), max(route_dt_time$Year), 1))

# test different requirements for route selection:

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
  
  # route_dt: df with columns RTENO (route ID), Year, Surveyed (1: route surveyed in this year, 0: route not surveyed)
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
                         version = "surveyed beginning & end, max. 5 years missing")


# selected routes as sf:
routes_sel_sf <- routes_sf %>% 
  filter(RTENO_BBS %in% routes_sel) %>% 
  group_by(RTENO_BBS) %>%  # merge lines if routes in shapefile consist of multiple adjacent lines
  summarise

# save this intermediate route selection step:
st_write(routes_sel_sf, file.path("data", "route_selection_1995_2019_surv_beg_end_max_5y_miss.shp"), append = FALSE) 


## 2) spatial thinning: ----

# require minimum distance of route centroids of 100 km

routes_sf_thinned <- thin(sf = routes_sel_sf, thin_dist = 50000, runs = 1, ncores = 1) # 100 km apart -> 50 km thinning distance (radius)
nrow(routes_sel_sf)
nrow(routes_sf_thinned) # 1991-2015: 533, 1995-2019: 632, 1991 - 2019: 491

# save intermediate route selection step:
st_write(routes_sf_thinned, file.path("data", "route_selection_1995_2019_surv_beg_end_max_5y_miss_v2_spat_thin_100km_centroids.shp"), append = FALSE) 


## 3) spatial coverage: ----
# remove routes to get approximately equal coverage of Bird Conservation Regions (BCRs):

# spatial subsampling as Jarzyna et al. 2017 & 2018:
# "We removed routes from BCRs with more than 30 routes (in order of
# proximity to remaining routes) until all BCRs had only 30 or fewer routes."

# routes left after checking temporal coverage and after spatial thinning:
route_dt_temp_ss <- route_dt %>% 
  filter(RTENO %in% routes_sf_thinned$RTENO_BBS)
  
# number of routes per BCR:
routes_per_BCR_dt <- route_dt_temp_ss %>% 
  select(RTENO, BCR, Latitude, Longitude) %>% 
  distinct() %>% 
  # convert to sf:
  left_join(routes_sel_sf[, c("RTENO_BBS")], by = c("RTENO" = "RTENO_BBS")) %>% 
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

# iterate over BCRs:

BCRs_to_subsample <- subsample_dt %>% 
  pull(BCR) %>% 
  unique

for (b in 1:length(BCRs_to_subsample)){
  
  subsample_dt_BCR <- subsample_dt %>% 
    filter(BCR == BCRs_to_subsample[b]) %>% 
    st_drop_geometry()

  subsample_BCR <- thin.max(x = subsample_dt_BCR, cols = c("Latitude", "Longitude"), npoints = 30)
  routes_keep_BCR <- subsample_BCR %>% pull(RTENO)
  
  # store:
  if(b == 1){
    routes_keep2 <- routes_keep_BCR
  } else {
    routes_keep2 <- c(routes_keep2, routes_keep_BCR)
  }
  
  # plots:
  subsample_dt_BCR %>%
    left_join(routes_sf_thinned[, "RTENO_BBS"], by = c("RTENO" = "RTENO_BBS")) %>%
    st_as_sf() %>%
    st_geometry() %>%
    plot(main = paste(BCRs_to_subsample[b]))

  subsample_dt_BCR %>%
    left_join(routes_sf_thinned[, "RTENO_BBS"], by = c("RTENO" = "RTENO_BBS")) %>%
    filter(RTENO %in% routes_keep_BCR) %>%
    st_as_sf() %>%
    st_geometry() %>%
    plot(add = TRUE, col = "red", pch = 20)
}

# subset data with selected routes:
route_dt_spat_temp_ss <- route_dt_temp_ss %>% 
  filter(RTENO %in% c(routes_keep1, routes_keep2))
length(unique(route_dt_spat_temp_ss$RTENO)) # 1995-2019: 539


# save output:
sel_routes_final <- route_dt_spat_temp_ss %>% 
  pull(RTENO) %>% 
  unique %>% 
  sort
save(sel_routes_final, 
     file = file.path("data", "route_selection_1995_2019_surv_beg_end_max_5y_miss_v2_spat_thin_100km_max_30_r_per_BCR.RData"))


# save selected routes as shapefile - centroids:
routes_sel_sf %>% 
  filter(RTENO_BBS %in% sel_routes_final) %>% 
  st_write(file.path("data", "route_selection_1995_2019_surv_beg_end_max_5y_miss_v2_spat_thin_100km_max_30_r_per_BCR_centroids.shp"),
           append = FALSE)


# map of selected routes: ------------------------------------------------------

# load data:

## selected routes:
routes_sel_sf <- read_sf(file.path("data", "route_selection_1995_2019_surv_beg_end_max_5y_miss_v2_spat_thin_100km_max_30_r_per_BCR_centroids.shp"))

## conterminous USA:
US_albers_sf <- read_sf(file.path("data", "US_outline_ESRI102003.shp")) # output of 1_0_dataprep_climate.R

## Bird Conservation Regions, clipped to conterminous US:
bcr_sf <- read_sf(file.path("data", "BCR_Terrestrial", "BCR_Terrestrial_master_International.shp")) %>% # from: https://www.birdscanada.org/bird-science/nabci-bird-conservation-regions
  st_transform(crs = "ESRI:102003") %>% 
  st_intersection(US_albers_sf) %>% 
  # add centroid coordinate for label placement:
  mutate(X = st_coordinates(st_centroid(.))[,1],
         Y = st_coordinates(st_centroid(.))[,2],
         Label2 = gsub(pattern = " ", replacement = "\n", x = Label))

# plot map:

route_map <- ggplot() +
  # BCRs as background map:
  geom_sf(data = bcr_sf, aes(fill = Label), show.legend = TRUE, alpha = 0.5) +
  geom_sf(data = routes_sel_sf, aes(colour = "black"), size = 2) +
  scale_fill_manual(values = pals::glasbey()) +
  scale_colour_identity(name = NULL, labels = c(black = "selected BBS route"), guide = "legend") +
  theme_light() +
  theme(legend.position = "bottom",
        legend.box = "vertical",
        text = element_text(size = 20),
        axis.title.x=element_blank(),
        axis.title.y=element_blank()) +
  ggspatial::annotation_scale(location = "bl") +
  ggspatial::annotation_north_arrow(location = "bl", 
                                    pad_x = unit(0.9, "in"), pad_y = unit(0.3, "in"),
                                    style = ggspatial::north_arrow_fancy_orienteering) +
  guides(fill = guide_legend(title = "Bird Conservation Region", position = "bottom", 
                        theme = theme(legend.text = element_text(face = "italic"),
                                      legend.title.position = "top"),
                        ncol = 3),
         colour = guide_legend(position = "bottom")
  )
route_map

# ggsave(filename = file.path("plots", "route_selection_BCRs_no_labels.svg"), 
#        plot = route_map,
#        device = "svg",
#        width = 32,
#        height = 29.7, # A4
#        units = "cm")
