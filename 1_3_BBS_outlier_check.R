# check whether for the selected species presences are recorded outside the breeding range / 
# far from the majority of other presences


# packages: ----

library(CoordinateCleaner)
library(sf)
library(dplyr)
library(ggplot2)


# load data: ----

# selected species:
load(file = file.path("data", "final_species_selection.RData")) # species_selection_final; output of 1_2_species_selection.R

# BBS route selection (centroids) to fit models:
routes_sel_sf <- st_read(file.path("data", "route_selection_1995_2019_surv_beg_end_max_5y_miss_v2_spat_thin_100km_max_30_r_per_BCR_centroids.shp")) # output of 1_1_route_selection.R

# route-year-species information (only surveyed)
load(file = file.path("data", "BBS_for_occ_spec_records.RData")) # bbs_dt_occ; output of 1_0_reformat_BBS_data.R


# find outliers: ----

RTENO_outliers <- vector(mode = "list", length = length(species_selection_final))
names(RTENO_outliers) <- species_selection_final

for(i in 1:length(species_selection_final)){
  
  print(i)
  
  spec <- species_selection_final[i]
  
  spec_pres_sf <- bbs_dt_occ %>% # also on non-selected routes !?
    filter(RTENO %in% routes_sel_sf$RTENO_BBS) %>% 
    filter(Year >= 1995 & Year <= 2019) %>% 
    filter(English_Common_Name == spec) %>% 
    left_join(routes_sel_sf, by = c("RTENO" = "RTENO_BBS")) %>% # species counts
    mutate(presence = rowSums(across(paste0("Count", seq(10, 50, 10))))) %>%
    filter(presence >= 1) %>% 
    st_as_sf() %>%
    select(species = English_Common_Name, RTENO) %>% 
    distinct()  # one row per species and route 
  
  # reformat for CoordinateCleaner:
  spec_pres_dt <- spec_pres_sf %>% 
    st_transform(crs = "EPSG:4326") %>% # requires geographic coordinates
    mutate(lon = st_coordinates(.)[,1],
           lat = st_coordinates(.)[,2]) %>% 
    as_tibble()
  
  # cc_outl: Identify Geographic Outliers in Species Distributions
  outlier_flags <- cc_outl(x = spec_pres_dt, lon = "lon", lat = "lat",
                                method =  "distance", # "quantile", #"mad",#"quantile",#
                                value = "flagged",
                                tdi = 1000
                                )
  
  RTENO_outliers[[i]] <- spec_pres_dt$RTENO[!outlier_flags]
  
  print(
    spec_pres_sf %>% 
    cbind(outlier_flags) %>% 
    ggplot() +
      geom_sf(aes(colour = outlier_flags)) +
      ggtitle(spec)
    )
}

# affected species:
RTENO_outliers[which(lengths(RTENO_outliers) != 0)]

# Marsh Wren: RTENO: 84083311, but this is close to breeding range patches -> keep
# Northern Bobwhite: RTENO: 84089079, but this is close to breeding range patches -> keep
# Red Crossbill: RTENO: 84007015
# Yellow-throated Warbler: RTENO: 84069052

# adapted 0_functions.R -> to exclude these routes
