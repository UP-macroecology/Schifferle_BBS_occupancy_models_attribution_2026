# reformat BBS data for occupancy modelling
# one row per route-year combination, routes either surveyed or not

library(dplyr)

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
