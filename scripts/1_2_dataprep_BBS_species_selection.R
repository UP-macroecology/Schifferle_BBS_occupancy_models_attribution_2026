# Script:   1_2_dataprep_BBS_species_selection.R
# Purpose:  Select BBS species for which to fit dynamic occupancy models
# Inputs:   data/BBS_data_merged.RData
#           data/BBS_for_occ_spec_records.RData
#           data/route_selection_1995_2019_surv_beg_end_max_5y_miss_v2_spat_thin_100km_max_30_r_per_BCR.RData
# Outputs:  data/final_species_selection.RData
# Runs on:  Local

# Steps:
# of all available species in BBS data, we excluded species which:
# 1) are nocturnal or water-related (not well captured by BBS method)
# 2) which are rare or very common (model fitting difficult)


source(file.path("scripts", "0_paths.R"))

# packages: --------------------------------------------------------------------

library(dplyr)
library(sf)
library(ggplot2)


# functions: -------------------------------------------------------------------

source(file.path("scripts", "0_functions.R"))


# load data: -------------------------------------------------------------------

# BBS cleaned, all species info:

load(file = file.path(dir, "data", "BBS_data_merged.RData")) # bbs_dt; output of 1_0_dataprep_BBS_bird_data.R
BBS_species_list <- bbs_dt  %>% 
  select(English_Common_Name, Scientific_Name, ORDER, Family) %>%
  distinct

# BBS data cleaned:
load(file = file.path(dir, "data", "BBS_for_occ_spec_records.RData")) # output of 1_0_dataprep_BBS_bird_data.R
bbs_dt_occ

# selected routes:
load(file = file.path(dir, "data", "route_selection_1995_2019_surv_beg_end_max_5y_miss_v2_spat_thin_100km_max_30_r_per_BCR.RData")) # output of 1_1_dataprep_BBS_route_selection.R
sel_routes_final


# species selection: -----------------------------------------------------------

# BBS data with only selected routes and focal time period and species info:

# time period selected in 1_1_dataprep_BBS_route_selection.R
start <- 1995
end <- 2019

bbs_dt_occ_sel <- bbs_dt_occ %>% 
  filter(RTENO %in% sel_routes_final) %>% 
  filter(Year >= start & Year <= end) %>% 
  left_join(BBS_species_list[, c("English_Common_Name", "Scientific_Name", "ORDER", "Family")], by = c("English_Common_Name"))
nrow(bbs_dt_occ_sel)

length(unique(bbs_dt_occ_sel$English_Common_Name)) # 506 species in total


## 1) exclude species based on order / AOU: ----

# exclude species not well sampled by BBS methods:

# exclude water-related birds (waterbirds & shorebirds: AOU <=2880, kingfishers: AOU >=3900 & <=3910 , dipper: AOU 7010),
# & nocturnal birds (owls: AOU >=3650 & <=3810 and nightjars: AOU >=4160 & <=4210): 

# (see https://ecologicaldata.org/wiki/breeding-bird-survey-north-america and Harris et al. 2018)

excl_orders <- bbs_dt_occ_sel %>% 
  filter(AOU <= 2880 | (AOU >=3900 & AOU <=3910)  | AOU == 7010  | (AOU >= 3650 & AOU <= 3810)| (AOU >=4160 & AOU <= 4210)) %>%
  pull(English_Common_Name) %>% 
  unique
excl_orders # 133 (112 water-related species, 21 nocturnal species)

# # how many species of which order excluded:
# bbs_dt_occ_sel %>% 
#   filter(English_Common_Name %in% excl_orders) %>% 
#   select(c(ORDER, English_Common_Name)) %>% 
#   distinct %>%
#   group_by(ORDER) %>% 
#   count()
# 
# # how many species of which order included:
# bbs_dt_occ_sel %>% 
#   filter(!English_Common_Name %in% excl_orders) %>% 
#   select(c(ORDER, English_Common_Name)) %>% 
#   distinct %>%
#   group_by(ORDER) %>% 
#   count()


## 2) exclude species based on data availability: -----

# at how many different routes is each species detected from 1995-2019:

n_routes_pres <- bbs_dt_occ_sel %>% # presences only
  #filter(!English_Common_Name %in% excl_orders) %>% 
  select(English_Common_Name, Family, Scientific_Name, ORDER, RTENO) %>%
  distinct %>% # one row per species ever recorded on a route
  group_by(English_Common_Name,  Family, Scientific_Name, ORDER) %>% 
  summarise(n_routes = n())


# # explorations:
# 
# # number of presence across all focal years for each species:
# spec_N_total <- bbs_dt_occ_sel %>% 
#   group_by(English_Common_Name, ORDER, Scientific_Name, Family) %>% 
#   summarise(spec_N = n())
# spec_N_total
# 
# # for how many species do we have a certain number of presences in the data from 1995-2019:
# N_species_presences <- spec_N_total %>% 
#   rename("N_presences" = spec_N) %>% 
#   group_by(N_presences) %>% 
#   summarise(N_species = n()) %>% 
#   ungroup() %>% 
#   mutate(spec_with_equal_or_less = cumsum(N_species)) %>% 
#   mutate(spec_with_more = length(unique(bbs_dt_occ_sel$English_Common_Name)) - spec_with_equal_or_less)
# N_species_presences


# exclude very rare and very common species:
# exclude species that are detected at less than 50 different routes across the whole time period
# and species that have less than 50 routes where they were never detected:

n_routes_pres %>% 
  filter(n_routes < 50) %>% 
  nrow
n_routes_pres %>% 
  filter(n_routes > (length(sel_routes_final)-50)) %>% 
  nrow

excl_data_av <- n_routes_pres %>% 
  filter(n_routes < 50 | n_routes > (length(sel_routes_final)-50)) %>%
  pull(English_Common_Name) %>% 
  sort
excl_data_av


## final species selection: ----

species_selection_final <- bbs_dt_occ_sel %>% 
  filter(!English_Common_Name %in% excl_orders) %>%
  filter(!English_Common_Name %in% excl_data_av) %>% 
  pull(English_Common_Name) %>% 
  unique %>% 
  sort

length(species_selection_final) # 192

# save:
save(species_selection_final, file = file.path(dir, "data", "final_species_selection.RData"))

# session info:
writeLines(capture.output(sessionInfo()), file.path(dir, "results", "sessionInfo", "1_2_dataprep_BBS_species_selection.txt"))
