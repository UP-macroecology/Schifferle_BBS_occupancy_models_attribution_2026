
# assigns BBS routes to 5 folds for spatially blocked cross validation:
# + plots and tests

# packages: ----

library(blockCV)
library(sf) # working with spatial vector data
library(terra) # working with spatial raster data
library(tmap) # plotting maps
library(dplyr)
library(ggplot2)

# can response be binary or continuous? (it said so somewhere)

# repeat this after having decided which buffer to use regarding what routes to include when fitting the model!

# load data: ----

# x - y - abs/pres:
load(file = file.path("data", "BBS_for_occ_spec_records.RData")) # output of 1_0_reformat_BBS_data.R

# BBS route selection (centroids) to fit models:
routes_sel_sf <- st_read(file.path("data", "route_selection_1991_2015_surv_beg_end_max_5y_miss_v2_spat_thin_100km_max_30_r_per_BCR_centroids.shp")) # output of 1_1_route_selection.R

# merged route, year, environment data:
load(file = file.path("data", "route_year_env_data.RData"))

# load env. data as raster stack to check environmental similarity between training and test folds:

# files:
bioclim_files <- list.files(file.path("data", "Env_data", "ISIMIP_CHELSA-W5E5v1.0", "bioclim"), full.names = TRUE)
lu_files <- list.files(file.path("data", "Env_data", "LUH2", "albers_proj"), full.names = TRUE)
sclim_files <- list.files(file.path("data", "Env_data", "ISIMIP_CHELSA-W5E5v1.0", "seasonal"), full.names = TRUE) #xx

bioclim_3yrs_sp <- rast(bioclim_files[which(grepl("1988_1991", bioclim_files))])
sclim_3yrs_sp <- rast(sclim_files[which(grepl("1988_1991", sclim_files))])
lu_3yrs_sp <- rast(lu_files[which(grepl("1988_1990", lu_files))])

# selected variables:
selvar <- c("bio1", "bio2", "bio3", "bio7", "bio14", "bio15", "pr_spring", "pr_summer", "pr_autumn", "pr_winter",
            "sum_annual_crops", "secdf", "pastr", "urban")

# reduce to selected variables:
bioclim_3yrs_sp_sel <- bioclim_3yrs_sp[[selvar[grepl(pattern = "bio", x = selvar)]]]
sclim_3yrs_sp_sel <- sclim_3yrs_sp[[selvar[grepl(pattern = "(spring|summer|autumn|winter)", x = selvar)]]]
lu_3yrs_sp_sel <- lu_3yrs_sp[[selvar[!grepl(pattern = "bio|pr_", x = selvar)]]] # xx

# combine:
env_rasters <- c(bioclim_3yrs_sp_sel, sclim_3yrs_sp_sel, lu_3yrs_sp_sel)


# assign routes to folds: ----

hexagon_size_m <- 500000

# # what is a suitable block size?
# 
# # for first choice look at existing autocorrelation in response or predictors:
# # = range over which observations are independent
# response data:
# sac2 <- cv_spatial_autocor(x = occ_spec_sf,
#                        column =  "presence_summarised")
# # 2 blocks


# since I use a buffer around the presence points of each species to determine the
# data used for each model, I also generate folds for each species separately:

# selected species:
load(file = file.path("data", "final_species_selection.RData")) # output of 1_2_species_selection.R

# iterate over species:
for(spec in species_selection_final){
  
  print(spec)
  
  # species data:
  
  ## presences:
  presences_spec <- bbs_dt_occ %>% 
    select(c(English_Common_Name, RTENO, Year, paste0("Count", seq(10, 50, 10)))) %>% 
    filter(English_Common_Name == spec)
  
  # match to routes-year-env to get presence-absence data:
  occ_dt_spec <- route_sel_env_dt_final %>% 
    # add observations:
    collapse::join(presences_spec, on = c("RTENO", "Year"), how = "left") %>% 
    # if route was surveyed but species not observed, replace NA with 0:
    mutate(across(Count10:Count50, ~ 
                    case_when(Surveyed == 1 & is.na(.) ~ 0,
                              .default = .))) %>%
    # convert bird counts to presence / absence:
    mutate(across(Count10:Count50, ~ 
                    case_when(. > 1 ~ 1,
                              .default = .)))
  
  # match with routes, summarise presences over time:
  occ_spec_sf <- routes_sel_sf %>%
    left_join(occ_dt_spec, by = c("RTENO_BBS" = "RTENO")) %>%
    # presence on route across all sections:
    mutate(presence = rowSums(across(paste0("Count", seq(10, 50, 10))))) %>%
    mutate(presence = ifelse(presence >= 1, 1, 0)) %>% 
    # presence on route across all years:
    group_by(RTENO_BBS) %>%
    summarise(presence_summarised = max(presence, na.rm=TRUE)) %>%
    mutate(presence_summarised = factor(presence_summarised, levels = c(1,0)))
  
  # buffer presences:
  pres_buffer <- occ_spec_sf %>% 
    filter(presence_summarised == 1) %>%
    st_buffer(dist = 750000) %>% # 750000
    st_union
  
  # routes within buffer:
  occ_spec_sf_buffered <- occ_spec_sf %>% 
    st_filter(., y = pres_buffer, join = st_within)
  
  # create blocks and assign routes to folds:
  sb_US <- cv_spatial(x = occ_spec_sf_buffered,
                      column = "presence_summarised", # the response column (binary or multi-class)
                      k = 5, # number of folds
                      size = hexagon_size_m, # size of the blocks in metres
                      selection = "random", # random blocks-to-fold; more even distribution of presence/absence instances between the train and test folds compared to ‘systematic’
                      iteration = 100, # find evenly dispersed folds; 50 no enough, 100 = default, 200 = no improvement for critical species
                      biomod2 = FALSE, # also create folds for biomod2
                      flat_top = TRUE,
                      seed = 5254,# 3456789: 39
                      )

  # write plot to file:
  jpeg(file = file.path("plots", "blockCV_folds", paste0("block_size_", hexagon_size_m/1000, "km"), paste0(spec, ".jpg")), 
       width = 1000, height = 800, quality = 100)
  print(cv_plot(cv = sb_US,
          x = occ_spec_sf_buffered, # sample points
          nrow = 2,
          points_alpha = 0.5) + 
    ggtitle(paste(spec, hexagon_size_m/1000, "km")) +
    theme(text = element_text(size = 18)))
  dev.off()
  
  # save route to block allocation:
  sb_US$blocks
  sb_US$records # number of presences in training data - number of absences in training data - pres. in test data - abs. in test data
  sb_US$folds_ids # which route in which fold
  sb_US$folds_list # for each fold, which routes are in training set, which routes are in test set
  
  save(sb_US, file = file.path("data", "CV_route_block_allocation", paste0("block_size_", hexagon_size_m/1000, "km"), paste0(spec, ".RData")))
  
  # environmental clustering
  
  # similarity measures to evaluate possible extrapolation in testing folds:
  # MESS = multivariate environmental similarity surface
  # how similar a point in a testing fold is to a training fold
  # negative values are sites where at least one variable of env. data has a value 
  # that is outside the range of environments over the reference set
  env_sim_US <- cv_similarity(cv = sb_US, # the environmental clustering
                              x = occ_spec_sf_buffered, 
                              r = env_rasters, 
                              progress = TRUE)
  
  # write plot to file:
  jpeg(file = file.path("plots", "blockCV_folds", paste0("block_size_", hexagon_size_m/1000, "km"), paste0(spec, "_env_similarity.jpg")), 
       width = 800, height = 600, quality = 100)
  print(env_sim_US +
    ggtitle(paste(spec, hexagon_size_m/1000, "km; neg. = novel env.")) +
    theme(text = element_text(size = 18)))
  dev.off()
}



# check blockCV results: ----

## number of presences and absences in training and test data ----

# write for each species number of presences and absences in training and test data of each fold to file:
sink(file.path("data", paste0("blockCV_", hexagon_size_m/1000, "km_pres_abs_per_fold.txt"))) # write console output here
for(spec in species_selection_final){
  load(file.path("data", "CV_route_block_allocation", paste0("block_size_", hexagon_size_m/1000, "km"), paste0(spec, ".RData")))
  print(spec)
  print(sb_US$records)
}
sink(file = NULL)

#st_write(sb_US$blocks, file.path("plots", "test_hex_blocks.shp"))

## presences for every year in the training data for each fold: ----

# do we have presences in every year of the training data for each fold (to be able to fit the model):

sink(file.path("data", paste0("blockCV_",  hexagon_size_m/1000, "km_fold_pres_per_year.txt"))) # write console output here

for(spec in species_selection_final){

  print(spec)
  
  n_species_without_enough_presences <- 0
  
  ## presences:
  presences_spec <- bbs_dt_occ %>% 
    select(c(English_Common_Name, RTENO, Year, paste0("Count", seq(10, 50, 10)))) %>% 
    filter(English_Common_Name == spec)
  
  # match to routes-year-env to get presence-absence data:
  occ_dt_spec <- route_sel_env_dt_final %>% 
    # add observations:
    collapse::join(presences_spec, on = c("RTENO", "Year"), how = "left", verbose = 0) %>% 
    # if route was surveyed but species not observed, replace NA with 0:
    mutate(across(Count10:Count50, ~ 
                    case_when(Surveyed == 1 & is.na(.) ~ 0,
                              .default = .))) %>%
    # convert bird counts to presence / absence:
    mutate(across(Count10:Count50, ~ 
                    case_when(. > 1 ~ 1,
                              .default = .))) %>% 
    # presence on route across all sections:
    mutate(presence = rowSums(across(paste0("Count", seq(10, 50, 10))))) %>%
    mutate(presence = ifelse(presence >= 1, 1, 0)) # save that somewhere!
  
  # match with routes, summarise presences over time:
  occ_spec_sf <- routes_sel_sf %>%
    left_join(occ_dt_spec, by = c("RTENO_BBS" = "RTENO")) %>%
    # presence on route across all sections:
    mutate(presence = rowSums(across(paste0("Count", seq(10, 50, 10))))) %>%
    mutate(presence = ifelse(presence >= 1, 1, 0)) %>% 
    # presence on route across all years:
    group_by(RTENO_BBS) %>%
    summarise(presence_summarised = max(presence, na.rm=TRUE)) %>%
    mutate(presence_summarised = factor(presence_summarised, levels = c(1,0)))
  
  # buffer presences:
  pres_buffer <- occ_spec_sf %>% 
    filter(presence_summarised == 1) %>%
    st_buffer(dist = 750000) %>% # 750000
    st_union
  
  # routes within buffer:
  occ_spec_sf_buffered <- occ_spec_sf %>% 
    st_filter(., y = pres_buffer, join = st_within)
  
  # data within buffer:
  occ_dt_spec_buff <- occ_dt_spec %>% 
    filter(RTENO %in% occ_spec_sf_buffered$RTENO_BBS)
  
  # assigned routes to folds:
  load(file.path("data", "CV_route_block_allocation", "block_size_500km", paste0(spec, ".RData")))
  
  for(f in 1:5){

    print(paste("fold", f))
    
    RTENO_training_fold <- occ_spec_sf_buffered$RTENO_BBS[which(sb_US$folds_ids != f)]
    
    presences_routes_year <- occ_dt_spec_buff %>% 
      filter(RTENO %in% RTENO_training_fold) %>% 
      filter(presence == 1) %>% 
      group_by(Year) %>%
      summarise(n_pres = n()) %>%
      print(n = 25) %>% 
      filter(n_pres == 0)
    
    if(nrow(presences_routes_year) != 0){
      print(paste("No routes with presences for", spec, "in fold", f))
      n_species_without_enough_presences <- n_species_without_enough_presences + 1
      }
  }

}
n_species_without_enough_presences
sink(file = NULL)
