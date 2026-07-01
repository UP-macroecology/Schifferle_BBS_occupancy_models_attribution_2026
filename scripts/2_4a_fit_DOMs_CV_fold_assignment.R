# Script:   2_4a_fit_DOMs_CV_fold_assignment.R
# Purpose:  Assign BBS routes to 5 folds for spatially blocked cross validation of DOMs (species specific: routes with 750 km buffer around species presences)
# Inputs:   data/BBS_for_occ_spec_records.RData
#           data/route_selection_1995_2019_surv_beg_end_max_5y_miss_v2_spat_thin_100km_max_30_r_per_BCR_centroids.shp
#           data/route_year_env_data.RData
#           data/final_species_selection.RData
#           data/selected_variables.RData
#           data/BBS_for_occ_selection.RData
# Outputs:  data/CV_route_block_allocation/block_size_500km/<species>.RData (one per species)
#           plots/blockCV_folds/block_size_500km/<species>.jpg (one per species)
#           plots/blockCV_folds/block_size_500km/<species>_env_similarity.jpg (one per species)
#           data/blockCV_500km_fold_pres_per_year.txt 
#           data/blockCV_500km_pres_abs_per_fold.txt 
# Runs on:  Local


source(file.path("scripts", "0_paths.R"))


# packages: --------------------------------------------------------------------

library(blockCV)
library(sf)
library(terra)
library(dplyr)
library(ggplot2)


# directories: -----------------------------------------------------------------

hexagon_size_m <- 500000 # tested different values

# directory to save plots regarding fold assignment:
if(!dir.exists(file.path(dir, "plots", "blockCV_folds", paste0("block_size_", hexagon_size_m/1000, "km")))){
  dir.create(file.path(dir, "plots", "blockCV_folds", paste0("block_size_", hexagon_size_m/1000, "km")), recursive = TRUE)
}
# directory to save block assignment:
if(!dir.exists(file.path(dir, "data", "CV_route_block_allocation", paste0("block_size_", hexagon_size_m/1000, "km")))){
  dir.create(file.path(dir, "data", "CV_route_block_allocation", paste0("block_size_", hexagon_size_m/1000, "km")), recursive = TRUE)
}


# functions: -------------------------------------------------------------------

source(file.path("scripts", "0_functions.R"))


# load data: -------------------------------------------------------------------

# route-year-species information (only surveyed)
load(file = file.path(dir, "data", "BBS_for_occ_spec_records.RData")) # bbs_dt_occ; output of 1_0_dataprep_BBS_bird_data.R

# selected routes spatial data (to buffer presences):
routes_sel_sf <- st_read(file.path(dir, "data", "route_selection_1995_2019_surv_beg_end_max_5y_miss_v2_spat_thin_100km_max_30_r_per_BCR_centroids.shp")) # output of 1_1_dataprep_BBS_route_selection.R

# selected routes and focal years matched to environmental data:
load(file = file.path(dir, "data", "route_year_env_data.RData")) # route_sel_env_dt_final; output 1_3_dataprep_match_BBS_routes_env_data.R

# selected species:
load(file = file.path(dir, "data", "final_species_selection.RData")) # species_selection_final; output of 1_2_dataprep_BBS_species_selection.R

# selected variables:
load(file = file.path(dir, "data", "selected_variables.RData")) # selvar_final; output of 1_2a_dataprep_env_variable_selection.R

# routes-years:
load(file = file.path(dir, "data", "BBS_for_occ_selection.RData")) # route_sel_dt; output of 1_3_dataprep_match_BBS_routes_env_data.R

# environmental data to check environmental similarity between training and test folds:

# files:
bioclim_files <- list.files(file.path(clim_path, "bioclim"), full.names = TRUE)
lu_files <- list.files(file.path(lu_path, "ISIMIP_LU_ESRI102003"), full.names = TRUE)
sclim_files <- list.files(file.path(clim_path, "seasonal"), full.names = TRUE)

# rasters:
bioclim_3yrs_sp <- rast(bioclim_files[which(grepl("1992_1995", bioclim_files))])
sclim_3yrs_sp <- rast(sclim_files[which(grepl("1992_1995", sclim_files))])
lu_3yrs_sp <- rast(lu_files[which(grepl("1992_1994", lu_files))])

# reduce to selected variables:
bioclim_3yrs_sp_sel <- bioclim_3yrs_sp[[selvar_final[grepl(pattern = "bio", x = selvar_final)]]]
sclim_3yrs_sp_sel <- sclim_3yrs_sp[[selvar_final[grepl(pattern = "(spring|summer|autumn|winter)", x = selvar_final)]]]
lu_3yrs_sp_sel <- lu_3yrs_sp[[selvar_final[!grepl(pattern = "bio|pr_", x = selvar_final)]]] 

# combine:
env_rasters <- c(bioclim_3yrs_sp_sel, sclim_3yrs_sp_sel, lu_3yrs_sp_sel)


# assign routes to folds: ------------------------------------------------------

# iterate over species:

# fold are generated for each species separately since relevant area (750 km buffer around presences)
# is species-specific

for(spec in species_selection_final){
  
  print(spec)
  
  # species data:
  occ_dt_spec <- BBS_pres_abs_spec(species = spec) # from 0_functions.R

  # routes within buffer:
  rel_routes <- training_routes(species = spec, output = "RTENO", buffer_km = 750)
  
  # match route spatial data to species observations, summarise presences over time:
  occ_spec_sf_buffered <- routes_sel_sf %>%
    left_join(occ_dt_spec, by = c("RTENO_BBS" = "RTENO")) %>%
    filter(RTENO_BBS %in% rel_routes) %>% 
    # presence on route across all years:
    group_by(RTENO_BBS) %>%
    summarise(presence_summarised = max(presence, na.rm=TRUE)) %>%
    mutate(presence_summarised = factor(presence_summarised, levels = c(1,0))) %>% 
    arrange(RTENO_BBS)
  
  # create blocks and assign routes to folds:
  sb_US <- cv_spatial(x = occ_spec_sf_buffered,
                      column = "presence_summarised", # the response column (binary or multi-class)
                      k = 5, # number of folds
                      size = hexagon_size_m, # size of the blocks in metres
                      selection = "random", # random blocks-to-fold; more even distribution of presence/absence instances between the train and test folds compared to ‘systematic’
                      iteration = 100, # find evenly dispersed folds; 50 not enough, 100 = default, 200 = no improvement for critical species
                      biomod2 = FALSE, # also create folds for biomod2
                      flat_top = TRUE,
                      seed = 3# 5254,# 3456789: 39
                      )

  # save plot:
  jpeg(file = file.path(dir, "plots", "blockCV_folds", paste0("block_size_", hexagon_size_m/1000, "km"), paste0(spec, ".jpg")), 
       width = 1000, height = 800, quality = 100)
  print(
    cv_plot(cv = sb_US, x = occ_spec_sf_buffered, nrow = 2, points_alpha = 0.5) + 
      ggtitle(paste(spec, hexagon_size_m/1000, "km")) +
      theme(text = element_text(size = 18))
    )
  dev.off()
  
  # save route to block allocation:
  sb_US$blocks
  sb_US$records # number of presences in training data - number of absences in training data - pres. in test data - abs. in test data
  sb_US$folds_ids # which route in which fold
  sb_US$folds_list # for each fold, which routes are in training set, which routes are in test set
  
  save(sb_US, file = file.path(dir, "data", "CV_route_block_allocation", paste0("block_size_", hexagon_size_m/1000, "km"), paste0(spec, ".RData")))
  
  
  # check environmental similarity to evaluate possible extrapolation in testing folds:

  # MESS = multivariate environmental similarity surface
  # how similar a point in a testing fold is to a training fold
  # negative values are sites where at least one variable of env. data has a value 
  # outside the range of environments over the reference set
  env_sim_US <- cv_similarity(cv = sb_US, # the environmental clustering
                              x = occ_spec_sf_buffered, 
                              r = env_rasters, 
                              progress = TRUE)
  
  # write plot to file:
  jpeg(file = file.path(dir, "plots", "blockCV_folds", paste0("block_size_", hexagon_size_m/1000, "km"), paste0(spec, "_env_similarity.jpg")), 
       width = 800, height = 600, quality = 100)
  print(
    env_sim_US +
      ggtitle(paste(spec, hexagon_size_m/1000, "km; neg. = novel env.")) +
      theme(text = element_text(size = 18))
    )
  dev.off()
}


# checks on fold assignment: ---------------------------------------------------

# do we have presences in every year of the training data for each fold (to be able to fit the model):

# write check results to text file:
sink(file.path(dir, "data", paste0("blockCV_",  hexagon_size_m/1000, "km_fold_pres_per_year.txt")))

for(spec in species_selection_final){

  print(spec)
  
  # data within buffer:
  
  rel_routes <- training_routes(species = spec, buffer_km = 750, output = "RTENOs")
  
  occ_dt_spec_buff <- occ_dt_spec %>% 
    filter(RTENO %in% rel_routes) %>% 
    arrange(RTENO)
  
  # assigned routes to folds:
  load(file.path(dir, "data", "CV_route_block_allocation", "block_size_500km", paste0(spec, ".RData")))
  
  # iterate over folds:
  
  for(f in 1:5){

    print(paste("fold", f))
    
    RTENO_training_fold <- occ_spec_sf_buffered$RTENO_BBS[which(sb_US$folds_ids != f)]
    
    # check whether the fold does contain observations for each year:
    presences_routes_year <- occ_dt_spec_buff %>% 
      filter(RTENO %in% RTENO_training_fold) %>% 
      filter(presence == 1) %>% 
      group_by(Year) %>%
      summarise(n_pres = n()) %>%
      print(n = 25) %>% 
      filter(n_pres == 0)
    
    if(nrow(presences_routes_year) != 0){
      print(paste("No routes with presences for", spec, "in fold", f))
      }
  }
}
sink(file = NULL)

## number of presences and absences in training and test data

# write for each species number of presences and absences in training and test data of each fold to file:
sink(file.path(dir, "data", paste0("blockCV_", hexagon_size_m/1000, "km_pres_abs_per_fold.txt")))
for(spec in species_selection_final){
  load(file.path(dir, "data", "CV_route_block_allocation", paste0("block_size_", hexagon_size_m/1000, "km"), paste0(spec, ".RData")))
  print(spec)
  print(sb_US$records)
}
sink(file = NULL)

# checked results, repeated fold assignment
# we want 80 % of presences in training data, 20 % in test data
# minimum 40 presences in training data and 10 presences in test data
# since the included species were detected on at least 50 routes

# session info:
writeLines(capture.output(sessionInfo()), file.path(dir, "results", "sessionInfo", "2_4a_fit_DOMs_CV_fold_assignment.txt"))
