# Script:   5_4a_attribution_maps_extract_species_ranges.R
# Purpose:  Extract species ranges from BirdLife data base (to then map spatial patterns)
# Inputs:   results/species_DOM_val_okay.RData
#           data/BBS_data_merged.RData
#           <datashare_Birdlife>/BOTW.gdb
# Outputs:  data/Birdlife_range_maps/<species> (one shapefile per species)
# Runs on:  HPC (NAS Potsdam)

source(file.path("scripts", "0_paths.R"))


# packages: --------------------------------------------------------------------

library(dplyr)
library(doParallel)
library(gdalUtilities)


# directories: -----------------------------------------------------------------

# directory to save extracted BirdLife ranges:
output_dir <- file.path(hpc_dir, "data", "Birdlife_range_maps")


# prepare data: ----------------------------------------------------------------

# species for attribution:
load(file = file.path(hpc_dir, "results", "species_DOM_val_okay.RData")) # output of 4_0_DOMs_predictions_y_routes_scenarios.R
spec_attr

# add scientific names:
load(file.path(hpc_dir, "data", "BBS_data_merged.RData")) # bbs_dt; output of 1_0_dataprep_BBS_bird_data.R
spec_names <- bbs_dt %>% 
  select(English_Common_Name, Scientific_Name) %>% 
  distinct %>% 
  filter(English_Common_Name %in% spec_attr)


## shapefile extraction from Birdlife gdb: -------------------------------------

# register cores for parallel computation:
registerDoParallel(cores = 20)

# only parts of range where species is considered extant (presence = 1) and which
# is used either throughout the whole year (seasonal = 1) or during the breeding season (seasonal = 2)

# # account for taxonomic changes between BBS and BL range maps:
# # (resources:
# # https://www.iucnredlist.org/
# # https://explorer.natureserve.org/Search)
spec_name_change_df <- data.frame("BBS_name" = "Dryobates villosus", "BL_name" = "Leuconotopicus villosus")
spec_name_change_df[2,] <- c("Dryocopus pileatus", "Hylatomus pileatus")

foreach(s = 1:nrow(spec_names),
        .packages = c("gdalUtilities"),
        .verbose = TRUE,
        .errorhandling = "remove",
        .inorder = FALSE) %dopar% {
          
          spec <- spec_names$Scientific_Name[s]
          
          # check if species was already processed:
          exists <- file.exists(file.path(hpc_dir, "data", "Birdlife_range_maps", paste0(spec_names$English_Common_Name[s], ".shp")))
          if(exists) {
            print(paste(spec, "ran already."))
            next
          }
          
          # adjust species name if BL uses other taxonomy than BBS:
          if(spec %in% spec_name_change_df$BBS_name){
            spec <- spec_name_change_df$BL_name[which(spec_name_change_df$BBS_name == spec)]
          }
          
          # extract shapefile:
          gdalUtilities::ogr2ogr(src_datasource_name = file.path(datashare_Birdlife, "BOTW.gdb"),
                                 layer = "All_Species",
                                 where = paste0("sci_name = '", spec, "' AND (SEASONAL = '1' OR SEASONAL = '2') AND PRESENCE = '1'"), # SEASONAL = 1: resident throughout the year, SEASONAL = 2: breeding season
                                 dst_datasource_name = file.path(output_dir, paste0(spec_names$English_Common_Name[s], ".shp")),
                                 overwrite = TRUE)
        }

# session info:
writeLines(capture.output(sessionInfo()), file.path(hpc_dir, "results", "sessionInfo", "5_4a_attribution_maps_extract_species_ranges.txt"))