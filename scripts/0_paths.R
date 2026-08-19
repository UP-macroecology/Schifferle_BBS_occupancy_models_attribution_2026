# Script: 0_paths.R
# Purpose: Defines paths and create folder structure used throughout the analyses
# Inputs:  -
# Outputs: -
# Runs on: Local


# project directory: -----------------------------------------------------------

# project directory from local computer:
dir <- file.path("//NAS-2-P-SN-01.ibb.uni-potsdam.de", "daten$", "AG26", "Transfer", "Schifferle_BBS_occupancy_models_2023")

# project directory from HPC, for computationally heavy steps:
hpc_dir <- file.path("/mnt", "ibb_share", "zurell_transfer", "Schifferle_BBS_occupancy_models_2023")


# create folder structure: -----------------------------------------------------

# dir before data and new plot structure? e.g. plots_manuscript/xx?

# data folder:
if(!dir.exists(file.path(dir, "data"))){
  dir.create(file.path(dir, "data"))
}

# results folder:
if(!dir.exists(file.path(dir, "results"))){
  dir.create(file.path(dir, "results"))
}

# plots folder:
if(!dir.exists(file.path(dir, "plots"))){
  dir.create(file.path(dir, "plots"))
}

# session info folder:
sessioninfo_dir <- file.path(dir, "results", "sessionInfo")
if(!dir.exists(sessioninfo_dir)){dir.create(sessioninfo_dir, recursive = TRUE)}


# directories downloaded data: --------------------------------------------

# directory BBS download:
datashare_BBS <- file.path("//NAS-2-P-SN-01.ibb.uni-potsdam.de/", "daten$", "AG26", "Arbeit", "datashare", "data", "biodat", "distribution", "BBS", "NABBS_2023")

# BBS routes download:
bbs_routes_file <- file.path("//NAS-2-P-SN-01.ibb.uni-potsdam.de/", "daten$", "AG26", "Arbeit", "datashare", "data", "biodat", "distribution", "BBS", "bbs_routes", "bbsrtsl020.shp")

# BirdLife range maps:
datashare_Birdlife <- file.path("/mnt", "ibb_share", "zurell", "biodat", "distribution", "Birdlife", "BOTW_2022", "BOTW.gdb")

# Bird Conservation Region:
bcr_dir <- file.path(dir, "data", "BCR_Terrestrial")

# climate and land use data:

# save data here as:
# gswp3-w5e5_obsclim_<var>_lat24.0to50.0lon-126.0to-66.0_daily_<year1>_<year2>.nc
clim_path <- file.path(dir, "data", "Env_data", "ISIMIP_GSWP3_W5E5")
if(!dir.exists(clim_path)){dir.create(clim_path, recursive = TRUE)}

# land use data:
lu_path <- file.path(dir, "data", "Env_data", "ISIMIP_land_use_and_irrigation")
if(!dir.exists(lu_path)){dir.create(lu_path, recursive = TRUE)}