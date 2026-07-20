# Script:   1_3_dataprep_match_BBS_routes_env_data.R
# Purpose:  Compile environmental data to fit dynamic occupancy models and to run counterfactual simulations
# Inputs:   data/BBS_for_occ.RData
#           data/route_selection_1995_2019_surv_beg_end_max_5y_miss_v2_spat_thin_100km_max_30_r_per_BCR_centroids.shp
#           data/selected_variables.RData
#           <clim_path>/bioclim/<var>_<year>.tif
#           <clim_path>/seasonal/<var>_<year>.tif
#           <lu_path>/ISIMIP_LU_ESRI102003/<var>_<year>_ESRI102003.tif
#           data/Counterfactual_env_data/ISIMIP_GSWP3_W5E5/attrici_detrending/output/ATTRICI_CLIM_ESRI102003_tifs/bioclim/<var>_<year>.tif
#           data/Counterfactual_env_data/ISIMIP_GSWP3_W5E5/attrici_detrending/output/ATTRICI_CLIM_ESRI102003_tifs/seasonal/<var>_<year>.tif
# Outputs:  data/BBS_for_occ_selection.RData
#           data/route_year_env_data.RData
#           data/route_year_env_data_cf.RData
#           plots/change_exposure_single_vars[1-3].svg (Fig. S7)
# Runs on:  Local
# Notes:    this script is run twice, once to compile factual data (set data <- "factual"),
#           once to compile counterfactual data (set data <- "counterfactual")

# Steps:
# extract values of selected climatic and land use variables at BBS route centroids;
# factual data to fit dynamic occupancy models
# counterfactual data to simulate scenarios without climate and land use change
## counterfactual climate: ISIMIP's climate data detrended with ATTRICI from 1995 onwards
## land use: fixed at 1995er level
## factual 3-year summaries as covariates for initial occupancy
# 1) generate subset of BBS data that covers selected routes and years
# 2) extract values of climatic and land use variables at each route for each year
# 3) extract 3-year summaries of before the focal period (as covariates of initial occupancy)


source(file.path("scripts", "0_paths.R"))


# packages: --------------------------------------------------------------------

library(sf)
library(dplyr)
library(terra)
library(ggplot2)
library(tidyterra)
library(patchwork)

# choose which environmental data set to process:

#data <- "factual"
data <- "counterfactual"

# load data: -------------------------------------------------------------------

# BBS data (which routes surveyed in which years) formatted for occupancy modelling:
load(file = file.path(dir, "data", "BBS_for_occ.RData")) # output of 1_0_dataprep_BBS_bird_data.R
route_dt
nrow(route_dt) # 224124

# BBS route selection (route centroids):
routes_sel_sf <- st_read(file.path(dir, "data", "route_selection_1995_2019_surv_beg_end_max_5y_miss_v2_spat_thin_100km_max_30_r_per_BCR_centroids.shp")) # output of 1_1_dataprep_BBS_route_selection.R
nrow(routes_sel_sf) # 539

# selected variables:
load(file = file.path(dir, "data", "selected_variables.RData")) # output of 1_2a_dataprep_env_variable_selection.R
selvar_final

# environmental data files:

if(data == "factual"){
  
  # bioclimatic variables:
  bioclim_files <- list.files(file.path(clim_path, "bioclim"), full.names = TRUE)
  # seasonal climatic variables:
  sclim_files <- list.files(file.path(clim_path, "seasonal"), full.names = TRUE)
  # land use variables:
  lu_files <- list.files(file.path(lu_path, "ISIMIP_LU_ESRI102003"), full.names = TRUE)
}

if(data == "counterfactual"){
  
  env_dir <- file.path(dir, "data", "Counterfactual_env_data", "ISIMIP_GSWP3_W5E5", "attrici_detrending", "output", "ATTRICI_CLIM_ESRI102003_tifs")
  
  # bioclimatic variables:
  bioclim_files <- list.files(file.path(env_dir, "bioclim"), full.names = TRUE)
  # seasonal climatic variables:
  sclim_files <- list.files(file.path(env_dir, "seasonal"), full.names = TRUE)
  # land use variables (factual of 1995):
  lu_files <- list.files(file.path(lu_path, "ISIMIP_LU_ESRI102003"), full.names = TRUE)
}


# 1) subset BBS routes based on selected routes and focal time period: ---------

route_sel_dt <- route_dt %>%
  filter(RTENO %in% routes_sel_sf$RTENO_BBS) %>%
  filter(Year >= 1995 & Year <= 2019) %>% 
  arrange(RTENO)
nrow(route_sel_dt) # 13475
#save(route_sel_dt, file = file.path(dir, "data", "BBS_for_occ_selection.RData"))


# 2) env. data for each route-year combination: --------------------------------


# iterate over years:

# extract BBS data matching year i, load bioclim and land use data of year i, extract values at route centroids:

years <- seq(min(route_sel_dt$Year), max(route_sel_dt$Year))

for(i in 1:length(years)){
  
  print(i)
  
  # routes surveyed in year i:
  route_IDs_year <- route_sel_dt %>% 
    filter(Year == years[i]) %>% 
    pull(RTENO)
  routes_sel_year_sf <- routes_sel_sf %>% 
    filter(RTENO_BBS %in% route_IDs_year) %>% 
    mutate(Year = years[i])
  
  # bioclimatic variables of year i:
  bioclim_year <- rast(bioclim_files[which(grepl(paste0("bio.{1,2}_", years[i], ".tif"), bioclim_files))])
  # reduce to selected variables:
  bioclim_year_sel <- bioclim_year[[selvar_final[grepl(pattern = "bio", x = selvar_final)]]]
  # extract values of each bioclimatic variable at each relevant route location:
  for(biovar in names(bioclim_year_sel)){
    routes_sel_year_sf[, biovar] <- bioclim_year_sel[[biovar]] %>% 
      terra::extract(y = routes_sel_year_sf) %>% 
      pull(biovar)
  }
  
  # seasonal climate variables of year i:
  sclim_year <- rast(sclim_files[which(grepl(paste0("(spring|summer|autumn|winter)", "_", years[i], ".tif"), sclim_files))])
  # reduce to selected variables:
  sclim_year_sel <- sclim_year[[selvar_final[grepl(pattern = "pr_", x = selvar_final)]]]
  # extract values of each variable at each relevant route location:
  for(sclimvar in names(sclim_year_sel)){
    routes_sel_year_sf[, sclimvar] <- sclim_year_sel[[sclimvar]] %>% 
      terra::extract(y = routes_sel_year_sf) %>% 
      pull(sclimvar)
  }
  
  # land use variables of year i:
  
  if(data == "factual"){
    lu_year <- rast(lu_files[which(grepl(paste0(years[i], "_ESRI102003.tif$"), lu_files))])
  }
  
  if(data == "counterfactual"){
    # constant values at level of 1995:
    lu_year <- rast(lu_files[which(grepl(paste0("1995", "_ESRI102003.tif$"), lu_files))])
  }
  
  # reduce to selected variables:
  lu_year_sel <- lu_year[[selvar_final[!grepl(pattern = "bio|pr_|mean", x = selvar_final)]]]
  # extract values of each land use class at each relevant route location:
  for(luvar in names(lu_year_sel)){
    routes_sel_year_sf[, luvar] <- lu_year_sel[[luvar]] %>% 
      terra::extract(y = routes_sel_year_sf) %>% 
      pull(luvar)
  }

  
  # merge data for all years:
  if(i == 1){
    routes_sel_all_sf <- routes_sel_year_sf
  } else{
    routes_sel_all_sf <- rbind(routes_sel_all_sf, routes_sel_year_sf)
  }
}
routes_sel_all_sf

# merge with BBS data:
route_sel_env_dt1 <- route_sel_dt %>% 
  left_join(routes_sel_all_sf, by = c(RTENO = "RTENO_BBS", Year = "Year")) %>% 
  select(-geometry)


# 3) variables summarising 3 years before focal period: ------------------------

# predictors for initial occupancy probability, same for factual and counterfactual simulations

# bioclimatic variables:
bioclim_files <- list.files(file.path("data", "Env_data", "ISIMIP_GSWP3_W5E5", "bioclim"), full.names = TRUE)
# seasonal climatic variables:
sclim_files <- list.files(file.path("data", "Env_data", "ISIMIP_GSWP3_W5E5", "seasonal"), full.names = TRUE)
# land use variables:
lu_files <- list.files(file.path("data", "Env_data", "ISIMIP_land_use_and_irrigation", "ISIMIP_LU_ESRI102003"), full.names = TRUE)


# bioclimatic variables:
bioclim_3yrs_sp <- rast(bioclim_files[which(grepl("1992_1995", bioclim_files))])
# reduce to selected variables:
bioclim_3yrs_sp_sel <- bioclim_3yrs_sp[[selvar_final[grepl(pattern = "bio", x = selvar_final)]]]
# extract value of each bioclimatic variable at each route centroid:
for(biovar in names(bioclim_3yrs_sp_sel)){
  routes_sel_sf[, paste0(biovar, "_3yrs")] <- bioclim_3yrs_sp_sel[[biovar]] %>% 
    terra::extract(y = routes_sel_sf) %>% 
    pull(biovar)
}

# seasonal climate variables:
sclim_3yrs_sp <- rast(sclim_files[which(grepl("1992_1995", sclim_files))])
# reduce to selected variables:
sclim_3yrs_sp_sel <- sclim_3yrs_sp[[selvar_final[grepl(pattern = "(spring|summer|autumn|winter)", x = selvar_final)]]]
# extract value of each bioclimatic variable at each route centroid:
for(sclimvar in names(sclim_3yrs_sp_sel)){
  routes_sel_sf[, paste0(sclimvar, "_3yrs")] <- sclim_3yrs_sp_sel[[sclimvar]] %>% 
    terra::extract(y = routes_sel_sf) %>% 
    pull(sclimvar)
}

# land use variables:
lu_3yrs_sp <- rast(lu_files[which(grepl("1992_1994", lu_files))])
# reduce to selected variables:
lu_3yrs_sp_sel <- lu_3yrs_sp[[selvar_final[!grepl(pattern = "bio|pr_|mean", x = selvar_final)]]]
# extract value of each lu variable at each route centroid:
for(luvar in names(lu_3yrs_sp_sel)){
  routes_sel_sf[, paste0(luvar, "_3yrs")] <- lu_3yrs_sp_sel[[luvar]] %>% 
    terra::extract(y = routes_sel_sf) %>% 
    pull(luvar)
}

# match to BBS data:
route_sel_env_dt_final <- route_sel_env_dt1 %>% 
  left_join(routes_sel_sf, by = c(RTENO = "RTENO_BBS")) %>% 
  select(-geometry) %>% 
  arrange(RTENO)

# write to file:

if(data == "factual") {
  save(route_sel_env_dt_final, file = file.path(dir, "data", "route_year_env_data.RData"))
  write.csv(route_sel_env_dt_final, file = file.path(dir, "data", "route_year_env_data.csv"),
            row.names = FALSE)
}
if(data == "counterfactual") {
  save(route_sel_env_dt_final, file = file.path(dir, "data", "route_year_env_data_cf.RData"))
  write.csv(route_sel_env_dt_final, file = file.path(dir, "data", "route_year_env_data_cf.csv"),
            row.names = FALSE)
}


# plot exposure to climate and land use change: --------------------------------

# compare mean of 1995-1997 to mean of 2017 - 2019?

# selected variables:
selvar_final


# bioclimatic variables:
bioclim_files <- list.files(file.path(clim_path, "bioclim"), 
                            pattern = paste0("(", paste0(grep("bio", selvar_final, value = TRUE), collapse = "_|"), ")"),
                            full.names = TRUE)
bioclim_1995 <- rast(bioclim_files[which(grepl(paste0("bio.{1,2}_", "1995", ".tif"), bioclim_files))])
bioclim_1996 <- rast(bioclim_files[which(grepl(paste0("bio.{1,2}_", "1996", ".tif"), bioclim_files))])
bioclim_1997 <- rast(bioclim_files[which(grepl(paste0("bio.{1,2}_", "1997", ".tif"), bioclim_files))])
bioclim_start <- mean(bioclim_1995, bioclim_1996, bioclim_1997)
bioclim_start$bio1 <- bioclim_start$bio1 - 273.15 # K to Celsius
bioclim_2017 <- rast(bioclim_files[which(grepl(paste0("bio.{1,2}_", "2017", ".tif"), bioclim_files))])
bioclim_2018 <- rast(bioclim_files[which(grepl(paste0("bio.{1,2}_", "2018", ".tif"), bioclim_files))])
bioclim_2019 <- rast(bioclim_files[which(grepl(paste0("bio.{1,2}_", "2019", ".tif"), bioclim_files))])
bioclim_end <- mean(bioclim_2017, bioclim_2018, bioclim_2019)
bioclim_end$bio1 <- bioclim_end$bio1 - 273.15 # K to Celsius

# seasonal climatic variables:
sclim_files <- list.files(file.path(clim_path, "seasonal"), full.names = TRUE)
sclim_1995 <- rast(sclim_files[which(grepl(paste0("pr_mean_", "(spring|summer|autumn|winter)", "_", "1995", ".tif"), sclim_files))])
sclim_1996 <- rast(sclim_files[which(grepl(paste0("pr_mean_", "(spring|summer|autumn|winter)", "_", "1996", ".tif"), sclim_files))])
sclim_1997 <- rast(sclim_files[which(grepl(paste0("pr_mean_", "(spring|summer|autumn|winter)", "_", "1997", ".tif"), sclim_files))])
sclim_start <- mean(sclim_1995, sclim_1996, sclim_1997)
sclim_2017 <- rast(sclim_files[which(grepl(paste0("pr_mean_", "(spring|summer|autumn|winter)", "_", "2017", ".tif"), sclim_files))])
sclim_2018 <- rast(sclim_files[which(grepl(paste0("pr_mean_", "(spring|summer|autumn|winter)", "_", "2018", ".tif"), sclim_files))])
sclim_2019 <- rast(sclim_files[which(grepl(paste0("pr_mean_", "(spring|summer|autumn|winter)", "_", "2019", ".tif"), sclim_files))])
sclim_end <- mean(sclim_2017, sclim_2018, sclim_2019)

# land use variables:
lu_files <- list.files(file.path(lu_path, "ISIMIP_LU_ESRI102003"), full.names = TRUE)

lu_1995 <- rast(grep(pattern = paste0("(", paste0(selvar_final[!grepl("(bio|pr_|mean)",  x = selvar_final)], collapse = ")|("), ")"),
                      x = lu_files[grep(pattern = "1995", x = lu_files)], value = TRUE))
lu_1996 <- rast(grep(pattern = paste0("(", paste0(selvar_final[!grepl("(bio|pr_|mean)",  x = selvar_final)], collapse = ")|("), ")"),
                     x = lu_files[grep(pattern = "1996", x = lu_files)], value = TRUE))
lu_1997 <- rast(grep(pattern = paste0("(", paste0(selvar_final[!grepl("(bio|pr_|mean)",  x = selvar_final)], collapse = ")|("), ")"),
                     x = lu_files[grep(pattern = "1997", x = lu_files)], value = TRUE))
lu_start <- mean(lu_1995, lu_1996, lu_1997)
lu_2017 <- rast(grep(pattern = paste0("(", paste0(selvar_final[!grepl("(bio|pr_|mean)",  x = selvar_final)], collapse = ")|("), ")"),
                     x = lu_files[grep(pattern = "2017", x = lu_files)], value = TRUE))
lu_2018 <- rast(grep(pattern = paste0("(", paste0(selvar_final[!grepl("(bio|pr_|mean)",  x = selvar_final)], collapse = ")|("), ")"),
                     x = lu_files[grep(pattern = "2018", x = lu_files)], value = TRUE))
lu_2019 <- rast(grep(pattern = paste0("(", paste0(selvar_final[!grepl("(bio|pr_|mean)",  x = selvar_final)], collapse = ")|("), ")"),
                     x = lu_files[grep(pattern = "2019", x = lu_files)], value = TRUE))
lu_end <- mean(lu_2017, lu_2018, lu_2019)

var_start <- c(bioclim_start, sclim_start, lu_start)
var_start$bio15 <- var_start$bio15*100 # precipitation seasonality in percent
var_start$urbanareas <- var_start$urbanareas*100 # in percent
var_start$managed_pastures <- var_start$managed_pastures*100 # in percent
var_start$primary_nonforests <- var_start$primary_nonforests*100 # in percent
var_start$secondary_nonforests <- var_start$secondary_nonforests*100 # in percent
var_start$sum_annual_crops <- var_start$sum_annual_crops*100 # in percent

# plot change in absolute values for each variable:

bioclim_diff <- bioclim_end - bioclim_start
sclim_diff <- sclim_end - sclim_start
lu_diff <- lu_end - lu_start
vars_diff <- c(bioclim_diff, sclim_diff, lu_diff)
vars_diff$bio15 <- vars_diff$bio15*100 # precipitation seasonality in percent
vars_diff$urbanareas <- vars_diff$urbanareas*100 # in percent
vars_diff$managed_pastures <- vars_diff$managed_pastures*100 # in percent
vars_diff$primary_nonforests <- vars_diff$primary_nonforests*100 # in percent
vars_diff$secondary_nonforests <- vars_diff$secondary_nonforests*100 # in percent
vars_diff$sum_annual_crops <- vars_diff$sum_annual_crops*100 # in percent


# plot mean of 1995-1997 next to change to show base on which to interpret change values:

# range of change for different sets of variables:
# monthly / seasonal precipitation:
prec_change <- subset(x = vars_diff, 
                      subset = c("bio14", "pr_mean_spring", "pr_mean_summer", "pr_mean_autumn", "pr_mean_winter"))
min_prec_change <- min(global(prec_change, "min", na.rm = TRUE))
max_prec_change <- max(global(prec_change, "max", na.rm = TRUE))

# # temperature change:
# temp_change <- subset(x = vars_diff, 
#                       subset = c("bio1", "bio2", "bio7"))
# min_temp_change <- min(global(temp_change, "min", na.rm = TRUE))
# max_temp_change <- max(global(temp_change, "max", na.rm = TRUE))

# land use change:
lu_change <- subset(x = vars_diff, 
                      subset = c("managed_pastures", "primary_nonforests", "secondary_nonforests",
                                 "sum_annual_crops", "urbanareas"))
min_lu_change <- min(global(lu_change, "min", na.rm = TRUE))
max_lu_change <- max(global(lu_change, "max", na.rm = TRUE))


plot_list <- vector(mode = "list", length = length(selvar_final))

names_long <- c("annual mean temperature", "diurnal temperature range", "isothermality", "annual temperature range",
                "precipitation driest month", "precipitation seasonality",
                "precipitation spring", "precipitation summer", "precipitation autumn", "precipitation winter",
                "urban", "managed pastures", "primary non-forests", "secondary non-forests",
                "annual crops")

unit <- c("°C", "°C", "%", "°C", 
          "kg/(m2*s)", "%",
          "kg/(m2*s)",  "kg/(m2*s)",  "kg/(m2*s)",  "kg/(m2*s)",
          "% cover", "% cover", "% cover", "% cover", "% cover")

for(i in 1:length(selvar_final)){
  
  var <- selvar_final[i]
  print(var)
  
  # limits of change colour scale:
  if(var %in% names(prec_change)){
    limits <- c(min_prec_change, max_prec_change)
 # } else if (var %in% names(temp_change)){
  #  limits <- c(min_temp_change, max_temp_change)
  } else if (var %in% names(lu_change)){
    limits <- c(min_lu_change, max_lu_change)
  } else {
    limits = NULL
  }
  
  # mean 1995-1997:
  plot_init <- ggplot() +
    geom_spatraster(data = var_start[[var]]) +
    # add selected routes:
    geom_sf(data = routes_sel_sf, colour = "grey30", shape = 19, size = 0.05) +
    scale_fill_viridis_c(na.value = "transparent",
                         name = paste0("mean\n1995-1997\n[", unit[i], "]"),
                         labels = scales::math_format(.x)) +
    theme_bw() +
    ggtitle(names_long[i]) +
    theme(text = element_text(size = 14),
          #legend.title = element_text(size = 14),
          plot.margin = unit(c(-1,0,-1,0), "cm"))
  
  # change to 2017-2019:
  plot_change <- ggplot() +
    geom_spatraster(data = vars_diff[[var]]) +
    # add selected routes:
    geom_sf(data = routes_sel_sf, colour = "grey30", shape = 19, size = 0.05) +
    scale_fill_gradient2(low = "dodgerblue", mid = "lightyellow", high = "red3", midpoint = 0, 
                         na.value = "transparent",
                         limits = limits,
                       name = paste0("change 1995-1997\nto 2017-2019\n[", unit[i], "]")) +
    theme_bw() +
    theme(text = element_text(size = 14),
          #legend.title = element_text(size = 12),
          plot.margin = unit(c(-1,0,-1,0), "cm"))
  
  plot_list[[i]] <- (plot_init | plot_change) + 
    plot_annotation(
      theme = theme(plot.title = element_text(size = 12))
    )
}

# save:
# pdf(file.path(dir, "plots", "change_exposure_single_vars.pdf"), 
#     onefile = TRUE,
#     width = 8, height = 11)
# wrap_plots(plot_list[1:5], ncol = 1, nrow = 5)
# wrap_plots(plot_list[6:10], ncol = 1, nrow = 5)
# wrap_plots(plot_list[11:15], ncol = 1, nrow = 5)
# dev.off()

svg(file.path(dir, "plots", "change_exposure_single_vars1.svg"), 
    onefile = TRUE,
    width = 10, height = 12)
wrap_plots(plot_list[1:5], ncol = 1, nrow = 5)
dev.off()
svg(file.path(dir, "plots", "change_exposure_single_vars2.svg"), 
    onefile = TRUE,
      width = 10, height = 12)
wrap_plots(plot_list[6:10], ncol = 1, nrow = 5)
dev.off()
svg(file.path(dir, "plots", "change_exposure_single_vars3.svg"), 
    onefile = TRUE,
    width = 10, height = 12)
wrap_plots(plot_list[11:15], ncol = 1, nrow = 5)
dev.off()

# session info:
writeLines(capture.output(sessionInfo()), file.path(dir, "results", "sessionInfo", "1_3_dataprep_match_BBS_routes_env_data.txt"))
