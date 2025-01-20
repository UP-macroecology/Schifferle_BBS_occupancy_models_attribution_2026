# temporal validation
# (model predictions generated with 2_1_2_DOM_flocker_temp_validation.R)

# prediction to last 10 years with model trained with first 15 years


# packages: ----

library(dplyr)
library(flocker)
library(sf)
library(ggplot2)
library(cmdstanr)
#set_cmdstan_path(path = NULL)#set_cmdstan_path("C:/Users/schifferle1/Documents/cmdstan-2.34.1") # xx
set_cmdstan_path("C:/Users/schifferle1/Documents/cmdstan-2.34.1")

# functions: ----

source("0_functions.R")

n_train_years <- 15
buffer_km <- 750

# directories: ----

print(tempdir())
#dir <- file.path("/import", "ecoc9z", "data-zurell", "schifferle", "BBS_occupancy_models_2023")
#dir <- getwd()
results_dir <- file.path("M:", "Documents", "DEBTs", "analysis", "Schifferle_BBS_occupancy_models_2023", 
                         "results", "temp_val_buffer_750_10yrs")
#results_dir <- file.path("//NAS-2-P-SN-01.ibb.uni-potsdam.de", "users$", "schifferle1", "Documents", 
# "DEBTs", "analysis", "Schifferle_BBS_occupancy_models_2023", "results", "temp_val") 


# load data: ----

# selected species:
load(file = file.path("data", "species_set_analysis.RData"))
final_species

# selected routes spatial data (to buffer presences):
routes_sel_sf <- st_read(file.path("data", "route_selection_1995_2019_surv_beg_end_max_5y_miss_v2_spat_thin_100km_max_30_r_per_BCR_centroids.shp")) # output of 1_1_route_selection.R

# routes-years:
load(file = file.path("data", "BBS_for_occ_selection.RData")) # route_sel_dt; output of 1_3_match_BBS_to_env_data.R 

# route-year-species information (only surveyed)
load(file = file.path("data", "BBS_for_occ_spec_records.RData")) # bbs_dt_occ; output of 1_0_reformat_BBS_data.R

years <- min(route_sel_dt$Year):max(route_sel_dt$Year) 


# data prep.: ----

# species:

C_temp_val_df <- data.frame("species" = final_species,
                            "C_ind_10yrs_preds" = NA)

for(i in 77:length(final_species)){
  
  spec <- final_species[i] # 1
  
  print(paste(i, spec))
  
  # check where to look for model output (did MCMC fitting work with less or only with more iterations?)
  if(file.exists(file.path(results_dir, "refit_2000_2000", paste0("out_", spec, "_temp_val_10yrs_buffer_750.RData")))){
    output_dir <- file.path(results_dir, "refit_2000_2000")
  } else {
    output_dir <- results_dir
  }

  print(output_dir)
  
  # observations:
  
  # relevant routes for the species, within distance of 750 km of presences:
  rel_routes <- training_routes(species = spec, buffer_km = buffer_km, output = "RTENOs")
  occ_dt_spec <- BBS_pres_abs_spec(species = spec) %>% 
    filter(RTENO %in% rel_routes)
  
  # sum all routes for each year (temporal trend)
  obs_temp_trend <- occ_dt_spec %>% 
    group_by(Year) %>% 
    summarise(pres_sum = sum(presence, na.rm = TRUE))
  
  # model predictions:
  load(file.path(output_dir, paste0("preds_", spec, "_temp_val_10yrs_buffer_750.RData")))

  #dim(res_list$y_preds) # routes - sections - years - draws
  # sum across route sections:
  preds_routes <- apply(res_list$y_preds, MAR = c(1,3,4), FUN = max)
  dim(preds_routes)
  # sum across routes for each year:
  preds_years <- apply(preds_routes, MAR = c(2,3), FUN = sum, na.rm = TRUE) # xx
  dim(preds_years)
  # mean:
  preds_years_mean <- apply(preds_years, MAR = 1, FUN = mean)
  # median:
  preds_years_median <- apply(preds_years, MAR = 1, FUN = median)
  # sd:
  preds_years_sd <- apply(preds_years, MAR = 1, FUN = sd)
  
  # C-index of last 10 years:
  C_temp_val <- Hmisc::rcorr.cens(x = preds_years_mean[(n_train_years + 1):length(years)], S = obs_temp_trend$pres_sum[(n_train_years + 1):length(years)])
  C_temp_val_df$C_ind_10yrs_preds[i] <- C_temp_val["C Index"]
  
  # plot predicted y against observations summed over years:
  
  # directory to save plots:
  if(!dir.exists(file.path(results_dir, "temp_eval", "10_years"))){
    dir.create(file.path(results_dir, "temp_eval", "10_years"), recursive = TRUE)
  }
  
  jpeg(file = file.path(results_dir, "temp_eval", "10_years", paste0("temp_val_10yrs_", spec,"_", buffer_km, "km.jpg")), 
       width = 1000, height = 700, quality = 100)
  print(
    ggplot(obs_temp_trend, aes(x = Year)) +
      geom_line(aes(y = pres_sum)) +
      geom_point(aes(y = pres_sum), size = 3) +
      geom_ribbon(aes(y = preds_years_mean, ymax = preds_years_mean + 1*preds_years_sd, ymin = preds_years_mean - 1*preds_years_sd), 
                  alpha = 0.2, fill = "cornflowerblue") +
      geom_line(aes(y = preds_years_mean), color = "cornflowerblue") +
      geom_point(aes(y = preds_years_mean), color = "cornflowerblue", size = 3) +
      ylab("N routes with presence") +
      theme_bw() +
      theme(text = element_text(size = 20)) +
      geom_vline(xintercept = years[1]+(n_train_years-1), linetype = "dashed") +
      ggtitle(paste0(spec, ": y preds (mean +/- 1 sd) vs. obs. (C index last 10 yrs: ", round(C_temp_val["C Index"], 2), ")"))
  )
  dev.off()
  
}


save(C_temp_val_df, file = file.path(results_dir, "temp_eval", "10_years", "C_temp_val_10yrs.RData"))

#load(file.path(results_dir, "temp_eval", "10_years", "C_temp_val_10yrs.RData"))
# decline was observed, this decline cannot be captured by the model, 
# this suggests that species climate and land use change in the breeding area is not the main reason for decline
# (or that climate and land use data we used are not properly capturing climate and land use change in the breeding area)

