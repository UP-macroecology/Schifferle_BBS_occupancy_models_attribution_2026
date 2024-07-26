# calculate evaluation measures for cross validation 
# (model predictions generated with 2_1_DOM_flocker_CV.R)

#  calculate Harrel's C indices
# - spatially, temporally


# packages: ----

library(dplyr)
library(flocker)
library(sf)
library(ggplot2)


# functions: ----

source("0_functions.R")


# directories: ----

print(tempdir())
#dir <- file.path("/import", "ecoc9z", "data-zurell", "schifferle", "BBS_occupancy_models_2023")
#dir <- getwd()
results_dir <- file.path("M:", "Documents", "DEBTs", "analysis", "Schifferle_BBS_occupancy_models_2023", "results") 


# load data: ----

# selected species:
load(file = file.path("data", "final_species_selection.RData")) # species_selection_final; output of 1_2_species_selection.R

# selected routes spatial data (to buffer presences):
routes_sel_sf <- st_read(file.path("data", "route_selection_1991_2015_surv_beg_end_max_5y_miss_v2_spat_thin_100km_max_30_r_per_BCR_centroids.shp")) # output of 1_1_route_selection.R

# routes-years:
load(file = file.path("data", "BBS_for_occ_selection.RData")) # route_sel_dt; output of 1_3_match_BBS_to_env_data.R 

# route-year-species information (only surveyed)
load(file = file.path("data", "BBS_for_occ_spec_records.RData")) # bbs_dt_occ; output of 1_0_reformat_BBS_data.R

years <- 1991:2015


# data prep.: ----

# species:
spec <- species_selection_final[1]


# load model predictions of each test fold:

res_lists_folds <- vector(mode = "list", length = 5)

for(fold in 1:5){
  
  print(fold)
  
  # assemble fold results:
  load(file = file.path(results_dir, paste0("out_", spec, "_CV_fold", fold, ".RData")))
  print(out)
  
  # assemble fold results:
  load(file = file.path(results_dir, paste0("test_preds_", spec, "_CV_fold", fold, ".RData")))
  res_lists_folds[[fold]] <- res_list
}


# find route IDs matching the predictions:

# relevant routes for the species, within distance of 750 km of presences:
rel_routes <- training_routes(species = spec, buffer_km = 750, output = "RTENOs")

# load fold assignment (sorting important!):
load(file.path("data", "CV_route_block_allocation", "block_size_500km", paste0(spec, ".RData")))

test_RTENOs <- c(sort(rel_routes[sb_US$folds_list[[1]][[2]]]), # test data fold 1
                 sort(rel_routes[sb_US$folds_list[[2]][[2]]]),
                 sort(rel_routes[sb_US$folds_list[[3]][[2]]]),
                 sort(rel_routes[sb_US$folds_list[[4]][[2]]]),
                 sort(rel_routes[sb_US$folds_list[[5]][[2]]]))


# observations can be compared either to the predicted occupancy probability or to predictions of y, in the latter 
# case, imperfect detection is taken into account


# # plot map comparing mean predictions and observations for single year: ---
# 
# # mean y predictions
# mean_y_preds <- apply(y_preds_routes, MAR = c(1,2), FUN = mean)
# 
# # observations:
# occ_dt_spec <- BBS_pres_abs_spec(species = spec)
# 
# obs_test <- occ_dt_spec %>% # RTENOs ordered
#   filter(RTENO %in% test_RTENOs) %>% 
#   filter(Year == years[2]) %>% 
#   arrange(RTENO) %>% 
#   mutate(mean_y_preds = mean_y_preds[,2]) %>% 
#   left_join(routes_sel_sf, by = c(RTENO = "RTENO_BBS")) %>% 
#   st_as_sf()
# 
# ggplot(obs_test) +
#   geom_sf(aes(fill = as.factor(presence)), pch = 21, size = 3) +
#   scale_fill_viridis_d(name = "observation") +
#   geom_sf(aes(color = mean_y_preds), size = 1) +
#   scale_color_viridis_c(name = "mean prediction") +
#   theme_bw() +
#   ggtitle(spec)



# spatial C index: ---

# do we get differences between sites right?

# for each year separately, calculate Harrel's C index for each MCMC jump (-> distribution of C values)

# C index near 1 -> good
# C index = 0.5 -> model as good as random guessing of which of two routes has a higher probability of being occupied


# 1) based on predictions of y (= occupancy prob. * detection prob.):

# predictions of y for each route section, sum across routes:

# bind predictions for all folds:
y_preds_all_routes <- abind::abind(res_lists_folds[[1]]$y_preds,
                                   res_lists_folds[[2]]$y_preds,
                                   res_lists_folds[[3]]$y_preds,
                                   res_lists_folds[[4]]$y_preds,
                                   res_lists_folds[[5]]$y_preds,
                                   along = 1)
dim(y_preds_all_routes) # routes - sections - years - draws

# sum over route sections:
y_preds_routes <- apply(y_preds_all_routes, MAR = c(1,3,4), FUN = sum)
y_preds_routes[which(y_preds_routes != 0)] <- 1 # convert detection sum across route sections to presence / absence
dim(y_preds_routes) # for each site and year 4000 draws


# for all years:
rcorr_ypreds_list <- vector(mode = "list", length = length(years))

for(i in 1:length(years)){
  
  print(i)
  
  # observations:
  obs_i <- occ_dt_spec %>% 
    filter(RTENO %in% test_RTENOs) %>% 
    # order same as RTENOs in test data:
    arrange(factor(RTENO, levels = test_RTENOs)) %>%
    filter(Year == years[i]) %>% 
    pull(presence) 
  
  # C index, predictions:
  C_distr_year_i <- apply(X = y_preds_routes[,i,], MARGIN = 2, FUN = Hmisc::rcorr.cens, x = obs_i, outx=FALSE) # TRUE? xx
  
  # reformat:
  rcorr_ypreds_list[[i]] <-  C_distr_year_i %>% 
    as_tibble(rownames = "metric") %>% 
    tidyr::pivot_longer(cols = V1:V4000, values_to = "draws", names_to = NULL) %>% 
    mutate(Year = years[i])
}

# bind all years:
rcorr_ypreds <- bind_rows(rcorr_ypreds_list) # C-index draws (and other output) for each year


# plot:
library(ggplot2)
rcorr_ypreds %>% 
  filter(metric == "C Index") %>% 
  ggplot() + 
  geom_density(aes(x = draws, group = Year, col = as.factor(Year))) +
  scale_color_viridis_d() +
  xlab("C-index") +
  theme_bw() +
  theme(text = element_text(size = 18)) +
  ggtitle(paste0(spec, ", spatial C-index (y preds based, each line = 1 year)")) +
  xlim(c(0, 1))
# NAs if there is no observation (?)


# 2) based on predictions of occupancy probability:

# bind predictions for all folds:
occ_all_routes <- abind::abind(res_lists_folds[[1]]$occ_posterior,
                               res_lists_folds[[2]]$occ_posterior,
                               res_lists_folds[[3]]$occ_posterior,
                               res_lists_folds[[4]]$occ_posterior,
                               res_lists_folds[[5]]$occ_posterior,
                               along = 1)
dim(occ_all_routes)


# for all years:
rcorr_occprob_list <- vector(mode = "list", length = length(years))

for(i in 1:length(years)){
  
  print(i)
  
  # observations:
  obs_i <- occ_dt_spec %>% 
    filter(RTENO %in% test_RTENOs) %>% 
    # order same as RTENOs in test data:
    arrange(factor(RTENO, levels = test_RTENOs)) %>% # IS ORDER RIGHT?
    filter(Year == years[i]) %>% 
    pull(presence)
  
  # C index, predictions:
  C_distr_year_i <- apply(X = occ_all_routes[,i,], MARGIN = 2, FUN = Hmisc::rcorr.cens, x = obs_i, outx=FALSE)
  
  # reformat:
  rcorr_occprob_list[[i]] <-  C_distr_year_i %>% 
    as_tibble(rownames = "metric") %>% 
    tidyr::pivot_longer(cols = V1:V4000, values_to = "draws", names_to = NULL) %>% 
    mutate(Year = years[i])
}

# bind all years:
rcorr_occprob <- bind_rows(rcorr_occprob_list)

# plot:
rcorr_occprob %>% 
  filter(metric == "C Index") %>% 
  ggplot() + 
  geom_density(aes(x = draws, group = Year, col = as.factor(Year))) +
  scale_color_viridis_d() +
  xlab("C-index") +
  theme_bw() +
  theme(text = element_text(size = 18), legend.position = "bottom") +
  ggtitle(paste0(spec, ", spatial C-index (occ. prob. based, each line = 1 year)")) +
  xlim(c(0, 1))


# temporal C indices: ---

# do we get trends right?

# sum predicted occupancy across all routes at each MCMC jump for each year
# vs.
# sum observations across all routes for each year

# observations:
sum_obs_year <- occ_dt_spec %>% 
  filter(RTENO %in% test_RTENOs) %>% 
  # order same as RTENOs in test data:
  arrange(factor(RTENO, levels = test_RTENOs)) %>%
  group_by(Year) %>% 
  summarise(sum_obs_year = sum(presence, na.rm = TRUE)) %>% 
  pull(sum_obs_year)

# predicted y:
sum_ypreds_year <- apply(X = y_preds_routes, MARGIN = c(2,3), FUN = sum)
dim(sum_ypreds_year)

# C index, y:
C_distr_temp_y <- apply(X = sum_ypreds_year, MARGIN = 2, FUN = Hmisc::rcorr.cens, x = sum_obs_year, outx=FALSE) # TRUE? xx

# reformat:
rcorr_temp_y_preds <-  C_distr_temp_y %>% 
  as_tibble(rownames = "metric") %>% 
  tidyr::pivot_longer(cols = V1:V4000, values_to = "draws", names_to = NULL)

# plot:
rcorr_temp_y_preds %>% 
  filter(metric == "C Index") %>% 
  ggplot() + 
  geom_density(aes(x = draws)) +
  xlab("C-index") +
  theme_bw() +
  theme(text = element_text(size = 18)) +
  ggtitle(paste0(spec, ", temp. C-index (y preds based)")) +
  xlim(c(0, 1))


# same for comparing observations with occuupancy probability:

# predicted occ. prob.:
sum_occprob_year <- apply(X = occ_all_routes, MARGIN = c(2,3), FUN = sum)

# C index, occ. prob.:
C_distr_temp_occ <- apply(X = sum_occprob_year, MARGIN = 2, FUN = Hmisc::rcorr.cens, x = sum_obs_year, outx=FALSE) # TRUE? xx

# reformat:
rcorr_temp_occprob <-  C_distr_temp_occ %>% 
  as_tibble(rownames = "metric") %>% 
  tidyr::pivot_longer(cols = V1:V4000, values_to = "draws", names_to = NULL)

# plot:
rcorr_temp_occprob %>% 
  filter(metric == "C Index") %>% 
  ggplot() + 
  geom_density(aes(x = draws)) +
  xlab("C-index") +
  theme_bw() +
  theme(text = element_text(size = 18)) +
  ggtitle(paste0(spec, ", temp. C-index (occ. based)")) +
  xlim(c(0, 1))

# save evaluation outputs:
save(list(rcorr_ypreds, rcorr_occprob, rcorr_temp_y_preds, rcorr_temp_occprob), 
     file = file.path(dir, "data", paste0("CV_eval_", spec, ".RData")))


# Briscoe metrics: ---