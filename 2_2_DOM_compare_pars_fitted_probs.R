
# compare output for same dynamic occupancy model but different underlying data (buffer sizes):

# packages: --------------------------------------------------------------------

library(dplyr)
library(sf)
library(flocker)
library(cmdstanr)
#set_cmdstan_path(path = NULL)
set_cmdstan_path("C:/Users/schifferle1/Documents/cmdstan-2.34.1")
library(ggplot2)
#library(bayesplot)


# function to compare two model outputs: ---------------------------------------

# load  function:
source("0_functions.R")

# compare fitted values for occ., col., ex. prob.:

compare_fl_fitted <- function(model1, model2, env_data = route_sel_env_dt_scaled, 
                              value = mean, spec,
                              route_subset_RTENOs = NA,
                              buffer = FALSE){ # only compare fitted values for routes with these RTENOs
  
  # compare fitted values:
  
  # fitted values model1:
  preds_m1 <- fitted_flocker(model1,
                             new_data = as.data.frame(env_data), # without providing new data it runs much longer ->?
                             components = c("occ", "col", "ex"),
                             summarise = FALSE)
  # fitted values model2:
  preds_m2 <- fitted_flocker(model2,
                             new_data = as.data.frame(env_data),
                             components = c("occ", "col", "ex"))
  
  # extract single value of posterior predictive distribution, given as argument value from distribution:
  preds_m1_aggr <- env_data %>%
    select(RTENO, Year) %>%
    cbind(psi_init1 = apply(preds_m1$linpred_occ, FUN = value, MAR = 1)) %>%
    cbind(col1 = apply(preds_m1$linpred_col, FUN = value, MAR = 1)) %>% 
    cbind(ex1 = apply(preds_m1$linpred_ex, FUN = value, MAR = 1)) %>% 
    filter(if(any(is.na(route_subset_RTENOs))){RTENO == RTENO} else {RTENO %in% route_subset_RTENOs})

  preds_m2_aggr <- env_data %>% 
    select(RTENO, Year) %>% 
    cbind(psi_init2 = apply(preds_m2$linpred_occ, FUN = value, MAR = 1)) %>% 
    cbind(col2 = apply(preds_m2$linpred_col, FUN = value, MAR = 1)) %>% 
    cbind(ex2 = apply(preds_m2$linpred_ex, FUN = value, MAR = 1)) %>% 
    filter(if(any(is.na(route_subset_RTENOs))){RTENO == RTENO} else {RTENO %in% route_subset_RTENOs})
  
  # merge with spatial data on routes:
  fitted_sf <- routes_sel_sf %>%
    right_join(preds_m1_aggr, by = c(RTENO_BBS = "RTENO")) %>%
    right_join(preds_m2_aggr, by = c(RTENO_BBS = "RTENO", Year = "Year"))
  
  if(buffer == TRUE){
    pres_buff_small <- training_routes(species = spec, buffer_km = 250, output = "buffer")
    pres_buff_large <- training_routes(species = spec, buffer_km = 750, output = "buffer")
  }
  
  # differences:
  
  print(fitted_sf %>% 
          filter(Year == 1991) %>% 
          mutate(diff = psi_init1 - psi_init2) %>%
          ggplot() +
          geom_sf(aes(color = diff), size = 2) +
          scale_colour_gradient2(name = "m1-m2", low = "cornflowerblue", high = "red3") +
          ggtitle(paste(spec, "initial occupancy")) +
          geom_sf(data = pres_buff_small, fill = "transparent") +
          geom_sf(data = pres_buff_large, fill = "transparent") 
      )
  
  print(fitted_sf %>% 
          mutate(diff = col1 - col2) %>% 
          ggplot() +
          facet_wrap(~Year) +
          geom_sf(aes(color = diff), size = 2) +
          scale_colour_gradient2(name = "m1-m2", low = "cornflowerblue", high = "red3") +
          ggtitle(paste(spec, "col. prob.")) +
          geom_sf(data = pres_buff_small, fill = "transparent") +
          geom_sf(data = pres_buff_large, fill = "transparent")
        )
  
  print(fitted_sf %>% 
          mutate(diff = ex1 - ex2) %>% 
          ggplot() +
          facet_wrap(~Year) +
          geom_sf(aes(color = diff), size = 2) +
          scale_colour_gradient2(name = "m1-m2", low = "cornflowerblue", high = "red3") +
          ggtitle(paste(spec, "ex. prob.")) +
          geom_sf(data = pres_buff_small, fill = "transparent") +
          geom_sf(data = pres_buff_large, fill = "transparent")
        )
}


# compare parameter estimates of two models:

compare_fl_par_estimates <- function(model1, model2, spec){
  
  posterior1 <- as.array(model1)
  posterior2 <- as.array(model2)
  
  # posterior interval estimates from MCMC draws:
  summary1 <- bayesplot::mcmc_intervals_data(posterior1, point_est = "median", prob_outer = 0.95)
  summary1$model <- factor(1, levels = c(1,2))
  summary2 <- bayesplot::mcmc_intervals_data(posterior2, point_est = "median", prob_outer = 0.95)
  summary2$model <- factor(2, levels = c(1,2))
  
  combined <- rbind(summary1, summary2) %>% 
    filter(grepl("b_", parameter))

  # plots:
  
  # initial occupancy:
  combined_occ <- combined %>% 
    filter(grepl("occ", parameter))
  
  pos <- position_nudge(y = case_when(
    combined_occ$model == "1" ~ 0,
    combined_occ$model == "2" ~ -0.1,
    TRUE ~ 0))
  
  print(ggplot(combined_occ, aes(x = m, y = parameter, color = model)) + 
    geom_linerange(aes(xmin = l, xmax = h), position = pos, linewidth = 0.8)+
    geom_linerange(aes(xmin = ll, xmax = hh), position = pos, linewidth = 0.3) +
    geom_point(position = pos, color="black", size = 0.8) +
    ggtitle(spec) +
    geom_vline(xintercept = 0, linetype = "dashed", colour = "grey60")
  )
  
  # colonisation probability:
  
  combined_colo <- combined %>% 
    filter(grepl("colo", parameter))
  
  print(ggplot(combined_colo, aes(x = m, y = parameter, color = model)) + 
    geom_linerange(aes(xmin = l, xmax = h), position = pos, linewidth = 0.8)+
    geom_linerange(aes(xmin = ll, xmax = hh), position = pos, linewidth = 0.3)+
    geom_point(position = pos, color="black") +
    ggtitle(spec) +
    geom_vline(xintercept = 0, linetype = "dashed", colour = "grey60"))
  
  # extinction probability:
  
  combined_ex <- combined %>% 
    filter(grepl("ex", parameter))
  
  print(ggplot(combined_ex, aes(x = m, y = parameter, color = model)) + 
    geom_linerange(aes(xmin = l, xmax = h), position = pos, linewidth=0.8)+
    geom_linerange(aes(xmin = ll, xmax = hh), position = pos, linewidth = 0.3)+
    geom_point(position = pos, color="black") +
    ggtitle(spec) +
    geom_vline(xintercept = 0, linetype = "dashed", colour = "grey60"))
  
}


# function to compare posterior densities between two models (with same parameters):

compare_post_dens <- function(model1, model2, regex_pars, spec) {
  
  pars <- names(model1$fit)[grepl(pattern = regex_pars, names(model1$fit))]
  
  m1_draws <- data.frame(model = "m1", par = rep(pars, each = 4000), draws = NA)
  m2_draws <- data.frame(model = "m2", par = rep(pars, each = 4000), draws = NA)
  for(p in pars){
    m1_draws[which(m1_draws$par == p), 3] <- unlist(lapply(model1$fit@sim$samples, "[", p))
    m2_draws[which(m2_draws$par == p), 3] <- unlist(lapply(model2$fit@sim$samples, "[", p))
  }
  
  combined <- rbind(m1_draws, m2_draws) %>% 
    mutate(model = factor(model))
  
  print(ggplot(combined) + 
          facet_wrap(~par) +
          geom_density(aes(x = draws, color = model))  +
          theme_bw() +
          labs(x = "estimate") +
          ggtitle(spec) +
          theme(legend.position="bottom") +
          #geom_vline(aes(xintercept = median(draws)), linetype = "dashed") +
          #geom_vline(aes(xintercept = median(draws)), linetype = "dashed") +
          geom_vline(xintercept = 0, colour = "grey60")
  )
}



# compare two models: ----

## load data: ----

# merged route, year, environment data:
load(file = file.path("data", "route_year_env_data.RData"))
# scale covariates:
route_sel_env_dt_scaled <- route_sel_env_dt_final %>% 
  mutate(across(bio2:pr_winter_3yrs, ~ (scale(.)) %>% as.vector()))

# BBS route data:
load(file = file.path("data", "BBS_for_occ_selection.RData")) # route_sel_dt; output of 1_3_match_BBS_to_env_data.R

# selected routes spatial data:
routes_sel_sf <- st_read(file.path("data", "route_selection_1991_2015_surv_beg_end_max_5y_miss_v2_spat_thin_100km_max_30_r_per_BCR_centroids.shp")) # output of 1_1_route_selection.R

# BBS species data:
# route-year-species information (only surveyed routes)
load(file = file.path("data", "BBS_for_occ_spec_records.RData")) # bbs_dt_occ; output of 1_0_reformat_BBS_data.R

# selected species:
load(file = file.path("data", "final_species_selection.RData")) # output of 1_2_species_selection.R


testspecs <- c("Acadian Flycatcher", "Alder Flycatcher", "American Goldfinch", "American Kestrel", 
            "American Redstart", "Ash-throated Flycatcher", "Baltimore Oriole",
           "Bank Swallow", "Bell's Vireo", "Bewick's Wren", "Black-and-white Warbler", "Black-billed Magpie")


## run functions: ----

# load model output:
load(file = file.path("results", paste0("out_fl_fm_buffer250_", testspecs[s], "_update_preds.RData")))
model1 <- out
load(file = file.path("M:/Documents/DEBTs/analysis/Schifferle_BBS_occupancy_models_2023", "results", paste0("out_fl_fm_buffer750_", testspecs[s], "_update_preds.RData")))
model2 <- out

# check convergence:
#bayesplot::mcmc_combo(model2, pars = "b_occ_Ibio15_3yrsE2")
model1
model2

# compare fitted values:
# # all routes:
# compare_fl_fitted(model1, model2, value = median, spec = testspecs[s])
# # only routes in training data of both models:
# route_subset_RTENOs <- training_routes(species = testspecs[s], buffer_km = 250)
# compare_fl_fitted(model1, model2, value = median, spec = testspecs[s], route_subset_RTENOs = route_subset_RTENOs)

# all routes:
compare_fl_fitted(model1, model2, value = median, spec = testspecs[s], buffer = TRUE)

# compare parameter estimates (intervals):
compare_fl_par_estimates(model1 = model1, model2 = model2, spec = testspecs[s])

# compare parameter estimates (posterior densities):
compare_post_dens(model1 = model1, model2 = model2, regex_pars = "b_occ", spec = testspecs[s])
compare_post_dens(model1 = model1, model2 = model2, regex_pars = "b_col", spec = testspecs[s])
compare_post_dens(model1 = model1, model2 = model2, regex_pars = "b_ex", spec = testspecs[s])


#
results_dir <- file.path("//NAS-2-P-SN-01.ibb.uni-potsdam.de", "users$", "schifferle1", "Documents", 
                         "DEBTs", "analysis", "Schifferle_BBS_occupancy_models_2023", "results", "full_model",
                         "final_run_750km") 

BMIP_sel_specs <- read.csv("C:/Users/schifferle1/Documents/BMIP_data_subset/BMIP chosen birds.csv")$common.name
spec <- BMIP_sel_specs[4]
flocker_fitted <- load(file.path(results_dir, paste0("postproc_", spec, "_fm_buffer750.RData")))
str(flocker_fitted)
dim(res_list$fitted$linpred_det) # route, surveys, years, samples

# detection probability is modelled only based on route section, so it doesn't differ between routes and years:
mean_det_surv <- apply(res_list$fitted$linpred_det, c(2), mean, na.rm = TRUE)

det_prob_along_transect_df <- data.frame("species" = BMIP_sel_specs, "section1"= NA,
                                         "section2"= NA,
                                         "section3"= NA,
                                         "section4"= NA,
                                         "section5"= NA)
for(spec in BMIP_sel_specs[4:18]) {
  
  print(spec)
  
  skip_to_next <- FALSE
  tryCatch(print(load(file.path(results_dir, paste0("postproc_", spec, "_fm_buffer750.RData")))),
           error = function(e) { skip_to_next <<- TRUE})
  if(skip_to_next) { next }
  
  mean_det_surv <- apply(res_list$fitted$linpred_det, c(2), mean, na.rm = TRUE)
  
  print(mean_det_surv)
  
  det_prob_along_transect_df[which(det_prob_along_transect_df$species == spec), 2:6] <- mean_det_surv
  
  print(det_prob_along_transect_df)
  
  rm(res_list)
}

write.csv(det_prob_along_transect_df, file = file.path("results", "detection_probs_DOM.csv"), row.names = FALSE)

det_prob_along_transect_df %>% 
  rowwise() %>% 
  mutate(mean_det = mean(c_across(starts_with("section")), na.rm = TRUE)) %>% 
  mutate(max_det = max(c_across(starts_with("section")))) %>%
  mutate(min_det = min(c_across(starts_with("section")))) %>%
  arrange(mean_det) %>% 
  filter(!is.na(section1)) %>% 
  tidyr::pivot_longer(cols = section1:section5, values_to = "det_prob", names_to = "section") %>% 
  ggplot() +
  geom_line(aes(x = section, y = det_prob, group = species, colour = species))

# same across sites, years,  