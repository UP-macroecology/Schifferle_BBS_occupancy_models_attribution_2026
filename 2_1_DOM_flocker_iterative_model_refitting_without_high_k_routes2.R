# fit DOM with
# - quadratic effects of bioclim covariates: bio1, bio2, bio3, spring, summer, autumn, winter prec., bio15
# - quadratic effect of land use (without perennial crops): sum annual crops, primn, secdn, pastr, urban
# - detection probability p: different intercepts for route sections
# with flocker, normal priors

# but: refit model for each species without the routes for which pareto k values are too high
# based on loo-CV of the previous round of model fitting

# k should be below 1 (for now)

# note that in each model fitting run the routes used in the respective model are just numbered consecutively,
# matching them to the correct routes (route identifiers, RTENOs) is a bit tedious


# packages: ----

library(dplyr)
library(sf)
library(doParallel)
library(flocker)
#install.packages("cmdstanr", repos = c("https://mc-stan.org/r-packages/", getOption("repos")))
library(cmdstanr)
set_cmdstan_path(path = NULL)#set_cmdstan_path("C:/Users/schifferle1/Documents/cmdstan-2.34.1") # xx


# directories: ----

print(tempdir())
dir <- file.path("/import", "ecoc9", "data-zurell", "schifferle", "BBS_occupancy_models_2023") # ecoc9z
#dir <- getwd()


# load data: ----

# selected routes and focal years matched to environmental data:
load(file = file.path("data", "route_year_env_data.RData"))

# scale covariates:
route_sel_env_dt_scaled <- route_sel_env_dt_final %>% 
  mutate(across(bio2:pr_winter_3yrs, ~ (scale(.)) %>% as.vector()))

# route-year-species information (only surveyed)
load(file = file.path("data", "BBS_for_occ_spec_records.RData")) # output of 1_0_reformat_BBS_data.R

# selected species:
load(file = file.path("data", "final_species_selection.RData")) # output of 1_2_species_selection.R

# selected routes spatial data (to buffer presences):
routes_sel_sf <- st_read(file.path("data", "route_selection_1991_2015_surv_beg_end_max_5y_miss_v2_spat_thin_100km_max_30_r_per_BCR_centroids.shp")) # output of 1_1_route_selection.R


# assemble data: ----

nyears <- length(unique(route_sel_env_dt_final$Year)) # 25
nsurveys <- 5

# route RTENOs:
route_RTENOs_all <- matrix(route_sel_env_dt_final$RTENO, nrow = 476, ncol = nyears, byrow = TRUE)[,1]


# register cores for parallel computation:
ncores <- 12 # models * 4 chains? 
cl <- makeCluster(ncores, setup_timeout = 0.5)
registerDoParallel(cl)


testspecs <- species_selection_final

# make blocks:
# 3 species per block:
spec_blocks_list <- split(testspecs, rep(1:ceiling(length(testspecs)/3), each = 3))
names(spec_blocks_list) <- NULL

for(i in 1:length(spec_blocks_list)){
  
  print(paste("block", i, "of", length(spec_blocks_list)))
  
  foreach(spec = spec_blocks_list[[i]],
          .packages = c("dplyr", "collapse", "flocker", "cmdstanr", "brms", "sf"), # xx
          .errorhandling = "pass", #"remove",
          .verbose = TRUE) %dopar% {
            
            # assemble data: ----
            
            # species presences:
            
            presences_spec <- bbs_dt_occ %>% 
              select(c(English_Common_Name, RTENO, Year, paste0("Count", seq(10, 50, 10)))) %>% 
              filter(English_Common_Name == spec)
            
            # match to routes-year-env:
            
            occ_dt_spec <- route_sel_env_dt_scaled %>% 
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
            
            # routes with presences (sf):
            
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
              st_buffer(dist = 750000) %>% 
              st_union
            
            # routes within buffer:
            routes_within <- occ_spec_sf %>% 
              st_filter(., y = pres_buffer, join = st_within)
            
            # library(ggplot2)
            # ggplot(occ_spec_sf) +
            #   geom_sf(data = pres_buffer) +
            #   geom_sf(aes(colour = presence_summarised), size = 0.7) +
            #   geom_sf(data = routes_within, colour = "yellow", size = 0.5)
            
            # RTENOs within the buffer = considered for model at the beginning:
            route_RTENOs_buff <- unique(routes_within$RTENO_BBS)
            
            
            # reformat obs. as array sites x surveys x years:
            years <- seq(min(occ_dt_spec$Year), max(occ_dt_spec$Year))
            nsites <- length(unique(routes_within$RTENO_BBS))
            
            y_array <- array(NA, dim = c(nsites, nsurveys, nyears))
            for (t in 1:nyears){
              y_array[1:nsites, 1:nsurveys, t] <- as.matrix(occ_dt_spec[which(occ_dt_spec$Year == years[t] & occ_dt_spec$RTENO %in% routes_within$RTENO_BBS), 
                                                                        c(paste0("Count", seq(10, 50, 10)))])
            }
            
            # reformat environmental covariates:
            
            route_sel_env_dt_scaled
            env_cov <- vector("list", length = nyears)
            for (t in 1:nyears){
              env_cov[[t]] <- route_sel_env_dt_scaled[which(route_sel_env_dt_scaled$Year == years[t] & occ_dt_spec$RTENO %in% routes_within$RTENO_BBS), 
                                                      c("bio1", "bio2", "bio3", "bio7", "bio14", "bio15", 
                                                        "pr_spring", "pr_summer","pr_autumn", "pr_winter",
                                                        "bio1_3yrs", "bio2_3yrs", "bio3_3yrs", "bio7_3yrs", "bio14_3yrs", "bio15_3yrs",
                                                        "pr_spring_3yrs", "pr_summer_3yrs", "pr_autumn_3yrs", "pr_winter_3yrs", 
                                                        "sum_annual_crops", "secdf","pastr", "urban",
                                                        "sum_annual_crops_3yrs", "secdf_3yrs", "pastr_3yrs", "urban_3yrs")]
            }
            
            # covariate for detection probability:
            
            det_cov <- vector("list", length = 1)
            names(det_cov) <- "route_section"
            det_cov$route_section <- array(NA, dim = c(nsites, nsurveys, nyears))
            det_cov$route_section[ , , 1:nyears] <- matrix(rep(c("Sect1", "Sect2", "Sect3", "Sect4", "Sect5"), nsites), nsites, byrow = TRUE)
            
            
            # 1st model run: -----------------------------------------------------
            
            round <- 1
            
            sink(paste0("out_fl_fm_buffer750_", spec, "_round_", round, ".txt")) # write console output here
            sink(type = "message")
            
            print(spec)
            
            print(paste("iteration round", round))
            
            # make flocker data:
            
            fd <- make_flocker_data_dynamic(
              obs = y_array,
              unit_covs = env_cov, 
              event_covs = det_cov, 
              quiet = TRUE
            )
            
            # fit model: 
            
            start.time <- Sys.time()
            
            out <- flock(
              f_occ = ~ bio1_3yrs + bio2_3yrs + bio3_3yrs + bio7_3yrs + bio14_3yrs + bio15_3yrs +
                pr_spring_3yrs + pr_summer_3yrs + pr_autumn_3yrs + pr_winter_3yrs  + 
                I(bio1_3yrs^2) + I(bio2_3yrs^2) + I(bio3_3yrs^2) + I(bio7_3yrs^2) + I(bio14_3yrs^2) + I(bio15_3yrs^2) + 
                I(pr_spring_3yrs^2) + I(pr_summer_3yrs^2) + I(pr_autumn_3yrs^2) + I(pr_winter_3yrs^2)  +
                sum_annual_crops_3yrs + secdf_3yrs + pastr_3yrs + urban_3yrs +
                I(sum_annual_crops_3yrs^2) + I(secdf_3yrs^2) + I(pastr_3yrs^2) + I(urban_3yrs^2),
              f_det = ~ route_section,
              f_col = ~ bio1 + bio2 + bio3 + bio7 + bio14 + bio15 +
                pr_spring + pr_summer + pr_autumn + pr_winter +
                I(bio1^2) + I(bio2^2) + I(bio3^2) + I(bio7^2)+ I(bio14^2)+ I(bio15^2) +
                I(pr_spring^2) + I(pr_summer^2) + I(pr_autumn^2) + I(pr_winter^2)  + 
                sum_annual_crops + secdf + pastr + urban +
                I(sum_annual_crops^2) + I(secdf^2) + I(pastr^2) + I(urban^2),
              f_ex = ~ bio1 + bio2 + bio3 + bio7 + bio14 + bio15 +
                pr_spring + pr_summer + pr_autumn + pr_winter +
                I(bio1^2) + I(bio2^2) + I(bio3^2) + I(bio7^2)+ I(bio14^2)+ I(bio15^2) +
                I(pr_spring^2) + I(pr_summer^2) + I(pr_autumn^2) + I(pr_winter^2)  + 
                sum_annual_crops + secdf + pastr + urban +
                I(sum_annual_crops^2) + I(secdf^2) + I(pastr^2) + I(urban^2),
              flocker_data = fd,
              prior = c(brms::set_prior("logistic(0,1)", class = "Intercept") + # flat on probability scale (https://cran.r-project.org/web/packages/flocker/vignettes/flocker_tutorial.html)
                          brms::set_prior("logistic(0,1)", class = "Intercept", dpar = "occ"),
                        brms::set_prior("logistic(0,1)", class = "Intercept", dpar = "colo"),
                        brms::set_prior("logistic(0,1)", class = "Intercept", dpar = "ex"),
                        brms::set_prior("normal(0,2)", class = "b"),
                        brms::set_prior("normal(0,2)", dpar = "occ", class = "b"),
                        brms::set_prior("normal(0,2)", dpar = "colo", class = "b"),
                        brms::set_prior("normal(0,2)", dpar = "ex", class = "b")),
              multiseason = "colex",
              multi_init = "explicit",
              backend = "cmdstanr",
              cores = 4,
              chains = 4,
              warmup = 250,
              iter = 250 + 1000
            )
            
            rm(fd)
            
            print(out)
            
            save(out, file = file.path(dir, "data", paste0("out_fl_fm_buffer750_", spec, "_round_", round, ".RData")))
            
            end.time <- Sys.time()
            print(round(end.time - start.time, 2))
            
            print("postprocessing")
            
            start.time <- Sys.time()
            
            # calculate further results:
            fitted_occ_col_ex <- fitted_flocker(out)
            occupancy_uncond <- get_Z(out, history_condition = FALSE) # default
            prediction_sites_uncond <- predict_flocker(out, history_condition = FALSE) # default, necessary for validation?
            loo_cv <- loo_flocker(out, thin = NULL)
            
            # route RTENOs with too high pareto k values:
            
            #high_k_routes <- loo::pareto_k_ids(loo_cv)
            high_k_routes <- which(loo::pareto_k_values(loo_cv) >= 1) # xx
            high_k_RTENOs <- route_RTENOs_buff[high_k_routes]
            
            res_list <- list("fitted" = fitted_occ_col_ex,
                             "Zs" = occupancy_uncond,
                             "preds_occ_uncond" = prediction_sites_uncond,
                             "loo_cv" = loo_cv,
                             "high_k_RTENOs" = high_k_RTENOs)
            
            print(res_list$loo_cv)
            
            rm(fitted_occ_col_ex, occupancy_uncond, prediction_sites_uncond, loo_cv)
            
            save(res_list, file = file.path(dir, "data", paste0("out_fl_fm_", spec, "_postproc_buffer750_round_", round, ".RData")))
            
            end.time <- Sys.time()
            print(round(end.time - start.time, 2))
            
            sink(file = NULL)
            
            
            # further model runs: ------------------------------------------------
            
            for(round in 2:30){
              
              sink(paste0("out_fl_fm_buffer750_", spec, "_round_", round, ".txt")) # write console output here
              sink(type = "message")
              
              print(spec)
              
              print(paste("iteration round", round))
              
              load(file.path(dir, "data", paste0("out_fl_fm_", spec, "_postproc_buffer750_round_", round-1, ".RData")))
              #load(file.path(dir, "results", paste0("out_flocker_", spec, "_fitted_preds_loo", ".RData")))
              #res_list$high_k_RTENOs <- route_RTENOs_all[loo::pareto_k_ids(res_list$loo_cv)]
              
              # discard routes that have been identified as having too high pareto k values in previous rounds:
              #RTENOs_still_in <- route_RTENOs_all[which(!route_RTENOs_all %in% res_list$high_k_RTENOs)]
              RTENOs_still_in <- route_RTENOs_buff[which(!route_RTENOs_buff %in% res_list$high_k_RTENOs)] # which routes that are in the buffer have not been previously identified as having to high k values
              
              
              # check whether any routes are left:
              if(length(RTENOs_still_in) == 0){
                print("no routes left")
                break
              }
              
              print(paste("routes still in:", length(RTENOs_still_in), "of", length(route_RTENOs_buff)))
              
              # now we don't need RTENO (route identifier), but consecutive number again:
              low_k_routes <- which(route_RTENOs_buff %in% RTENOs_still_in)
              
              print(low_k_routes)
              
              y_array2 <- y_array[low_k_routes, ,]
              env_cov2 <- lapply(env_cov, function(x) x[low_k_routes,])
              det_cov2 <- lapply(det_cov, function(x) x[low_k_routes,,])
              
              # check whether enough routes with presences are left:
              n_routes_presence <- length(unique(which(y_array2 == 1, arr.ind = TRUE)[,1]))
              if(n_routes_presence < 5){
                print("less than 5 routes with presences left, stop here")
                break
              }
              
              # make flocker data:
              
              fd <- make_flocker_data_dynamic(
                obs = y_array2,
                unit_covs = env_cov2, 
                event_covs = det_cov2, 
                quiet = TRUE
              )
              
              
              # fit model: 
              
              start.time <- Sys.time()
              
              out <- flock(
                f_occ = ~ bio1_3yrs + bio2_3yrs + bio3_3yrs + bio7_3yrs + bio14_3yrs + bio15_3yrs +
                  pr_spring_3yrs + pr_summer_3yrs + pr_autumn_3yrs + pr_winter_3yrs  + 
                  I(bio1_3yrs^2) + I(bio2_3yrs^2) + I(bio3_3yrs^2) + I(bio7_3yrs^2) + I(bio14_3yrs^2) + I(bio15_3yrs^2) + 
                  I(pr_spring_3yrs^2) + I(pr_summer_3yrs^2) + I(pr_autumn_3yrs^2) + I(pr_winter_3yrs^2)  +
                  sum_annual_crops_3yrs + secdf_3yrs + pastr_3yrs + urban_3yrs +
                  I(sum_annual_crops_3yrs^2) + I(secdf_3yrs^2) + I(pastr_3yrs^2) + I(urban_3yrs^2),
                f_det = ~ route_section,
                f_col = ~ bio1 + bio2 + bio3 + bio7 + bio14 + bio15 +
                  pr_spring + pr_summer + pr_autumn + pr_winter +
                  I(bio1^2) + I(bio2^2) + I(bio3^2) + I(bio7^2)+ I(bio14^2)+ I(bio15^2) +
                  I(pr_spring^2) + I(pr_summer^2) + I(pr_autumn^2) + I(pr_winter^2)  + 
                  sum_annual_crops + secdf + pastr + urban +
                  I(sum_annual_crops^2) + I(secdf^2) + I(pastr^2) + I(urban^2),
                f_ex = ~ bio1 + bio2 + bio3 + bio7 + bio14 + bio15 +
                  pr_spring + pr_summer + pr_autumn + pr_winter +
                  I(bio1^2) + I(bio2^2) + I(bio3^2) + I(bio7^2)+ I(bio14^2)+ I(bio15^2) +
                  I(pr_spring^2) + I(pr_summer^2) + I(pr_autumn^2) + I(pr_winter^2)  + 
                  sum_annual_crops + secdf + pastr + urban +
                  I(sum_annual_crops^2) + I(secdf^2) + I(pastr^2) + I(urban^2),
                flocker_data = fd,
                prior = c(brms::set_prior("logistic(0,1)", class = "Intercept") + # flat on probability scale (https://cran.r-project.org/web/packages/flocker/vignettes/flocker_tutorial.html)
                          brms::set_prior("logistic(0,1)", class = "Intercept", dpar = "occ"),
                          brms::set_prior("logistic(0,1)", class = "Intercept", dpar = "colo"),
                          brms::set_prior("logistic(0,1)", class = "Intercept", dpar = "ex"),
                          brms::set_prior("normal(0,2)", class = "b"),
                          brms::set_prior("normal(0,2)", dpar = "occ", class = "b"),
                          brms::set_prior("normal(0,2)", dpar = "colo", class = "b"),
                          brms::set_prior("normal(0,2)", dpar = "ex", class = "b")),
                multiseason = "colex",
                multi_init = "explicit",
                backend = "cmdstanr",
                cores = 4,
                chains = 4,
                warmup = 250,
                iter = 250 + 1000
              )
              
              rm(fd)
              print(out)
              
              save(out, file = file.path(dir, "data", paste0("out_fl_fm_buffer750_", spec, "_round_", round, ".RData")))
              
              end.time <- Sys.time()
              print(round(end.time - start.time, 2))
              
              print("postprocessing")
              
              start.time <- Sys.time()
              
              # calculate further results:
              fitted_occ_col_ex <- fitted_flocker(out)
              occupancy_uncond <- get_Z(out, history_condition = FALSE) # default
              prediction_sites_uncond <- predict_flocker(out, history_condition = FALSE) # default, necessary for validation?
              loo_cv <- loo_flocker(out, thin = NULL)
              
              # route RTENOs with too high pareto k values:
              
              #high_k_routes <- loo::pareto_k_ids(loo_cv)
              high_k_routes <- which(loo::pareto_k_values(loo_cv) >= 1) # xx
              high_k_RTENOs <- RTENOs_still_in[high_k_routes]
              
              # check whether there are still too high k values:
              if(length(high_k_routes) == 0){
                print("no routes with pareto k values > 1 left")
                break
              }
              
              print(paste("routes with k >= 1 in current round:", paste(high_k_RTENOs, collapse = ",")))
              
              res_list <- list("fitted" = fitted_occ_col_ex,
                               "Zs" = occupancy_uncond,
                               "preds_occ_uncond" = prediction_sites_uncond,
                               "loo_cv" = loo_cv,
                               "high_k_RTENOs" = subset(route_RTENOs_buff, (!route_RTENOs_buff %in% RTENOs_still_in) | (route_RTENOs_buff %in% high_k_RTENOs))) # high_k_RTENOs now: all that have been discarded in previous rounds and in current round:
              
              print(res_list$loo_cv)
              print("routes now out:")
              print(res_list$high_k_RTENOs)
              
              rm(fitted_occ_col_ex, occupancy_uncond, prediction_sites_uncond, loo_cv)
              rm(out)
              
              save(res_list, file = file.path(dir, "data", paste0("out_fl_fm_", spec, "_postproc_buffer750_round_", round, ".RData")))
              rm(res_list)
              
              end.time <- Sys.time()
              print(round(end.time - start.time, 2))
              
              sink(file = NULL)
            }
            
          }
  
}


stopCluster(cl)

rm(list=ls())
gc()