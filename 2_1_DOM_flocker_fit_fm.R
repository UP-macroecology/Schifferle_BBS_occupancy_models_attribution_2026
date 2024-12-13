# fit DOM with
# - quadratic effects of bioclim covariates: bio1, bio2, bio3, bio7, bio14, bio15, spring, summer, autumn, winter prec.
# - quadratic effect of land use: urban areas, managed_pastures, primary_nonforests, secondary_nonforests, sum_annual_crops
# - p: different intercepts for route sections
# with flocker, normal priors

# with 750 km buffer
# compare with 250 km buffer

buffer_km <- 750 # 250


# packages: ----

library(dplyr)
library(sf)
library(doParallel)
library(flocker)
#install.packages("cmdstanr", repos = c("https://mc-stan.org/r-packages/", getOption("repos")))
library(cmdstanr)
set_cmdstan_path(path = NULL)#set_cmdstan_path("C:/Users/schifferle1/Documents/cmdstan-2.34.1") # xx

# register cores for parallel computation:
ncores <- 40 # models * 4 chains? 
cl <- makeCluster(ncores, setup_timeout = 0.5)
registerDoParallel(cl)

# directories: ----

print(tempdir())

dir <- file.path("/import", "ecoc9z", "data-zurell", "schifferle", "BBS_occupancy_models_2023")
#dir <- getwd()

# directory for logfiles:
log_dir <- file.path("logfiles", "fm_buffer750km")
# directory for results:
res_dir <- file.path(dir, "results", "fm_buffer750km")

# functions: -----

source("0_functions.R")


# load data: ----

# selected routes and focal years matched to environmental data:
# merged route, year, environment data:
load(file = file.path("data", "route_year_env_data.RData")) # route_sel_env_dt_final; output of 1_3_match_BBS_to_env_data.R 

# selected variables:
load(file = file.path("data", "selected_variables.RData")) # selvar_final; output of 1_2_variable_selection.R

# routes-years:
load(file = file.path("data", "BBS_for_occ_selection.RData")) # route_sel_dt; output of 1_3_match_BBS_to_env_data.R 

# selected routes spatial data (to buffer presences):
routes_sel_sf <- st_read(file.path("data", "route_selection_1995_2019_surv_beg_end_max_5y_miss_v2_spat_thin_100km_max_30_r_per_BCR_centroids.shp")) # output of 1_1_route_selection.R

# route-year-species information (only surveyed)
load(file = file.path("data", "BBS_for_occ_spec_records.RData")) # bbs_dt_occ; output of 1_0_reformat_BBS_data.R

# selected species:
#load(file = file.path("data", "final_species_selection.RData")) # species_selection_final; output of 1_2_species_selection.R
# sorted by ecoregion:
load(file = file.path("data", "final_species_selection_eco_sorted.RData")) # final_species_eco_sorted; output of 1_2_species_selection.R


# assemble overall data: ----


# scale covariates:

# save center and scale used:

env_means <- route_sel_env_dt_final %>% 
  select(-c(Latitude, Longitude, BCR, ObsN, doy)) %>% 
  summarise(across(!c(RTENO, Year, Surveyed), mean))

env_sds <- route_sel_env_dt_final %>% 
  select(-c(Latitude, Longitude, BCR, ObsN, doy)) %>% 
  summarise(across(!c(RTENO, Year, Surveyed), sd))

env_scale_pars <- list("center" = env_means, "scale" =  env_sds)

save(env_scale_pars, file = file.path("data", "route_env_dt_scale_pars.RData"))

# scale variables:

route_sel_env_dt_scaled <- route_sel_env_dt_final %>% 
  select(-c(Latitude, Longitude, BCR, ObsN, doy)) %>% 
  mutate(across(!c(RTENO, Year, Surveyed), ~ as.numeric(scale(., center=mean(.), scale = sd(.)))))

# # check:
# route_sel_env_dt_scaled %>% 
#   summarise(across(!c(RTENO, Year, Surveyed), ~ round(mean(.x), 5))) %>% 
#   as.numeric()
# route_sel_env_dt_scaled %>% 
#   summarise(across(!c(RTENO, Year, Surveyed), ~ round(sd(.x), 5))) %>% 
#   as.numeric()

rm(route_sel_env_dt_final)

nyears <- length(unique(route_sel_env_dt_scaled$Year)) # 25
nsurveys <- 5


# make blocks:

# 10 species per block
spec_blocks_list <- split(final_species_eco_sorted, c(rep(1:(floor(length(final_species_eco_sorted)/10)), each = 10), 
                                                      rep(ceiling(length(final_species_eco_sorted)/10), length(final_species_eco_sorted) %% 10))) # xx change once it runs
names(spec_blocks_list) <- NULL

# log overall progress:
prog_log_file <- file(file.path(log_dir, "fm_buffer_750_progress.txt"), open = "wt") # write console output here
sink(prog_log_file, type = "message")
sink(prog_log_file, type = "output")

for(i in 1:3){#length(spec_blocks_list)){

  print(paste("block", i, "of", length(spec_blocks_list)))
  
  foreach(spec = spec_blocks_list[[i]],
          .packages = c("dplyr", "collapse", "flocker", "cmdstanr", "brms", "sf"), # xx
          .errorhandling = "pass", #"remove",
          .verbose = TRUE) %dopar% {
            
            
            # check whether species has run already:
            model_run <- file.exists(file.path(res_dir, paste0("out_", spec, "_fm_buffer", buffer_km, ".RData")))
            post_proc_run <- file.exists(file.path(res_dir, paste0("postproc_", spec, "_fm_buffer", buffer_km, ".RData")))
            
            if(model_run & post_proc_run) {
              print(paste(spec, "ran already."))
              next
            }
            
            # assemble data: ----

            spec_fm_log_file <- file(file.path(log_dir, paste0("out_", spec, "_fm_buffer", buffer_km, ".txt")), open = "wt") # write console output here
            sink(spec_fm_log_file, type = "message")
            sink(spec_fm_log_file, type = "output")
            
            print(spec)
            
            # species presences-absences:
            
            occ_dt_spec <- BBS_pres_abs_spec(species = spec)
            
            # relevant routes, within distance of 750 km of species records:
            rel_routes <- training_routes(species = spec, buffer_km = buffer_km, output = "RTENOs")
            
            print(paste("training routes:", length(rel_routes)))
            
            # reformat obs. as array sites x surveys x years:
            
            years <- seq(min(occ_dt_spec$Year), max(occ_dt_spec$Year))
            nsites <- length(rel_routes)
            
            y_array <- array(NA, dim = c(nsites, nsurveys, nyears))
            for (t in 1:nyears){
              y_array[1:nsites, 1:nsurveys, t] <- as.matrix(occ_dt_spec[which(occ_dt_spec$Year == years[t] & occ_dt_spec$RTENO %in% rel_routes), 
                                                                        c(paste0("Count", seq(10, 50, 10)))])
            }
            
            # reformat environmental covariates:
            
            env_cov <- vector("list", length = nyears)
            for (t in 1:nyears){
              env_cov[[t]] <- route_sel_env_dt_scaled[which(route_sel_env_dt_scaled$Year == years[t] & route_sel_env_dt_scaled$RTENO %in% rel_routes), 
                                                      c(selvar_final, paste0(selvar_final, "_3yrs"))]
            }
            
            # covariate for detection probability:
            
            det_cov <- vector("list", length = 1)
            names(det_cov) <- "route_section"
            det_cov$route_section <- array(NA, dim = c(nsites, nsurveys, nyears))
            det_cov$route_section[ , , 1:nyears] <- matrix(rep(c("Sect1", "Sect2", "Sect3", "Sect4", "Sect5"), nsites), nsites, byrow = TRUE)
            
            # make flocker data:
            
            fd <- make_flocker_data_dynamic(
              obs = y_array,
              unit_covs = env_cov, 
              event_covs = det_cov, 
              quiet = TRUE
            )
            
            
            # fit model: ----
            

            start.time <- Sys.time()
            
            out <- flock(
              f_occ = ~ bio1_3yrs + bio2_3yrs + bio3_3yrs + bio7_3yrs + bio14_3yrs + bio15_3yrs +
                pr_mean_spring_3yrs + pr_mean_summer_3yrs + pr_mean_autumn_3yrs + pr_mean_winter_3yrs  + 
                I(bio1_3yrs^2) + I(bio2_3yrs^2) + I(bio3_3yrs^2) + I(bio7_3yrs^2) + I(bio14_3yrs^2) + I(bio15_3yrs^2) + 
                I(pr_mean_spring_3yrs^2) + I(pr_mean_summer_3yrs^2) + I(pr_mean_autumn_3yrs^2) + I(pr_mean_winter_3yrs^2)  +
                urbanareas_3yrs + managed_pastures_3yrs + primary_nonforests_3yrs + secondary_nonforests_3yrs + sum_annual_crops_3yrs +
                I(urbanareas_3yrs^2) + I(managed_pastures_3yrs^2) + I(primary_nonforests_3yrs^2) + I(secondary_nonforests_3yrs^2) + I(sum_annual_crops_3yrs^2),
              f_det = ~ route_section,
              f_col = ~ bio1 + bio2 + bio3 + bio7 + bio14 + bio15 +
                pr_mean_spring + pr_mean_summer + pr_mean_autumn + pr_mean_winter +
                I(bio1^2) + I(bio2^2) + I(bio3^2) + I(bio7^2)+ I(bio14^2)+ I(bio15^2) +
                I(pr_mean_spring^2) + I(pr_mean_summer^2) + I(pr_mean_autumn^2) + I(pr_mean_winter^2)  + 
                urbanareas + managed_pastures + primary_nonforests + secondary_nonforests + sum_annual_crops +
                I(urbanareas^2) + I(managed_pastures^2) + I(primary_nonforests^2) + I(secondary_nonforests^2) + I(sum_annual_crops^2),
              f_ex = ~ bio1 + bio2 + bio3 + bio7 + bio14 + bio15 +
                pr_mean_spring + pr_mean_summer + pr_mean_autumn + pr_mean_winter +
                I(bio1^2) + I(bio2^2) + I(bio3^2) + I(bio7^2)+ I(bio14^2)+ I(bio15^2) +
                I(pr_mean_spring^2) + I(pr_mean_summer^2) + I(pr_mean_autumn^2) + I(pr_mean_winter^2)  + 
                urbanareas + managed_pastures + primary_nonforests + secondary_nonforests + sum_annual_crops +
                I(urbanareas^2) + I(managed_pastures^2) + I(primary_nonforests^2) + I(secondary_nonforests^2) + I(sum_annual_crops^2),
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
              warmup = 1000,
              iter = 1000 + 1000
            )
            
            print(out)
            
            save(out, file = file.path(res_dir, paste0("out_", spec, "_fm_buffer", buffer_km, ".RData")))
            
            end.time <- Sys.time()
            print(round(end.time - start.time, 2))
            
            
            # postprocessing: ----
            
            print("postprocessing")
            
            start.time <- Sys.time()
            
            # calculate further results:
            fitted_initocc_col_ex <- fitted_flocker(out) 
            occ_posterior <- get_Z(out, history_condition = FALSE) # default
            y_predictions <- predict_flocker(out, history_condition = FALSE) # default, necessary for validation?
            loo_cv <- loo_flocker(out, thin = NULL)

            res_list <- list("fitted" = fitted_initocc_col_ex,
                             "occ_posterior" = occ_posterior,
                             "y_preds" = y_predictions,
                             "loo_cv" = loo_cv)
            
            print(res_list$loo_cv)
            
            high_k_routes <- loo::pareto_k_ids(loo_cv)
            high_k_RTENO <- rel_routes[high_k_routes]
            
            print(paste("RTENO with too high pareto k:", paste(high_k_RTENO, collapse = ", ")))
            
            rm(fitted_initocc_col_ex, occ_posterior, y_predictions, loo_cv)
            
            save(res_list, file = file.path(res_dir, paste0("postproc_", spec, "_fm_buffer", buffer_km, ".RData")))
            
            rm(res_list) 
            
            end.time <- Sys.time()
            print(round(end.time - start.time, 2))
            
            sink(type="message")
            sink(type="output")
            
          }
}
 
sink(type="message")
sink(type="output")

stopCluster(cl)

rm(list=ls())
gc()