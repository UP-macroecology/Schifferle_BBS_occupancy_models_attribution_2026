# fit DOM with
# - quadratic effects of bioclim covariates: bio1, bio2, bio3, spring, summer, autumn, winter prec., bio15
# - quadratic effect of land use (without perennial crops): sum annual crops, primn, secdn, pastr, urban
# - p: different intercepts for route sections
# with flocker, normal priors

# with 750 km buffer


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
dir <- file.path("/import", "ecoc9z", "data-zurell", "schifferle", "BBS_occupancy_models_2023")
#dir <- getwd()


# load data: ----

# selected routes and focal years matched to environmental data:
# merged route, year, environment data:
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

# register cores for parallel computation:
ncores <- 16 # models * 4 chains? 
cl <- makeCluster(ncores, setup_timeout = 0.5)
registerDoParallel(cl)

testspecs0 <- c("Black-billed Cuckoo",
                "Great-tailed Grackle",
               "Scissor-tailed Flycatcher",
                "Common Raven",
                "Wild Turkey",
                "Bald Eagle",
                "Black-chinned Hummingbird",
                "Broad-winged Hawk",
                "Greater Roadrunner",
                "Eurasian Collared-Dove"
)

testspecs1 <- c("Mourning Dove",
               "Red-winged Blackbird",
               "Brown-headed Cowbird",
               "Barn Swallow",
               "American Robin",
               "European Starling",
               "American Crow",
               "House Sparrow",
               "Common Grackle",
               "Common Yellowthroat")

testspecs2 <- c("Northern Flicker",
               "Blue Jay",
               "Eastern Meadowlark",
               "Indigo Bunting",
               "American Kestrel",
               "American Goldfinch",
               "Cedar Waxwing",
               "Northern Mockingbird",
               "Yellow Warbler",
               "Song Sparrow")

testspecs3 <- c("Black-capped Chickadee",
               "House Wren",
               "Eastern Towhee",
               "Eastern Kingbird",
               "Gray Catbird",
               "Northern Cardinal",
               "Horned Lark",
               "Dickcissel",
               "White-breasted Nuthatch",
               "Tufted Titmouse")

testspecs <- species_selection_final[which(!species_selection_final %in% c(testspecs0, testspecs1, testspecs2, testspecs3))]


foreach(spec = testspecs,
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
          
          
          # reformat obs. as array sites x surveys x years:
          
          years <- seq(min(occ_dt_spec$Year), max(occ_dt_spec$Year))
          nsites <- length(unique(routes_within$RTENO_BBS))
          
          y_array <- array(NA, dim = c(nsites, nsurveys, nyears))
          for (t in 1:nyears){
            y_array[1:nsites, 1:nsurveys, t] <- as.matrix(occ_dt_spec[which(occ_dt_spec$Year == years[t] & occ_dt_spec$RTENO %in% routes_within$RTENO_BBS), 
                                                                      c(paste0("Count", seq(10, 50, 10)))])
          }
          
          # reformat environmental covariates:
          
          env_cov <- vector("list", length = nyears)
          for (t in 1:nyears){
            env_cov[[t]] <- route_sel_env_dt_scaled[which(route_sel_env_dt_scaled$Year == years[t] & occ_dt_spec$RTENO %in% routes_within$RTENO_BBS), 
                                                    c("bio1", "bio2", "bio3", "pr_spring", "pr_summer","pr_autumn", 
                                                      "pr_winter", "bio15","bio1_3yrs", "bio2_3yrs", "bio3_3yrs", 
                                                      "pr_spring_3yrs", "pr_summer_3yrs", "pr_autumn_3yrs", "pr_winter_3yrs", 
                                                      "bio15_3yrs",
                                                      "sum_annual_crops", "secdn","pastr", "urban", "primn",
                                                      "sum_annual_crops_3yrs", "secdn_3yrs", "pastr_3yrs", "urban_3yrs", "primn_3yrs")]
            env_cov[[t]]$year <- as.character(years[t]) # factor
          }
          
          # covariate for detection probability:
          
          det_cov <- vector("list", length = 1)
          names(det_cov) <- "route_section"
          det_cov$route_section <- array(NA, dim = c(nsites, nsurveys, nyears))
          det_cov$route_section[ , , 1:nyears] <- matrix(rep(c("Sect1", "Sect2", "Sect3", "Sect4", "Sect5"), nsites), nsites, byrow = TRUE)
          
          det_cov$year <- array(NA, dim = c(nsites, nsurveys, nyears))
          for (t in 1:nyears){
            det_cov$year[ , , t] <- matrix(as.character(years[t]), nrow = nsites, ncol = nsurveys)
          }
          # make flocker data:
          
          fd <- make_flocker_data_dynamic(
            obs = y_array,
            unit_covs = env_cov, 
            event_covs = det_cov, 
            quiet = TRUE
          )
          
          
          # fit model: ----
          
          sink(paste0("out_fl_fm_buffer750_", spec, ".txt")) # write console output here
          
          sink(type = "message")
          
          print(spec)
          
          start.time <- Sys.time()
          
          out <- flock(
            f_occ = ~ bio1_3yrs + bio2_3yrs + bio3_3yrs + pr_spring_3yrs + pr_summer_3yrs + pr_autumn_3yrs + 
              pr_winter_3yrs + bio15_3yrs + I(bio1_3yrs^2) + I(bio2_3yrs^2) + I(bio3_3yrs^2) + 
              I(pr_spring_3yrs^2) + I(pr_summer_3yrs^2) + I(pr_autumn_3yrs^2) + I(pr_winter_3yrs^2) + I(bio15_3yrs^2) +
              sum_annual_crops_3yrs + secdn_3yrs + pastr_3yrs + urban_3yrs + primn_3yrs +
              I(sum_annual_crops_3yrs^2) + I(secdn_3yrs^2) + I(pastr_3yrs^2) + I(urban_3yrs^2) + I(primn_3yrs^2),
            f_det = ~ route_section,
            f_col = ~ bio1 + bio2 + bio3 + pr_spring + pr_summer + pr_autumn + pr_winter + bio15 +
              I(bio1^2) + I(bio2^2) + I(bio3^2) + I(pr_spring^2) + I(pr_summer^2) + I(pr_autumn^2) + 
              I(pr_winter^2) + I(bio15^2) + sum_annual_crops + secdn + pastr + urban + primn +
              I(sum_annual_crops^2) + I(secdn^2) + I(pastr^2) + I(urban^2) + I(primn^2),
            f_ex = ~ bio1 + bio2 + bio3 + pr_spring + pr_summer + pr_autumn + pr_winter + bio15 +
              I(bio1^2) + I(bio2^2) + I(bio3^2) + I(pr_spring^2) + I(pr_summer^2) + I(pr_autumn^2) + 
              I(pr_winter^2) + I(bio15^2) + sum_annual_crops + secdn + pastr + urban + primn +
              I(sum_annual_crops^2) + I(secdn^2) + I(pastr^2) + I(urban^2) + I(primn^2),
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
          
          print(out)
          
          save(out, file = file.path(dir, "data", paste0("out_fl_fm_buffer750_", spec, ".RData")))
          
          end.time <- Sys.time()
          time.taken <- round(end.time - start.time,2)
          print(time.taken)
          
          print("postprocessing")
          
          start.time <- Sys.time()
          
          # calculate further results:
          fitted_occ_col_ex <- fitted_flocker(out)
          occupancy_uncond <- get_Z(out, history_condition = FALSE) # default
          prediction_sites_uncond <- predict_flocker(out, history_condition = FALSE) # default, necessary for validation?
          loo_cv <- loo_flocker(out, thin = NULL)
          
          res_list <- list("fitted" = fitted_occ_col_ex,
                           "Zs" = occupancy_uncond,
                           "preds_occ_uncond" = prediction_sites_uncond,
                           "loo_cv" = loo_cv)
          
          print(res_list$loo_cv)
          
          rm(fitted_occ_col_ex, occupancy_uncond, prediction_sites_uncond, loo_cv)
          
          save(res_list, file = file.path(dir, "data", paste0("out_fl_fm_", spec, "_postproc_buffer750.RData")))
          
          end.time <- Sys.time()
          time.taken <- round(end.time - start.time,2)
          print(time.taken)
          
          sink(file = NULL)
          
        }


stopCluster(cl)

rm(list=ls())
gc()