# extract results from flocker DOM:
# - traceplots
# - fitted values for initial occupancy, colonisation and extinction probability
# - predictions for new (old?) data
# - Z values = occupancy over time?
# - loo validation


# packages: ----

library(dplyr)
library(doParallel)
library(flocker)
#install.packages("cmdstanr", repos = c("https://mc-stan.org/r-packages/", getOption("repos")))
library(cmdstanr)
#set_cmdstan_path(path = NULL)
set_cmdstan_path("C:/Users/schifferle1/Documents/cmdstan-2.34.1") # xx
library(sf)
library(ggplot2)

# load model outputs: ----

load(file.path("results", "American_Goldfinch_cl_lu_p_Stan_no_hs.RData"))
#model <- "AG_cl_lu_p_flocker_normprior"
out <- multi_colex_cl_lu_p_prior
# load(file.path("results", "AG_cl_lu_p_flock_horseshoe(df = 1, par_ratio = 0.01, df_global = 1, scale_slab = 2, df_slab = 4).RData"))
# model <- "AG_cl_lu_p_flocker_hsprior"
# out <- multi_colex_cl_lu_p_prior

# warmup tests:
# load(file.path("results", "AG_cl_lu_p_flock_warmup250.RData"))
# load(file.path("results", "AG_cl_lu_p_flock_warmup500.RData"))
# load(file.path("results", "AG_cl_lu_p_flock_warmup750.RData"))
# load(file.path("results", "AG_cl_lu_p_flock_warmup1250.RData"))




# check MCMC: ------------------------------------------------------------------

print(out)
# Bulk_ESS and Tail_ESS looks fine, Rhat as well

# traceplots:

plot(out)
brms::mcmc_plot(out, type = "trace")
brms::mcmc_plot(out, type = "trace", variable = "^b_" , regex = TRUE) # same
brms::mcmc_plot(out, type = "trace", variable = "^b_occ" , regex = TRUE)
brms::mcmc_plot(out, type = "trace", variable = "^b_colo" , regex = TRUE)
brms::mcmc_plot(out, type = "trace", variable = "^b_ex" , regex = TRUE)
# look fine

# effective sampling size:

brms::mcmc_plot(out, type = "neff") # documentation says between 0.1 and 0.5 is good, larger is high
# for 3 parameters neff is smaller than half the number of sampling iterations
# more problems with horseshoe
# Rhat:

brms::mcmc_plot(out, type = "rhat")
# looks fine

# -> MCMC seems to have worked okay


# fitted values for initial occ, col., ext., occ. prob: ------------------------

# get values from model:

# all here: ===
# fitted_occ_col_ex <- fitted_flocker(out)
# occupancy_uncond <- get_Z(out, history_condition = FALSE) # default
# prediction_sites_uncond <- predict_flocker(out, history_condition = FALSE) # default, necessary for validation?
# loo_cv <- loo_flocker(out, thin = NULL)
# 
# res_list <- list("fitted" = fitted_occ_col_ex,
#                  "Zs" = occupancy_uncond,
#                  "preds_occ_uncond" = prediction_sites_uncond,
#                  "loo_cv" = loo_cv)
# save(res_list, file = file.path("results", paste0(model, "_fitted_preds_loo.RData")))
load(file.path("results", "AG_cl_lu_p_flocker_normpriorfitted_preds_loo.RData"))

fitted_occ_col_ex <- res_list$fitted
occupancy_uncond <- res_list$Zs
prediction_sites_uncond <- res_list$preds_od_uncond
loo_cv <- res_list$loo_cv

#===



## fitted values for initial occ, col., ext.:
#fitted_occ_col_ex <- fitted_flocker(out)
str(fitted_occ_col_ex) # list with 4 elements: occ, col, ext, det
fitted_occ_col_ex[[1]][1,1,1,] # [route, survey, year, draw] 4000 draws
## predictions for occupancy in each year:
#occupancy_uncond <- get_Z(out, history_condition = FALSE) # default
str(occupancy_uncond) # occ. probability for each site and year

# predictions for old data:----
#prediction_sites_uncond <- predict_flocker(out, history_condition = FALSE) # default, necessary for validation?
str(prediction_sites_uncond) # 0 and 1s, for each site, survey and year
# predictions for new data:


# get observations:

# routes:
routes_sel_sf <- sf::st_read(file.path("data", "route_selection_1991_2015_surv_beg_end_max_5y_miss_v2_spat_thin_100km_max_30_r_per_BCR_centroids.shp")) # output of 1_1_route_selection.R

# route-year-species information (only surveyed)
load(file = file.path("data", "BBS_for_occ_spec_records.RData")) # output of 1_0_reformat_BBS_data.R

spec <- "American Goldfinch"

presences_spec <- bbs_dt_occ %>% 
  select(c(English_Common_Name, RTENO, Year, paste0("Count", seq(10, 50, 10)))) %>% 
  filter(English_Common_Name == spec)

# match to routes-year-env:
occ_dt_spec <- route_sel_env_dt_final %>% 
  # add observations:
  collapse::join(presences_spec, on = c("RTENO", "Year"), how = "left") %>% 
  # if route was surveyed but species not observed, replace NA with 0:
  mutate(across(Count10:Count50, ~ 
                  case_when(Surveyed == 1 & is.na(.) ~ 0,
                            .default = .))) %>%
  # convert bird counts to presence / absence:
  mutate(across(Count10:Count50, ~ 
                  case_when(. > 1 ~ 1,
                            .default = .))) %>% 
  # sum presence/absence on route:
  mutate(occ_route = rowSums(across(Count10:Count50))) %>%
  # convert presence absence:
  mutate(occ_route = ifelse(occ_route > 0, 1, 0))


# initial occupancy:

# match site number to route:
load(file = file.path("data", "route_year_env_data.RData")) # merged route, year, environment data
route_sel_env_dt_final
nsites <- length(unique(route_sel_env_dt_final$RTENO))
nyears <- length(unique(route_sel_env_dt_final$Year))
route_nrs <- matrix(route_sel_env_dt_final$RTENO, nrow = nsites, ncol = nyears, byrow = TRUE)[,1]

# fitted values:
route_preds <- data.frame(RTENO = route_nrs)
# mean:
route_preds$mean_psi_init <- apply(fitted_occ_col_ex$linpred_occ[, 1, 1, ], FUN = mean, MAR = 1) # all routes, first survey, first year, all draws
# standard deviation:
route_preds$sd_psi_init <- apply(fitted_occ_col_ex$linpred_occ[, 1, 1, ], FUN = sd, MAR = 1)
# merge with spatial data on routes:
psi_init_sf <- routes_sel_sf %>% 
  left_join(route_preds, by = c(RTENO_BBS = "RTENO"))

# observations in year 1:
psi_init_obs <- occ_dt_spec %>% 
  filter(Year == 1991) %>% 
  select(RTENO, occ_route)
psi_init_obs_sf <- routes_sel_sf %>% 
  left_join(psi_init_obs, by = c(RTENO_BBS = "RTENO"))


# mean estimated initial occupancy:
occ_colors <- c("1" = "black", "0" = "transparent")
ggplot() +
  geom_sf(data = psi_init_obs_sf[which(psi_init_obs_sf$occ_route == 1),], # only routes where species was observed
          aes(fill = as.character(occ_route)), size = 3, shape = 21) +
  geom_sf(data = psi_init_sf, aes(color = mean_psi_init), size = 2) +
  scale_fill_manual(values = occ_colors, name = "observation") +
  scale_color_viridis_c(name = "estimate") +
  ggtitle("mean initial occupancy")

# sd of estimated initial occupancy:
ggplot() +
  geom_sf(data = psi_init_obs_sf[which(psi_init_obs_sf$occ_route == 1),], # only routes where species was observed
          aes(fill = as.character(occ_route)), size = 3, shape = 21) +
  geom_sf(data = psi_init_sf, aes(color = sd_psi_init), size = 2) +
  scale_fill_manual(values = occ_colors, name = "observation") +
  scale_color_viridis_c(name = "estimate") +
  ggtitle("sd initial occupancy")


# maps of colonisation, extinction and occupancy probability for each year:

# prepare data:

## observations:
routes_years_obs <- occ_dt_spec %>% 
  select(RTENO, Year, occ_route)
routes_years_obs_sf <- routes_sel_sf %>% 
  left_join(routes_years_obs, by = c(RTENO_BBS = "RTENO")) %>% 
  # calculate raw dynamics:
  group_by(RTENO_BBS) %>% 
  mutate(raw_dyn = occ_route - lag(occ_route)) %>% # -1: raw extinction, 1 = raw colonisation
  ungroup

years <- unique(preds_routes_years$Year)

## model results:
# mean and sd for col, ext, occ:
preds_routes_years <- route_sel_env_dt_final[, c("RTENO", "Year")]
preds_routes_years$mean_col <- NA
preds_routes_years$sd_col <- NA
preds_routes_years$mean_ex <- NA
preds_routes_years$sd_ex <- NA
preds_routes_years$mean_occ_prob <- NA
preds_routes_years$sd_occ_prob <- NA

for(i in 1:nyears){
  print(i)
  preds_routes_years$mean_col[which(preds_routes_years$Year == years[i])] <- apply(fitted_occ_col_ex$linpred_col[,1,i,], FUN = mean, MARGIN = 1) # not for all sites was colonisation probability estimated
  preds_routes_years$sd_col[which(preds_routes_years$Year == years[i])] <- apply(fitted_occ_col_ex$linpred_col[,1,i,], FUN = sd, MARGIN = 1)
  preds_routes_years$mean_ex[which(preds_routes_years$Year == years[i])] <- apply(fitted_occ_col_ex$linpred_ex[,1,i,], FUN = mean, MARGIN = 1)
  preds_routes_years$sd_ex[which(preds_routes_years$Year == years[i])] <- apply(fitted_occ_col_ex$linpred_ex[,1,i,], FUN = sd, MARGIN = 1)
  preds_routes_years$mean_occ_prob[which(preds_routes_years$Year == years[i])] <- apply(occupancy_uncond[, i,], FUN = mean, MARGIN = 1)
  preds_routes_years$sd_occ_prob[which(preds_routes_years$Year == years[i])] <- apply(occupancy_uncond[, i,], FUN = sd, MARGIN = 1)
  
}
# join to spatial data:
col_ex_occ_sf <- routes_sel_sf %>% 
  right_join(preds_routes_years, by = c(RTENO_BBS = "RTENO"))


# maps:

## colonisation probability:

ggplot() +
  # geom_sf(data = routes_years_obs_sf[which(routes_years_obs_sf$raw_dyn == 1),], # raw colonisations
  #         aes(fill = as.character(occ_route)), size = 2, shape = 21) +
  # or: observation of species instead of raw colonisations to examine potential colonisation credits:
  # geom_sf(data = routes_years_obs_sf[which(routes_years_obs_sf$occ_route == 1),], # raw observations
  #         aes(fill = as.character(occ_route)), size = 2, shape = 21) +
  # colonisation prob.:
  geom_sf(data = col_ex_occ_sf, aes(fill = mean_col), size = 3, shape = 21) +
  # or better occupancy predictions to account for imperfect detection:
  geom_sf(data = col_ex_occ_sf, aes(color = mean_occ_prob), size = 1) +
  scale_color_viridis_c(name = "mean occ. prob.") +
  scale_fill_gradient(low = "white", high = "black", name = "mean col. prob.") +
  #scale_fill_manual(values = occ_colors, name = "observation") +
  facet_wrap(~Year) +
  ggtitle("mean col. prob.")
# highest modelled colonisation prob. not where colonisations were observed

ggplot() +
  geom_sf(data = routes_years_obs_sf[which(routes_years_obs_sf$raw_dyn == 1),], # raw colonisations
          aes(fill = as.character(occ_route)), size = 2, shape = 21) +
  geom_sf(data = col_ex_occ_sf, aes(color = sd_col), size = 1) +
  scale_color_viridis_c(name = "estimate") +
  scale_fill_manual(values = occ_colors, name = "observation") +
  facet_wrap(~Year) +
  ggtitle("sd col. prob.")

## extinction probability:

ggplot() +
  # geom_sf(data = routes_years_obs_sf[which(routes_years_obs_sf$raw_dyn == -1),], # raw extinctions
  #         aes(fill = as.character(occ_route)), size = 2, shape = 21) +
  # or: observation of species instead of raw extinctions to examine potential extinction debts:
  # geom_sf(data = routes_years_obs_sf[which(routes_years_obs_sf$occ_route == 1),], # raw observations
  #         aes(fill = as.character(occ_route)), size = 2, shape = 21) +
  # extinction prob.:
  geom_sf(data = col_ex_occ_sf, aes(fill = mean_ex), size = 3, shape = 21) +
  # or better occupancy predictions to account for imperfect detection:
  geom_sf(data = col_ex_occ_sf, aes(color = mean_occ_prob), size = 1) +
  scale_color_viridis_c(name = "mean occ. prob.") +
  scale_fill_gradient(low = "white", high = "black", name = "mean ex. prob.") +
  #scale_fill_manual(values = occ_colors, name = "observation") +
  facet_wrap(~Year) +
  ggtitle("mean ex. prob.")
# highest modelled extinction prob. not where extinctions were observed

ggplot() +
  geom_sf(data = routes_years_obs_sf[which(routes_years_obs_sf$raw_dyn == -1),], # raw extinctions
          aes(fill = as.character(occ_route)), size = 2, shape = 21) +
  geom_sf(data = col_ex_occ_sf, aes(color = sd_ex), size = 1) +
  scale_color_viridis_c(name = "estimate") +
  scale_fill_manual(values = occ_colors, name = "observation") +
  facet_wrap(~Year) +
  ggtitle("sd ex. prob.")

## occupancy:

ggplot() +
  geom_sf(data = routes_years_obs_sf[which(routes_years_obs_sf$occ_route == 1),], # raw observations
          aes(fill = as.character(occ_route)), size = 2, shape = 21) +
  geom_sf(data = col_ex_occ_sf, aes(color = mean_occ_prob), size = 1) +
  scale_color_viridis_c(name = "estimate") +
  scale_fill_manual(values = occ_colors, name = "observation") +
  facet_wrap(~Year) +
  ggtitle("mean occ. prob.")

ggplot() +
  geom_sf(data = routes_years_obs_sf[which(routes_years_obs_sf$occ_route == 1),], # raw observations
          aes(fill = as.character(occ_route)), size = 2, shape = 21) +
  geom_sf(data = col_ex_occ_sf, aes(color = sd_occ_prob), size = 1) +
  scale_color_viridis_c(name = "estimate") +
  scale_fill_manual(values = occ_colors, name = "observation") +
  facet_wrap(~Year) +
  ggtitle("sd occ. prob.")





# loo cross validation: --------------------------------------------------------

#loo_cv <- loo_flocker(out, thin = NULL) # here I can change thinning
# for which is it problematic:
loo::pareto_k_table(loo_cv)
loo::pareto_k_ids(loo_cv)
loo::pareto_k_influence_values(loo_cv)
plot(loo_cv)
plot(loo_cv, diagnostic = "n_eff")

# which routes are these:
plot(st_geometry(routes_sel_sf))
plot(st_geometry(routes_sel_sf %>% filter(RTENO_BBS %in% route_nrs[loo::pareto_k_ids(loo_cv)])), 
     add = TRUE, col = "red", pch = 19)
# these are influential observations regarding the posterior distribution, cause overfitting risk
# pareto-smoothing: approximation for these is unreliable
# If leaving out an observation changes the posterior too much then importance sampling is not able to give a reliable estimate

str(loo_cv) # p_loo = "effective number of parameters"
loo_cv$pointwise[loo::pareto_k_ids(loo_cv),]
# https://mc-stan.org/loo/reference/loo-glossary.html
loo_cv$pointwise %>%  View
brms::variables(out)

loo_cv$diagnostics$n_eff[loo::pareto_k_ids(loo_cv)] # effective sample size is very small for most values with high k values
summary(loo_cv$diagnostics$n_eff)
order(loo_cv$diagnostics$n_eff)
# will increasing the number of samples help?

# -> leave-one-series-out cross validation approximated with importance sampling had issues!

# save processed flocker output: -----------------------------------------------

# results as list:
res_list <- list("fitted" = fitted_occ_col_ex,
                 "Zs" = occupancy_uncond,
                 "preds_od_uncond" = prediction_sites_uncond,
                 "loo_cv" = loo_cv)
save(res_list, file = file.path("results", paste0(model, "_fitted_preds_loo.RData")))
                                                  
                                                  
# compare different priors with loo: ----

load(file.path("results", "AG_cl_lu_p_flocker_logistic(0,2)_normal(0,2)_fitted_preds_loo.RData"))
res_list$loo_cv

loo_compare_flocker()

load(file.path("results", "AG_cl_lu_p_flocker_logistic(0,1)_normal(0,3)_fitted_preds_loo.RData"))
res_list$loo_cv
load(file.path("results", "AG_cl_lu_p_flocker_logistic(0,3)_normal(0,2)_fitted_preds_loo.RData"))
res_list$loo_cv
load(file.path("results", "AG_cl_lu_p_flocker_logistic(0,1)_normal(0,1)_fitted_preds_loo.RData"))
res_list$loo_cv
load(file.path("results", "AG_cl_lu_p_flocker_logistic(0,2)_horseshoe(df = 3, par_ratio = 0.3, df_global = 1, scale_slab = 2, df_slab = 4)_fitted_preds_loo.RData"))
res_list$loo_cv
load(file.path("results", "AG_cl_lu_p_flocker_normpriorfitted_preds_loo.RData"))
res_list$loo_cv
load(file.path("results", "AG_cl_lu_p_flocker_hsprior_fitted_preds_loo.RData"))
res_list$loo_cv

load(file.path("results", "AG_cl_lu_p_flock_logistic(0,2)_normal(0,2).RData"))
out_norm_prior <- out
load(file.path("results", "AG_cl_lu_p_flock_logistic(0,2)_horseshoe(df = 3, par_ratio = 0.3, df_global = 1, scale_slab = 2, df_slab = 4).RData"))
out_hs_prior <- out

test <- loo_compare_flocker(list(out_norm_prior, out_hs_prior))



# misc: ----

# # divergent transitions:
# 
# posterior_cp <- as.array(fit_cp) # extract posterior draws for later use
# np_cp <- nuts_params(fit_cp)
# mcmc_parcoord(posterior_cp, np = np_cp)


## posterior predictive checks:
# plot observed data vs. replicated data based on model:
# doesn't work for repeated surveys as in Gabry et al. 2019 (otherwise it would probably have been in flocker documentation)
# library(bayesplot)
# ppc_dens_overlay(y = y_array, # xx must be vector, doesn't work for repeated surveys!
#                  yrep = prediction_sites_uncond) # use output of predict_flocker for posterior predictive checking
