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
library(terra)

# load model outputs: ----

load(file.path("results", "AG_cl_lu_p_flock_logistic(0,1)_normal(0,3).RData"))
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

# different species:
testspecs <- c("Black-billed Cuckoo",
               "Great-tailed Grackle",
               "Scissor-tailed Flycatcher",
               "Common Raven",
               "Wild Turkey",
               "Bald Eagle",
               "Black-chinned Hummingbird",
               "Broad-winged Hawk",
               "Greater Roadrunner",
               "Eurasian Collared-Dove")

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
spec <- testspecs[10]

# full model:
load(file.path("results", paste0("out_flocker_", spec, "_buffer.RData")))
out_buff1000 <- out
load(file.path("results", paste0("out_flocker_", spec, "_buffer750.RData")))
out_buff750 <- out
load(file.path("results", paste0("out_flocker_", spec, ".RData")))
out_no_buff <- out

load(file.path("results", paste0("out_flocker_", spec, "_fitted_preds_loo.RData")))
res_list$loo_cv

load(file.path("results", paste0("out_flocker_", spec, "_fitted_preds_loo_buffer750.RData")))
res_list$loo_cv


load(file.path("results", paste0("out_flocker_", spec, ".RData")))
out_1000 <- out
load(file.path("results", paste0("out_flocker_", spec, "_it_10000_2.RData")))
out_10k <- out
load(file.path("results", paste0("out_flocker_", spec, "_fitted_preds_loo_it_10000_2.RData")))
res_list$loo_cv

# compare predictor sets:
load(file.path("results", paste0("out_flocker_", spec, "_buffer750.RData")))
out_old_preds <- out
load(file.path("results", paste0("out_fl_fm_buffer750_", spec, "_update_preds.RData")))
out_new_preds <- out

load(file.path("results", paste0("out_flocker_", spec, ".RData")))
out1 <- out
load(file.path("results", paste0("out_flocker_", spec, "_round2.RData")))
out2 <- out
load(file.path("results", paste0("out_flocker_", spec, "_round3.RData")))
out3 <- out
load(file.path("results", paste0("out_flocker_", spec, "_round4.RData")))
out4 <- out
load(file.path("results", paste0("out_flocker_", spec, "_round5.RData")))
out5 <- out
load(file.path("results", paste0("out_flocker_", spec, "_round6.RData")))
out6 <- out
load(file.path("results", paste0("out_flocker_", spec, "_round7.RData")))
out7 <- out
load(file.path("results", paste0("out_flocker_", spec, "_round8.RData")))
out8 <- out
load(file.path("results", paste0("out_flocker_", spec, "_round9.RData")))
out9 <- out
load(file.path("results", paste0("out_flocker_", spec, "_round10.RData")))
out10 <- out
load(file.path("results", paste0("out_flocker_", spec, "_round11.RData")))
out11 <- out
load(file.path("results", paste0("out_flocker_", spec, "_round12.RData")))
out12 <- out
load(file.path("results", paste0("out_flocker_", spec, "_round13.RData")))
out13 <- out
load(file.path("results", paste0("out_flocker_", spec, "_round14.RData")))
out14 <- out
load(file.path("results", paste0("out_flocker_", spec, "_round15.RData")))
out15 <- out


load(file.path("results", paste0("out_flocker_", spec, "_fitted_preds_loo_round15.RData")))
res_list$loo_cv

r_high_k <- c(13,12,23,6,8,22,11,21,8,29,6,15,7,13,6,8,7,16,5,13, 29,13,8,14,21,9,16,6,12,11)
summary(r_high_k)

# detection intercept only:
load(file.path("results", paste0("out_flocker_det_int", spec, ".RData")))
out_det_int <- out

load(file.path("results", paste0("out_flocker_", spec, "_fitted_preds_loo_det_int.RData")))
cv_det_int <- res_list$loo_cv

# without land use:
load(file.path("results", paste0("out_flocker_no_lu", spec, ".RData")))
out_no_lu <- out

load(file.path("results", paste0("out_flocker_", spec, "_fitted_preds_loo_no_lu.RData")))
cv_no_lu <- res_list$loo_cv

# without climate:
load(file.path("results", paste0("out_flocker_no_cl", spec, ".RData")))
out_no_cl <- out

load(file.path("results", paste0("out_flocker_", spec, "_fitted_preds_loo_no_cl.RData")))
cv_no_cl <- res_list$loo_cv

# PSIS-loo model comparison:
load(file.path("results", paste0("loo_comp_", spec, ".RData")))
loo_comp


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


# effective sampling size:

brms::mcmc_plot(out, type = "neff") # documentation says between 0.1 and 0.5 is good, larger is high
# for 3 parameters neff is smaller than half the number of sampling iterations

# Rhat:
brms::mcmc_plot(out, type = "rhat")


# -> MCMC seems to have worked okay

# plot posteriors: ----
library(ggplot2)
brms::variables(out)
posterior <- as.array(out)


bayesplot::mcmc_intervals(posterior, pars = c("b_occ_bio1_3yrs"
                                              ))
bayesplot::mcmc_intervals(posterior, regex_pars = c("^b_occ"))
bayesplot::mcmc_intervals(posterior, regex_pars = c("^b_colo")) +
  scale_y_discrete(labels = "test") +
  ggtitle("colonisation probability")

bayesplot::mcmc_areas(posterior, regex_pars = c("^b_colo")) +
  scale_y_discrete(labels = c("Intercept", "annual temp.", 
                              "diurnal range", "isothermality", "spring prec.",
                              "summer prec.", "autumn prec.", 
                              "winter prec.", "prec. seasonality",
                              "(annual temp.)^2", "(diurnal range)^2",
                              "(isothermality)^2",
                              "(spring prec.)^2", "(summer prec.)^2", 
                              "(autumn prec.)^2", "(winter prec.)^2", 
                              "(prec. seasonality)^2",
                              "annual crops", "sec. non-forest", "pasture",
                              "urban", "prim. non-forest", "(annual crops)^2",
                              "(sec. non-forest)^2", "(pasture)^2", "(urban)^2",
                              "(prim. non-forest)^2")) +
  ggtitle("colonization probability")

bayesplot::mcmc_intervals(posterior, regex_pars = c("^b_colo"),
                          point_size = 2) +
  scale_y_discrete(labels = c("b_colo_Intercept" = "Intercept",
                              "b_colo_bio1" = "annual temp.", 
                              "b_colo_bio2" = "diurnal range",
                              "b_colo_bio3" = "isothermality", 
                              "b_colo_pr_spring" = "spring prec.",
                              "b_colo_pr_summer" = "summer prec.", 
                              "b_colo_pr_autumn" = "autumn prec.", 
                              "b_colo_pr_winter" = "winter prec.", 
                              "b_colo_bio15" = "prec. seasonality",
                              "b_colo_Ibio1E2" = "(annual temp.)^2",
                              "b_colo_Ibio2E2" = "(diurnal range)^2",
                              "b_colo_Ibio3E2" = "(isothermality)^2",
                              "b_colo_Ipr_springE2" = "(spring prec.)^2",
                              "b_colo_Ipr_summerE2" = "(summer prec.)^2",
                              "b_colo_Ipr_autumnE2" = "(autumn prec.)^2",
                              "b_colo_Ipr_winterE2" = "(winter prec.)^2",
                              "b_colo_Ibio15E2" = "(prec. seasonality)^2",
                              "b_colo_sum_annual_crops" = "annual crops",
                              "b_colo_secdn" = "sec. non-forest",
                              "b_colo_pastr" = "pasture",
                              "b_colo_urban" = "urban",
                              "b_colo_primn" = "prim. non-forest",
                              "b_colo_Isum_annual_cropsE2" = "(annual crops)^2",
                              "b_colo_IsecdnE2" = "(sec. non-forest)^2",
                              "b_colo_IpastrE2" = "(pasture)^2",
                              "b_colo_IurbanE2" = "(urban)^2",
                              "b_colo_IprimnE2" = "(prim. non-forest)^2"),
                   
                   limits = c("b_colo_IprimnE2",
                              "b_colo_primn",
                              "b_colo_IurbanE2",
                              "b_colo_urban",
                              "b_colo_IpastrE2",
                              "b_colo_pastr",
                              "b_colo_IsecdnE2",
                              "b_colo_secdn",
                              "b_colo_Isum_annual_cropsE2",
                              "b_colo_sum_annual_crops",
                              "b_colo_Ibio15E2",
                              "b_colo_bio15",
                              "b_colo_Ipr_winterE2",
                              "b_colo_pr_winter",
                              "b_colo_Ipr_autumnE2",
                              "b_colo_pr_autumn",
                              "b_colo_Ipr_summerE2",
                              "b_colo_pr_summer",
                              "b_colo_Ipr_springE2",
                              "b_colo_pr_spring",
                              "b_colo_Ibio3E2",
                              "b_colo_bio3",
                              "b_colo_Ibio2E2",
                              "b_colo_bio2",
                              "b_colo_Ibio1E2",
                              "b_colo_bio1")
                              ) +
  xlim(c(-0.9, 0.9)) + 
  ggtitle("colonization probability")


bayesplot::mcmc_dens(posterior, regex_pars = c("^b_colo"))

bayesplot::mcmc_areas(posterior, pars = c("b_colo_bio1")) +
  scale_y_discrete(labels = "test")
bayesplot::mcmc_dens(posterior, pars = c("b_colo_bio1"))

## compare posteriors of two models: ----

library(bayesplot)
posterior1 <- as.array(out1)
posterior2 <- as.array(out2)
posterior3 <- as.array(out3)
posterior4 <- as.array(out4)
posterior5 <- as.array(out5)
posterior6 <- as.array(out6)
posterior7 <- as.array(out7)
posterior8 <- as.array(out8)
posterior9 <- as.array(out9)
posterior10 <- as.array(out10)
posterior11 <- as.array(out11)
posterior12 <- as.array(out12)
posterior13 <- as.array(out13)
posterior14 <- as.array(out14)
posterior15 <- as.array(out15)

posterior1 <- as.array(out_no_buff)
posterior2 <- as.array(out_buff1000)
posterior3 <- as.array(out_buff750)

posterior1 <- as.array(out_old_preds)
posterior2 <- as.array(out_new_preds)

summary1 <- mcmc_intervals_data(posterior1)
summary1$model <- factor(1, levels = c(1,2))
summary2 <- mcmc_intervals_data(posterior2)
summary2$model <- factor(2, levels = c(1,2))

combined <- rbind(summary1, summary2)

# combined <- rbind(mcmc_intervals_data(posterior1), mcmc_intervals_data(posterior2)#,
#                   # mcmc_intervals_data(posterior3), mcmc_intervals_data(posterior4), 
#                   # mcmc_intervals_data(posterior5), mcmc_intervals_data(posterior6),
#                   # mcmc_intervals_data(posterior7), mcmc_intervals_data(posterior8),
#                   # mcmc_intervals_data(posterior9), mcmc_intervals_data(posterior10),
#                   # mcmc_intervals_data(posterior11), mcmc_intervals_data(posterior12),
#                   # mcmc_intervals_data(posterior13), mcmc_intervals_data(posterior14),
#                   # mcmc_intervals_data(posterior15)
#                   )
# combined$model <- rep(as.character(1:15), each = (nrow(combined)/15))
# #combined$model <- rep(as.character(1:8), each = (nrow(combined)/8))
# #combined$model <- rep(as.character(1:3), each = (nrow(combined)/3))
# #combined$model <- rep(as.character(1:2), each = (nrow(combined)/2))
combined <- combined %>% 
  filter(!parameter %in% c("lp__", "lprior"))

# make the plot using ggplot 
library(ggplot2)
theme_set(bayesplot::theme_default())

combined_occ <- combined %>% 
  filter(grepl("occ", parameter))
# pos <- position_nudge(y = case_when(
#   combined_occ$model == "1" ~ 0,
#   combined_occ$model == "2" ~ -0.06,
#   combined_occ$model == "3" ~ -0.06*2,
#   combined_occ$model == "4" ~ -0.06*3,
#   combined_occ$model == "5" ~ -0.06*4,
#   combined_occ$model == "6" ~ -0.06*5,
#   combined_occ$model == "7" ~ -0.06*6,
#   combined_occ$model == "8" ~ -0.06*7,
#   combined_occ$model == "9" ~ -0.06*8,
#   combined_occ$model == "10" ~ -0.06*9,
#   combined_occ$model == "11" ~ -0.06*10,
#   combined_occ$model == "12" ~ -0.06*11,
#   combined_occ$model == "13" ~ -0.06*12,
#   combined_occ$model == "14" ~ -0.06*13,
#   combined_occ$model == "15" ~ -0.06*14,
#   TRUE ~ 0))

pos <- position_nudge(y = case_when(
  combined_occ$model == "1" ~ 0,
  combined_occ$model == "2" ~ -0.1,
  TRUE ~ 0))
                        
ggplot(combined_occ, aes(x = m, y = parameter, color = model)) + 
  geom_linerange(aes(xmin = l, xmax = h), position = pos, linewidth = 0.8)+
  geom_linerange(aes(xmin = ll, xmax = hh), position = pos, linewidth = 0.3) +
  geom_point(position = pos, color="black", size = 0.8) +
  ggtitle(spec) +
  geom_vline(xintercept = 0, linetype = "dashed", colour = "grey60")

combined_colo <- combined %>% 
  filter(grepl("colo", parameter))
ggplot(combined_colo, aes(x = m, y = parameter, color = model)) + 
  geom_linerange(aes(xmin = l, xmax = h), position = pos, linewidth = 0.8)+
  geom_linerange(aes(xmin = ll, xmax = hh), position = pos, linewidth = 0.3)+
  geom_point(position = pos, color="black") +
  ggtitle(spec) +
  geom_vline(xintercept = 0, linetype = "dashed", colour = "grey60")

combined_ex <- combined %>% 
  filter(grepl("ex", parameter))
ggplot(combined_ex, aes(x = m, y = parameter, color = model)) + 
  geom_linerange(aes(xmin = l, xmax = h), position = pos, linewidth=0.8)+
  geom_linerange(aes(xmin = ll, xmax = hh), position = pos, linewidth = 0.3)+
  geom_point(position = pos, color="black") +
  ggtitle(spec) +
  geom_vline(xintercept = 0, linetype = "dashed", colour = "grey60")


# check loo cross validation / outliers: ---------------------------------------

# read routes to see which are problematic:
routes_sel_sf <- sf::st_read(file.path("data", "route_selection_1991_2015_surv_beg_end_max_5y_miss_v2_spat_thin_100km_max_30_r_per_BCR_centroids.shp")) # output of 1_1_route_selection.R
# match site number in model to route:
load(file = file.path("data", "route_year_env_data.RData")) # merged route, year, environment data
route_sel_env_dt_final
nsites <- length(unique(route_sel_env_dt_final$RTENO))
nyears <- length(unique(route_sel_env_dt_final$Year))
route_nrs <- matrix(route_sel_env_dt_final$RTENO, nrow = nsites, ncol = nyears, byrow = TRUE)[,1]

res_list$loo_cv

# for which is it problematic:
loo::pareto_k_table(res_list$loo_cv)
loo::pareto_k_ids(res_list$loo_cv) # routes with Pareto k estimates above threshold

loo::pareto_k_influence_values(res_list$loo_cv)
plot(res_list$loo_cv)
plot(res_list$loo_cv, diagnostic = "n_eff")




## routes with too high pareto k values, when iteratively refitting model, each time discarding routes with too high values: ----

load(file.path("results", paste0("out_flocker_", spec, "_fitted_preds_loo.RData")))
probl_routes1 <- loo::pareto_k_ids(res_list$loo_cv)
load(file.path("results", paste0("out_flocker_", spec, "_fitted_preds_loo_round2.RData")))
probl_routes2 <- loo::pareto_k_ids(res_list$loo_cv)
load(file.path("results", paste0("out_flocker_", spec, "_fitted_preds_loo_round3.RData")))
probl_routes3 <- loo::pareto_k_ids(res_list$loo_cv)
load(file.path("results", paste0("out_flocker_", spec, "_fitted_preds_loo_round4.RData")))
probl_routes4 <- loo::pareto_k_ids(res_list$loo_cv)
load(file.path("results", paste0("out_flocker_", spec, "_fitted_preds_loo_round5.RData")))
probl_routes5 <- loo::pareto_k_ids(res_list$loo_cv)
load(file.path("results", paste0("out_flocker_", spec, "_fitted_preds_loo_round6.RData")))
probl_routes6 <- loo::pareto_k_ids(res_list$loo_cv)
load(file.path("results", paste0("out_flocker_", spec, "_fitted_preds_loo_round7.RData")))
probl_routes7 <- loo::pareto_k_ids(res_list$loo_cv)
load(file.path("results", paste0("out_flocker_", spec, "_fitted_preds_loo_round8.RData")))
probl_routes8 <- loo::pareto_k_ids(res_list$loo_cv)
load(file.path("results", paste0("out_flocker_", spec, "_fitted_preds_loo_round9.RData")))
probl_routes9 <- loo::pareto_k_ids(res_list$loo_cv)
load(file.path("results", paste0("out_flocker_", spec, "_fitted_preds_loo_round10.RData")))
probl_routes10 <- loo::pareto_k_ids(res_list$loo_cv)
load(file.path("results", paste0("out_flocker_", spec, "_fitted_preds_loo_round11.RData")))
probl_routes11 <- loo::pareto_k_ids(res_list$loo_cv)
load(file.path("results", paste0("out_flocker_", spec, "_fitted_preds_loo_round12.RData")))
probl_routes12 <- loo::pareto_k_ids(res_list$loo_cv)
load(file.path("results", paste0("out_flocker_", spec, "_fitted_preds_loo_round13.RData")))
probl_routes13 <- loo::pareto_k_ids(res_list$loo_cv)
load(file.path("results", paste0("out_flocker_", spec, "_fitted_preds_loo_round14.RData")))
probl_routes14 <- loo::pareto_k_ids(res_list$loo_cv)
load(file.path("results", paste0("out_flocker_", spec, "_fitted_preds_loo_round15.RData")))
probl_routes15 <- loo::pareto_k_ids(res_list$loo_cv)

plot(st_geometry(routes_sel_sf), main = paste0(spec, " round 1"))
routes_sel_sf %>% 
  filter(RTENO_BBS %in% route_nrs[probl_routes1]) %>% 
  st_geometry() %>% 
  plot(., add = TRUE, col = "red", pch = 19)

plot(st_geometry(routes_sel_sf), main = paste0(spec, " round 2"))
routes_sel_sf %>% 
  filter(RTENO_BBS %in% route_nrs[probl_routes1]) %>% 
  st_geometry() %>% 
  plot(., add = TRUE, col = "yellow", pch = 19)
routes_sel_sf %>% 
  filter(RTENO_BBS %in% route_nrs[-probl_routes1][probl_routes2]) %>% 
  st_geometry() %>% 
  plot(., add = TRUE, col = "red", pch = 19)


plot(st_geometry(routes_sel_sf), main = paste0(spec, " round 3"))
routes_sel_sf %>% 
  filter(RTENO_BBS %in% route_nrs[probl_routes1]) %>% 
  st_geometry() %>% 
  plot(., add = TRUE, col = "yellow", pch = 19)
routes_sel_sf %>% 
  filter(RTENO_BBS %in% route_nrs[-probl_routes1][probl_routes2]) %>% 
  st_geometry() %>% 
  plot(., add = TRUE, col = "yellow", pch = 19)
routes_sel_sf %>% 
  filter(RTENO_BBS %in% route_nrs[-probl_routes1][-probl_routes2][probl_routes3]) %>% 
  st_geometry() %>% 
  plot(., add = TRUE, col = "red", pch = 19)

plot(st_geometry(routes_sel_sf), main = paste0(spec, " round 4"))
routes_sel_sf %>% 
  filter(RTENO_BBS %in% route_nrs[probl_routes1]) %>% 
  st_geometry() %>% 
  plot(., add = TRUE, col = "yellow", pch = 19)
routes_sel_sf %>% 
  filter(RTENO_BBS %in% route_nrs[-probl_routes1][probl_routes2]) %>% 
  st_geometry() %>% 
  plot(., add = TRUE, col = "yellow", pch = 19)
routes_sel_sf %>% 
  filter(RTENO_BBS %in% route_nrs[-probl_routes1][-probl_routes2][probl_routes3]) %>% 
  st_geometry() %>% 
  plot(., add = TRUE, col = "yellow", pch = 19)
routes_sel_sf %>% 
  filter(RTENO_BBS %in% route_nrs[-probl_routes1][-probl_routes2][-probl_routes3][probl_routes4]) %>% 
  st_geometry() %>% 
  plot(., add = TRUE, col = "red", pch = 19)

plot(st_geometry(routes_sel_sf), main = paste0(spec, " round 5"))
routes_sel_sf %>% 
  filter(RTENO_BBS %in% route_nrs[probl_routes1]) %>% 
  st_geometry() %>% 
  plot(., add = TRUE, col = "yellow", pch = 19)
routes_sel_sf %>% 
  filter(RTENO_BBS %in% route_nrs[-probl_routes1][probl_routes2]) %>% 
  st_geometry() %>% 
  plot(., add = TRUE, col = "yellow", pch = 19)
routes_sel_sf %>% 
  filter(RTENO_BBS %in% route_nrs[-probl_routes1][-probl_routes2][probl_routes3]) %>% 
  st_geometry() %>% 
  plot(., add = TRUE, col = "yellow", pch = 19)
routes_sel_sf %>% 
  filter(RTENO_BBS %in% route_nrs[-probl_routes1][-probl_routes2][-probl_routes3][probl_routes4]) %>% 
  st_geometry() %>% 
  plot(., add = TRUE, col = "yellow", pch = 19)
routes_sel_sf %>% 
  filter(RTENO_BBS %in% route_nrs[-probl_routes1][-probl_routes2][-probl_routes3][-probl_routes4][probl_routes5]) %>% 
  st_geometry() %>% 
  plot(., add = TRUE, col = "red", pch = 19)


plot(st_geometry(routes_sel_sf), main = paste0(spec, " round 6"))
routes_sel_sf %>% 
  filter(RTENO_BBS %in% route_nrs[probl_routes1]) %>% 
  st_geometry() %>% 
  plot(., add = TRUE, col = "yellow", pch = 19)
routes_sel_sf %>% 
  filter(RTENO_BBS %in% route_nrs[-probl_routes1][probl_routes2]) %>% 
  st_geometry() %>% 
  plot(., add = TRUE, col = "yellow", pch = 19)
routes_sel_sf %>% 
  filter(RTENO_BBS %in% route_nrs[-probl_routes1][-probl_routes2][probl_routes3]) %>% 
  st_geometry() %>% 
  plot(., add = TRUE, col = "yellow", pch = 19)
routes_sel_sf %>% 
  filter(RTENO_BBS %in% route_nrs[-probl_routes1][-probl_routes2][-probl_routes3][probl_routes4]) %>% 
  st_geometry() %>% 
  plot(., add = TRUE, col = "yellow", pch = 19)
routes_sel_sf %>% 
  filter(RTENO_BBS %in% route_nrs[-probl_routes1][-probl_routes2][-probl_routes3][-probl_routes4][probl_routes5]) %>% 
  st_geometry() %>% 
  plot(., add = TRUE, col = "yellow", pch = 19)
routes_sel_sf %>% 
  filter(RTENO_BBS %in% route_nrs[-probl_routes1][-probl_routes2][-probl_routes3][-probl_routes4][-probl_routes5][probl_routes6]) %>% 
  st_geometry() %>% 
  plot(., add = TRUE, col = "red", pch = 19)

plot(st_geometry(routes_sel_sf), main = paste0(spec, " round 7"))
routes_sel_sf %>% 
  filter(RTENO_BBS %in% route_nrs[probl_routes1]) %>% 
  st_geometry() %>% 
  plot(., add = TRUE, col = "yellow", pch = 19)
routes_sel_sf %>% 
  filter(RTENO_BBS %in% route_nrs[-probl_routes1][probl_routes2]) %>% 
  st_geometry() %>% 
  plot(., add = TRUE, col = "yellow", pch = 19)
routes_sel_sf %>% 
  filter(RTENO_BBS %in% route_nrs[-probl_routes1][-probl_routes2][probl_routes3]) %>% 
  st_geometry() %>% 
  plot(., add = TRUE, col = "yellow", pch = 19)
routes_sel_sf %>% 
  filter(RTENO_BBS %in% route_nrs[-probl_routes1][-probl_routes2][-probl_routes3][probl_routes4]) %>% 
  st_geometry() %>% 
  plot(., add = TRUE, col = "yellow", pch = 19)
routes_sel_sf %>% 
  filter(RTENO_BBS %in% route_nrs[-probl_routes1][-probl_routes2][-probl_routes3][-probl_routes4][probl_routes5]) %>% 
  st_geometry() %>% 
  plot(., add = TRUE, col = "yellow", pch = 19)
routes_sel_sf %>% 
  filter(RTENO_BBS %in% route_nrs[-probl_routes1][-probl_routes2][-probl_routes3][-probl_routes4][-probl_routes5][probl_routes6]) %>% 
  st_geometry() %>% 
  plot(., add = TRUE, col = "yellow", pch = 19)
routes_sel_sf %>% 
  filter(RTENO_BBS %in% route_nrs[-probl_routes1][-probl_routes2][-probl_routes3][-probl_routes4][-probl_routes5][-probl_routes6][probl_routes7]) %>% 
  st_geometry() %>% 
  plot(., add = TRUE, col = "red", pch = 19)

plot(st_geometry(routes_sel_sf), main = paste0(spec, " round 8"))
routes_sel_sf %>% 
  filter(RTENO_BBS %in% route_nrs[probl_routes1]) %>% 
  st_geometry() %>% 
  plot(., add = TRUE, col = "yellow", pch = 19)
routes_sel_sf %>% 
  filter(RTENO_BBS %in% route_nrs[-probl_routes1][probl_routes2]) %>% 
  st_geometry() %>% 
  plot(., add = TRUE, col = "yellow", pch = 19)
routes_sel_sf %>% 
  filter(RTENO_BBS %in% route_nrs[-probl_routes1][-probl_routes2][probl_routes3]) %>% 
  st_geometry() %>% 
  plot(., add = TRUE, col = "yellow", pch = 19)
routes_sel_sf %>% 
  filter(RTENO_BBS %in% route_nrs[-probl_routes1][-probl_routes2][-probl_routes3][probl_routes4]) %>% 
  st_geometry() %>% 
  plot(., add = TRUE, col = "yellow", pch = 19)
routes_sel_sf %>% 
  filter(RTENO_BBS %in% route_nrs[-probl_routes1][-probl_routes2][-probl_routes3][-probl_routes4][probl_routes5]) %>% 
  st_geometry() %>% 
  plot(., add = TRUE, col = "yellow", pch = 19)
routes_sel_sf %>% 
  filter(RTENO_BBS %in% route_nrs[-probl_routes1][-probl_routes2][-probl_routes3][-probl_routes4][-probl_routes5][probl_routes6]) %>% 
  st_geometry() %>% 
  plot(., add = TRUE, col = "yellow", pch = 19)
routes_sel_sf %>% 
  filter(RTENO_BBS %in% route_nrs[-probl_routes1][-probl_routes2][-probl_routes3][-probl_routes4][-probl_routes5][-probl_routes6][probl_routes7]) %>% 
  st_geometry() %>% 
  plot(., add = TRUE, col = "yellow", pch = 19)
routes_sel_sf %>% 
  filter(RTENO_BBS %in% route_nrs[-probl_routes1][-probl_routes2][-probl_routes3][-probl_routes4][-probl_routes5][-probl_routes6][-probl_routes7][probl_routes8]) %>% 
  st_geometry() %>% 
  plot(., add = TRUE, col = "red", pch = 19)

plot(st_geometry(routes_sel_sf), main = paste0(spec, " round 9"))
routes_sel_sf %>% 
  filter(RTENO_BBS %in% route_nrs[probl_routes1]) %>% 
  st_geometry() %>% 
  plot(., add = TRUE, col = "yellow", pch = 19)
routes_sel_sf %>% 
  filter(RTENO_BBS %in% route_nrs[-probl_routes1][probl_routes2]) %>% 
  st_geometry() %>% 
  plot(., add = TRUE, col = "yellow", pch = 19)
routes_sel_sf %>% 
  filter(RTENO_BBS %in% route_nrs[-probl_routes1][-probl_routes2][probl_routes3]) %>% 
  st_geometry() %>% 
  plot(., add = TRUE, col = "yellow", pch = 19)
routes_sel_sf %>% 
  filter(RTENO_BBS %in% route_nrs[-probl_routes1][-probl_routes2][-probl_routes3][probl_routes4]) %>% 
  st_geometry() %>% 
  plot(., add = TRUE, col = "yellow", pch = 19)
routes_sel_sf %>% 
  filter(RTENO_BBS %in% route_nrs[-probl_routes1][-probl_routes2][-probl_routes3][-probl_routes4][probl_routes5]) %>% 
  st_geometry() %>% 
  plot(., add = TRUE, col = "yellow", pch = 19)
routes_sel_sf %>% 
  filter(RTENO_BBS %in% route_nrs[-probl_routes1][-probl_routes2][-probl_routes3][-probl_routes4][-probl_routes5][probl_routes6]) %>% 
  st_geometry() %>% 
  plot(., add = TRUE, col = "yellow", pch = 19)
routes_sel_sf %>% 
  filter(RTENO_BBS %in% route_nrs[-probl_routes1][-probl_routes2][-probl_routes3][-probl_routes4][-probl_routes5][-probl_routes6][probl_routes7]) %>% 
  st_geometry() %>% 
  plot(., add = TRUE, col = "yellow", pch = 19)
routes_sel_sf %>% 
  filter(RTENO_BBS %in% route_nrs[-probl_routes1][-probl_routes2][-probl_routes3][-probl_routes4][-probl_routes5][-probl_routes6][-probl_routes7][probl_routes8]) %>% 
  st_geometry() %>% 
  plot(., add = TRUE, col = "yellow", pch = 19)
routes_sel_sf %>% 
  filter(RTENO_BBS %in% route_nrs[-probl_routes1][-probl_routes2][-probl_routes3][-probl_routes4][-probl_routes5][-probl_routes6][-probl_routes7][-probl_routes8][probl_routes9]) %>% 
  st_geometry() %>% 
  plot(., add = TRUE, col = "red", pch = 19)


plot(st_geometry(routes_sel_sf), main = paste0(spec, " round 10"))
routes_sel_sf %>% 
  filter(RTENO_BBS %in% route_nrs[probl_routes1]) %>% 
  st_geometry() %>% 
  plot(., add = TRUE, col = "yellow", pch = 19)
routes_sel_sf %>% 
  filter(RTENO_BBS %in% route_nrs[-probl_routes1][probl_routes2]) %>% 
  st_geometry() %>% 
  plot(., add = TRUE, col = "yellow", pch = 19)
routes_sel_sf %>% 
  filter(RTENO_BBS %in% route_nrs[-probl_routes1][-probl_routes2][probl_routes3]) %>% 
  st_geometry() %>% 
  plot(., add = TRUE, col = "yellow", pch = 19)
routes_sel_sf %>% 
  filter(RTENO_BBS %in% route_nrs[-probl_routes1][-probl_routes2][-probl_routes3][probl_routes4]) %>% 
  st_geometry() %>% 
  plot(., add = TRUE, col = "yellow", pch = 19)
routes_sel_sf %>% 
  filter(RTENO_BBS %in% route_nrs[-probl_routes1][-probl_routes2][-probl_routes3][-probl_routes4][probl_routes5]) %>% 
  st_geometry() %>% 
  plot(., add = TRUE, col = "yellow", pch = 19)
routes_sel_sf %>% 
  filter(RTENO_BBS %in% route_nrs[-probl_routes1][-probl_routes2][-probl_routes3][-probl_routes4][-probl_routes5][probl_routes6]) %>% 
  st_geometry() %>% 
  plot(., add = TRUE, col = "yellow", pch = 19)
routes_sel_sf %>% 
  filter(RTENO_BBS %in% route_nrs[-probl_routes1][-probl_routes2][-probl_routes3][-probl_routes4][-probl_routes5][-probl_routes6][probl_routes7]) %>% 
  st_geometry() %>% 
  plot(., add = TRUE, col = "yellow", pch = 19)
routes_sel_sf %>% 
  filter(RTENO_BBS %in% route_nrs[-probl_routes1][-probl_routes2][-probl_routes3][-probl_routes4][-probl_routes5][-probl_routes6][-probl_routes7][probl_routes8]) %>% 
  st_geometry() %>% 
  plot(., add = TRUE, col = "yellow", pch = 19)
routes_sel_sf %>% 
  filter(RTENO_BBS %in% route_nrs[-probl_routes1][-probl_routes2][-probl_routes3][-probl_routes4][-probl_routes5][-probl_routes6][-probl_routes7][-probl_routes8][probl_routes9]) %>% 
  st_geometry() %>% 
  plot(., add = TRUE, col = "yellow", pch = 19)
routes_sel_sf %>% 
  filter(RTENO_BBS %in% route_nrs[-probl_routes1][-probl_routes2][-probl_routes3][-probl_routes4][-probl_routes5][-probl_routes6][-probl_routes7][-probl_routes8][-probl_routes9][probl_routes10]) %>% 
  st_geometry() %>% 
  plot(., add = TRUE, col = "red", pch = 19)

plot(st_geometry(routes_sel_sf), main = paste0(spec, " round 15"))
routes_sel_sf %>% 
  filter(RTENO_BBS %in% route_nrs[probl_routes1]) %>% 
  st_geometry() %>% 
  plot(., add = TRUE, col = "yellow", pch = 19)
routes_sel_sf %>% 
  filter(RTENO_BBS %in% route_nrs[-probl_routes1][probl_routes2]) %>% 
  st_geometry() %>% 
  plot(., add = TRUE, col = "yellow", pch = 19)
routes_sel_sf %>% 
  filter(RTENO_BBS %in% route_nrs[-probl_routes1][-probl_routes2][probl_routes3]) %>% 
  st_geometry() %>% 
  plot(., add = TRUE, col = "yellow", pch = 19)
routes_sel_sf %>% 
  filter(RTENO_BBS %in% route_nrs[-probl_routes1][-probl_routes2][-probl_routes3][probl_routes4]) %>% 
  st_geometry() %>% 
  plot(., add = TRUE, col = "yellow", pch = 19)
routes_sel_sf %>% 
  filter(RTENO_BBS %in% route_nrs[-probl_routes1][-probl_routes2][-probl_routes3][-probl_routes4][probl_routes5]) %>% 
  st_geometry() %>% 
  plot(., add = TRUE, col = "yellow", pch = 19)
routes_sel_sf %>% 
  filter(RTENO_BBS %in% route_nrs[-probl_routes1][-probl_routes2][-probl_routes3][-probl_routes4][-probl_routes5][probl_routes6]) %>% 
  st_geometry() %>% 
  plot(., add = TRUE, col = "yellow", pch = 19)
routes_sel_sf %>% 
  filter(RTENO_BBS %in% route_nrs[-probl_routes1][-probl_routes2][-probl_routes3][-probl_routes4][-probl_routes5][-probl_routes6][probl_routes7]) %>% 
  st_geometry() %>% 
  plot(., add = TRUE, col = "yellow", pch = 19)
routes_sel_sf %>% 
  filter(RTENO_BBS %in% route_nrs[-probl_routes1][-probl_routes2][-probl_routes3][-probl_routes4][-probl_routes5][-probl_routes6][-probl_routes7][probl_routes8]) %>% 
  st_geometry() %>% 
  plot(., add = TRUE, col = "yellow", pch = 19)
routes_sel_sf %>% 
  filter(RTENO_BBS %in% route_nrs[-probl_routes1][-probl_routes2][-probl_routes3][-probl_routes4][-probl_routes5][-probl_routes6][-probl_routes7][-probl_routes8][probl_routes9]) %>% 
  st_geometry() %>% 
  plot(., add = TRUE, col = "yellow", pch = 19)
routes_sel_sf %>% 
  filter(RTENO_BBS %in% route_nrs[-probl_routes1][-probl_routes2][-probl_routes3][-probl_routes4][-probl_routes5][-probl_routes6][-probl_routes7][-probl_routes8][-probl_routes9]
         [probl_routes10]) %>% 
  st_geometry() %>% 
  plot(., add = TRUE, col = "yellow", pch = 19)
routes_sel_sf %>% 
  filter(RTENO_BBS %in% route_nrs[-probl_routes1][-probl_routes2][-probl_routes3][-probl_routes4][-probl_routes5][-probl_routes6][-probl_routes7][-probl_routes8][-probl_routes9]
         [-probl_routes10][probl_routes11]) %>% 
  st_geometry() %>% 
  plot(., add = TRUE, col = "yellow", pch = 19)
routes_sel_sf %>% 
  filter(RTENO_BBS %in% route_nrs[-probl_routes1][-probl_routes2][-probl_routes3][-probl_routes4][-probl_routes5][-probl_routes6][-probl_routes7][-probl_routes8][-probl_routes9]
         [-probl_routes10][-probl_routes11][probl_routes12]) %>% 
  st_geometry() %>% 
  plot(., add = TRUE, col = "yellow", pch = 19)
routes_sel_sf %>% 
  filter(RTENO_BBS %in% route_nrs[-probl_routes1][-probl_routes2][-probl_routes3][-probl_routes4][-probl_routes5][-probl_routes6][-probl_routes7][-probl_routes8][-probl_routes9]
         [-probl_routes10][-probl_routes11][-probl_routes12][probl_routes13]) %>% 
  st_geometry() %>% 
  plot(., add = TRUE, col = "yellow", pch = 19)
routes_sel_sf %>% 
  filter(RTENO_BBS %in% route_nrs[-probl_routes1][-probl_routes2][-probl_routes3][-probl_routes4][-probl_routes5][-probl_routes6][-probl_routes7][-probl_routes8][-probl_routes9]
         [-probl_routes10][-probl_routes11][-probl_routes12][-probl_routes13][probl_routes14]) %>% 
  st_geometry() %>% 
  plot(., add = TRUE, col = "yellow", pch = 19)
routes_sel_sf %>% 
  filter(RTENO_BBS %in% route_nrs[-probl_routes1][-probl_routes2][-probl_routes3][-probl_routes4][-probl_routes5][-probl_routes6][-probl_routes7][-probl_routes8][-probl_routes9]
         [-probl_routes10][-probl_routes11][-probl_routes12][-probl_routes13][-probl_routes14][probl_routes15]) %>% 
  st_geometry() %>% 
  plot(., add = TRUE, col = "red", pch = 19)


# these are influential observations regarding the posterior distribution, cause overfitting risk
# pareto-smoothing: approximation for these is unreliable
# If leaving out an observation changes the posterior too much then importance sampling is not able to give a reliable estimate

str(res_list$loo_cv) # p_loo = "effective number of parameters"
res_list$loo_cv$pointwise[loo::pareto_k_ids(res_list$loo_cv),]
# https://mc-stan.org/loo/reference/loo-glossary.html
res_list$loo_cv$pointwise %>%  View

res_list$loo_cv$diagnostics$n_eff[loo::pareto_k_ids(res_list$loo_cv)] # effective sample size is very small for most values with high k values
summary(res_list$loo_cv$diagnostics$n_eff)
order(res_list$loo_cv$diagnostics$n_eff)

# -> leave-one-series-out cross validation approximated with importance sampling had issues!

# elpd:
res_list$loo_cv$estimates

# gather problematic routes per species:
spec_probl_routes <- vector(mode = "list", length = 10)
names(spec_probl_routes) <- testspecs
for(i in 1:length(testspecs)){
  print(testspecs[i])
  load(file.path("results", paste0("out_flocker_", testspecs[i], "_fitted_preds_loo.RData")))
  spec_probl_routes[[testspecs[i]]] <- loo::pareto_k_ids(res_list$loo_cv)
}
save(spec_probl_routes, file = file.path("data", "problematic_routes_10_testspecs.RData"))


low_n_eff <- which(cv_it2000$diagnostics$n_eff < 2200)
high_k <- which(cv_it2000$diagnostics$pareto_k > 0.7)
which(low_n_eff %in% high_k)
which(high_k %in% low_n_eff)
plot(cv_it2000, diagnostic = "k")



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
load(file.path("results", "out_flocker_Eurasian Collared-Dove_fitted_preds_loo.RData"))

fitted_occ_col_ex <- res_list$fitted
occupancy_uncond <- res_list$Zs
prediction_sites_uncond <- res_list$preds_occ_uncond
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

#test <- loo_compare_flocker(list(out_norm_prior, out_hs_prior))

# test comparing two models: ----

testspecs <- c("Black-billed Cuckoo",
               "Great-tailed Grackle",
               "Scissor-tailed Flycatcher",
               "Common Raven",
               "Wild Turkey",
               "Bald Eagle",
               "Black-chinned Hummingbird",
               "Broad-winged Hawk",
               "Greater Roadrunner",
               "Eurasian Collared-Dove")
spec <- testspecs[1]

load(file.path("results", paste0("out_flocker_", spec, ".RData")))
out1 <- out
load(file.path("results", paste0("out_flocker_", spec, "_round2.RData")))
out2 <- out
rm(out)
load(file.path("results", paste0("out_flocker_", spec, "_fitted_preds_loo.RData")))
res_list1 <- res_list
load(file.path("results", paste0("out_flocker_", spec, "_fitted_preds_loo_round2.RData")))
res_list2 <- res_list
rm(res_list)

res_list1$loo_cv
res_list2$loo_cv



# which routes are these:
route_nrs <- matrix(route_sel_env_dt_final$RTENO, nrow = nsites, ncol = nyears, byrow = TRUE)[,1]

route_nrs2 <- route_nrs[-c(loo::pareto_k_ids(res_list1$loo_cv))]

plot(st_geometry(routes_sel_sf), main = spec)
routes_sel_sf %>% 
  filter(RTENO_BBS %in% route_nrs[loo::pareto_k_ids(res_list1$loo_cv)]) %>% 
  st_geometry() %>% 
  plot(., add = TRUE, col = "red", pch = 19)
routes_sel_sf %>% 
  filter(RTENO_BBS %in% route_nrs2[loo::pareto_k_ids(res_list2$loo_cv)]) %>% 
  st_geometry() %>% 
  plot(., add = TRUE, col = "blue", pch = 22)


plot(res_list1$loo_cv)
plot(res_list2$loo_cv)

#loo_comp <- loo_compare_flocker(list(out1, out2)) # can't compare models with different underlying data / different number of data points (full model, model with subset of data)

# test loglik flocker function
# 4000 samples of loglik for each route / time series
ll1 <- log_lik_flocker(out1)
ll1
dim(ll1) # 4000, 476

# as_tibble(ll1) %>% 
#   tidyr::pivot_longer(cols = V1:V476, names_to = "route", values_to = "draw") %>% 
#   ggplot() +
#   geom_density(aes(x = draw, colour = route)) +
#   theme(legend.position = "none")





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
