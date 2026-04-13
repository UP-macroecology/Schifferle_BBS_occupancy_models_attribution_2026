# calculate time series of occupancy change between 1995 and 2019 aggregated 
# across the conterminous USA based on: 
# - observations
# - dynamic occupancy model predictions for factual data (climate + land use change)
# - dynamic occupancy model predictions for counterfactual scenarios: 
#   - no climate change
#   - no land use change
#   - no climate + no land use change
# as basis for attribution
# ( + time series plots)


# packages: --------------------------------------------------------------------

library(dplyr)
library(sf)
library(doParallel)
library(bayestestR)
library(ggplot2)
library(ggrepel)


# functions: -------------------------------------------------------------------

source(file.path("scripts", "0_functions.R"))


# directories: -----------------------------------------------------------------

# project directory:
# dir <- file.path("//NAS-2-P-SN-01.ibb.uni-potsdam.de", "daten$", "AG26", "Transfer", "Schifferle_BBS_occupancy_models_2023")
dir <- file.path("/mnt", "ibb_share", "zurell_transfer", "Schifferle_BBS_occupancy_models_2023")

# directory with fitted models:
res_dir <- file.path(dir, "results", "fm_buffer750km")

# save observations time series:
obs_dir <- file.path(dir, "data", "observed_time_series_1995_2019")

# save predicted time series for factual data:
fact_dir <- file.path(dir, "results", "fm_buffer750km", "fact_pred_time_series_1995_2019")

# save predicted time series for counterfactual data:
cfact_dir <- file.path(dir, "results", "attribution", "cfact_pred_time_series_1995_2019")

if(!dir.exists(obs_dir)){dir.create(obs_dir)}
if(!dir.exists(fact_dir)){dir.create(fact_dir)}
if(!dir.exists(cfact_dir)){dir.create(cfact_dir)}

# directory to save plots:
plot_dir <- file.path("plots", "attribution", "fm_y_preds_routes_cf_1995_all")
if(!dir.exists(plot_dir)){dir.create(plot_dir, recursive = TRUE)}


# load data: -------------------------------------------------------------------

# species for attribution:
load(file = file.path("data", "species_DOM_val_okay.RData")) # output of 4_0_DOMs_predictions_y_routes_scenarios.R
spec_attr

# routes-years:
load(file = file.path("data", "BBS_for_occ_selection.RData")) # route_sel_dt; output of 1_3_dataprep_match_BBS_routes_env_data.R

# selected routes spatial data (to buffer presences):
routes_sel_sf <- st_read(file.path("data", "route_selection_1995_2019_surv_beg_end_max_5y_miss_v2_spat_thin_100km_max_30_r_per_BCR_centroids.shp")) # output of 1_1_dataprep_BBS_route_selection.R

# route-year-species information (only surveyed)
load(file = file.path("data", "BBS_for_occ_spec_records.RData")) # bbs_dt_occ; output of 1_0_dataprep_BBS_bird_data.R


# calculate time series: -------------------------------------------------------

# register cores for parallel computation:
ncores <- 2
cl <- makeCluster(ncores, setup_timeout = 0.5)
registerDoParallel(cl)


# iterate over species:

foreach(s = 1:length(spec_attr),
  .packages = c("dplyr", "collapse", "sf", "bayestestR"),
  .errorhandling = "pass", #"remove",
  .verbose = TRUE) %dopar% {
  
  spec <- spec_attr[s]
  print(paste(s, spec))
  
  # observations: -----------------------------
  
  print("observations")
  
  rel_routes <- training_routes(species = spec, buffer_km = 750, output = "RTENOs")
  occ_dt_spec <- BBS_pres_abs_spec(species = spec) %>%
    filter(RTENO %in% rel_routes)
  
  # sum presences across all routes within buffer for each year:
  ts_obs <- occ_dt_spec %>%
    rename("year" = Year) %>% 
    group_by(year) %>%
    summarise(Npres = sum(presence, na.rm = TRUE))
  
  save(ts_obs, file = file.path(obs_dir, paste0(spec, "_obs_ts_sum_occ_routes.RData")))
  
  
  # predictions for factual environmental data: ---------
  
  print("factual predictions")
  
  # check where to look for model output (did MCMC fitting work with less or only with more iterations?)
  if(file.exists(file.path(res_dir, "refit_2000_2000", paste0("out_", spec, "_fm_buffer_750.RData")))){ # output of 2_1_fit_DOMs_full_model.R
    output_dir <- file.path(res_dir, "refit_2000_2000")
  } else {
    output_dir <- res_dir
  }
  
  # aggregate predicted y for route sections at route level:
  
  # load(file.path(output_dir, paste0("postproc_", spec, "_fm_buffer750.RData")))
  # preds_routes <- apply(res_list$y_preds, MAR = c(1,3,4), FUN = max)
  # save(preds_routes, file = file.path(res_dir, "y_preds_route_level_section_sum", paste0(spec, "_y_preds_route_level_section_sum.RData")))
  load(file.path(res_dir, "y_preds_route_level_section_sum", paste0(spec, "_y_preds_route_level_section_sum.RData"))) 
  preds_routes # sites, years, draws
  
  # aggregate across the conterminous USA:
  # sum across routes for each year:
  preds_years <- apply(preds_routes, MAR = c(2,3), FUN = sum, na.rm = TRUE) 
  
  # median:
  ts_median <- apply(preds_years, MAR = 1, FUN = median)
  
  # 100 draws of posterior distribution for each year:
  n_draws <- 100
  draws <- t(apply(preds_years, MAR = 1, FUN = function(x) sample(x = x, size = n_draws, replace = FALSE))) 
  colnames(draws) <- paste0("draw", 1:n_draws)
  draws <- draws %>% 
    as_tibble() %>% 
    mutate(year = 1995:2019) %>% 
    select(year, everything())
  
  # 90% credible interval
  ts_ci90 <- apply(preds_years, MAR = 1, FUN = bayestestR::ci, ci = 0.9, method = "ETI")
  ts_ci90_low <- unlist(lapply(ts_ci90, FUN = function(x) x$CI_low))
  ts_ci90_high <- unlist(lapply(ts_ci90, FUN = function(x) x$CI_high))
  # 80% credible interval
  ts_ci80 <- apply(preds_years, MAR = 1, FUN = bayestestR::ci, ci = 0.8, method = "ETI")
  ts_ci80_low <- unlist(lapply(ts_ci80, FUN = function(x) x$CI_low))
  ts_ci80_high <- unlist(lapply(ts_ci80, FUN = function(x) x$CI_high))
  
  # assemble df:
  ts_preds_fact <- tibble(year = 1995:2019, 
                          median_Nocc = ts_median, 
                          CI90low = ts_ci90_low, 
                          CI90high = ts_ci90_high,
                          CI80low = ts_ci80_low, 
                          CI80high = ts_ci80_high) %>% 
    left_join(draws)
  
  # # plot:
  # ggplot(ts_preds_fact) +
  #   geom_line(aes(x = year, y = median_Nocc)) +
  #   geom_ribbon(aes(x = year, ymax = CI90high, ymin = CI90low),
  #               alpha = 0.2, fill = "cornflowerblue") +
  #   geom_point(data = ts_preds_fact %>% tidyr::pivot_longer(starts_with("draw"), names_to = "draw", values_to = "value"),
  #              aes(x = year, y = value)) +
  #   ggtitle(spec) +
  #   ylab("N routes") +
  #   theme_bw()
  
  save(ts_preds_fact, file = file.path(fact_dir, paste0(spec, "_ts_sum_occ_routes_f_preds.RData")))
  
  
  # predictions for counterfactual environmental data: ---------
  
  print("counterfactual predictions")
  
  # species counterfactual predictions:
  cf_files <- list.files(file.path(dir, "results", "attribution", "fm_y_preds_routes_cf_1995_all"), 
                         pattern = spec,
                         full.names = TRUE) # output of 4_0_DOMs_predictions_y_routes_scenarios.R
  
  # iterate over counterfactual scenarios:
  
  for(v in 1:length(cf_files)){
    
    print(v)
    
    # extract scenario from file name:
    scen <- unlist(stringr::str_split(cf_files[v], pattern = "(y_preds_cf_)|(.RData)"))[2]
    
    # aggregate predictions across the conterminous USA:
    # sum across routes for each year:
    load(cf_files[v]) # y_preds_route_cf
    preds_years <- apply(y_preds_route_cf, MAR = c(2,3), FUN = sum, na.rm = TRUE)
    
    # median:
    ts_median <- apply(preds_years, MAR = 1, FUN = median)
    
    # 100 draws of posterior distribution for each year:
    n_draws <- 100
    draws <- t(apply(preds_years, MAR = 1, FUN = function(x) sample(x = x, size = n_draws, replace = FALSE))) 
    colnames(draws) <- paste0("draw", 1:n_draws)
    draws <- draws %>% 
      as_tibble() %>% 
      mutate(year = 1995:2019) %>% 
      select(year, everything())
    
    # 90% credible interval
    ts_ci90 <- apply(preds_years, MAR = 1, FUN = bayestestR::ci, ci = 0.9, method = "ETI")
    ts_ci90_low <- unlist(lapply(ts_ci90, FUN = function(x) x$CI_low))
    ts_ci90_high <- unlist(lapply(ts_ci90, FUN = function(x) x$CI_high))
    # 80% credible interval
    ts_ci80 <- apply(preds_years, MAR = 1, FUN = bayestestR::ci, ci = 0.8, method = "ETI")
    ts_ci80_low <- unlist(lapply(ts_ci80, FUN = function(x) x$CI_low))
    ts_ci80_high <- unlist(lapply(ts_ci80, FUN = function(x) x$CI_high))
    
    # assemble df:
    ts_preds_cfact_df <- tibble(year = 1995:2019, 
                                median_Nocc = ts_median, 
                                CI90low = ts_ci90_low, 
                                CI90high = ts_ci90_high,
                                CI80low = ts_ci80_low, 
                                CI80high = ts_ci80_high) %>% 
      left_join(draws)
    
    if(scen == "counterclim"){
      ts_preds_cfact_clim <- ts_preds_cfact_df
    } else if(scen == "1995soc"){
      ts_preds_cfact_1995soc <- ts_preds_cfact_df
    } else {
      ts_preds_cfact_clim_1995soc <- ts_preds_cfact_df
    }
  }
  
  # save as list:
  ts_preds_cfact <- list(ts_preds_cfact_clim, ts_preds_cfact_1995soc, ts_preds_cfact_clim_1995soc)
  names(ts_preds_cfact) <- c("cf_clim", "cf_1995soc", "cf_clim_1995soc")
  
  save(ts_preds_cfact, file = file.path(cfact_dir, paste0(spec, "_ts_sum_occ_routes_cf_preds.RData")))
  
  }


# # plot time series: ------------------------------------------------------------
# 
# start <- 1995
# end <- 2019
# 
# for(s in 1:length(spec_attr)){
#   
#   spec <- spec_attr[s]
#   
#   print(paste(s, spec))
#   
#   # observations time series:
#   load(file.path(dir, "data", "observed_time_series_1995_2019", paste0(spec, "_obs_ts_sum_occ_routes.RData"))) 
#   ts_obs
#   
#   # factual predictions time series:
#   load(file.path(dir, "results", "fm_buffer750km", "fact_pred_time_series_1995_2019", 
#                  paste0(spec, "_ts_sum_occ_routes_f_preds.RData")))
#   ts_preds_f <- ts_preds_fact
#   
#   # counterfactual predictions time series:
#   load(file.path(dir, "results", "attribution", "cfact_pred_time_series_1995_2019", 
#                  paste0(spec, "_ts_sum_occ_routes_cf_preds.RData")))
#   ts_preds_cfact
#   
#   # assemble df:
#   
#   time_series_df <- tibble(spec = spec, year = start:end) %>% 
#     left_join(ts_obs, by = "year")
#   
#   # factual:
#   time_series_df$fact <- ts_preds_f$median_Nocc
#   time_series_df$fact_CIlow <- ts_preds_f$CI80low
#   time_series_df$fact_CIhigh <- ts_preds_f$CI80high
#   
#   # counterfactual climate:
#   time_series_df$cfclim <- ts_preds_cfact$cf_clim$median_Nocc
#   time_series_df$cfclim_CIlow <- ts_preds_cfact$cf_clim$CI80low
#   time_series_df$cfclim_CIhigh <- ts_preds_cfact$cf_clim$CI80high
#   
#   # counterfactual land use:
#   time_series_df$cflu <- ts_preds_cfact$cf_1995soc$median_Nocc
#   time_series_df$cflu_CIlow <- ts_preds_cfact$cf_1995soc$CI80low
#   time_series_df$cflu_CIhigh <- ts_preds_cfact$cf_1995soc$CI80high
#   
#   # counterfactual climate + land use:
#   time_series_df$cfclimlu <- ts_preds_cfact$cf_clim_1995soc$median_Nocc
#   time_series_df$cfclimlu_CIlow <- ts_preds_cfact$cf_clim_1995soc$CI80low
#   time_series_df$cfclimlu_CIhigh <- ts_preds_cfact$cf_clim_1995soc$CI80high
#   
#   # plot:
#   # https://r-graph-gallery.com/web-line-chart-with-labels-at-end-of-line.html
#   
#   # line plot:
#   
#   # data:
#   time_series_df_lf <- time_series_df %>%
#     select(-matches("CI")) %>% 
#     tidyr::pivot_longer(cols = -c(spec, year), names_to = "scenario", values_to = "N_routes_median") %>%
#     mutate(scenario = factor(scenario),
#            name_lab = if_else(year == 2019, scenario, NA_character_),
#            name_lab = case_match(name_lab,
#                                  "Npres" ~ "observations",
#                                  "fact" ~ "all factual",
#                                  "cfclim" ~ "all climate counterfact.",
#                                  "cflu" ~ "all land use counterfact.",
#                                  "cfclimlu" ~ "all counterfact."))
#   
#   # colours:
#   cols <- c("cfclim" = "#0D98BA",
#             "cflu" = "#B7410E",
#             "cfclimlu" = "#046865",
#             "Npres" = "black",
#             "fact" = "#85CB33")
#   
#   # linetype:
#   line_type <- c("Npres" = 1,
#                  "fact" = 1,
#                  "cfclim" = 5,
#                  "cflu" = 5,
#                  "cfclimlu" = 5)
#   # plot:
#   p <- ggplot(data = time_series_df_lf,
#               aes(x = year, y = N_routes_median, group = scenario)) +
#     # geometric annotations that play the role of grid lines (to avoud grid lines on the right where labels are)
#     geom_vline(xintercept = seq(1995, 2020, by = 5), color = "grey90", linewidth = .6) +
#     geom_segment(data = tibble(y = seq(round(min(time_series_df_lf$N_routes_median), digits = -1),
#                                        round(max(time_series_df_lf$N_routes_median), digits = -1), length = 5),
#                                x1 = 1995, x2 = 2020),
#                  aes(x = x1, xend = x2, y = y, yend = y), inherit.aes = FALSE, color = "grey90", linewidth = .6) +
#     # data lines:
#     geom_line(aes(colour = scenario, linetype = scenario), linewidth = 1) +
#     geom_point(data = time_series_df_lf %>%  filter(scenario == "Npres"),
#                aes(colour = scenario), size = 3) +
#     ylab("N routes with presences") +
#     theme_bw() +
#     ggtitle(spec) +
#     geom_text_repel(
#       aes(color = scenario, label = name_lab),
#       fontface = "bold",
#       size = 6,
#       direction = "y", xlim = c(2020, NA), hjust = 0,
#       segment.size = .7, segment.alpha = .5, segment.linetype = 1,
#       #box.padding = .4,
#       segment.curvature = -0.1, segment.ncp = 3, segment.angle = 20,
#       max.overlaps = 30
#     ) +
#     scale_x_continuous(expand = c(0, 0), limits = c(1995, 2027)) +
#     scale_linetype_manual(values = line_type) +
#     scale_colour_manual(values = cols) +
#     theme(legend.position = "none", panel.grid = element_blank(),
#           text = element_text(size = 20))
#   p
#   
#   # jpeg(file = file.path(plot_dir, paste0(spec, "_ts_y_preds_cf_1995_all_line.jpg")), 
#   #      width = 1000, height = 700, quality = 100)
#   # print(p)
#   # dev.off()
#   # 
#   
#   # # version with credible intervals:
#   # 
#   # ts_median_df_lf <- time_series_df %>%
#   #   select(!matches("CI")) %>%
#   #   tidyr::pivot_longer(cols = c(Npres, fact, cfclim, cflu, cfclimlu), names_to = "scenario", values_to = "N_routes_median") %>%
#   #   mutate(scenario = factor(scenario, levels = c("Npres", "cfclim", "fact", "cflu", "cfclimlu")),
#   #          name_lab = if_else(year == 2019, scenario, NA_character_), # to get only one label
#   #          name_lab = case_match(name_lab,
#   #                                "Npres" ~ "observations",
#   #                                "fact" ~ "obsclim_histsoc",
#   #                                "cfclim" ~ "climate counterfact.\n1995",
#   #                                "cflu" ~ "land use counterfact.\n1995",
#   #                                "cfclimlu" ~ "climate and land use counterfactual 1995"))
#   # 
#   # ts_CI_df_lf <- time_series_df %>%
#   #   select(c(spec, year, matches("CI"))) %>%
#   #   tidyr::pivot_longer(cols = matches("CI"), names_to = "scenario", values_to = "N_routes") %>%
#   #   mutate(CI = gsub(".*_", "", scenario),
#   #          scenario = gsub("_.*", "", scenario)) %>%
#   #   tidyr::pivot_wider(names_from = CI, values_from = N_routes)
#   # 
#   # time_series_df_lf <- ts_median_df_lf %>%
#   #   left_join(ts_CI_df_lf, by = c("spec", "year", "scenario"))
#   # 
#   # 
#   # # colours:
#   # cols <- c("cfclim" = "#0D98BA",
#   #           "obs" = "grey40",
#   #           "fact" = "#FFCC17", #85CB33
#   #           "cflu" = "#C7653A",
#   #           "cfclimlu" = "#331832")#"#B7410E"
#   # 
#   #  # plot:
#   #  p <- time_series_df_lf %>%
#   #    filter(!scenario %in% c("Npres", "cfclimlu")) %>%
#   #    ggplot(aes(x = year, y = N_routes_median, group = scenario)) +
#   #    # geometric annotations that play the role of grid lines (to avoid grid lines on the right where labels are)
#   #    geom_vline(xintercept = seq(1995, 2020, by = 5), color = "grey90", linewidth = .6) +
#   #    geom_segment(data = tibble(y = seq(round(min(time_series_df_lf$CIlow, na.rm = TRUE), digits = -1),
#   #                                       round(max(time_series_df_lf$CIhigh, na.rm = TRUE), digits = -1), length = 5),
#   #                               x1 = 1995, 
#   #                               x2 = 2020),
#   #                 aes(x = x1, xend = x2, y = y, yend = y), 
#   #                 inherit.aes = FALSE, color = "grey90", linewidth = .6) +
#   #    # data lines:
#   #    geom_ribbon(aes(ymin = CIlow, ymax = CIhigh, fill = scenario), alpha = 0.3) +
#   #    geom_line(aes(colour = scenario), linewidth = 1) +
#   #    # geom_point(data = time_series_df_lf %>%  filter(scenario == "Npres"),
#   #    #            aes(colour = scenario), size = 2) +
#   #    labs(title = spec,  subtitle = "median + 80 % CI", y = "Number of occupied routes", x = "") +
#   #    geom_text_repel(
#   #      aes(color = scenario, label = name_lab),
#   #      fontface = "bold",
#   #      size = 6,
#   #      direction = "y", xlim = c(2020, NA), hjust = 0,
#   #      segment.size = .7, segment.alpha = .5, segment.linetype = 1,
#   #      #box.padding = .4,
#   #      segment.curvature = -0.1, segment.ncp = 3, segment.angle = 20,
#   #      max.overlaps = 30
#   #    ) +
#   #    scale_x_continuous(expand = c(0, 0), limits = c(1995, 2027)) +
#   #    scale_colour_manual(values = cols) +
#   #    scale_fill_manual(values = cols) +
#   #    theme_bw() +
#   #    theme(legend.position = "none", panel.grid = element_blank(),
#   #          text = element_text(size = 25))
#   #  
#   #  p
#   #  
#   #  jpeg(file = file.path(plot_dir, paste0(spec, "_ts_y_preds_cf_1995_all_unc.jpg")), 
#   #       width = 1000, height = 700, quality = 100)
#   #  print(p)
#   #  dev.off()
#   # 
# }
