# calculate time series of sum of routes with observed or predicted species presence:
# median per year, CI per year, 100 draws of posterior per year

# packages: --------------------------------------------------------------------

library(dplyr)
library(sf)
library(ggplot2)
library(doParallel)
library(ggrepel)

# register cores for parallel computation:
ncores <- 2
cl <- makeCluster(ncores, setup_timeout = 0.5)
registerDoParallel(cl)

# functions: -------------------------------------------------------------------

source("0_functions.R")


# directories: -----------------------------------------------------------------

# main directory:
# main_dir <- file.path("//NAS-2-P-SN-01.ibb.uni-potsdam.de", "daten$", "AG26", "Transfer",
#                      "Schifferle_BBS_occupancy_models_2023")

main_dir <- file.path("/mnt", "ibb_share", "zurell_transfer", "Schifferle_BBS_occupancy_models_2023")

start <- 1904
end <- 2019


# load data: -------------------------------------------------------------------

# selected species:
load(file = file.path("data", "species_set_analysis.RData")) # output of 3_1_DOM_CV_evaluation_metrics.R
final_species


# calculate time series: -------------------------------------------------------

scenarios <- c("obsclim_histsoc", "counterclim_histsoc", "obsclim_1901soc")

# iterate over species:
foreach(s = 1:length(final_species),
        .packages = c("dplyr", "collapse", "sf", "bayestestR"),
        .errorhandling = "pass", #"remove",
        .verbose = TRUE) %dopar% {
          
          spec <- final_species[s]
          print(paste(s, spec))
          
          # iterate over scenarios:
          for(i in 1:length(scenarios)){
            
            print(scenarios[i])
            
            # check if species ran already:
            skip_to_next <- FALSE
            tryCatch(print(load(file.path(main_dir, "results", "fm_preds_ISIMIP", scenarios[i], 
                                          paste0(spec, "_y_preds_", scenarios[i],".RData")))), # output of 2_1_DOM_flocker_single_model.R
                     error = function(e) { skip_to_next <<- TRUE})
            if(skip_to_next) {
              print("predictions not found")
              next 
            }
            
            # load predictions:
            load(file.path(main_dir, "results", "fm_preds_ISIMIP", scenarios[i], 
                           paste0(spec, "_y_preds_", scenarios[i],".RData")))
            
            dim(y_preds_route) # sites, years, draws
            
            # sum across routes for each year:
            preds_years <- apply(y_preds_route, MAR = c(2,3), FUN = sum, na.rm = TRUE) 
            
            # median and 90% credible interval:
            ts_median <- apply(preds_years, MAR = 1, FUN = median)
            
            ts_ci90 <- apply(preds_years, MAR = 1, FUN = bayestestR::ci, ci = 0.9, method = "ETI")
            ts_ci90_low <- unlist(lapply(ts_ci90, FUN = function(x) x$CI_low))
            ts_ci90_high <- unlist(lapply(ts_ci90, FUN = function(x) x$CI_high))
            
            # 80 % CI:
            ts_ci80 <- apply(preds_years, MAR = 1, FUN = bayestestR::ci, ci = 0.8, method = "ETI")
            ts_ci80_low <- unlist(lapply(ts_ci80, FUN = function(x) x$CI_low))
            ts_ci80_high <- unlist(lapply(ts_ci80, FUN = function(x) x$CI_high))
            
            # add 100 draws of posterior distribution for each year:
            n_draws <- 100
            draws <- t(apply(preds_years, MAR = 1, FUN = function(x) sample(x = x, size = n_draws, replace = FALSE))) 
            colnames(draws) <- paste0("draw", 1:n_draws)
            draws <- draws %>% 
              as_tibble() %>% 
              mutate(year = start:end) %>% 
              select(year, everything())
            
            # draws %>%
            #   tidyr::pivot_longer(starts_with("draw"), names_to = "draw", values_to = "value") %>%
            #   ggplot(aes(x = year, y = value)) +
            #   geom_point() +
            #   geom_smooth(method = "lm")
            
            # assemble df:
            ts_preds <- tibble(year = start:end, 
                               median_Nocc = ts_median, 
                               CIlow = ts_ci90_low, 
                               CI90high = ts_ci90_high,
                               CI80low = ts_ci80_low, 
                               CI80high = ts_ci80_high) %>% 
              left_join(draws)
            
            # # plot:
            # ggplot(ts_preds) +
            #   geom_line(aes(x = year, y = median_Nocc_f)) +
            #   geom_ribbon(aes(x = year, ymax = CI90high_f, ymin = CIlow_f),
            #               alpha = 0.2, fill = "cornflowerblue") +
            #   ggtitle(spec) +
            #   theme_bw()
            
            # save predicted time series:
            out_dir <- file.path(main_dir, "results", "fm_preds_ISIMIP", scenarios[i],
                                 "time_series")
            if(!dir.exists(out_dir)){dir.create(out_dir)}
            
            save(ts_preds, file = file.path(out_dir, paste0(spec, "_ts_sum_occ_routes_", scenarios[i], ".RData")))
          }
        }

# 
# plots: -----------------

plot_dir <- file.path("plots", "attribution", "time_series_ISIMIP")
if(!dir.exists(plot_dir)){dir.create(plot_dir, recursive = TRUE)}

# load data to include observations:
# routes-years:
load(file = file.path("data", "BBS_for_occ_selection.RData")) # route_sel_dt; output of 1_3_match_BBS_to_env_data.R

# selected routes spatial data (to buffer presences):
routes_sel_sf <- st_read(file.path("data", "route_selection_1995_2019_surv_beg_end_max_5y_miss_v2_spat_thin_100km_max_30_r_per_BCR_centroids.shp")) # output of 1_1_route_selection.R

# route-year-species information (only surveyed)
load(file = file.path("data", "BBS_for_occ_spec_records.RData")) # bbs_dt_occ; output of 1_0_reformat_BBS_data.R



for(s in 104:length(final_species)){ # 84 Marsh Wren; 92-93; 100, 102, 106, 108, 113

  spec <- final_species[s]

  print(paste(s, spec))

  # observations:
  rel_routes <- training_routes(species = spec, buffer_km = 750, output = "RTENOs")
  occ_dt_spec <- BBS_pres_abs_spec(species = spec) %>%
    filter(RTENO %in% rel_routes)
  # route-level presence:
  # sum all routes for each year (temporal trend)
  obs_temp_trend <- occ_dt_spec %>% #xx
    group_by(Year) %>% #xx
    summarise(pres_sum = sum(presence, na.rm = TRUE)) #xx


  # factual predictions:
  load(file.path(main_dir, "results", "fm_preds_ISIMIP", "obsclim_histsoc", "time_series",
                 paste0(spec, "_ts_sum_occ_routes_obsclim_histsoc.RData")))
  ts_preds_f <- ts_preds

  # counterfactual predictions climate:
  load(file.path(main_dir, "results", "fm_preds_ISIMIP", "counterclim_histsoc", "time_series",
                 paste0(spec, "_ts_sum_occ_routes_counterclim_histsoc.RData")))
  ts_preds_cfclim <- ts_preds

  # counterfactual predictions land use:
  load(file.path(main_dir, "results", "fm_preds_ISIMIP", "obsclim_1901soc", "time_series",
                 paste0(spec, "_ts_sum_occ_routes_obsclim_1901soc.RData")))
  ts_preds_cflu <- ts_preds


  time_series_df <- data.frame(spec = spec, year = start:end)
  time_series_df <- time_series_df %>%
    left_join(obs_temp_trend, by = c(year = "Year")) %>%
    rename(obs = "pres_sum")

  time_series_df$fact <- ts_preds_f$median_Nocc
  time_series_df$fact_CIlow <- ts_preds_f$CI80low
  time_series_df$fact_CIhigh <- ts_preds_f$CI80high

  time_series_df$cfclim <- ts_preds_cfclim$median_Nocc
  time_series_df$cfclim_CIlow <- ts_preds_cfclim$CI80low
  time_series_df$cfclim_CIhigh <- ts_preds_cfclim$CI80high

  time_series_df$cflu <- ts_preds_cflu$median_Nocc
  time_series_df$cflu_CIlow <- ts_preds_cflu$CI80low
  time_series_df$cflu_CIhigh <- ts_preds_cflu$CI80high

  # compare with using predictions for initial occupancy based on factual data
  # ggplot, mark lines:
  # https://r-graph-gallery.com/web-line-chart-with-labels-at-end-of-line.html


  ts_median_df_lf <- time_series_df %>%
    select(!matches("CI")) %>%
    tidyr::pivot_longer(cols = c(obs, fact, cfclim, cflu), names_to = "scenario", values_to = "N_routes_median") %>%
    mutate(scenario = factor(scenario, levels = c("obs", "cfclim", "fact", "cflu")),
           name_lab = if_else(year == 2019, scenario, NA_character_),
           name_lab = case_match(name_lab,
                                 "obs" ~ "observations",
                                 "fact" ~ "obsclim_histsoc",
                                 "cfclim" ~ "counterclim_histsoc",
                                 "cflu" ~ "obsclim_1901soc"))

  ts_CI_df_lf <- time_series_df %>%
    select(c(spec, year, matches("CI"))) %>%
    tidyr::pivot_longer(cols = matches("CI"), names_to = "scenario", values_to = "N_routes") %>%
    mutate(CI = gsub(".*_", "", scenario),
           scenario = gsub("_.*", "", scenario)) %>%
    tidyr::pivot_wider(names_from = CI, values_from = N_routes)

  time_series_df_lf <- ts_median_df_lf %>%
    left_join(ts_CI_df_lf, by = c("spec", "year", "scenario"))

  ## linetype:
  line_type1 <- rep(1, 2)
  names(line_type1) <- c("obs", "fact")
  line_type2 <- rep(5, 2)
  names(line_type2) <- c("cfclim", "cflu")
  line_type <- c(line_type1, line_type2)

  ## colours:
  cols <- c("cfclim" = "#0D98BA",
            "obs" = "grey40",
            "fact" = "#FFCC17", #85CB33
            "cflu" = "#C7653A")#"#B7410E"

  p <- time_series_df_lf %>%
    filter(scenario != "obs") %>%
    ggplot(aes(x = year, y = N_routes_median, group = scenario)) +
    # geometric annotations that play the role of grid lines (to avoid grid lines on the right where labels are)
    geom_vline(xintercept = seq(1900, 2020, by = 10), color = "grey90", linewidth = .6) +
    geom_segment(data = tibble(y = seq(round(min(time_series_df_lf$CIlow, na.rm = TRUE), digits = -1),
                                       round(max(time_series_df_lf$CIhigh, na.rm = TRUE), digits = -1), length = 5), # xx
                               x1 = 1900, x2 = 2020),
                 aes(x = x1, xend = x2, y = y, yend = y), inherit.aes = FALSE, color = "grey90", linewidth = .6) +
    # data lines:
    geom_ribbon(aes(ymin = CIlow, ymax = CIhigh, fill = scenario), alpha = 0.3) +
    geom_line(aes(colour = scenario
                  #, linetype = scenario
                  ), linewidth = 1) +
    # geom_point(data = time_series_df_lf %>%  filter(scenario == "obs"),
    #            aes(colour = scenario), size = 2) +
    labs(title = spec,  subtitle = "median + 80 % CI", y = "Number of occupied routes", x = "") +
    geom_text_repel(
      aes(color = scenario, label = name_lab),
      fontface = "bold",
      size = 6,
      direction = "y", xlim = c(2020, NA), hjust = 0,
      segment.size = .7, segment.alpha = .5, segment.linetype = 1,
      #box.padding = .4,
      segment.curvature = -0.1, segment.ncp = 3, segment.angle = 20,
      max.overlaps = 30
    ) +
    scale_x_continuous(expand = c(0, 0), limits = c(1900, 2040), breaks = seq(1900, 2020, 10)) +
    #scale_linetype_manual(values = line_type) +
    scale_colour_manual(values = cols) +
    scale_fill_manual(values = cols) +
    theme_bw() +
    theme(legend.position = "none", panel.grid = element_blank(),
          text = element_text(size = 25))

  p


  jpeg(file = file.path(plot_dir, paste0(spec, "_ts_unc.jpg")),
       width = 1400, height = 700, quality = 100)
  print(p)
  dev.off()

}

