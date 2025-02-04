# explore routes for which approximate leave-one-out cross validation (PSIS-LOO)
# (flocker::loo_flocker) yields too high Pareto k values (default threshold):

# packages:

library(flocker)
library(dplyr)
library(sf)
library(ggplot2)


# functions: -----

source("0_functions.R")


# directories: ----

results_dir <- file.path("M:", "Documents", "DEBTs", "analysis", "Schifferle_BBS_occupancy_models_2023",
                         "results", "fm_buffer750km")

buffer_km <- 750


# load data: ----

# route-year-species information (only surveyed)
load(file = file.path("data", "BBS_for_occ_spec_records.RData")) # bbs_dt_occ; output of 1_0_reformat_BBS_data.R

# routes-years:
load(file = file.path("data", "BBS_for_occ_selection.RData")) # route_sel_dt; output of 1_3_match_BBS_to_env_data.R 

# selected routes spatial data (to buffer presences):
routes_sel_sf <- st_read(file.path("data", "route_selection_1995_2019_surv_beg_end_max_5y_miss_v2_spat_thin_100km_max_30_r_per_BCR_centroids.shp")) # output of 1_1_route_selection.R

# selected species:
load(file = file.path("data", "species_set_analysis.RData")) # output of 3_1_DOM_CV_evaluation_metrics.R

# elevation:
US_contours_subset_simpl <- st_read(file.path("data", "US_elev_contours_simplified.shp")) # output of 5_0_DOM_predictions_USA.R

# outline US (to buffer presences):
US_sf <- st_read(file.path("data", "US_outline_ESRI102003.shp"))


# extract routes flagged as having too high pareto k values for each species: ----

high_k_list <- vector(mode = "list", length = length(final_species))

for(i in 1:length(final_species)){
  
  spec <- final_species[i]
  
  print(paste0(i, spec))
  
  names(high_k_list)[i] <- spec
  
  # check where to look for model output (did MCMC fitting work with less or only with more iterations?)
  if(file.exists(file.path(results_dir, "refit_2000_2000", paste0("out_", spec, "_fm_buffer750.RData")))){
    output_dir <- file.path(results_dir, "refit_2000_2000")
  } else {
    output_dir <- results_dir
  }
  
  print(output_dir)
  
  # load model post-processing results:
  skip_to_next <- FALSE
  tryCatch(print(load(file.path(output_dir, paste0("postproc_", spec, "_fm_buffer750.RData")))),
           error = function(e) { skip_to_next <<- TRUE})
  if(skip_to_next) { next }
  
  # high k routes indices:
  high_k_routes <- loo::pareto_k_ids(res_list$loo_cv)
  
  # species-specific relevant routes, within distance of 750 km of species records:
  rel_routes <- training_routes(species = spec, buffer_km = buffer_km, output = "RTENOs")
  
  # high k routes RTENO number:
  high_k_RTENO <- rel_routes[high_k_routes]
  print(high_k_RTENO)

  high_k_list[[i]] <- high_k_RTENO
  
  rm(res_list)
}

#save(high_k_list, file = file.path("data", "high_k_routes_per_species_update2.RData"))


# save data frame: species - total number of routes - routes with presences - ----
# routes with too high k values - % of routes with too high k values: 

high_k_overview_df <- data.frame("species" = final_species,
                                 "N_routes" = NA,
                                 "N_routes_pres" = NA,
                                 "N_routes_high_k" = NA,
                                 "prop_routes_high_k" = NA,
                                 "prop_pres_routes_high_k" = NA)

for(i in 1:nrow(high_k_overview_df)){
  
  print(i)
  
  rel_routes <- training_routes(species = high_k_overview_df$species[i], buffer_km = buffer_km, output = "RTENOs")
  
  high_k_overview_df$N_routes[i] <- length(rel_routes)
  
  pres_abs <- BBS_pres_abs_spec(species = high_k_overview_df$species[i])
  
  high_k_overview_df$N_routes_pres[i] <- pres_abs %>% 
    filter(RTENO %in% rel_routes) %>% 
    filter(presence == 1) %>% 
    pull(RTENO) %>% 
    unique() %>% 
    length()
    
  high_k_overview_df$N_routes_high_k[i] <- length(high_k_list[[high_k_overview_df$species[i]]])
      
  high_k_overview_df$prop_routes_high_k[i] <- round(high_k_overview_df$N_routes_high_k[i]/high_k_overview_df$N_routes[i], 3)

  # number of routes with presence that have too high k values:
  pres_high_k_routes <- pres_abs %>% 
    filter(RTENO %in% rel_routes) %>% 
    filter(presence == 1) %>% 
    filter(RTENO %in% high_k_list[[high_k_overview_df$species[i]]]) %>% 
    pull(RTENO) %>% 
    unique()
  
  # proportion of routes with presence that have high k values:
  high_k_overview_df$prop_pres_routes_high_k[i] <- round(length(pres_high_k_routes)/high_k_overview_df$N_routes_pres[i], 3)
  }

save(high_k_overview_df, file = file.path("data", "high_k_routes_overview.RData"))
summary(high_k_overview_df)


# explorations: ----------------------------------------------------------------

load(file.path("data", "high_k_routes_per_species_update2.RData"))

# how many routes are flagged as highly influential across species:
high_k_routes <- unique(unlist(high_k_list)) # 467 of 539

# number of species for which each route has too high k values:
sort(table(unlist(high_k_list)), decreasing = TRUE)

# number of highly influential routes per species:
sort(lengths(high_k_list)) # sort species by number of high k routes
summary(lengths(high_k_list))


## high k routes for many species: ----

n_specs_high_k <- sort(table(unlist(high_k_list)), decreasing = TRUE) # number of species for which each route has too high k values

jpeg(file = file.path("plots", "high_k_maps", "high_k_20_and_more_specs.jpg"), 
     width = 900, height = 600, quality = 100)
ggplot() +
  geom_sf(data = US_sf, fill = "floralwhite") +
  geom_sf(data = US_contours_subset_simpl, aes(colour = Contour), linewidth = 0.2) +
  geom_sf(data = routes_sel_sf, size = 3, shape = 20, colour = "grey70") +
  geom_sf(data = routes_sel_sf %>% 
            filter(RTENO_BBS %in% names(which(n_specs_high_k >= 20))), colour = "#FC4E07", size = 5) +
  scale_colour_gradientn(colours = c(terrain.colors(8)[-9], "grey80"),
                         transform = "sqrt", 
                         breaks = c(100, 200, 500, 700, 1000, 1500, 2100, 2500, 2900)) +
  labs(col = "elevation \n[m]") +
  ggtitle("high k for more than 20 species") +
  theme_bw() +
  theme(text = element_text(size = 20), legend.key.height = unit(3, "lines"))
dev.off()


## number of species with high k: ----

# map number of species for which each route is flagged as highly influential:
route_n_spec_df <- tibble(as.data.frame((table(unlist(high_k_list))))) %>% 
  rename("RTENO_BBS" = "Var1", "n_specs" = "Freq") %>% 
  mutate(RTENO_BBS = as.numeric(as.character((RTENO_BBS))))

jpeg(file = file.path("plots", "high_k_maps", "high_k_n_species.jpg"), 
     width = 900, height = 600, quality = 100)
ggplot() +
  geom_sf(data = US_sf, fill = "floralwhite") +
  geom_sf(data = US_contours_subset_simpl, aes(colour = Contour), linewidth = 0.2) +
  geom_sf(data = routes_sel_sf %>% 
            left_join(route_n_spec_df), aes(fill = n_specs), size = 5, shape = 21) +
  scale_colour_gradientn(colours = c(terrain.colors(8)[-9], "grey80"),
                         transform = "sqrt", 
                         breaks = c(100, 200, 500, 700, 1000, 1500, 2100, 2500, 2900)) +
  scale_fill_viridis_c(na.value = "white", transform = "sqrt") +
  labs(col = "elevation \n[m]", fill = "N species") +
  ggtitle("number of species for which route has too high pareto k values") +
  theme_bw() +
  theme(text = element_text(size = 20), legend.key.height = unit(3, "lines"))
dev.off()

## species specific high k routes: ----

# how many routes are highly influential only for one species:

unique_high_k_df <- data.frame("species" = final_species, "unique_high_k" = NA)

for(i in 1:length(final_species)){
  
  spec <- final_species[i]
  
  print(paste(i, spec))
  
  unique_k_routes_ind <- which(!high_k_list[[spec]] %in% unique(unlist(high_k_list[-which(final_species == spec)])))
  unique_high_k_df$unique_high_k[i] <- length(unique_k_routes_ind)
}

unique_high_k_df %>% 
  filter(unique_high_k != 0) %>% 
  arrange(-unique_high_k)
# for 52 of 159 species there are routes for which pareto k is high on this route only for this species
# most for Blackburnian Warbler


## N high k ~ N presences: ----
# is number of highly influential routes related to number of presences?

n_high_k <- data.frame("English_Common_Name" = names(sort(lengths(high_k_list))),
                       "n_high_k" = sort(lengths(high_k_list)))

n_presences <- bbs_dt_occ %>% 
  filter(Year >= 1995 & Year <= 2019) %>% 
  filter(RTENO %in% route_sel_dt$RTENO) %>% 
  filter(English_Common_Name %in% final_species) %>% 
  group_by(English_Common_Name) %>% 
  summarise(n_presences = n())

n_pres_high_k <- n_high_k %>% 
  left_join(n_presences)

cor.test(n_pres_high_k$n_high_k, n_pres_high_k$n_presences, method = "p")
# sign. negative correlation -> the more presences, the fewer highly influential routes


## map high k routes for each species: ----

for(i in 1:length(final_species)){
  
  spec <- final_species[i]
  print(paste(i, spec))
  
  rel_routes <- training_routes(species = spec, buffer_km = buffer_km, output = "RTENOs")
  
  pres_RTENO <- BBS_pres_abs_spec(species = spec) %>% 
    filter(presence == 1) %>% 
    pull(RTENO) %>% 
    unique
  
  jpeg(file = file.path("plots", "high_k_maps", "buffer_750km", 
                        paste0(spec, ".jpg")), 
       width = 900, height = 600, quality = 100)
  
  print(
    ggplot() +
      geom_sf(data = US_sf, fill = "floralwhite") +
      geom_sf(data = US_contours_subset_simpl, aes(colour = Contour), linewidth = 0.2) +
      geom_sf(data = routes_sel_sf, size = 3, shape = 20, colour = "grey70") +
      geom_sf(data = routes_sel_sf %>% 
                filter(RTENO_BBS %in% high_k_list[[spec]]), colour = "#FC4E07", size = 5) +
      geom_sf(data = routes_sel_sf %>%
                filter(RTENO_BBS %in% rel_routes), colour = "#E7B800", size = 3) +
      geom_sf(data = routes_sel_sf %>%
                filter(RTENO_BBS %in% pres_RTENO), colour = "#00AFBB", size = 3) +
      scale_colour_gradientn(colours = c(terrain.colors(8)[-9], "grey80"),
                             transform = "sqrt", 
                             breaks = c(100, 200, 500, 700, 1000, 1500, 2100, 2500, 2900)) +
      labs(col = "elevation \n[m]") +
      ggtitle(spec) +
      theme_bw() +
      theme(text = element_text(size = 20), legend.key.height = unit(3, "lines"))
    )
  dev.off()
}

## compare observations on high k routes: ----

# compare observations between routes that have high k values for one but not for another species:

for(i in 1:length(final_species)){
  
  spec <- final_species[i]
  
  print(paste(i, spec))
  
  obs <- BBS_pres_abs_spec(spec) %>% 
    select(Year, RTENO, presence) %>% 
    group_by(RTENO) %>% 
    summarise(n_pres = sum(presence, na.rm = TRUE), 
              n_abs = sum(presence == 0, na.rm = TRUE),
              n_NA = sum(is.na(presence))) %>% 
    mutate(high_k = as.factor(ifelse(RTENO %in% high_k_list[[spec]], 1, 0)))
  
  obs %>% 
    filter(n_pres != 0) %>% 
    ggplot() +
    geom_boxplot(aes(y = n_abs, x = high_k, group = high_k), fill = c("limegreen", "red3")) +
    ggtitle(spec) +
    theme_bw() +
    ylab("N absences")
  
  obs %>% 
    filter(n_pres != 0) %>% 
    ggplot() +
    geom_boxplot(aes(y = n_pres, x = high_k, group = high_k), fill = c("limegreen", "red3")) +
    ggtitle(spec) +
    theme_bw() +
    ylab("N presences")
}


## N raw col. events and high k: ----

# do difficult routes have a higher number of colonisation events in raw data?

col_ev <- data.frame("species" = final_species,
                     "n_col_high_k" = NA,
                     "n_col_ok_k" = NA)

# only look at routes with at least one presence:

for(i in 1:length(final_species)){
  
  spec <- final_species[i]
  print(paste(i, spec))
  
  rel_routes <- training_routes(species = spec, buffer_km = buffer_km, output = "RTENOs")
  
  # routes with at least one detection:
  pres_routes <- BBS_pres_abs_spec(spec) %>% 
    select(Year, RTENO, presence) %>% 
    group_by(RTENO) %>% 
    summarise(n_pres = sum(presence, na.rm = TRUE)) %>%
    filter(n_pres != 0) %>% 
    pull(RTENO)
  
  # mean number of raw colonisation events across high k routes with at least one detection:
  n_col_high_k <- BBS_pres_abs_spec(species = spec) %>% 
    select(RTENO, Year, presence) %>% 
    filter(RTENO %in% pres_routes) %>% 
    filter(RTENO %in% high_k_list[[spec]]) %>% 
    group_by(RTENO) %>% 
    mutate(diff = presence - lag(presence)) %>% # diff = 1 = colonisation event
    mutate(diff = ifelse(diff == 1, 1, 0)) %>% 
    mutate(n_col = sum(diff, na.rm = TRUE)) %>%
    select(RTENO, n_col) %>% 
    distinct() %>% 
    pull(n_col) %>% 
    mean
  
  col_ev$n_col_high_k[i] <- n_col_high_k
  
  n_col_ok_k <- BBS_pres_abs_spec(species = spec) %>% 
    filter(RTENO %in% pres_routes) %>%
    select(RTENO, Year, presence) %>% 
    filter(RTENO %in% rel_routes) %>% 
    filter(!RTENO %in% high_k_list[[spec]]) %>% 
    group_by(RTENO) %>% 
    mutate(diff = presence - lag(presence)) %>% # diff = 1 = colonisation event
    mutate(diff = ifelse(diff == 1, 1, 0)) %>%
    mutate(n_col = sum(diff, na.rm = TRUE)) %>%
    select(RTENO, n_col) %>% 
    distinct() %>% 
    pull(n_col) %>% 
    mean
  
  col_ev$n_col_ok_k[i] <- n_col_ok_k
  
}

#save(col_ev, file = file.path("data", "high_k_routes_raw_col_events.RData"))

col_ev %>% 
  tidyr::pivot_longer(cols = n_col_high_k:n_col_ok_k, names_to = "col_k", 
                      values_to = "mean_col_ev") %>% 
  mutate(col_k = ifelse(col_k == "n_col_high_k", "high pareto k", "okay pareto k")) %>% 
  ggplot() +
  geom_boxplot(aes(y = mean_col_ev, x = col_k, group = as.factor(col_k))) +
  theme_bw() +
  theme(text = element_text(size = 20)) +
  xlab("") +
  ylab("mean number of raw colonisations")

# routes with high pareto k values have on average a larger number of raw colonisation events
# than routes with okay pareto k values


## high k and env. data: ----

# does environment differ between routes with high k values and other routes:

# prepare env. data:

# merged route, year, environment data:
load(file = file.path("data", "route_year_env_data.RData"))

# scale variables:
route_sel_env_dt_scaled <- route_sel_env_dt_final %>% 
  select(-c(Latitude, Longitude, BCR, ObsN, doy)) %>% 
  mutate(across(!c(RTENO, Year, Surveyed), ~ as.numeric(scale(., center=mean(.), scale = sd(.)))))

env_dt_lf <- route_sel_env_dt_scaled %>%
  tidyr::pivot_longer(cols = c(bio1:sum_annual_crops), names_to = "variable", values_to = "value") %>% 
  select(RTENO, Year, variable, value)

# plot time series for env. variables, highlight certain routes:

high_k_RTENO <- "84047011"

ggplot(data = env_dt_lf)  +
  facet_wrap(facets = ~variable, scales = "free_y") +
  geom_line(data = env_dt_lf,
            aes(x = Year, y = value, group = RTENO), colour = "gray60") +
  geom_line(data = env_dt_lf %>% filter(RTENO %in% high_k_RTENO),
            aes(x = Year, y = value, colour = as.factor(RTENO))) +
  #scale_color_brewer(palette = "Paired") +
  theme_bw() +
  theme(legend.position = "none")

# which route is it:
ggplot() +
  geom_sf(data = US_sf, fill = "floralwhite") +
  geom_sf(data = routes_sel_sf, size = 3, shape = 20, colour = "grey70") +
  geom_sf(data = routes_sel_sf %>% 
            filter(RTENO_BBS %in% high_k_RTENO), colour = "#FC4E07", size = 5)


# some notes:
# 84034014: Chicago: highest percentage of urban land cover; probl. for 126 species
# 84027005: largest percentage of secondary nonforest
# 84083311: near Houston at the south coast, peak in summer precipitation and prec. seasonality in one year 
# 84014108: Californian coast: quite high isothermality
# 84017047: Rocky Mountains, Colorado, among coldest routes and highest diurnal temp range, among lowest percentage of urban lc
# 84014155: West Coast, rain quite low, pastures low
# 84089009: along NW coast, temp. annual range and mean diurnal range lowest, pr spring, autumn, winter quite high
# 84083060: high pasture around Dallas
# 84047011: North-East Boston, 2nd highest urban areas
# 84046110: Washington area second highest percentage of urban land cover
# 84089111: quite urban, NW corner of US
# 84089001: NW corner of US, quite extreme in several variables, quite wet, not urban at all, low diurnal temp. rang
# 84042026: highest secondary forest, low other land use
# 84017003-84017055, 84060010 (6 routes): highest parts of Rocky Mountains
# 84054012: among lowest secondary forest 
# 84050008: relative low isothermality (annual temp. range high compared to diurnal temp. range)#
# 84063022: relative wet



### boxplot high k vs. okay k routes for each env. variable: ----

spec <- "Northern Rough-winged Swallow" #"Hairy Woodpecker"
rel_routes <- training_routes(species = spec, buffer_km = buffer_km, output = "RTENOs")

env_pareto_k <- env_dt_lf %>% 
  filter(RTENO %in% rel_routes) %>% 
  mutate(pareto = ifelse(RTENO %in% high_k_list[[spec]], "high", "okay"))

# consider species presence-absence?!
occ_dt_spec <- BBS_pres_abs_spec(species = spec)

env_pareto_k_occ <- env_pareto_k %>% 
  left_join(occ_dt_spec %>% select(RTENO, Year, presence))

env_pareto_k_occ %>% 
  #filter(pareto == "okay") %>%
  ggplot()+
  facet_wrap(facets = ~variable, scales = "free_y") +
  geom_boxplot(aes(y = value, fill = as.factor(pareto)))


### plot for each species the env. variables for routes with high k values compared to others: ----

# plot for each species the env. variables for routes with high k values compared to others:

for(i in 1:length(final_species)){
  
  spec <- final_species[i]
  print(paste(i, spec))
  
  rel_routes <- training_routes(species = spec, buffer_km = buffer_km, output = "RTENOs")
  
  pres_RTENO <- BBS_pres_abs_spec(species = spec) %>% 
    filter(presence == 1) %>% 
    pull(RTENO) %>% 
    unique
  
  jpeg(file = file.path("plots", "high_k_maps", "env_data_comp_750km", 
                        paste0(spec, ".jpg")), 
       width = 1400, height = 1000, quality = 100)
  
  print(
    ggplot(data = env_dt_lf)  +
      facet_wrap(facets = ~variable, scales = "free_y") +
      geom_line(data = env_dt_lf %>% filter(RTENO %in% rel_routes),
                aes(x = Year, y = value, group = RTENO), colour = "gray80") +
      geom_line(data = env_dt_lf %>% filter(RTENO %in% pres_RTENO),
                aes(x = Year, y = value, group = RTENO), colour = "grey20", size = 1) +
      geom_line(data = env_dt_lf %>% filter(RTENO %in% high_k_list[[spec]]),
                aes(x = Year, y = value, colour = as.factor(RTENO))) +
      #scale_color_brewer(palette = "Paired") +
      theme_bw() +
      ggtitle(spec) +
      theme(legend.position = "none", text = element_text(size = 20))
  )
  dev.off()
}

### PCA: where are high k routes within 2-dim env. space?: ----

library(ade4)

high_k_20_specs <- names(n_specs_high_k)[which(n_specs_high_k >= 20)]

for(i in 1:length(final_species)){
  
  spec <- final_species[i]
  print(paste(i, spec))
  
  rel_routes <- training_routes(species = spec, buffer_km = buffer_km, output = "RTENOs")
  occ_dt_spec <- BBS_pres_abs_spec(species = spec)
  occ_dt_spec_buff <- occ_dt_spec %>% 
    select(RTENO, Year, presence) %>% 
    filter(RTENO %in% rel_routes)
  
  
  env_data_spec <-  route_sel_env_dt_scaled %>% 
    filter(Year == 2005) %>% # choose arbitrary year to have one data point per route
    filter(RTENO %in% rel_routes) %>% 
    select(RTENO, bio1:sum_annual_crops)
  
  pca.env <- dudi.pca(env_data_spec[,2:ncol(env_data_spec)],
                      scannf = FALSE,
                      nf = ncol(env_data_spec)) # 2
  
  k_RTENO <- high_k_list[[spec]]
  
  PC_dt <- tibble(pca.env$li) %>% 
    rename("PC1" = Axis1, "PC2" = Axis2) %>% 
    cbind(RTENO = env_data_spec$RTENO) %>% 
    mutate(high_k = factor(ifelse(RTENO %in% k_RTENO, "yes", "no"), levels = c("no", "yes"))) %>% 
    arrange(high_k)
  
  jpeg(file = file.path("plots", "high_k_maps", "buffer_750km", "env_PCA2",
                        paste0(spec, ".jpg")), 
       width = 900, height = 600, quality = 100)
  print(
    PC_dt %>% 
      ggplot() +
      # mark points that have high k for more than 20 species:
      geom_point(aes(x = PC1, y = PC2), 
                 size = 5, colour = "black", data = PC_dt %>% filter(RTENO %in% high_k_20_specs)) +
      geom_point(aes(x = PC1, y = PC2, colour = high_k), size = 3) +
      theme_bw() +
      ggtitle(spec) +
      theme(text = element_text(size = 20))
  )
  dev.off()
  
}

final_species
i = 133

BBS_pres_abs_spec("Tufted Titmouse") %>% 
  filter(RTENO == 84025081) %>% 
  pull(presence)
BBS_pres_abs_spec("Tufted Titmouse") %>% 
  filter(RTENO == 84025087) %>% 
  pull(presence)
BBS_pres_abs_spec("Tufted Titmouse") %>% 
  filter(RTENO == 84025073) %>% 
  pull(presence)


BBS_pres_abs_spec("Carolina Wren") %>% 
  filter(RTENO == 84025087) %>% 
  pull(presence)
BBS_pres_abs_spec("Carolina Wren") %>% 
  filter(RTENO == 84025073) %>% 
  pull(presence)


### map frequency of env.: -----

# may differ when looking at route cells only xx
library(terra)

# selected variables:
load(file = file.path("data", "selected_variables.RData")) # selvar_final; output of 1_2_variable_selection.R
selvar_final

# load data of arbitrary year:

# bioclimatic vars:
bioclim_path <- file.path("data", "Env_data", "ISIMIP_GSWP3_W5E5", "bioclim")
# seasonal climatic vars:
sclim_path <- file.path("data", "Env_data", "ISIMIP_GSWP3_W5E5", "seasonal")
# land use:
lu_path <- file.path("data", "Env_data", "ISIMIP_land_use_and_irrigation", "ISIMIP_LU_ESRI102003")

freq_layers <- vector(mode = "list", length = length(selvar_final))
for(v in 1:length(selvar_final)){
  
  print(v)
  
  var <- selvar_final[v]
  
  if(grepl(pattern = "bio", x = var)){
    v_tif <- rast(file.path(bioclim_path, paste0(var, "_2000.tif")))
  }
  else if(grepl(pattern = "pr_mean", x = var)){
    v_tif <- rast(file.path(sclim_path, paste0(var, "_2000.tif")))
  } else {
    v_tif <- rast(file.path(lu_path, paste0(var, "_2000_ESRI102003.tif")))
  }
  
  v_tif_scaled <- scale(v_tif)
  freq_cell_values <- freq(v_tif_scaled, digits = 1)
  freq_layers[[v]] <- test <- subst(round(v_tif_scaled, 1), from = freq_cell_values$value, to = freq_cell_values$count)
  
  
}

freq_layers

plot(freq_layers[[1]], alpha = 0.2, col = grey.colors(10, start=0, end=1))
plot(freq_layers[[2]], alpha = 0.2, col = grey.colors(10, start=0, end=1), add = TRUE)
plot(freq_layers[[3]], alpha = 0.2, col = grey.colors(10, start=0, end=1), add = TRUE)
plot(freq_layers[[4]], alpha = 0.2, col = grey.colors(10, start=0, end=1), add = TRUE)
plot(freq_layers[[5]], alpha = 0.2, col = grey.colors(10, start=0, end=1), add = TRUE)
plot(freq_layers[[6]], alpha = 0.2, col = grey.colors(10, start=0, end=1), add = TRUE)
plot(freq_layers[[7]], alpha = 0.2, col = grey.colors(10, start=0, end=1), add = TRUE)
plot(freq_layers[[8]], alpha = 0.2, col = grey.colors(10, start=0, end=1), add = TRUE)
plot(freq_layers[[9]], alpha = 0.2, col = grey.colors(10, start=0, end=1), add = TRUE)
plot(freq_layers[[10]], alpha = 0.2, col = grey.colors(10, start=0, end=1), add = TRUE)
plot(freq_layers[[11]], alpha = 0.2, col = grey.colors(10, start=0, end=1), add = TRUE)
plot(freq_layers[[12]], alpha = 0.2, col = grey.colors(10, start=0, end=1), add = TRUE)
plot(freq_layers[[13]], alpha = 0.2, col = grey.colors(10, start=0, end=1), add = TRUE)
plot(freq_layers[[14]], alpha = 0.2, col = grey.colors(10, start=0, end=1), add = TRUE)
plot(freq_layers[[15]], alpha = 0.2, col = grey.colors(10, start=0, end=1), add = TRUE)
plot(st_geometry(routes_sel_sf), add = TRUE)
