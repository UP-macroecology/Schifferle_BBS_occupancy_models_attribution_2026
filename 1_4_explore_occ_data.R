# exploration regarding data selected for occupancy modelling:

# packages: ----

library(dplyr)
library(sf)
library(ggplot2)


# load data: ----

source("0_functions.R") # BBS_pres_abs_spec()

# BBS route selection (centroids):
routes_sel_sf <- st_read(file.path("data", "route_selection_1995_2019_surv_beg_end_max_5y_miss_v2_spat_thin_100km_max_30_r_per_BCR_centroids.shp")) # output of 1_1_route_selection.R

# selected species:
load(file = file.path("data", "final_species_selection.RData")) # output of 1_2_species_selection.R

load(file = file.path("data", "final_species_selection_eco_sorted.RData")) # final_species_eco_sorted; output of 1_2_species_selection.R

# route-year-species information (only surveyed)
load(file = file.path("data", "BBS_for_occ_spec_records.RData")) # bbs_dt_occ; output of 1_0_reformat_BBS_data.R

# selected routes and focal years matched to environmental data: # route_sel_env_dt_final; output of 1_3_match_BBS_to_env_data.R
# merged route, year, environment data:
load(file = file.path("data", "route_year_env_data.RData"))

# routes-years:
load(file = file.path("data", "BBS_for_occ_selection.RData")) # route_sel_dt; output of 1_3_match_BBS_to_env_data.R 

# selected variables:
load(file = file.path("data", "selected_variables.RData")) # selvar_final; output of 1_2_variable_selection.R


# overview map species detections: ------------

plot_list <- vector(mode = "list", length = length(species_selection_final))

for(i in 1:length(species_selection_final)){
  
  print(species_selection_final[i])
  
  spec_pres_abs_cum <- BBS_pres_abs_spec(species_selection_final[i]) %>%
    right_join(routes_sel_sf, by = c("RTENO" = "RTENO_BBS")) %>% 
    st_as_sf() %>% 
    group_by(RTENO) %>%
    # presence = at least one record between 1995-2019:
    summarise(presence_summarised = max(presence, na.rm=TRUE)) %>%
    mutate(presence_summarised = factor(presence_summarised, levels = c(1,0)))

  plot_list[[i]] <- ggplot(spec_pres_abs_cum) +
    geom_sf(aes(colour = presence_summarised), size = 0.5) +
    ggtitle(species_selection_final[i]) +
    theme_bw() +
    theme(legend.position="bottom")
  
}
names(plot_list) <- species_selection_final

ggsave(
  filename = "plots/detection_maps_sel_routes_sel_species.pdf", 
  plot = gridExtra::marrangeGrob(plot_list, nrow=4, ncol=4), 
  width = 15, height = 9
)

# export shapefile with presences: ----

# presences only:
presences <- bbs_dt_occ %>% 
  select(c(English_Common_Name, RTENO, Year)) %>% 
  arrange(RTENO)

# merge with route spatial data:
pres_sf <- routes_sel_sf %>% 
  right_join(presences, by = c("RTENO_BBS" = "RTENO")) %>% 
  select(c("RTENO" = RTENO_BBS, Year, "species" = English_Common_Name))

st_write(pres_sf, file.path("data", "spec_presences.shp"), append = FALSE)



# trend in raw number of occupied routes over time: ----

# quantify trend as slope of linear model (sum occupied routes ~ year):

lm_trend_df <- data.frame("species" = final_species_eco_sorted,
                          "trend" = NA,
                          "p_trend" = NA)

for(i in 1:length(final_species_eco_sorted)){
  
  spec <- final_species_eco_sorted[i]
  
  print(paste(i, spec))
  
  # relevant routes for the species, within distance of 750 km of presences:
  rel_routes <- training_routes(species = spec, buffer_km = 750, output = "RTENOs")
  # observations:
  occ_dt_spec <- BBS_pres_abs_spec(species = spec) %>% 
    filter(RTENO %in% rel_routes)
  
  obs_temp_trend <- occ_dt_spec %>% 
    group_by(Year) %>% 
    summarise(pres_sum = sum(presence, na.rm = TRUE))
  
  trend_lm <- lm(pres_sum ~ Year, data = obs_temp_trend)
  summary_lm <- summary(trend_lm)
  trend <- summary_lm$coefficients[2,1]
  print(trend)
  pvalue_trend <- summary_lm$coefficients[2,4]
  print(pvalue_trend)
  
  lm_trend_df$trend[i] <- trend
  lm_trend_df$p_trend[i] <- pvalue_trend
}
lm_trend_df

write.csv(lm_trend_df, file = file.path("results", "raw_occ_trends_species.csv"), row.names = FALSE)

lm_trend_df %>% 
  filter(p_trend < 0.05) %>% 
  mutate(trend_abs = abs(trend)) %>% 
  View

lm_trend_df %>% 
  filter(p_trend > 0.05) %>% 
  mutate(trend_abs = abs(trend)) %>% 
  View


# trend in environmental variables over time: ----

# quantify trend as slope of linear model (env. variable ~ year):
# species specific: consider only data from routes within 750 km buffer around presences

# # these species have been discarded in a later step because MCMC fitting didn't work for them:
# species_discard <- c("Golden Eagle", "Prairie Falcon", "Sharp-shinned Hawk", "Broad-winged Hawk",
#                      "Cooper's Hawk", "Osprey", "Evening Grosbeak", "Golden-crowned Kinglet", "Olive-sided Flycatcher",
#                      "Wilson's Warbler", "Bullock's Oriole", "Western Wood-Pewee", "Black-throated Gray Warbler",
#                      "MacGillivray's Warbler", "Ruffed Grouse", "Northern Bobwhite",
#                      "Common Raven", "Common Yellowthroat", "Hairy Woodpecker", "Indigo Bunting",
#                      "Louisiana Waterthrush", "Northern Mockingbird", "Northern Parula", "Prairie Warbler",
#                      "Bewick's Wren", "Brewer's Blackbird")

#species_set <- final_species_eco_sorted[!final_species_eco_sorted %in% species_discard]
species_set <- species_selection_final

# selected variables (output of 1_2_variable_selection.R):
load(file = file.path("data", "selected_variables.RData"))
selvar_final

route_sel_env_dt <- route_sel_env_dt_final %>% 
  # only selected variables:
  select(c("RTENO", "Year", all_of(selvar_final)))

# data frame to store results:

colnames <- c("species", sort(paste0(rep(selvar_final, 2), c("_trend", "_pval"))))
lm_env_trend_df <- as.data.frame(matrix(NA, nrow = length(species_set), ncol = length(colnames)))
colnames(lm_env_trend_df) <- colnames
lm_env_trend_df$species <- species_set

for(i in 1:length(species_set)){
  
  spec <- species_set[i]
  
  print(paste(i, spec)) 
  
  # relevant routes for the species, within distance of 750 km of presences:
  rel_routes <- training_routes(species = spec, buffer_km = 750, output = "RTENOs")
  
  env_dt_spec <- route_sel_env_dt %>% 
    filter(RTENO %in% rel_routes)
  
  for(v in selvar_final){
    
    print(v)
    
    summary_lm <- summary(lm(unlist(env_dt_spec[,v]) ~ unlist(env_dt_spec[,"Year"])))
    lm_env_trend_df[i, paste0(v, "_trend")] <- summary_lm$coefficients[2,1]
    lm_env_trend_df[i, paste0(v, "_pval")] <- summary_lm$coefficients[2,4]

  }

}

write.csv(lm_env_trend_df, file = file.path("results", "env_vars_trends_species.csv"), row.names = FALSE)

# plot
route_sel_env_dt %>% 
  filter(RTENO %in% rel_routes) %>% 
  ggplot(aes(x = Year, y = secondary_nonforests)) +
  geom_line(aes(group = RTENO)) +
  geom_smooth(method = 'lm')




# counterfactual land use data: ----

# across US:

# factual environmental data for each US grid cell:
load(file = file.path("data", "US_grid_env_data.RData")) # output of 1_3_env_data_df_contUS.R; clim_lu_cells_sf
fact_env_df <- sf::st_drop_geometry(clim_lu_cells_sf)

# counterfactual environmental data for each US grid cell:
#load(file = file.path("data", "US_grid_env_data_cf.RData")) # xx output of 1_3_env_data_df_contUS.R; clim_lu_cells_sf
#counterfact_env_df <- sf::st_drop_geometry(clim_lu_cells_sf)

load(file = file.path("data", "US_grid_env_data_cf_1950soc.RData")) # xx output of 1_3_env_data_df_contUS.R; clim_lu_cells_sf
counterfact_env_df <- sf::st_drop_geometry(clim_lu_cells_sf_cf_year)
# scale env. data:

# load mean and sd with which training data were scaled:
load(file.path("data", "route_env_dt_scale_pars.RData")) # output of 2_1_DOM_flocker_fit_fm.R; env_scale_pars

# factual data:
fact_env_df_scaled <- fact_env_df
# iterate over columns:
for(c in colnames(fact_env_df[, 3:ncol(fact_env_df)])){
  print(c)
  fact_env_df_scaled[, c] <- as.numeric(scale(fact_env_df[, c], 
                                              center = as.numeric(env_scale_pars$center[c]),
                                              scale = as.numeric(env_scale_pars$scale[c])))
}

# counterfactual data:
counterfact_env_df_scaled <- counterfact_env_df
# iterate over columns:
for(c in colnames(counterfact_env_df[, 3:ncol(counterfact_env_df)])){
  print(c)
  counterfact_env_df_scaled[, c] <- as.numeric(scale(counterfact_env_df[, c], 
                                                     center = as.numeric(env_scale_pars$center[c]),
                                                     scale = as.numeric(env_scale_pars$scale[c])))
}


# factual all US:
fact_env_df_scaled %>% 
  select(c(year, grep(x = selvar_final, pattern = "(bio)|(pr_mean)", value = TRUE, invert = TRUE)))  %>% 
  tidyr::pivot_longer(cols = -year, names_to = "variable", values_to = "value") %>% 
  group_by(year, variable) %>% 
  summarise(sum = sum(value)) %>% 
  ggplot() +
  geom_line(aes(x = year, y = sum, colour = variable), linewidth = 2) +
  theme_bw() 

# counterfactual all US:
counterfact_env_df_scaled %>% 
  select(c(year, grep(x = selvar_final, pattern = "(bio)|(pr_mean)", value = TRUE, invert = TRUE)))  %>% 
  tidyr::pivot_longer(cols = -year, names_to = "variable", values_to = "value") %>% 
  group_by(year, variable) %>% 
  summarise(sum = sum(value)) %>% 
  ggplot() +
  geom_line(aes(x = year, y = sum, colour = variable), linewidth = 2) +
  theme_bw() 

# -> counterfactual data static over 1995 - 2019

# compare factual (histsoc) and counterfcatual (1901soc) land use data from (1850) 1901 to 2019:
# (via sum across US)

# load tifs:
library(terra)

# counterfactual data:

lu_vars <- grep(x = selvar_final, pattern = "(bio)|(pr_mean)", value = TRUE, invert = TRUE)

for(i in 1:length(lu_vars)){
  
  lu_v <- lu_vars[i]
  
  print(lu_v)
  
  files <- list.files(path = file.path("data", "Env_data", "ISIMIP_land_use_and_irrigation", 
                                       "ISIMIP_LU_ESRI102003"),
                      pattern = paste0(lu_v, "_[0-9]{4}"),
                      full.names = TRUE)
  
  rast_stack_f <- rast(files)
  #names(rast_stack_f) <- paste0(lu_v, "_", 1850:2019) 
  
  US_sum <- global(rast_stack_f, "sum", na.rm = TRUE)$sum
  US_sum 
  # for few variables data from 1850 onwards, remove data unitil 1901:
  if(length(US_sum) == 170){
    US_sum <- US_sum[52:length(US_sum)]
  }
  
  df_var <- data.frame(var = lu_v, sum_US_fact = US_sum, year = 1901:2019)
  
  # counterfactual data:
  files <- list.files(path = file.path("data", "Counterfactual_env_data", "ISIMIP_land_use_and_irrigation",
                                       "ISIMIP_LU_ESRI102003"),
                      pattern = paste0(lu_v, "_[0-9]{4}"),
                      full.names = TRUE)
  
  # files <- list.files(path = file.path("data", "Env_data", "ISIMIP_land_use_and_irrigation",
  #                                      "ISIMIP_LU_ESRI102003"),
  #                     pattern = paste0(lu_v, "_1995"),
  #                     full.names = TRUE)
  
  rast_stack_cf <- rast(files)
  
  #names(rast_stack_cf) <- paste0(lu_v, "_", 1850:2019) 
  
  sum_cf <- global(rast_stack_cf, "sum", na.rm = TRUE)$sum
  # for few variables data from 1850 onwards, remove data unitil 1901:
  if(length(sum_cf) == 170){
    sum_cf <- sum_cf[52:length(sum_cf)]
  }
  
  df_var$sum_US_cf <- sum_cf
  
  if(i == 1){
    df_all <- df_var
  } else {
    df_all <- rbind(df_all, df_var)
  }
  
}

# plot changes in cell sums of land use classes, factual vs. counterfactual:
df_all %>% 
  tidyr::pivot_longer(cols = c(sum_US_fact , sum_US_cf), names_to = "scen", values_to = "sum_US") %>% 
  mutate(scen = case_match(scen, "sum_US_cf" ~ "counterfact.", "sum_US_fact" ~ "fact.")) %>% 
  ggplot() +
  geom_point(aes(x = year, y = sum_US, colour = as.factor(scen)), size = 3) +
  geom_vline(xintercept = 1995, linetype = 2) +
  facet_wrap(~var, scales = "free") +
  theme_bw() +
  ylab("cells summed across US") +
  theme(legend.title = element_blank(), text = element_text(size = 25))

# how many routes were surveyed each year: ----
occ_dt_spec %>% 
  select(RTENO, Year, presence) %>% 
  #filter(presence == 1) %>% 
  filter(!is.na(presence)) %>% 
  group_by(Year) %>% 
  summarise(n = n()) %>% View
  ggplot() +
  geom_line(aes(x= Year, y = n))
  
