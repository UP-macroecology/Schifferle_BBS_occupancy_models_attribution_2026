# exploration regarding data selected for occupancy modelling:

# packages: ----

library(dplyr)
library(sf)
library(ggplot2)


# load data: ----

# BBS route selection (centroids):
routes_sel_sf <- st_read(file.path("data", "route_selection_1991_2015_surv_beg_end_max_5y_miss_v2_spat_thin_100km_max_30_r_per_BCR_centroids.shp")) # output of 1_1_route_selection.R

# selected species:
load(file = file.path("data", "final_species_selection.RData")) # output of 1_2_species_selection.R
species_selection_final <- sort(species_selection_final)
load(file = file.path("data", "final_species_selection_eco_sorted.RData")) # final_species_eco_sorted; output of 1_2_species_selection.R

# route-year-species information (only surveyed)
load(file = file.path("data", "BBS_for_occ_spec_records.RData")) # bbs_dt_occ; output of 1_0_reformat_BBS_data.R

# selected routes and focal years matched to environmental data: # route_sel_env_dt_final; output of 1_3_match_BBS_to_env_data.R
# merged route, year, environment data:
load(file = file.path("data", "route_year_env_data.RData"))

# routes-years:
load(file = file.path("data", "BBS_for_occ_selection.RData")) # route_sel_dt; output of 1_3_match_BBS_to_env_data.R 



# overview map species detections: ------------

plot_list <- vector(mode = "list", length = length(species_selection_final))

for(i in 1:length(species_selection_final)){
  
  print(species_selection_final[i])
  
  presences_spec <- bbs_dt_occ %>% 
    select(c(English_Common_Name, RTENO, Year, paste0("Count", seq(10, 50, 10)))) %>% 
    filter(English_Common_Name == species_selection_final[i])
  
  # match to routes-year-env:
  occ_dt_spec <- route_sel_env_dt_final %>% 
    select(c(RTENO, Year, Surveyed)) %>% 
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
  
  # match to spatial data:
  occ_spec_sf <- routes_sel_sf %>%
    left_join(occ_dt_spec, by = c("RTENO_BBS" = "RTENO")) %>% # species counts
    mutate(presence = rowSums(across(paste0("Count", seq(10, 50, 10))))) %>%
    mutate(presence = ifelse(presence >= 1, 1, 0))
  
  # summarise presences across years:
  spec_pres_abs_cum <- occ_spec_sf %>%
    group_by(RTENO_BBS) %>%
    summarise(presence_summarised = max(presence, na.rm=TRUE)) %>%
    mutate(presence_summarised = factor(presence_summarised, levels = c(1,0))) %>% 
    select(presence_summarised)
  
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




# explore correlation between env. variables taking only route locations (and not whole US) into account: ----

M <- cor(route_sel_env_dt_final[, c(9:32)], method = "s")
corrplot::corrplot(M, method = "square", order = "hclust",
         addCoef.col = "black",
         diag = FALSE,
         tl.cex = 1,
         number.cex = 0.8,
         number.digits= 2)

M <- cor(route_sel_env_dt_final[, c("bio2", "bio3", "bio5", "bio6", "bio13", "bio14",
                                    "sum_annual_crops","urban", "pastr", "secdn", "primn")], method = "s")
corrplot::corrplot(M, method = "square", order = "hclust",
                   addCoef.col = "black",
                   diag = FALSE,
                   tl.cex = 1,
                   number.cex = 0.8,
                   number.digits= 2)
# secdn and primn are too highly correlated considering only route locations

M <- cor(route_sel_env_dt_final[, c("bio1", "bio2", "bio3", "bio10", "bio12", "bio16")], method = "s")
corrplot::corrplot(M, method = "square", order = "hclust",
                   addCoef.col = "black",
                   diag = FALSE,
                   tl.cex = 1,
                   number.cex = 0.8,
                   number.digits= 2)
# bio1 and bio10, bio12 and bio16 are highly correlated taking only route locations in account







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

# these species have been discarded in a later step because MCMC fitting didn't work for them:
species_discard <- c("Golden Eagle", "Prairie Falcon", "Sharp-shinned Hawk", "Broad-winged Hawk",
                     "Cooper's Hawk", "Osprey", "Evening Grosbeak", "Golden-crowned Kinglet", "Olive-sided Flycatcher",
                     "Wilson's Warbler", "Bullock's Oriole", "Western Wood-Pewee", "Black-throated Gray Warbler",
                     "MacGillivray's Warbler", "Ruffed Grouse", "Northern Bobwhite",
                     "Common Raven", "Common Yellowthroat", "Hairy Woodpecker", "Indigo Bunting",
                     "Louisiana Waterthrush", "Northern Mockingbird", "Northern Parula", "Prairie Warbler",
                     "Bewick's Wren", "Brewer's Blackbird")

species_set <- final_species_eco_sorted[!final_species_eco_sorted %in% species_discard]

env_vars <- c("bio1", "bio2", "bio3", "bio7", "bio14", "bio15", 
              "pr_spring", "pr_summer","pr_autumn", "pr_winter",
              "sum_annual_crops", "secdf","pastr", "urban")

route_sel_env_dt <- route_sel_env_dt_final %>% 
  # only selected variables:
  select(c("RTENO", "Year", "bio1", "bio2", "bio3", "bio7", "bio14", "bio15", 
           "pr_spring", "pr_summer","pr_autumn", "pr_winter",
           "sum_annual_crops", "secdf","pastr", "urban"))

lm_env_trend_df <- data.frame("species" = species_set,
                              "bio1_trend" = NA,
                              "bio1_p_value" = NA,
                              "bio2_trend" = NA,
                              "bio2_p_value" = NA,
                              "bio3_trend" = NA,
                              "bio3_p_value" = NA,
                              "bio7_trend" = NA,
                              "bio7_p_value" = NA,
                              "bio14_trend" = NA,
                              "bio14_p_value" = NA,
                              "bio15_trend" = NA,
                              "bio15_p_value" = NA,
                              "pr_spring_trend" = NA,
                              "pr_spring_p_value" = NA,
                              "pr_summer_trend" = NA,
                              "pr_summer_p_value" = NA,
                              "pr_autumn_trend" = NA,
                              "pr_autumn_p_value" = NA,
                              "pr_winter_trend" = NA,
                              "pr_winter_p_value" = NA,
                              "sum_annual_crops_trend" = NA,
                              "sum_annual_crops_p_value" = NA,
                              "secdf_trend" = NA,
                              "secdf_p_value" = NA,
                              "pastr_trend" = NA,
                              "pastr_p_value" = NA,
                              "urban_trend" = NA,
                              "urban_p_value" = NA)

for(i in 1:length(species_set)){
  
  spec <- species_set[i]
  
  print(paste(i, spec)) 
  
  # relevant routes for the species, within distance of 750 km of presences:
  rel_routes <- training_routes(species = spec, buffer_km = 750, output = "RTENOs")
  
  env_dt_spec <- route_sel_env_dt %>% 
    filter(RTENO %in% rel_routes)
  
  for(v in env_vars){
    
    print(v)
    
    summary_lm <- summary(lm(unlist(env_dt_spec[,v]) ~ unlist(env_dt_spec[,"Year"])))
    lm_env_trend_df[i, paste0(v, "_trend")] <- summary_lm$coefficients[2,1]
    lm_env_trend_df[i, paste0(v, "_p_value")] <- summary_lm$coefficients[2,4]

  }

}

write.csv(lm_env_trend_df, file = file.path("results", "env_vars_trends_species.csv"), row.names = FALSE)

# plot
route_sel_env_dt %>% 
  filter(RTENO %in% rel_routes) %>% 
  ggplot(aes(x = Year, y = urban)) +
  geom_line(aes(group = RTENO)) +
  geom_smooth(method = 'lm')
