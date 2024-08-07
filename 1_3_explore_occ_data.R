# exploration regarding data selected for occupancy modelling:

# packages: 

library(dplyr)
library(sf)
library(ggplot2)


# load data: ----

# BBS route selection (centroids) to fit models:
routes_sel_sf <- st_read(file.path("data", "route_selection_1991_2015_surv_beg_end_max_5y_miss_v2_spat_thin_100km_max_30_r_per_BCR_centroids.shp")) # output of 1_1_route_selection.R

# selected species:
load(file = file.path("data", "final_species_selection.RData")) # output of 1_2_species_selection.R
species_selection_final <- sort(species_selection_final)

# route-year-species information (only surveyed)
load(file = file.path("data", "BBS_for_occ_spec_records.RData")) # output of 1_0_reformat_BBS_data.R

# selected routes and focal years matched to environmental data:
# merged route, year, environment data:
load(file = file.path("data", "route_year_env_data.RData"))


# plot routes on which species was detected at least once or never: ------------

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

# how variable are covariates across selected routes: --------------------------

# scale covariates:
route_sel_env_dt_scaled <- route_sel_env_dt_final %>% 
  mutate(across(bio14:primn_3yrs, ~ (scale(.)) %>% as.vector()))

hist(route_sel_env_dt_scaled$bio2)
hist(route_sel_env_dt_scaled$bio3)
hist(route_sel_env_dt_scaled$bio5)
hist(route_sel_env_dt_scaled$bio6)
hist(route_sel_env_dt_scaled$bio13) # not that variable
hist(route_sel_env_dt_scaled$bio14) # not variable
hist(route_sel_env_dt_scaled$bio12)

length(unique(route_sel_env_dt_final$bio2))
length(unique(route_sel_env_dt_final$bio11))

hist(route_sel_env_dt_final$bio2)
hist(route_sel_env_dt_final$bio3)
hist(route_sel_env_dt_final$bio10)
hist(route_sel_env_dt_final$bio11)
hist(route_sel_env_dt_final$bio16)
hist(route_sel_env_dt_final$bio17)


hist(route_sel_env_dt_final$bio1)
hist(route_sel_env_dt_final$bio13)
hist(route_sel_env_dt_final$bio15)
hist(route_sel_env_dt_final$bio16)
hist(route_sel_env_dt_final$bio16)
hist(route_sel_env_dt_final$bio17)



# land use:
hist(route_sel_env_dt_final$sum_annual_crops)
hist(route_sel_env_dt_final$primn)
hist(route_sel_env_dt_final$secdn)
hist(route_sel_env_dt_final$pastr)
hist(route_sel_env_dt_final$urban)
hist(route_sel_env_dt_final$c4per)
hist(route_sel_env_dt_final$c3per)

length(unique(route_sel_env_dt_final$primf))
load(file = file.path("data", "selected_variables_seasonal.RData"))
selvar_seasonal
# compare values for species presences and absences:
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
  mutate(presence = rowSums(across(paste0("Count", seq(10, 50, 10))))) %>%
  mutate(presence = ifelse(presence >= 1, 1, 0))

ggplot(occ_dt_spec) +
  facet_wrap(~factor(presence)) +
  geom_histogram(aes(x = bio4)) # bio2, bio3, bio5, bio6, bio13, bio14

# correlation looking at route locations only: ----
M <- cor(route_sel_env_dt_final[, c(8:32)], method = "s")
corrplot::corrplot(M, method = "square", order = "hclust",
         addCoef.col = "black",
         diag = FALSE,
         tl.cex = 1,
         number.cex = 0.8,
         number.digits= 2)

which(M["bio17",] > 0.7)
which(M["sum_annual_crops",] > 0.7)
which(M["c4per",] > 0.7)
which(M["c3per",] > 0.7)
which(M["urban",] > 0.7)
which(M["secdn",] > 0.7)
which(M["pastr",] > 0.7)
which(M["primn",] > 0.7)
which(M["secdf",] > 0.7)

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


# shapefile with occupancy data: ----

# routes-years:
load(file = file.path("data", "BBS_for_occ_selection.RData")) # route_sel_dt; output of 1_3_match_BBS_to_env_data.R 

# route-year-species information (only surveyed)
load(file = file.path("data", "BBS_for_occ_spec_records.RData")) # bbs_dt_occ; output of 1_0_reformat_BBS_data.R


# presences only:
presences <- bbs_dt_occ %>% 
  select(c(English_Common_Name, RTENO, Year)) %>% 
  arrange(RTENO)

# merge with route spatial data:
pres_sf <- routes_sel_sf %>% 
  #collapse::join(presences, on = c("RTENO_BBS" = "RTENO"), how = "right") %>% # something goes wrong in writing shapefile when using this
  right_join(presences, by = c("RTENO_BBS" = "RTENO")) %>% 
  select(c("RTENO" = RTENO_BBS, Year, "species" = English_Common_Name))

st_write(pres_sf, file.path("data", "spec_presences.shp"), append = FALSE)

