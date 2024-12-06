# Select BBS species for which occupancy models are calculated:

# Jarzyna et al. 2018: remove nocturnal, crepuscular, pelagic species (poorly captured by survey methodology) (no details, 494 species left)
# Harris et al. 2018: exclude all nocturnal, crepuscular, and aquatic species (not well sampled by BBS methods) (385 species left, in supplementals is code, exclude:
# same as in https://ecologicaldata.org/wiki/breeding-bird-survey-north-america: "For community analyses it is generally best to exclude 
# nocturnal, crepuscular, and aquatic species as they are not well sampled. 
# (That is, exclude AOU species codes <=2880 [waterbirds, shorebirds, etc], 
# (>=3650 & <=3810) [owls], (>=3900 & <=3910) [kingfishers], (>=4160 & <=4210) [nightjars], 7010 [dipper].)")
# Hurlbert & White 2005: excluding nocturnal, crepuscular, and otherwise difficult to survey groups (e.g., raptors) (372 species left)


# data availability for each species after route selection:
# consider only species that have been detected at >= 50 (40) different BBS routes over the whole time period
# (following Briscoe et al. 2021)


# packages: --------------------------------------------------------------------

library(dplyr)
library(sf)
library(ggplot2)

# load data: -------------------------------------------------------------------

source("0_functions.R")

start <- 1995
end <- 2019

# BBS species records:
load(file = file.path("data", "BBS_data_merged.RData")) # output of DEBTs\analysis\Schifferle_BBS_explorations_2023\BBS_data_prep.R

# species info:
BBS_species_list <- bbs_dt %>% 
  select(English_Common_Name, Scientific_Name, ORDER, Family) %>% 
  distinct

load(file = file.path("data", "BBS_for_occ_spec_records.RData")) # output of 1_0_reformat_BBS_data.R
bbs_dt_occ
nrow(bbs_dt_occ)

# selected routes:
load(file = file.path("data", "route_selection_1995_2019_surv_beg_end_max_5y_miss_v2_spat_thin_100km_max_30_r_per_BCR.RData")) # output of 1_1_route_selection.R
sel_routes_final

# BBS data, only selected routes and focal time period, join species info:
bbs_dt_occ_sel <- bbs_dt_occ %>% 
  filter(RTENO %in% sel_routes_final) %>% 
  filter(Year >= start & Year <= end) %>% 
  left_join(BBS_species_list[, c("English_Common_Name", "Scientific_Name", "ORDER", "Family")], by = c("English_Common_Name"))
nrow(bbs_dt_occ_sel)

length(unique(bbs_dt_occ_sel$English_Common_Name)) # 506 species in total


# species selection: -----------------------------------------------------------

## exclude species based on order / AOU: ----

# exclude water-related birds (waterbirds & shorebirds: AOU <=2880, kingfishers: AOU >=3900 & <=3910 , dipper: AOU 7010),
# & nocturnal birds (owls: AOU >=3650 & <=3810 and nightjars: AOU >=4160 & <=4210): 
excl_orders <- bbs_dt_occ_sel %>% 
  filter(AOU <= 2880 | (AOU >=3900 & AOU <=3910)  | AOU == 7010  | (AOU >= 3650 & AOU <= 3810)| (AOU >=4160 & AOU <= 4210)) %>%
  pull(English_Common_Name) %>% 
  unique
excl_orders # 133 (112 water-related species, 21 nocturnal species)

# how many species of which order excluded:
bbs_dt_occ_sel %>% 
  filter(English_Common_Name %in% excl_orders) %>% 
  select(c(ORDER, English_Common_Name)) %>% 
  distinct %>%
  group_by(ORDER) %>% 
  count()

# how many species of which order included:
bbs_dt_occ_sel %>% 
  filter(!English_Common_Name %in% excl_orders) %>% 
  select(c(ORDER, English_Common_Name)) %>% 
  distinct %>%
  group_by(ORDER) %>% 
  count()


## pelagic specialists and crepuscular/ nocturnal species based on EltonTraits (Wilman et al. 2014): ----

# load data:
EltonTraits <- read.csv2(file.path("data", "BirdFuncDat.txt"), header = TRUE, sep = "\t") # BirdFuncDat from the EltonTraits database: https://figshare.com/articles/dataset/Data_Paper_Data_Paper/3559887?backTo=/collections/EltonTraits_1_0_Species-level_foraging_attributes_of_the_world_s_birds_and_mammals/3306933

# pelagic specialists:
seabirds_dt <- subset(EltonTraits, PelagicSpecialist == 1)[, c("Scientific", "English")]
# which of these included in BBS:
seabirdsET_BBS <- bbs_dt_occ_sel %>% 
  filter(Scientific_Name %in% seabirds_dt$Scientific | English_Common_Name %in% seabirds_dt$English) %>% 
  pull(English_Common_Name) %>% 
  unique # 17 pelagic specialists

# nocturnal species:
nocturnal_dt <- subset(EltonTraits, Nocturnal ==1)[, c("Scientific", "English")]
# which of these included in BBS:
nocturnalET_BBS <- bbs_dt_occ_sel %>% 
  filter(Scientific_Name %in% nocturnal_dt$Scientific | English_Common_Name %in% nocturnal_dt$English) %>% 
  pull(English_Common_Name) %>% 
  unique # 19 nocturnal species (owls and nightjars)
# 2 nightjar species didn't match with EltonTraits

setdiff(seabirdsET_BBS, excl_orders) # all pelagic specialists are water-related birds
# notes:
# according to EltonTraits all BBS Alcidae, all BBS Fregatidae, most BBS Laridae except for 8 species and 2 out of 4 BBS cormorant species are considered pelagic specialists


## exclude species based on data availability: ----------------------------------

# at how many different routes is each species detected from 1991-2015:

n_routes_pres <- bbs_dt_occ_sel %>% # presences only
  #filter(!English_Common_Name %in% excl_orders) %>% 
  select(English_Common_Name, Family, Scientific_Name, ORDER, RTENO) %>%
  distinct %>% # one row per species ever recorded on a route
  group_by(English_Common_Name,  Family, Scientific_Name, ORDER) %>% 
  summarise(n_routes = n())

# exploration: number of species with N presences / N routes with presences:

# for each species number of presence across all focal years:
spec_N_total <- bbs_dt_occ_sel %>% 
  group_by(English_Common_Name, ORDER, Scientific_Name, Family) %>% 
  summarise(spec_N = n())
spec_N_total

# number of total presences from 1995-2019:
N_species_presences <- spec_N_total %>% 
  rename("N_presences" = spec_N) %>% 
  group_by(N_presences) %>% 
  summarise(N_species = n()) %>% 
  ungroup() %>% 
  mutate(spec_with_equal_or_less = cumsum(N_species)) %>% 
  mutate(spec_with_more = length(unique(bbs_dt_occ_sel$English_Common_Name)) - spec_with_equal_or_less)
N_species_presences

# number of total routes with from 1995-2019:
N_species_routes <- n_routes_pres %>% 
  group_by(n_routes) %>% 
  summarise(N_species = n()) %>% 
  ungroup() %>% 
  mutate(spec_with_equal_or_less = cumsum(N_species)) %>% 
  mutate(spec_with_more = length(unique(bbs_dt_occ_sel$English_Common_Name)) - spec_with_equal_or_less)
N_species_routes

# exclude very rare and very common species:
# exclude species that are detected at less than 50 different routes across the whole time period
# and species that have less than 50 routes where they were never detected:

excl_data_av <- n_routes_pres %>% 
  filter(n_routes < 50 | n_routes > (length(sel_routes_final)-50)) %>%
  pull(English_Common_Name)
sort(excl_data_av) # 291


## final species selection: ----

species_selection_final <- bbs_dt_occ_sel %>% 
  filter(!English_Common_Name %in% excl_orders) %>% # exclude nocturnal and water-related species (136)
  filter(!English_Common_Name %in% excl_data_av) %>% # 280 in total, 188 that are neither nocturnal nor water-related
  pull(English_Common_Name) %>% 
  unique %>% 
  sort

length(species_selection_final) # 192

# save:
save(species_selection_final, file = file.path("data", "final_species_selection.RData"))



# assign ecoregions to species: ----

# load data:

# ecoregions level 1:
# https://www.epa.gov/eco-research/ecoregions-north-america
eco_sf <- read_sf(file.path("data", "na_cec_eco_l1", "NA_CEC_Eco_Level1.shp")) %>% 
  st_transform(crs = "ESRI:102003") %>% 
  select(NA_L1NAME)

# selected routes spatial data (to buffer presences):
routes_sel_sf <- st_read(file.path("data", "route_selection_1995_2019_surv_beg_end_max_5y_miss_v2_spat_thin_100km_max_30_r_per_BCR_centroids.shp")) # output of 1_1_route_selection.R

# BBS data (which routes surveyed in which years) formatted for occupancy modelling:
load(file = file.path("data", "BBS_for_occ.RData")) # route_dt; output of 1_0_reformat_BBS_data.R
# only selected routes and focal time period:
route_sel_dt <- route_dt %>%
  filter(RTENO %in% routes_sel_sf$RTENO_BBS) %>%
  filter(Year >= 1995 & Year <= 2019) %>% 
  arrange(RTENO)
nrow(route_sel_dt) # 13475

# route-year-species information (only surveyed)
load(file = file.path("data", "BBS_for_occ_spec_records.RData")) # bbs_dt_occ; output of 1_0_reformat_BBS_data.R


# iterate over species:

spec_eco_df <- data.frame(species = species_selection_final, eco = NA)

for(i in 1:length(species_selection_final)){
  
  spec <- species_selection_final[i]
  
  print(paste(i, spec))
  
  # species presences:
  
  pres_spec <- BBS_pres_abs_spec(species = spec)
  
  pres_routes <- pres_spec %>% 
    filter(presence == 1) %>% 
    pull(RTENO) %>% 
    unique
  
  pres_routes_sf <- routes_sel_sf %>% 
    filter(RTENO_BBS %in% pres_routes) # 159
  
  # count presence per ecoregion:
  inters <- st_intersection(pres_routes_sf, eco_sf) %>% 
    st_drop_geometry() %>% 
    group_by(NA_L1NAME) %>% 
    summarise(n_routes = n())
  
  print(inters)
  
  # in which ecoregion are most presences?
  
  # are > 50% of presences in one ecoregion?
  eco_spec <- inters %>% 
    filter(n_routes > nrow(pres_routes_sf) * 0.5) %>% 
    pull(NA_L1NAME)
  
  print(eco_spec)
  
  spec_eco_df$eco[i] <- ifelse(length(eco_spec) == 0, NA, eco_spec)
  
}

spec_eco_df
table(spec_eco_df$eco)
save(spec_eco_df, file = file.path("data", "species_ecoregions.RData"))

# sort selected species by ecoregion: 
## Eastern temperate forests - Northern forests - great plains - Northwestern forested mountains - North American Deserts:

final_species_eco_sorted <- c(spec_eco_df %>% filter(eco == "EASTERN TEMPERATE FORESTS") %>% pull(species),
                              spec_eco_df %>% filter(eco == "NORTHERN FORESTS") %>% pull(species),
                              spec_eco_df %>% filter(eco == "GREAT PLAINS") %>% pull(species),
                              spec_eco_df %>% filter(eco == "NORTHWESTERN FORESTED MOUNTAINS") %>% pull(species),
                              spec_eco_df %>% filter(eco == "NORTH AMERICAN DESERTS") %>% pull(species),
                              spec_eco_df %>% filter(is.na(eco)) %>% pull(species))

save(final_species_eco_sorted, file = file.path("data", "final_species_selection_eco_sorted.RData"))


# plots to explore data availability for different species groups: -------------

# folder to store plots:
plot_dir <- file.path("plots", "spec_sel2")
if(!dir.exists(plot_dir)){dir.create(plot_dir, recursive = TRUE)}

## number of routes on which each species was detected: ---- 
jpeg(file = file.path(plot_dir, "excl_total_routes.jpg"), 
     width = 6000, height = 1000, quality = 100)
n_routes_pres %>% 
  mutate(exclude = ifelse(!English_Common_Name %in% excl_orders, 0, 1)) %>% 
  mutate(exclude = ifelse(n_routes < 50, 2, exclude)) %>% 
  mutate(exclude = ifelse(n_routes > (476-50), 2, exclude)) %>% 
  mutate(exclude = factor(exclude, levels = c(0, 1, 2))) %>% 
  ggplot(aes(x = reorder(English_Common_Name, n_routes, mean, decreasing = TRUE), y = n_routes, colour = exclude)) +
  geom_linerange(aes(ymin = 0, ymax = n_routes, colour = exclude), 
                 linewidth = 1,
                 position = position_dodge(0.3)
  ) +
  geom_point() +
  geom_hline(yintercept = 40, lty = 2) +
  geom_hline(yintercept = 50, lty = 2) +
  theme_bw() + 
  xlab("species") + ylab("N routes") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  ggtitle("0 = included, 1 = excluded water/nocturnal spec., 2 = excluded too few data")
dev.off()

## excluded vs. included based on orders, total presences: ----
jpeg(file = file.path(plot_dir, "excl_orders_total_presences.jpg"), 
     width = 6000, height = 1000, quality = 100)
spec_N_total %>% 
  mutate(exclude = ifelse(!English_Common_Name %in% excl_orders, 0, 1)) %>% 
  mutate(exclude = factor(exclude, levels = c(0, 1))) %>% 
  ggplot(aes(x = reorder(English_Common_Name, spec_N, mean, decreasing = TRUE), y = spec_N, colour = exclude)) +
  geom_linerange(aes(ymin = 0, ymax = spec_N, colour = exclude), 
                 linewidth = 1,
                 position = position_dodge(0.3)
  ) +
  geom_point() +
  geom_hline(yintercept = 100, lty = 2) +
  geom_hline(yintercept = 500, lty = 2) +
  theme_bw() + 
  xlab("species") + ylab("N presences") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  ggtitle("include vs. exclude based on orders")
dev.off()

## same for total routes on which species were detected: ----
jpeg(file = file.path(plot_dir, "excl_orders_total_routes.jpg"), 
     width = 6000, height = 1000, quality = 100)
n_routes_pres %>% 
  mutate(exclude = ifelse(!English_Common_Name %in% excl_orders, 0, 1)) %>% 
  mutate(exclude = factor(exclude, levels = c(0, 1))) %>% 
  ggplot(aes(x = reorder(English_Common_Name, n_routes, mean, decreasing = TRUE), y = n_routes, colour = exclude)) +
  geom_linerange(aes(ymin = 0, ymax = n_routes, colour = exclude), 
                 linewidth = 1,
                 position = position_dodge(0.3)
  ) +
  geom_point() +
  geom_hline(yintercept = 40, lty = 2) +
  geom_hline(yintercept = 50, lty = 2) +
  theme_bw() + 
  xlab("species") + ylab("N routes") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  ggtitle("include vs. exclude based on orders")
dev.off()

## same, but excluding only nocturnal species and pelagic specialists based on EltonTraits: ----
jpeg(file = file.path(plot_dir, "excl_noct_pel_total_presences.jpg"), 
     width = 6000, height = 1000, quality = 100)
spec_N_total %>% 
  mutate(exclude = ifelse(English_Common_Name %in% nocturnalET_BBS, "nocturnal", "include")) %>% 
  mutate(exclude = ifelse(English_Common_Name %in% seabirdsET_BBS, "pelagic", exclude)) %>% 
  ggplot(aes(x = reorder(English_Common_Name, spec_N, mean, decreasing = TRUE), y = spec_N, colour = exclude)) +
  geom_linerange(aes(ymin = 0, ymax = spec_N, colour = exclude), 
                 linewidth = 1,
                 position = position_dodge(0.3)
  ) +
  geom_point() +
  geom_hline(yintercept = 100, lty = 2) +
  geom_hline(yintercept = 500, lty = 2) +
  theme_bw() + 
  xlab("species") + ylab("N presences") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  ggtitle("include vs. exclude nocturnal species / pelagic specialists")
dev.off()

## same for total routes on which species were detected: ----
jpeg(file = file.path(plot_dir, "excl_noct_pel_total_routes.jpg"), 
     width = 6000, height = 1000, quality = 100)
n_routes_pres %>% 
  mutate(exclude = ifelse(English_Common_Name %in% nocturnalET_BBS, "nocturnal", "include")) %>% 
  mutate(exclude = ifelse(English_Common_Name %in% seabirdsET_BBS, "pelagic", exclude)) %>% 
  ggplot(aes(x = reorder(English_Common_Name, n_routes, mean, decreasing = TRUE), y = n_routes, colour = exclude)) +
  geom_linerange(aes(ymin = 0, ymax = n_routes, colour = exclude), 
                 linewidth = 1,
                 position = position_dodge(0.3)
  ) +
  geom_point() +
  geom_hline(yintercept = 40, lty = 2) +
  geom_hline(yintercept = 50, lty = 2) +
  theme_bw() + 
  xlab("species") + ylab("N routes") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  ggtitle("include vs. exclude nocturnal species / pelagic specialists")
dev.off()


## only pelagic specialists (number of routes): ----
jpeg(file = file.path(plot_dir, "pelagic_specialists_total_routes_1991-2015.jpg"), 
     width = 500, height = 500, quality = 100)
n_routes_pres %>% 
  filter(English_Common_Name %in% seabirdsET_BBS) %>%
  ggplot(aes(x = reorder(English_Common_Name, n_routes, mean, decreasing = TRUE), y = n_routes, colour = ORDER)) +
  geom_linerange(aes(ymin = 0, ymax = n_routes, colour = ORDER), 
                 linewidth = 0.8,
                 position = position_dodge(0.3)
  ) +
  geom_point() +
  geom_hline(yintercept = 40, lty = 2) +
  geom_hline(yintercept = 50, lty = 2) +
  theme_bw() + 
  xlab("species") + ylab("N routes") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  ggtitle("pelagic specialists")
dev.off()

## only nocturnal species (number of routes): ----
jpeg(file = file.path(plot_dir, "nocturnal_species_total_routes_1991-2015.jpg"), 
     width = 500, height = 500, quality = 100)
n_routes_pres %>% 
  filter(ORDER %in% c("Caprimulgiformes", "Strigiformes")) %>% # 
  ggplot(aes(x = reorder(English_Common_Name, n_routes, mean, decreasing = TRUE), y = n_routes, colour = ORDER)) +
  geom_linerange(aes(ymin = 0, ymax = n_routes, colour = ORDER), 
                 size = 0.8,
                 position = position_dodge(0.3)
  ) +
  geom_point() +
  geom_hline(yintercept = 40, lty = 2) +
  geom_hline(yintercept = 50, lty = 2) +
  theme_bw() + 
  xlab("species") + ylab("N routes") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  ggtitle("nocturnal species")
dev.off()



## only water-related species (waterbirds/shorebirds, kingfishers, dipper) (number of routes): ----
jpeg(file = file.path(plot_dir, "water_related_birds_total_routes_1991-2015.jpg"), 
     width = 1500, height = 1000, quality = 100)
n_routes_pres %>% 
  dplyr::left_join(y = distinct(bbs_dt_occ_sel[, c("English_Common_Name", "AOU")]), by = "English_Common_Name") %>%
  filter(AOU <= 2880 | (AOU >=3900 & AOU <=3910) | AOU == 7010) %>%  
  ggplot(aes(x = reorder(English_Common_Name, n_routes, mean, decreasing = TRUE), y = n_routes, colour = ORDER)) +
  geom_linerange(aes(ymin = 0, ymax = n_routes, colour = ORDER), 
                 size = 0.8,
                 position = position_dodge(0.3)
  ) +
  geom_point() +
  geom_hline(yintercept = 40, lty = 2) +
  geom_hline(yintercept = 50, lty = 2) +
  theme_bw() + 
  xlab("species") + ylab("N routes") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  ggtitle("waterbirds, shorebirds, kingfisher, dipper")
dev.off()

## pelagic specialists vs. waterbirds, shorebirds, kingfisher, dipper (number of routes): ----
jpeg(file = file.path(plot_dir, "water_related_vs_pelagic_total_routes_1991-2015.jpg"), 
     width = 1500, height = 800, quality = 100)
n_routes_pres %>% 
  dplyr::left_join(y = distinct(bbs_dt_occ_sel[, c("English_Common_Name", "AOU")]), by = "English_Common_Name") %>%
  filter(AOU <= 2880 | (AOU >=3900 & AOU <=3910) | AOU == 7010) %>%
  mutate(pelagic_specialist = factor(ifelse(English_Common_Name %in% seabirdsET_BBS, 1, 0))) %>% 
  ggplot(aes(x = reorder(English_Common_Name, n_routes, mean, decreasing = TRUE), y = n_routes, colour = pelagic_specialist)) +
  geom_linerange(aes(ymin = 0, ymax = n_routes, colour = pelagic_specialist), 
                 size = 0.8,
                 position = position_dodge(0.3)
  ) +
  geom_point() +
  geom_hline(yintercept = 40, lty = 2) +
  geom_hline(yintercept = 50, lty = 2) +
  theme_bw() + 
  xlab("species") + ylab("N routes") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  ggtitle("pelagic specialists vs. waterbirds, shorebirds, kingfisher, dipper")
dev.off()

## number of routes with presences over time for single species / families...: -----
#spec <- "Wood Duck"
spec <- "Killdeer"
#spec <- "Gray Kingbird" # only in Florida and around Cuba and other Caribbean islands
#spec <- "Olive Warbler" # only in very small area with ~ 4 routes

bbs_dt_occ_sel %>% 
  filter(English_Common_Name == spec) %>% 
  group_by(Year) %>% 
  summarise(count_year = n()) %>% 
  ggplot(aes(x = Year, y = count_year)) +
  geom_point() + 
  geom_line() +
  ggtitle(spec) + 
  theme_bw()

# Cormorants:
bbs_dt_occ_sel %>% 
  #filter(AOU <= 2880 | (AOU >=3900 & AOU <=3910) | AOU == 7010) %>% # 
  filter(Family == "Phalacrocoracidae") %>%
  #filter(English_Common_Name %in% seabirdsET_BBS) %>% 
  group_by(Year, English_Common_Name) %>% 
  summarise(count_year = n()) %>% 
  ggplot(aes(x = Year, y = count_year, color = English_Common_Name, group = English_Common_Name)) +
  geom_point() + 
  geom_line() +
  ggtitle("") + 
  theme_bw()

# Kingfishers, dipper:
bbs_dt_occ_sel %>% 
  filter((AOU >=3900 & AOU <=3910) | AOU == 7010) %>% # 
  group_by(Year, English_Common_Name) %>% 
  summarise(count_year = n()) %>% 
  ggplot(aes(x = Year, y = count_year, color = English_Common_Name, group = English_Common_Name)) +
  geom_point() + 
  geom_line() +
  ggtitle("") + 
  theme_bw()
# American dipper breeding range only in parts of the western US
# Belted Kingfisher widely distributed
# Ringed Kingfisher mainly distributed in South America