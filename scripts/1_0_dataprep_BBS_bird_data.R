# Script:   1_0_dataprep_BBS_bird_data.R
# Purpose:  Preprocess data of North American Breeding Bird Survey (BBS) used to fit dynamic occupancy models
# Inputs:   downloaded BBS data, downloaded BBS routes shapefile
# Outputs:  data/BBS_species_list.csv
#           data/BBS_data_merged.RData
#           data/BBS_for_occ_spec_records.RData
#           data/BBS_for_occ.RData
#           data/route_starting_points.shp
#           data/BBS_routes_all_centroids.shp
# Runs on:  Local

# Steps:
# 1) import and merge BBS files, filter data, merge records of subspecies
# 2) reformat data for occupancy models
# 3) write shapefiles of spatial data (routes starting points and centroids)


source(file.path("scripts", "0_paths.R"))

# packages: --------------------------------------------------------------------

#devtools::install_github("trashbirdecology/bbsAssistant")
library(bbsAssistant)
library(dplyr)
library(sf)


# functions: -------------------------------------------------------------------

# modified version of bbsAssistant::import_bbs_data to import 10-stop summaries which are available
# also for years before 1997 (replaced "50-StopData.zip" with "States.zip")
# and adjusted to slight data structure changes in 2023er data release:

import_species_list_2023 <- function (bbs_dir){
  fn <- list.files(bbs_dir, full.names = TRUE, pattern = "SpeciesList")
  species_list <- readr::read_fwf(fn, skip = c(10))
  temp <- lapply(species_list[1, ], as.character)
  species_list <- species_list[-c(1:2), ]
  colnames(species_list) <- temp
  species_list$AOU <- as.integer(as.character(species_list$AOU))
  species_list <- dplyr::left_join(species_list, bbsAssistant::species_list[, c("AOU", "Scientific_Name")], by = "AOU") %>% # Seq column (Phylogenetic sequence number) has changed between releases, don't join using this column
    mutate(Scientific_Name = ifelse(is.na(Scientific_Name), paste(Genus, Species), Scientific_Name)) # checked that for each species
  return(species_list)
}

import_bbs_data_states_2023 <- function(bbs_dir) {
  
  ObsN <- RTENO <- Date <- TotalSpp  <- NULL 
  
  # where to save the unzipped files
  tempdir = tempdir()
  
  # create a vector of desired file locations.
  zipF <-
    list.files(path = bbs_dir,
               pattern = "States.zip",
               full.names = TRUE)
  utils::unzip(zipF, exdir = tempdir) # unzip the top directory
  fns.states <-
    paste0(tempdir,
           "/",
           unzip(
             zipfile = zipF,
             list = TRUE,
             exdir = tempdir()
           )$Name)

  fns.routes <- list.files(path = paste0(bbs_dir),
                           pattern = "routes.csv",
                           full.names = TRUE)
  fns.vehicle <- list.files(path = paste0(bbs_dir),
                            pattern = "ehicle",
                            full.names = TRUE)
  fns.weather <- list.files(path = paste0(bbs_dir),
                            pattern = "eather.csv",
                            full.names = TRUE)
  
  # define potential columns and desired types to ensure consistency across data files
  col_types <- readr::cols(
    AOU = readr::col_integer(),
    CountryNum = readr::col_integer(),
    Route = readr::col_character(),
    RouteDataID = readr::col_integer(),
    RPID = readr::col_integer(),
    StateNum = readr::col_integer(),
    Year = readr::col_integer()
  )
  
  # Get observations and routes ---
  
  observations <- list()
  for (i in seq_along(fns.states)) {
    f <- fns.states[i]
    observations[[i]]  <- readr::read_csv(f, col_types = col_types)
  }
  observations <- dplyr::bind_rows(observations)

  # Get species list ---
  species_list <- import_species_list_2023(bbs_dir)
  
  # Get route metadata ---
  routes <- suppressWarnings(readr::read_csv(fns.routes, col_types = col_types))
  weather <- suppressWarnings(readr::read_csv(fns.weather, col_types = col_types))
  vehicle_data <- suppressWarnings(readr::read_csv(fns.vehicle, col_types = col_types)) 
  
  observers <- weather %>%
    make.dates() %>% 
    make.rteno() %>%
    dplyr::select(ObsN, RTENO, Date, TotalSpp) %>%
    ##create binary for if observer's first year on the BBS and on the route
    dplyr::group_by(ObsN) %>% #observation identifier (number)
    dplyr::mutate(ObsFirstYearOnBBS = ifelse(Date==min(Date), 1, 0)) %>%
    dplyr::group_by(ObsN, RTENO) %>%
    dplyr::mutate(ObsFirstYearOnRTENO = ifelse(Date==min(Date), 1, 0)) %>%
    dplyr::ungroup() # to be safe
  
  # Create a list of data and information to export or return
  list.elements <-
    list("observations",
         "routes",
         "observers",
         "weather",
         "species_list",
         "vehicle_data"
    )
  bbs <- lapply(
    list.elements,
    FUN = function(x) {
      eval(parse(text = paste(x))) %>%
        make.rteno()
    }
  )
  names(bbs) <- list.elements
  
  # END FUNCTION ---
  return(bbs)
}


# BBS data prep: ---------------------------------------------------------------

# import BBS counts (aggregated to 5 sections - States.zip):
bbs_agg <- import_bbs_data_states_2023(bbs_dir = file.path(datashare_BBS))

# save species list (later used to merge species ids to species names):
bbs_agg$species_list$English_Common_Name[which(bbs_agg$species_list$Scientific_Name == "Zosterops simplex")] <- "Swinhoe's White-eye" # typo / non-UTF-8 character in SpeciesList.txt

species_list <- bbs_agg$species_list[, c("AOU", "English_Common_Name", "ORDER", "Family", "Genus", "Species", "Scientific_Name")]

# merge BBS datasets:
bbs_dt <- bbs_agg$observations %>%
  left_join(bbs_agg$routes) %>%
  left_join(bbs_agg$weather) %>%
  # only conterminous United States routes (no detailed spatial info for other BBS countries)
  filter(CountryNum == 840) %>%
  filter(StateNum != 3) %>%  # exclude Alaska
  filter(RunType == 1) %>% # Quality Control
  mutate(Stratum = factor(Stratum)) %>% 
  mutate(BCR = factor(BCR))
  
rm(bbs_agg)

# taxonomic adjustments:
bbs_dt <- bbs_dt %>% 
  left_join(species_list[, c("AOU", "English_Common_Name", "Scientific_Name", "ORDER", "Family", "Genus")],
            by = "AOU") %>% # AOU 4882 and 4890 don't appear in SpeciesList.txt -> may be typo in BBS data (80 cases)? -> deleted
  # merge records for subspecies by replacing English and scientific name:
  # based on https://www.pwrc.usgs.gov/bbs/bbsnews/aousplt1.htm (and https://www.allaboutbirds.org/guide)
  mutate(English_Common_Name = case_when(
    grepl("Northern Flicker", English_Common_Name) ~ "Northern Flicker", # includes subspecies as well as unidentified and hybrid Northern Flickers ("The red-shafted and yellow-shafted forms of the Northern Flicker 
    # formerly were considered different species." (https://www.allaboutbirds.org/guide/Northern_Flicker/overview#))
    grepl("Dark-eyed Junco", English_Common_Name) ~ "Dark-eyed Junco", # "There is a huge range of geographic variation in the Dark-eyed Junco. Among the 15 
    # described races, six forms are easily recognizable in the field and five used to be 
    # considered separate species until the 1980s" (https://www.allaboutbirds.org/guide/Dark-eyed_Junco/id)
    grepl("Yellow-rumped Warbler", English_Common_Name) ~ "Yellow-rumped Warbler", # Yellow-rumped Warbler: "The Yellow-rumped Warbler has two distinct subspecies that 
    # used to be considered separate species: the "Myrtle" Warbler of the eastern U.S. and Canada's 
    # boreal forest, and "Audubon’s" Warbler of the mountainous West." (https://www.allaboutbirds.org/guide/Yellow-rumped_Warbler/id)
    grepl("Great Blue Heron", English_Common_Name) ~ "Great Blue Heron", #(Great White Heron) Great Blue Heron: "An all-white subspecies [of Great Blue Heron], the Great White Heron, is found in coastal areas 
    # of southern Florida" (https://www.allaboutbirds.org/guide/Great_Blue_Heron/id)
    # "Eastern Towhee is assigned AOU number 05870, Spotted Towhee is AOU number 05880, and unidentified towhees are AOU number 05871. [...] 
    # These two species have largely allopatric distributions where the retroactive assignment of data was routine 
    # following the taxonomic change in 1995. Their breeding ranges potentially or actually overlap on the north-central 
    # Great Plains along a zone extending from central Nebraska north to southern Saskatchewan. 
    # Within this range of sympatry, all BBS data for this species complex obtained before 1996 is currently assigned to 
    # the unidentified towhee category. (https://www.pwrc.usgs.gov/bbs/bbsnews/aousplt1.htm) -> only few selected routes (~12) in this area of overlap -> keep species
    grepl("Red-tailed Hawk", English_Common_Name) ~ "Red-tailed Hawk", # (Harlan's Hawk) Red-tailed Hawk: subspecies of Red-tailed Hawk
    grepl("American Crow", English_Common_Name) ~ "American Crow", # (Northwestern Crow) American Crow: subspecies of American Crow
    TRUE ~ English_Common_Name)) %>%  
  mutate(Scientific_Name = case_when(
    grepl("Northern Flicker", English_Common_Name) ~ "Colaptes auratus", # includes subspecies as well as unidentified and hybrid Northern Flickers
    grepl("Dark-eyed Junco", English_Common_Name) ~ "Junco hyemalis",
    grepl("Yellow-rumped Warbler", English_Common_Name) ~ "Setophaga coronata",
    grepl("Great Blue Heron", English_Common_Name) ~ "Ardea herodias",
    grepl("Red-tailed Hawk", English_Common_Name) ~ "Buteo jamaicensis",
    grepl("American Crow", English_Common_Name) ~ "Corvus brachyrhynchos",
    TRUE ~ Scientific_Name)) %>% 
  # remove records not identified at species level and hybrids between species:
  filter(!grepl(pattern = "unid|Unid|hybrid| x ", English_Common_Name)) %>% 
  filter(!is.na(English_Common_Name))

# save merged BBS data:
save(bbs_dt, file = file.path(dir, "data", "BBS_data_merged.RData"))


# reformat BBS data for occupancy models: --------------------------------------

# from only surveyed routes to one row per route-year combination, routes either surveyed or not

# further clean data set:
bbs_dt_occ <- bbs_dt %>% 
  select(English_Common_Name, AOU, RTENO, Latitude, Longitude, BCR, Year, paste0("Count", seq(10, 50, 10)), Month, Day, ObsN) %>% 
  # convert month and date to day of year:
  mutate(date = lubridate::ymd(paste(Year, Month, Day, sep = "/"))) %>% 
  mutate(doy = lubridate::yday(date)) %>% 
  select(-c(Month, Day, date)) %>% 
  # add column on whether route was surveyed (needed later):
  mutate(Surveyed = 1) %>% 
  # site needs to be numeric:
  mutate(RTENO = as.numeric(RTENO)) %>% 
  arrange(RTENO)

save(bbs_dt_occ, file = file.path(dir, "data", "BBS_for_occ_spec_records.RData"))

# expand data to have one row per route and year:
route_dt <- tidyr::expand_grid(RTENO = unique(bbs_dt_occ$RTENO),
                               Year = min(bbs_dt_occ$Year):max(bbs_dt_occ$Year)) %>% # 224'124
  # join route data:
  collapse::join(bbs_dt_occ[, c("RTENO", "Latitude", "Longitude", "BCR")], on = c("RTENO"), how = "left") %>% 
  
  # add observer and date when route was surveyed:
  collapse::join(bbs_dt_occ[, c("RTENO", "Year", "ObsN", "doy")], on = c("RTENO", "Year"), how = "left") %>%
  # all route-year combinations without date / observer haven't been surveyed:
  mutate(Surveyed = if_else(is.na(doy), 0, 1)) %>% 
  arrange(RTENO)

save(route_dt, file = file.path(dir, "data", "BBS_for_occ.RData"))


# BBS routes - spatial data: ---------------------------------------------------

routes_sf <- read_sf(bbs_routes_file) %>% # downloaded from https://purl.stanford.edu/vy474dv5024
  st_transform(crs = "ESRI:102003") %>% # Albers Equal Area projection
  mutate(RTENO_BBS = as.integer(paste0("840", stringr::str_pad(RTENO, width = 5, side = "left", pad = "0")))) # reformat RTENO to match RTENO from BBS data imported with the bbsAssistant package

# not for all routes complete spatial route information available in this shapefile;

# 1) generate shapefile with starting points, available for all routes:

route_startpoints_sf <- bbs_dt %>%
  select(c(Latitude, Longitude, RTENO, RouteName, Active, Stratum, BCR, RouteTypeID, RouteTypeDetailID)) %>% # for ~ 40 entries different information for same starting points -> skip for now
  distinct() %>% 
  rename(Name = RouteName, RTID = RouteTypeID, RTDetID = RouteTypeDetailID) %>% # shorten names
  st_as_sf(coords = c("Longitude", "Latitude"), crs = "4269") # NAD83

# write to shapefile:
st_write(route_startpoints_sf, file.path(dir, "data", "route_starting_points.shp"), append = FALSE)

# 2) generate shapefile with centroids of routes for which full spatial information is available:

routes_sf2 <- routes_sf %>% 
  group_by(RTENO_BBS) %>% summarise %>% # merge lines if routes in shapefile consist of multiple adjacent lines
  mutate(centroid_X = NA) %>% 
  mutate(centroid_Y = NA)

# coordinates of route centroids, if route geometry is available:
for(j in 1:nrow(routes_sf2)){
  
  print(j)
  
  centroid_coords <- routes_sf2[j,] %>% 
    st_bbox %>% 
    st_as_sfc %>% 
    st_centroid %>% 
    st_coordinates
  
  routes_sf2$centroid_X[j] <- centroid_coords[1]
  routes_sf2$centroid_Y[j] <- centroid_coords[2]
}
# change geometry from route line to centroid:
routes_sf_centr <- routes_sf2 %>% 
  st_drop_geometry %>% 
  st_as_sf(coords = c("centroid_X", "centroid_Y"), crs = "ESRI:102003")

# write to shapefile:
st_write(routes_sf_centr, file.path(dir, "data", "BBS_routes_all_centroids.shp"), append = FALSE)

# session info:
writeLines(capture.output(sessionInfo()), file.path(dir, "results", "sessionInfo", "1_0_dataprep_BBS_bird_data.txt"))
