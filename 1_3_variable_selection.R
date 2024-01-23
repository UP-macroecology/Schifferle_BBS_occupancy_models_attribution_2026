# pre-filtering of variables for occupancy models:

# following König et al. 2021 and Briscoe et al. 2021

# common set of predictors for all species
# -> identify predictor subset that best captures the dimensions of environmental variation 
# and is not substantially correlated (|r| < 0.7)

# 1.) PCA of all candidate variables (data of all years rbinded) to quantify how much 
# environmental variation each variable captures
# 2.) find for each principal axis the variable that loads most on it
# 3.) rank variables according to axis (most important variable = variable that 
# loads most on first axis, second most important = variable that loads most on second axis etc.)
# 4.) find which variables are correlated >= 0.7
# 5.) of highly correlated variables choose the one with the highest rank based on PCA

# potential static and dynamic site-year covariates:
# static: bioclimatic variables based on 3 years before survey started
# dynamic: bioclimatic variables and land use class

# land use classes: sum_annual_crops, c3per, c4per, pastr, primf, primn, range, secdf, secdn, urban
# secma and secmb?

# packages:
library(ade4)
library(terra)
library(tidyverse)

# look at environment at route locations or across US? here: across US (do both?)

# mask:
mask <- rast(file.path("data", "Env_data", "ISIMIP_CHELSA-W5E5v1.0", "monthly_albers_proj", "pr_199001_ESRI102003.tif"))
values(mask)[!is.na(values(mask))] <- 1
plot(mask)

# load rasters:
# files:
bioclim_files <- list.files(file.path("data", "Env_data", "ISIMIP_CHELSA-W5E5v1.0", "bioclim"), 
                            full.names = TRUE)
# exclude static bioclims xx?
lu_files <- list.files(file.path("data", "Env_data", "LUH2", "albers_proj"), full.names = TRUE)

years <- seq(1995, 2015) # Chelsa data for historic period only until 2016, LUH2 until 2015 xx max(route_sel_dt$Year))
# store data for each year:
bioclim_lst <- vector(mode = "list", length = length(years))
landuse_lst <- vector(mode = "list", length = length(years))

for(i in 1:length(years)){
  
  print(years[i])
  
  # bioclims:
  bioclim_year_files <- bioclim_files[which(grepl(paste0("bio.{1,2}_", years[i], ".tif"), bioclim_files))]
  bioclim_year_rast <- rast(bioclim_year_files) %>% mask(mask)
  bioclim_lst[[i]] <- values(bioclim_year_rast, dataframe = TRUE) # each row = one cell
  
  # land use:
  lu_year_rast <- rast(lu_files[which(grepl(paste0(years[i], "_ESRI102003_ave.tif$"), lu_files))]) %>% mask(mask)
  landuse_lst[[i]] <- values(lu_year_rast, dataframe = TRUE) # each row = one cell, keep only rows with data (non NA cells)
}

bioclim_dt <- bind_rows(bioclim_lst)

# # add 3yrs bioclims as columns:
# bioclim_3yr_files <- bioclim_files[which(grepl("1992_1995.tif", bioclim_files))]
# bioclim_3yr_rast <- rast(bioclim_3yr_files) %>% mask(mask)
# names(bioclim_3yr_rast) <- paste0(names(bioclim_3yr_rast), "_3yrs")
# bioclim_3yr_df <- values(bioclim_3yr_rast, dataframe = TRUE)
# bioclim_dt <- cbind(bioclim_dt, bioclim_3yr_df)

lu_dt <- bind_rows(landuse_lst)

# combine bioclim and land use:
bioclim_lu_dt <- cbind(bioclim_dt, lu_dt)
# keep only complete rows (no NAs allowed in data for dudi.pca)
bioclim_lu_dt_cc <- bioclim_lu_dt[complete.cases(bioclim_lu_dt), ]
nrow(bioclim_lu_dt_cc)

# discard as predictors:
## bio8, bio9, bio18, bio19, bio8_3yrs, bio9_3yrs, bio18_3yrs, bio19_3yrs: statistical artifacts
## c3ann, c4ann, c3nfx: summarised in sum_annual_crops
bioclim_lu_dt_cc_ss <- bioclim_lu_dt_cc %>% 
  select(-c("bio8", "bio9", "bio18", "bio19", "c3ann", "c4ann", "c3nfx")) #%>% 
  #select(-c("bio8_3yrs", "bio9_3yrs", "bio18_3yrs", "bio19_3yrs"))

# PCA of bioclim and lu data (all years together):
var.pca <- dudi.pca(df = scale(bioclim_lu_dt_cc_ss),
                    nf = ncol(bioclim_lu_dt_cc_ss),# 27 number of kept axes xx
                    scannf = FALSE) # should screeplot be displayed?

# explore PCA results:
library(factoextra)
fviz_eig(var.pca)
fviz_pca_var(var.pca,
             axes = c(1, 2),
             col.var = "contrib", # Color by contributions to the PC
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
             repel = TRUE     # Avoid text overlapping
)
View(var.pca$c1) # principal axes, column normed scores -> loadings of variables, variables’ coordinates, normed to 1. 
# (l1 gives the individuals’ coordinates, normed to 1. It is also called the loadings of individuals; indiviuduals here = raster cells)
dim(var.pca$c1) # one row per original column, nf columns (1 less, why?)


# for each principal axis, select the variable with the highest loading:
maxloads <- var.pca$c1 %>% # loadings of variables, variables’ coordinates, normed to 1
  mutate(across(paste0("CS", 1:26), abs)) %>% # convert to absolute values
  summarise(across(paste0("CS", 1:26), which.max)) %>% # variable with max. values
  unlist
loads <- rownames(var.pca$c1)[maxloads[!duplicated(maxloads)]]

# rank variables: best variable gets highest rank/importance
var.loading <- rev(seq_len(length(loads)))
names(var.loading) <- loads
var.loading

# load select07 function:
source("0_functions.R")
# select variables (of correlated variables choose the one with the highest rank):
selvar <- select07(imp = var.loading[colnames(bioclim_lu_dt_cc_ss)], # importance of variables based on PCA
                   X = bioclim_lu_dt_cc_ss)
selvar # 14 variables left

ranks <- rev(seq_len(length(selvar)))
names(ranks) <- selvar
ranks

# explorations:

# 14 variables selected (of 27)
# 6 x bioclim
# 8 x land use
# most important = most environmental variation covered:

## 1.) bio3 isothermality (mean diurnal temp range / annual range):
bio3 <- rast(file.path("data", "Env_data", "ISIMIP_CHELSA-W5E5v1.0", "bioclim", "bio3_2010.tif"))
plot(bio3)  
## 2.) bio7 temperature annual range:
bio7 <- rast(file.path("data", "Env_data", "ISIMIP_CHELSA-W5E5v1.0", "bioclim", "bio7_2010.tif"))
plot(bio7) 
## 3.) sum annual crops:
sum_annual_crops <- rast(file.path("data", "Env_data", "LUH2", "albers_proj", "annual_crops_2010_ESRI102003_ave.tif"))
plot(sum_annual_crops) 
## 4.) bio15 precipitation seasonality:
bio15 <- rast(file.path("data", "Env_data", "ISIMIP_CHELSA-W5E5v1.0", "bioclim", "bio15_2010.tif"))
plot(bio15) 
## 5.) urban:
urban <- rast(file.path("data", "Env_data", "LUH2", "albers_proj", "urban_2010_ESRI102003_ave.tif"))
plot(urban)
## 6.) c4per:
c4per <- rast(file.path("data", "Env_data", "LUH2", "albers_proj", "c4per_2010_ESRI102003_ave.tif"))
plot(c4per)
## 7.) c3per:
c3per <- rast(file.path("data", "Env_data", "LUH2", "albers_proj", "c3per_2010_ESRI102003_ave.tif"))
plot(c3per)
# bio10 = mean temp. of warmest quarter
# bio14 = precipitation of driest month
# bio2 = mean diurnal range

# land use: discarded are primary forest, secondary forest, rangeland, secondary mean biomass
# rangeland: 
cor(bioclim_lu_dt_cc_ss$range, bioclim_lu_dt_cc_ss$bio3, method = "s")
# correlation 0.5


#----
route_env_dt <- read.csv(file = file.path("data", "route_year_env_data.csv"))
colnames(route_env_dt)