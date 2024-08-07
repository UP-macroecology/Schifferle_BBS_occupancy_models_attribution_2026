# Pre-filtering of variables to use as predictors in occupancy models:

# common set of predictors for all species
# -> identify predictor subset that best captures the dimensions of environmental variation 
# and is not substantially correlated (|r| < 0.7)

# roughly following König et al. 2021 and Briscoe et al. 2021 (difference in determining variable importance)

# 1.) PCA of all candidate variables (bioclimatic variables and land use data of all years rbinded) 
# to quantify how much environmental variation each variable captures
# 2.) rank variables according to importance: importance defined as the variance explained by the principal axis on 
# which the variable loads most times its loading on this axis
# 3.) find which variables are correlated >= 0.7, of highly correlated variables choose the one with the highest rank

# candidate variables:
# yearly bioclimatic variables except for bio8, bio9, bio18 and bio19
# land use classes: sum_annual_crops (sum of c3ann, c4ann, c3nfx), c3per, c4per, pastr, primf, primn, range, secdf, secdn, urban

# packages: --------------------------------------------------------------------

library(ade4)
library(terra)
library(tidyverse)
library(factoextra) # for PCA related plots

# data preparation: ------------------------------------------------------------

# mask:

# outline conterminous US:
# library(spData)
# if (requireNamespace("sf", quietly = TRUE)) {
#   data(us_states)
# }
# US_albers_sf <- us_states %>%
#   st_union() %>%
#   st_transform(crs = "ESRI:102003")
# # save as shp:
# write_sf(US_albers_sf, file.path("data", "US_outline_ESRI102003.shp"))
US_albers_sf <- read_sf(file.path("data", "US_outline_ESRI102003.shp"))


# load bioclimatic and land use rasters:

bioclim_files <- list.files(file.path("data", "Env_data", "ISIMIP_CHELSA-W5E5v1.0", "bioclim"), 
                            full.names = TRUE)
lu_files <- list.files(file.path("data", "Env_data", "LUH2", "albers_proj"), full.names = TRUE)

years <- seq(1991, 2015) # historic period of LUH2 data until 2015 (for Chelsa until 2016)

# store data of each year:
bioclim_lst <- vector(mode = "list", length = length(years))
landuse_lst <- vector(mode = "list", length = length(years))

for(i in 1:length(years)){
  
  print(years[i])
  
  # bioclims:
  bioclim_year_files <- bioclim_files[which(grepl(paste0("bio.{1,2}_", years[i], ".tif"), bioclim_files))]
  bioclim_year_rast <- rast(bioclim_year_files) %>% 
    terra::mask(US_albers_sf)
  bioclim_lst[[i]] <- values(bioclim_year_rast, dataframe = TRUE) # each row = one cell
  
  # land use:
  lu_year_rast <- rast(lu_files[which(grepl(paste0(years[i], "_ESRI102003_ave.tif$"), lu_files))]) %>% 
    terra::mask(US_albers_sf)
  landuse_lst[[i]] <- values(lu_year_rast, dataframe = TRUE) 
}

bioclim_dt <- bind_rows(bioclim_lst)
lu_dt <- bind_rows(landuse_lst)

# bioclim and land use in one df for PCA:

bioclim_lu_dt <- cbind(bioclim_dt, lu_dt)
# keep only complete rows (no NAs allowed in data for dudi.pca)
bioclim_lu_dt_cc <- bioclim_lu_dt[complete.cases(bioclim_lu_dt), ]
nrow(bioclim_lu_dt_cc)

# discard as predictors:
## bio8, bio9, bio18, bio19 due to statistical artifacts in US
## c3ann, c4ann, c3nfx: summarised in sum_annual_crops
## secma, secmb: would require interaction between these and secdf / secdn
bioclim_lu_dt_cc_ss <- bioclim_lu_dt_cc %>% 
  select(-c("bio8", "bio9", "bio18", "bio19", "c3ann", "c4ann", "c3nfx", "secma", "secmb"))

ncol(bioclim_lu_dt_cc_ss) # 25: 15 bioclim, 10 land use variables


# 1.) PCA of all candidate variables: -------------------------------------------

# PCA of bioclim and lu data (all years together):
var.pca <- dudi.pca(df = scale(bioclim_lu_dt_cc_ss),
                    nf = ncol(bioclim_lu_dt_cc_ss), # number of kept axes
                    scannf = FALSE) # display screeplot?

## plots to explore PCA results:

dir.create("plots/PCA_var_sel")

fviz_eig(var.pca) # eigenvalues (corresponds to variance explained by principal components)

#jpeg(file = file.path("plots", "PCA_var_sel", "contributions_PC1_PC2.jpg"), 
#     width = 800, height = 500, quality = 100)
fviz_pca_var(var.pca, axes = c(1, 2), col.var = "contrib",
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),repel = TRUE)
#dev.off()

#jpeg(file = file.path("plots", "PCA_var_sel", "contributions_PC3_PC4.jpg"), 
#     width = 800, height = 500, quality = 100)
fviz_pca_var(var.pca, axes = c(3, 4), col.var = "contrib",
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),repel = TRUE)
#dev.off()

#jpeg(file = file.path("plots", "PCA_var_sel", "contributions_PC5_PC6.jpg"), 
#     width = 800, height = 500, quality = 100)
fviz_pca_var(var.pca, axes = c(5, 6), col.var = "contrib",
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),repel = TRUE)
#dev.off()

View(var.pca$c1) # principal axes, column normed scores -> loadings of variables, variables’ coordinates, normed to 1. 
# (l1 gives the individuals’ coordinates, normed to 1. It is also called the loadings of individuals; indiviuduals here = raster cells)

# contributions of variables to PCs
fviz_contrib(var.pca, choice = "var", axes = 1, top = 25)
fviz_contrib(var.pca, choice = "var", axes = 2, top = 25)
fviz_contrib(var.pca, choice = "var", axes = 3, top = 25)
fviz_contrib(var.pca, choice = "var", axes = 4, top = 25)


# 2.) rank variables according to importance: ----------------------------------------------

# importance = variance explained by principal axis on which variables loads most * loading on this axis

pc_maxloads <- apply(abs(var.pca$c1), MARGIN = 1, which.max) # axis on which each variable loads most
maxloads <- apply(abs(var.pca$c1), MARGIN = 1, max) # loading on this axis
var_expl <- var.pca$eig / sum(var.pca$eig) # variance explained by each axis, same as get_eigenvalue(var.pca)

variable_importance <- maxloads * var_expl[pc_maxloads]

# rank variables: most important variable gets highest rank:
variable_rank <- rev(seq_len(length(variable_importance)))
names(variable_rank) <- names(maxloads)[order(variable_importance, decreasing = TRUE)]
variable_rank


# 3.) select variables: --------------------------------------------------------
# find which variables are correlated >= 0.7 and of highly correlated variables 
# choose the one with the highest rank

# load select07 function:
source("0_functions.R")

selvar <- select07(imp = variable_rank[colnames(bioclim_lu_dt_cc_ss)], # importance of variables
                    X = bioclim_lu_dt_cc_ss)
selvar 
# 13 variables, 7 land use + 6 bioclim:
# bioclim:
# mean diurnal temperature range (bio2), isothermality (bio3), max. temp. of warmest month (bio5),
# min. temperature of coldest month (bio6), precipitation of wettest month (bio13), precipitation of driest month (bio14)
# land use:
# annual crops, perennial C4 and C3 crops, urban, pasture, secondary non-forest, primary non-forest

# save variable selection:
save(selvar, file = file.path("data", "selected_variables.RData"))


## plot selected variables: ----------------------------------------------------

## bio2 mean diurnal temperature range:
rast(file.path("data", "Env_data", "ISIMIP_CHELSA-W5E5v1.0", "bioclim", "bio2_2010.tif")) %>% 
  plot
## bio3 isothermality:
rast(file.path("data", "Env_data", "ISIMIP_CHELSA-W5E5v1.0", "bioclim", "bio3_2010.tif")) %>% 
  plot
## bio5 max. temp. of warmest month:
rast(file.path("data", "Env_data", "ISIMIP_CHELSA-W5E5v1.0", "bioclim", "bio5_2010.tif")) %>% 
  plot
## bio6 min. temperature of coldest month:
rast(file.path("data", "Env_data", "ISIMIP_CHELSA-W5E5v1.0", "bioclim", "bio6_2010.tif")) %>% 
  plot
## bio13 precipitation of wettest month:
rast(file.path("data", "Env_data", "ISIMIP_CHELSA-W5E5v1.0", "bioclim", "bio13_2010.tif")) %>% 
  plot
## bio14 precipitation of driest month:
rast(file.path("data", "Env_data", "ISIMIP_CHELSA-W5E5v1.0", "bioclim", "bio14_2010.tif")) %>% 
  plot
## sum annual crops:
rast(file.path("data", "Env_data", "LUH2", "albers_proj", "annual_crops_2010_ESRI102003_ave.tif")) %>% 
  plot() 
## perennial C4:
rast(file.path("data", "Env_data", "LUH2", "albers_proj", "c4per_2010_ESRI102003_ave.tif")) %>% 
  plot() 
## perennial C3:
rast(file.path("data", "Env_data", "LUH2", "albers_proj", "c3per_2010_ESRI102003_ave.tif")) %>% 
  plot() 
## urban:
rast(file.path("data", "Env_data", "LUH2", "albers_proj", "urban_2010_ESRI102003_ave.tif")) %>% 
  plot() 
## pasture:
rast(file.path("data", "Env_data", "LUH2", "albers_proj", "pastr_2010_ESRI102003_ave.tif")) %>% 
  plot()
## secondary non-forest:
rast(file.path("data", "Env_data", "LUH2", "albers_proj", "secdn_2010_ESRI102003_ave.tif")) %>% 
  plot() 
## primary non-forest:
rast(file.path("data", "Env_data", "LUH2", "albers_proj", "primn_2010_ESRI102003_ave.tif")) %>% 
  plot() 


# (Misc. explorations:) --------------------------------------------------------

## other option to define variable importance following König et. al. and Briscoe et al.: ----

# importance based on principal axis with max. loading:
## most important variable = variable that loads most on PC1, then variable that loads most on PC2 etc.

# for each principal axis, select the variable with the highest loading:
maxloads2 <- var.pca$c1 %>%
  mutate(across(paste0("CS", 1:24), abs)) %>% # convert to absolute values
  summarise(across(paste0("CS", 1:24), which.max)) %>% # variable with max. values
  unlist

variable_importance2 <- rownames(var.pca$c1)[maxloads2[!duplicated(maxloads2)]] # remove duplicates (same variable may load most on multiple axis)

# rank variables: best variable gets highest rank/importance
variable_rank2 <- rev(seq_len(length(variable_importance2)))
names(variable_rank2) <- variable_importance2
variable_rank2

# select variables (of correlated variables choose the one with the highest rank):
selvar2 <- select07(imp = variable_rank2[colnames(bioclim_lu_dt_cc_ss)],
                   X = bioclim_lu_dt_cc_ss)
selvar2 # 13 variables left (7 land use, 6 bioclim)

setdiff(selvar, selvar2)
setdiff(selvar2, selvar)
# same selection except that bio6 is replaced by bio7 and bio13 is replaced by bio15 (each highly correlated)


## explore correlations: --------------------------------------------------------

# correlation matrix: 

library(corrplot)

M <- cor(bioclim_lu_dt_cc_ss, method = "s")
#jpeg(file = file.path("plots", "PCA_var_sel", "spearmans_corr_M.jpg"), 
#     width = 800, height = 800, quality = 100)
corrplot(M, method = "square", order = "hclust",
         addCoef.col = "black",
         diag = FALSE,
         tl.cex = 1,
         number.cex = 0.8,
         number.digits= 2)
#dev.off()

#jpeg(file = file.path("plots", "PCA_var_sel", "spearmans_corr_M_2.jpg"), 
#     width = 800, height = 800, quality = 100)
corrplot(M, method = "square", order = "hclust", type = "lower",
         addCoef.col = "black",
         diag = FALSE,
         tl.cex = 1,
         number.cex = 0.8,
         number.digits = 2,
         col = c(COL2('RdBu', 20)[1:3], rep("white", 14), COL2('RdBu', 20)[18:20]))
#dev.off()

# cluster Analysis:

library(Hmisc)

par(mar = c(2,5,5,2))
#jpeg(file = file.path("plots", "PCA_var_sel", "cluster_complete.jpg"), 
#     width = 700, height = 600, quality = 100)
plot(varclus(as.matrix(bioclim_lu_dt_cc_ss), similarity = "spearman",
             method = "complete")) 
abline(h = 1-(0.7^2), lty = 2, col = "grey") # |r|^2 = 0.7^2 = 0.5
abline(h = 0.1, lty = 2, col = "grey") 
#dev.off()
test <- varclus(as.matrix(bioclim_lu_dt_cc_ss), similarity = "spearman",
                method = "complete")
test


# seasonal clim. variables: ----

sclim_files <- list.files(file.path("data", "Env_data", "ISIMIP_CHELSA-W5E5v1.0", "seasonal"), 
                            full.names = TRUE)


years <- seq(1991, 2015) # historic period of LUH2 data until 2015 (for Chelsa until 2016)

# store data of each year:
sclim_lst <- vector(mode = "list", length = length(years))

for(i in 1:length(years)){
  
  print(years[i])
  
  # seasonal clims:
  sclim_year_files <- sclim_files[which(grepl(paste0("(spring|summer|autumn|winter)", "_", years[i], ".tif"), sclim_files))]
  sclim_year_rast <- rast(sclim_year_files) %>% 
    terra::mask(US_albers_sf)
  sclim_lst[[i]] <- values(sclim_year_rast, dataframe = TRUE) # each row = one cell
  
}

sclim_dt <- dplyr::bind_rows(sclim_lst)

# which variables to inspect:
# all bioclims except bio8, 9, 18, 19
# + seasonal
# + land use:

# seasonal clim., bioclim and land use in one df for PCA:

bioclim_sclim_lu_dt <- cbind(bioclim_dt, lu_dt, sclim_dt)
# keep only complete rows (no NAs allowed in data for dudi.pca)
bioclim_lu_sclim_dt_cc <- bioclim_sclim_lu_dt[complete.cases(bioclim_sclim_lu_dt), ]
nrow(bioclim_lu_sclim_dt_cc)

bioclim_lu_sclim_dt_cc_ss <- bioclim_lu_sclim_dt_cc %>% 
  select(-c("bio8", "bio9", "bio18", "bio19", "c3ann", "c4ann", "c3nfx", "secma", "secmb"))

ncol(bioclim_lu_sclim_dt_cc_ss) # 41


## 1.) PCA of all candidate variables: -------------------------------------------

var.pca <- dudi.pca(df = scale(bioclim_lu_sclim_dt_cc_ss),
                    nf = ncol(bioclim_lu_sclim_dt_cc_ss), # number of kept axes
                    scannf = FALSE) # display screeplot?

## 2.) rank variables according to importance: ----------------------------------------------

# importance = variance explained by principal axis on which variables loads most * loading on this axis

pc_maxloads <- apply(abs(var.pca$c1), MARGIN = 1, which.max) # axis on which each variable loads most
maxloads <- apply(abs(var.pca$c1), MARGIN = 1, max) # loading on this axis
var_expl <- var.pca$eig / sum(var.pca$eig) # variance explained by each axis, same as get_eigenvalue(var.pca)

variable_importance <- maxloads * var_expl[pc_maxloads]

# rank variables: most important variable gets highest rank:
variable_rank <- rev(seq_len(length(variable_importance)))
names(variable_rank) <- names(maxloads)[order(variable_importance, decreasing = TRUE)]
variable_rank


## 3.) select variables: --------------------------------------------------------
# find which variables are correlated >= 0.7 and of highly correlated variables 
# choose the one with the highest rank

# load select07 function:
source("0_functions.R")

selvar_seasonal <- select07(imp = variable_rank[colnames(bioclim_lu_sclim_dt_cc_ss)], # importance of variables
                            X = bioclim_lu_sclim_dt_cc_ss)
# save variable selection:
save(selvar_seasonal, file = file.path("data", "selected_variables_seasonal.RData"))

###

M <- cor(bioclim_lu_sclim_dt_cc_ss, method = "s")

corrplot::corrplot(M, method = "square", order = "hclust",
         addCoef.col = "black",
         diag = FALSE,
         tl.cex = 0.6,#1
         number.cex = 0.4, # 0.8
         number.digits= 2)
