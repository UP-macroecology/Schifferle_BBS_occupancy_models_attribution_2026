# Script:   1_2a_dataprep_env_variable_selection.R
# Purpose:  Filtering of variables to use as covariates in dynamic occupancy models
# Inputs:   data/US_outline_ESRI102003.shp 
#           <clim_path>/bioclim/<var>_<year>.tif
#           <clim_path>/seasonal/<var>_<year>.tif
#           <lu_path>/ISIMIP_LU_ESRI102003/<var>_<year>_ESRI102003.tif
# Outputs:  data/selected_variables.RData"
# Runs on:  Local

# Steps:
# common set of covariates for all species
# -> identify subset that best captures environmental variation across the conterminous USA
# and is not substantially correlated (|r| < 0.7)
# (roughly following König et al. 2021 and Briscoe et al. 2021)
# 1) PCA of all candidate variables (bioclimatic variables and land use data of all years) 
# to quantify how much environmental variation each variable captures
# 2) rank variables according to importance: importance defined as the variance explained by the principal axis on 
# which the variable loads most times its loading on this axis
# 3) find which variables are correlated >= 0.7, of highly correlated variables choose the one with the highest rank


source(file.path("scripts", "0_paths.R"))


# packages: --------------------------------------------------------------------

library(sf)
library(terra)
library(ade4)
library(dplyr)
library(factoextra) # for PCA related plots


# functions: -----------------------------------------------------------------

source(file.path("scripts", "0_functions.R")) # to get select07 function


# load data: -------------------------------------------------------------------

# outline conterminous US:
US_albers_sf <- read_sf(file.path(dir, "data", "US_outline_ESRI102003.shp")) # output of 1_0_dataprep_climate.R

# bioclimatic vars:
bioclim_files <- list.files(file.path(clim_path, "bioclim"), full.names = TRUE) # output of 1_0_dataprep_climate.R

# seasonal climatic vars:
sclim_files <- list.files(file.path(clim_path, "seasonal"), full.names = TRUE) # output of 1_0_dataprep_climate.R

# land use:
lu_files <- list.files(file.path(lu_path, "ISIMIP_LU_ESRI102003"), full.names = TRUE) # output of 1_0_dataprep_landuse.R


# 1) PCA of all candidate variables: -------------------------------------------

# restructure the data:

# iterate over years: 
years <- seq(1995, 2019)

# store data of each year:
bioclim_lst <- sclim_lst <- landuse_lst <- vector(mode = "list", length = length(years))

for(i in 1:length(years)){
  
  print(years[i])
  
  # bioclims:
  bioclim_year_files <- bioclim_files[which(grepl(paste0("bio.{1,2}_", years[i], ".tif"), bioclim_files))]
  bioclim_year_rast <- rast(bioclim_year_files)
  bioclim_lst[[i]] <- values(bioclim_year_rast, dataframe = TRUE) # each row = one cell
  
  # seasonal climate:
  sclim_year_files <- sclim_files[which(grepl(paste0("(spring|summer|autumn|winter)", "_", years[i], ".tif"), sclim_files))]
  sclim_year_rast <- rast(sclim_year_files)
  sclim_lst[[i]] <- values(sclim_year_rast, dataframe = TRUE) # each row = one cell
  
  # land use:
  lu_year_rast <- rast(lu_files[which(grepl(paste0(years[i], "_ESRI102003.tif$"), lu_files))])
  landuse_lst[[i]] <- values(lu_year_rast, dataframe = TRUE) 
}

bioclim_dt <- bind_rows(bioclim_lst)
sclim_dt <- bind_rows(sclim_lst)
lu_dt <- bind_rows(landuse_lst)

# bioclimatic variables and land use in one dataframe for PCA:

env_dt <- cbind(bioclim_dt, sclim_dt, lu_dt)

# keep only complete rows (no NAs allowed in data for dudi.pca)
env_dt_cc <- env_dt[complete.cases(env_dt), ]
nrow(env_dt_cc)

# discard as predictors:
## bio8, bio9, bio18, bio19 due to statistical artifacts in the conterminous USA
## secondary mean age and secondary mean biomass
## primary forests (mainly only at border to Canada)
env_dt_cc_ss <- env_dt_cc %>% 
  select(-c("bio8", "bio9", "bio18", "bio19", 
            "secondary_mean_age", "secondary_mean_biomass", "primary_forests"))

ncol(env_dt_cc_ss) # 38: 15 bioclim, 16 seasonal clim, 7 land use variables


# scale variables:

env_means <- colMeans(env_dt_cc_ss)
env_sds <- apply(env_dt_cc_ss, 2, sd)
env_scale_pars <- list("center" = env_means, "scale" =  env_sds)

env_dt_cc_ss_scaled <- scale(env_dt_cc_ss, center=env_means, scale=env_sds)


# PCA (all years together):
var.pca <- dudi.pca(df = env_dt_cc_ss_scaled,
                    nf = ncol(env_dt_cc_ss_scaled), # number of kept axes
                    scannf = FALSE)

# ## explore PCA results:
# 
# fviz_eig(var.pca) # eigenvalues (corresponds to variance explained by principal components)
# 
# fviz_pca_var(var.pca, axes = c(1, 2), col.var = "contrib",
#              gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),repel = TRUE)
# 
# 
# fviz_pca_var(var.pca, axes = c(3, 4), col.var = "contrib",
#              gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),repel = TRUE)
# 
# 
# fviz_pca_var(var.pca, axes = c(5, 6), col.var = "contrib",
#              gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),repel = TRUE)
# 
# # View(var.pca$c1) # principal axes, column normed scores -> loadings of variables, variables’ coordinates, normed to 1. 
# # (l1 gives the individuals’ coordinates, normed to 1, called the loadings of individuals; individuals here = raster cells)
# 
# # contributions of variables to PCs
# fviz_contrib(var.pca, choice = "var", axes = 1, top = 25)
# fviz_contrib(var.pca, choice = "var", axes = 2, top = 25)
# fviz_contrib(var.pca, choice = "var", axes = 3, top = 25)
# fviz_contrib(var.pca, choice = "var", axes = 4, top = 25)



# 2) rank variables according to importance: ----------------------------------------------

# importance = variance explained by principal axis on which variables loads most * loading on this axis

pc_maxloads <- apply(abs(var.pca$c1), MARGIN = 1, which.max) # axis on which each variable loads most
maxloads <- apply(abs(var.pca$c1), MARGIN = 1, max) # loading on this axis
var_expl <- var.pca$eig / sum(var.pca$eig) # variance explained by each axis, same as get_eigenvalue(var.pca)

variable_importance <- maxloads * var_expl[pc_maxloads]

sort(variable_importance, decreasing = TRUE)

# rank variables: most important variable gets highest rank:
variable_rank <- rev(seq_len(length(variable_importance)))
names(variable_rank) <- names(maxloads)[order(variable_importance, decreasing = TRUE)]
variable_rank


# 3) select variables: --------------------------------------------------------
# find which variables are correlated >= 0.7 and of highly correlated variables 
# choose the one with the highest rank

selvar <- select07(imp = variable_rank[colnames(env_dt_cc_ss_scaled)], # importance of variables
                    X = env_dt_cc_ss_scaled)
selvar 

# we prioritise bio1, setting this as most important variable yields:
variable_rank_bio1 <- variable_rank
variable_rank_bio1["bio1"] <- 40
selvar_bio1 <- select07(imp = variable_rank_bio1[colnames(env_dt_cc_ss_scaled)], # importance of variables
                   X = env_dt_cc_ss_scaled)
sort(selvar_bio1)

# use this as final selection, reorder:
selvar_final <- c("bio1", "bio2", "bio3", "bio7", "bio14",  "bio15",
                  "pr_mean_spring", "pr_mean_summer", "pr_mean_autumn", "pr_mean_winter",
                  "urbanareas", "managed_pastures", "primary_nonforests","secondary_nonforests", "sum_annual_crops")

# save variable selection:
save(selvar_final, file = file.path(dir, "data", "selected_variables.RData"))

# session info:
writeLines(capture.output(sessionInfo()), file.path(dir, "results", "sessionInfo", "1_2a_dataprep_env_variable_selection.txt"))
