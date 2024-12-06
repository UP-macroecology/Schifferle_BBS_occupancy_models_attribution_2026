# Pre-filtering of variables to use as predictors in occupancy models:

# common set of predictors for all species
# -> identify predictor subset that best captures the dimensions of environmental variation 
# and is not substantially correlated (|r| < 0.7)

# roughly following König et al. 2021 and Briscoe et al. 2021 (difference in determining variable importance)

# 1.) PCA of all candidate variables (bioclimatic variables and land use data of all years) 
# to quantify how much environmental variation each variable captures
# 2.) rank variables according to importance: importance defined as the variance explained by the principal axis on 
# which the variable loads most times its loading on this axis
# 3.) find which variables are correlated >= 0.7, of highly correlated variables choose the one with the highest rank


# packages: --------------------------------------------------------------------

library(ade4)
library(terra)
library(tidyverse)
library(factoextra) # for PCA related plots
library(sf)

# data preparation: ------------------------------------------------------------

# mask:

# outline conterminous US:
US_albers_sf <- read_sf(file.path("data", "US_outline_ESRI102003.shp"))


# env. rasters:

# bioclimatic vars:
bioclim_files <- list.files(file.path("data", "Env_data", "ISIMIP_GSWP3_W5E5", "bioclim"), 
                            full.names = TRUE)
# seasonal climatic vars:
sclim_files <- list.files(file.path("data", "Env_data", "ISIMIP_GSWP3_W5E5", "seasonal"), 
                          full.names = TRUE)
# land use:
lu_files <- list.files(file.path("data", "Env_data", "ISIMIP_land_use_and_irrigation", "ISIMIP_LU_ESRI102003"), full.names = TRUE)

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

# bioclim and land use in one df for PCA:

env_dt <- cbind(bioclim_dt, sclim_dt, lu_dt)
# keep only complete rows (no NAs allowed in data for dudi.pca)
env_dt_cc <- env_dt[complete.cases(env_dt), ]
nrow(env_dt_cc)

# discard as predictors:
## bio8, bio9, bio18, bio19 due to statistical artifacts in US
## secondary mean age and secondary mean biomass
## primary forests (mainly at border to Canada)
env_dt_cc_ss <- env_dt_cc %>% 
  select(-c("bio8", "bio9", "bio18", "bio19", 
            "secondary_mean_age", "secondary_mean_biomass", "primary_forests"))

ncol(env_dt_cc_ss) # 39: 15 bioclim, 16 seasonal clim, 8 land use variables


# 1.) PCA of all candidate variables: -------------------------------------------

# scale variables:

env_means <- colMeans(env_dt_cc_ss)
env_sds <- apply(env_dt_cc_ss, 2, sd)
env_scale_pars <- list("center" = env_means, "scale" =  env_sds)
# save:
save(env_scale_pars, file = file.path("data", "all_US_env_dt_scale_pars.RData"))

# # same as:
# scale_pars <- scale(env_dt_cc_ss)
# scale_pars_att <- attributes(scale_pars)
# scale_center <- scale_pars_att$`scaled:center`
# scale_scale <- scale_pars_att$`scaled:scale`

env_dt_cc_ss_scaled <- scale(env_dt_cc_ss, center=env_means, scale=env_sds)


# PCA of bioclim and lu data (all years together):
var.pca <- dudi.pca(df = env_dt_cc_ss_scaled,
                    nf = ncol(env_dt_cc_ss_scaled), # number of kept axes
                    scannf = FALSE) # display screeplot?

## plots to explore PCA results:

# folder to store plots:
plot_dir <- file.path("plots", "PC_var_sel2")
if(!dir.exists(plot_dir)){dir.create(plot_dir, recursive = TRUE)}

fviz_eig(var.pca) # eigenvalues (corresponds to variance explained by principal components)

#jpeg(file = file.path(plot_dir, "contributions_PC1_PC2.jpg"), 
     #width = 800, height = 500, quality = 100)
fviz_pca_var(var.pca, axes = c(1, 2), col.var = "contrib",
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),repel = TRUE)
#dev.off()

#jpeg(file = file.path(plot_dir, "contributions_PC3_PC4.jpg"), 
#     width = 800, height = 500, quality = 100)
fviz_pca_var(var.pca, axes = c(3, 4), col.var = "contrib",
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),repel = TRUE)
#dev.off()

#jpeg(file = file.path(plot_dir, "contributions_PC5_PC6.jpg"), 
#     width = 800, height = 500, quality = 100)
fviz_pca_var(var.pca, axes = c(5, 6), col.var = "contrib",
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),repel = TRUE)
#dev.off()

# View(var.pca$c1) # principal axes, column normed scores -> loadings of variables, variables’ coordinates, normed to 1. 
# (l1 gives the individuals’ coordinates, normed to 1. It is also called the loadings of individuals; individuals here = raster cells)

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

sort(variable_importance, decreasing = TRUE)

# rank variables: most important variable gets highest rank:
variable_rank <- rev(seq_len(length(variable_importance)))
names(variable_rank) <- names(maxloads)[order(variable_importance, decreasing = TRUE)]
variable_rank


# 3.) select variables: --------------------------------------------------------
# find which variables are correlated >= 0.7 and of highly correlated variables 
# choose the one with the highest rank

# load select07 function:
source("0_functions.R")

selvar <- select07(imp = variable_rank[colnames(env_dt_cc_ss_scaled)], # importance of variables
                    X = env_dt_cc_ss_scaled)
selvar 

# we want to keep bio1, then we cannot use bio5, bio10, bio6, bio11 and all seasonal tas vars
# if we don't use bio5 and bio6 we can use bio7
# -> we'll use selvar, but replace bio6, bio5 and tasmin_mean_summer with bio1 and bio7
# test: we want to keep bio1, so we set this as most important variable:
variable_rank_bio1 <- variable_rank
variable_rank_bio1["bio1"] <- 40
selvar_bio1 <- select07(imp = variable_rank_bio1[colnames(env_dt_cc_ss_scaled)], # importance of variables
                   X = env_dt_cc_ss_scaled)

# use selvar_bio1 as final selection, reorder:
selvar_final <- c("bio1", "bio2", "bio3", "bio7", "bio14",  "bio15",
                  "pr_mean_spring", "pr_mean_summer", "pr_mean_autumn", "pr_mean_winter",
                  "urbanareas", "managed_pastures", "primary_nonforests","secondary_nonforests", "sum_annual_crops")

# save variable selection:
save(selvar_final, file = file.path("data", "selected_variables.RData"))


# Misc. explorations: --------------------------------------------------------

## other option to define variable importance following König et. al. and Briscoe et al.: ----

# importance based on principal axis with max. loading:
## most important variable = variable that loads most on PC1, then variable that loads most on PC2 etc.

# for each principal axis, select the variable with the highest loading:
maxloads2 <- var.pca$c1 %>%
  mutate(across(paste0("CS", 1:33), abs)) %>% # convert to absolute values
  summarise(across(paste0("CS", 1:33), which.max)) %>% # variable with max. values
  unlist

variable_importance2 <- rownames(var.pca$c1)[maxloads2[!duplicated(maxloads2)]] # remove duplicates (same variable may load most on multiple axis)

# rank variables: best variable gets highest rank/importance
variable_rank2 <- rev(seq_len(length(variable_importance2)))
names(variable_rank2) <- variable_importance2
variable_rank2

# select variables (of correlated variables choose the one with the highest rank):
selvar2 <- select07(imp = variable_rank2[colnames(env_dt_cc_ss)],
                   X = env_dt_cc_ss)
selvar2

setdiff(selvar, selvar2)
setdiff(selvar2, selvar)

# selvar: 
# - pr_spring, summer, autumn, winter instead of annual precipitation (bio12)
# - will use bio1 and bio7 instead of bio6 and tasmin_mean_summer -> bio7 (temp. annual range) instead of bio4 (temp seasonality)
# - secondary non-forest and primary non-forest instead of only secondary nonforests
# - bio6 and tasmin_mean_summer instead of bio1


## explore correlations: --------------------------------------------------------

# correlation matrix: 

library(corrplot)

M <- cor(env_dt_cc_ss, method = "s")
jpeg(file = file.path(plot_dir, "spearmans_corr_M.jpg"), 
     width = 1200, height = 1200, quality = 100)
corrplot(M, method = "square", order = "hclust",
         addCoef.col = "black",
         diag = FALSE,
         tl.cex = 1,#1
         number.cex = 0.8, # 0.8
         number.digits= 2)
dev.off()

M <- cor(env_dt_cc_ss[, c( "primary_nonforests", "secondary_forests",
                          "secondary_nonforests", "rangeland")], method = "s")
corrplot(M, method = "square", order = "hclust",
         addCoef.col = "black",
         diag = FALSE,
         tl.cex = 1,#1
         number.cex = 0.8, # 0.8
         number.digits= 2)

#jpeg(file = file.path(plot_dir, "spearmans_corr_M_2.jpg"), 
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
jpeg(file = file.path(plot_dir, "cluster_complete.jpg"), 
     width = 1000, height = 800, quality = 100)
plot(varclus(as.matrix(env_dt_cc_ss), similarity = "spearman",
             method = "complete")) 
abline(h = 1-(0.7^2), lty = 2, col = "grey") # |r|^2 = 0.7^2 = 0.5
abline(h = 0.1, lty = 2, col = "grey") 
dev.off()
