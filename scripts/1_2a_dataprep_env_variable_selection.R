# Pre-filtering of variables to use as covariates in occupancy models:

# common set of covariates for all species
# -> identify subset that best captures environmental variation across the conterminous USA
# and is not substantially correlated (|r| < 0.7)
# (roughly following König et al. 2021 and Briscoe et al. 2021)

# 1) PCA of all candidate variables (bioclimatic variables and land use data of all years) 
# to quantify how much environmental variation each variable captures
# 2) rank variables according to importance: importance defined as the variance explained by the principal axis on 
# which the variable loads most times its loading on this axis
# 3) find which variables are correlated >= 0.7, of highly correlated variables choose the one with the highest rank


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
US_albers_sf <- read_sf(file.path("data", "US_outline_ESRI102003.shp")) # output of 1_0_dataprep_climate.R

# bioclimatic vars:
bioclim_files <- list.files(file.path("data", "Env_data", "ISIMIP_GSWP3_W5E5", "bioclim"), 
                            full.names = TRUE) # output of 1_0_dataprep_climate.R
# seasonal climatic vars:
sclim_files <- list.files(file.path("data", "Env_data", "ISIMIP_GSWP3_W5E5", "seasonal"), 
                          full.names = TRUE) # output of 1_0_dataprep_climate.R
# land use:
lu_files <- list.files(file.path("data", "Env_data", "ISIMIP_land_use_and_irrigation", "ISIMIP_LU_ESRI102003"), full.names = TRUE) # output of 1_0_dataprep_landuse.R


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
:
#save(env_scale_pars, file = file.path("data", "all_US_env_dt_scale_pars.RData"))

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
save(selvar_final, file = file.path("data", "selected_variables.RData"))


# explorations: ----------------------------------------------------------------

# # correlation matrix: 
# 
# library(corrplot)
# 
# M <- cor(env_dt_cc_ss, method = "s")
# 
# #jpeg(file = file.path(plot_dir, "spearmans_corr_M_2.jpg"), 
# #     width = 800, height = 800, quality = 100)
# corrplot(M, method = "square", order = "hclust", type = "lower",
#          addCoef.col = "black",
#          diag = FALSE,
#          tl.cex = 1,
#          number.cex = 0.8,
#          number.digits = 2,
#          col = c(COL2('RdBu', 20)[1:3], rep("white", 14), COL2('RdBu', 20)[18:20]))
# #dev.off()
# 
# # cluster analysis:
# 
# library(Hmisc)
# 
# par(mar = c(2,5,5,2))
# jpeg(file = file.path(plot_dir, "cluster_complete.jpg"), 
#      width = 1000, height = 800, quality = 100)
# plot(varclus(as.matrix(env_dt_cc_ss), similarity = "spearman",
#              method = "complete")) 
# abline(h = 1-(0.7^2), lty = 2, col = "grey") # |r|^2 = 0.7^2 = 0.5
# abline(h = 0.1, lty = 2, col = "grey") 
# dev.off()
# 
# 
# # plot selected climatic variables of each year:
# 
# files1 <- list.files(file.path("data", "Env_data", "ISIMIP_GSWP3_W5E5", "bioclim"), full.names = TRUE)
# files2 <- list.files(file.path("data", "Env_data", "ISIMIP_GSWP3_W5E5", "seasonal"), full.names = TRUE)
# 
# #dir.create(file.path("plots", "clim_1995_2019"))
# 
# # selected climatic variables:
# sel_clim_var <- grep(pattern = "(bio)|(pr_mean)", x = selvar_final, value = TRUE) 
# 
# sel_clim_files1 <- grep(pattern = paste0(paste0(sel_clim_var, "_"), collapse = "|"), x = files1, value = TRUE)
# sel_clim_files1_years <- grep(pattern = paste0(1995:2019, collapse = "|"), x = sel_clim_files1, value = TRUE)
# sel_clim_files1_years <- grep(pattern = "_1992_1995", x = sel_clim_files1_years, value = TRUE, invert = TRUE)
# 
# sel_clim_files2 <- grep(pattern = paste0(sel_clim_var, collapse = "|"), x = files2, value = TRUE)
# sel_clim_files2_years <- grep(pattern = paste0(1995:2019, collapse = "|"), x = sel_clim_files2, value = TRUE)
# sel_clim_files2_years <- grep(pattern = "_1992_1995", x = sel_clim_files2_years, value = TRUE, invert = TRUE)
# 
# sel_clim_files <- c(sel_clim_files1_years, sel_clim_files2_years)
# 
# # scale variables:
# load(file = file.path("data", "route_env_dt_scale_pars.RData")) # output of 2_1_DOM_flocker_fit_fm.R
# env_scale_pars
# 
# for(v in sel_clim_var){
#   
#   print(v)
#   
#   cl_rast <- terra::rast(grep(pattern = paste0(v, "_"), x = sel_clim_files, value = TRUE))
#   
#   # scale:
#   cl_rast_scaled <- scale(cl_rast, center = as.numeric(env_scale_pars$center[v]), scale = as.numeric(env_scale_pars$scale[v]))
#   
#   range <- range(values(cl_rast_scaled), na.rm = TRUE)
#   
#   for(y in 1:length(1995:2019)){
#     
#     print(y)
#     
#     jpeg(file = file.path("plots", "clim_1995_2019", paste0(names(cl_rast_scaled)[1], "_", (1995:2019)[y], ".jpg")), 
#          width = 800, height = 500, quality = 100)
#     terra::plot(cl_rast_scaled[[y]], main = paste(names(cl_rast_scaled)[1], (1995:2019)[y]), range = range)
#     dev.off()
#     
#   }
# }
# 
# 
# # plot selected land use variables:
# 
# files <- list.files(res_dir_proj, full.names = TRUE)
# #dir.create(file.path("plots", "land_use_1995_2010"))
# 
# # selected land use variables:
# sel_lu_var <- c("urbanareas", "managed_pastures", "primary_nonforests", "secondary_nonforests", "sum_annual_crops")
# 
# sel_lu_files <- grep(pattern = paste0(sel_lu_var, collapse = "|"), x = files, value = TRUE)
# sel_lu_files_years <- grep(pattern = paste0(1995:2019, collapse = "_|_"), x = sel_lu_files, value = TRUE)
# 
# # scale variables:
# load(file = file.path("data", "route_env_dt_scale_pars.RData")) # output of 2_1_DOM_flocker_fit_fm.R
# env_scale_pars
# 
# for(v in sel_lu_var){
#   
#   print(v)
#   
#   lu_rast <- terra::rast(grep(pattern = v, x = sel_lu_files_years, value = TRUE))
#   
#   # scale:
#   lu_rast_scaled <- scale(lu_rast, center = as.numeric(env_scale_pars$center[v]), scale = as.numeric(env_scale_pars$scale[v]))
#   
#   range <- range(values(lu_rast_scaled), na.rm = TRUE)
#   
#   for(y in 1:length(1995:2019)){
#     
#     print(y)
#     
#     jpeg(file = file.path("plots", "land_use_1995_2010", paste0(names(lu_rast_scaled)[1], "_", (1995:2019)[y], ".jpg")), 
#          width = 800, height = 500, quality = 100)
#     terra::plot(lu_rast_scaled[[y]], main = paste(names(lu_rast_scaled)[1], (1995:2019)[y]), range = range)
#     dev.off()
#     
#   }
# }