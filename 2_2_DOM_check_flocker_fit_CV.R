# write MCMC diagnostics for CV model runs to file

# packages: ----

library(dplyr)
library(flocker)
library(sf)
library(ggplot2)
library(cmdstanr)
library(bayesplot)
set_cmdstan_path("C:/Users/schifferle1/Documents/cmdstan-2.34.1")
library(brms)

# settings: ----
results_dir <- file.path("M:", "Documents", "DEBTs", "analysis", "Schifferle_BBS_occupancy_models_2023",
                         "results", "CV_cluster")
# results_dir <- file.path("//NAS-2-P-SN-01.ibb.uni-potsdam.de", "users$", "schifferle1", "Documents", "DEBTs", "analysis", 
#                          "Schifferle_BBS_occupancy_models_2023", "results", "CV_cluster")
buffer_km <- 750

# log output / store initial information as text file:
MCMC_check_file <- file(file.path(results_dir, "check_output", paste0("MCMC_check_CV_buffer_",  buffer_km, "_final_final.txt")), open = "wt") # write console output here
sink(MCMC_check_file, type = "output")

print(paste("buffer distance:", buffer_km))
print(paste("model outputs:", results_dir))

# iterate over species: ----

# selected species, sorted by ecoregion:
load(file = file.path("data", "final_species_selection_eco_sorted.RData")) # final_species_eco_sorted; output of 1_2_species_selection.R

species_set <- unique(gsub(pattern = "(out_)|(_CV_fold..RData)", x = list.files(results_dir, pattern = "out_"), replacement = ""))

# species and folds for which MCMC failed: xx
spec_folds_MCMC_fail <- vector(mode = "list", length = length(species_set))
names(spec_folds_MCMC_fail) <- species_set

for(i in 1:length(species_set)){
  
  spec <- species_set[i] # 1
  
  print(paste(i, spec))
  
  # load model fits of each fold:
  out_folds <- vector(mode = "list", length = 5)
  
  # test whether species data are there:
  skip_to_next <- FALSE
  tryCatch(print(load(file = file.path(results_dir, paste0("test_preds_", spec, "_CV_fold5.RData")))),
           error = function(e) { skip_to_next <<- TRUE})
  if(skip_to_next) { next }
  
  for(fold in 1:5){
    
    print(fold)
    
    # assemble fold results:
    load(file = file.path(results_dir, paste0("out_", spec, "_CV_fold", fold, ".RData")))
    out_folds[[fold]] <- out
  }
  
  
  # check MCMC:
  
  # check for divergent transitions:
  hmc_diagnostics <- lapply(out_folds, function(x) nuts_params(x))
  div_trans <- lapply(hmc_diagnostics, function(x) sum(subset(x, Parameter == "divergent__")$Value))
  
  print(paste("divergent transitions:", paste(div_trans, collapse = ",")))
  
  # is effective number of MCMC samples large enough: 
  n_eff_ratios <- lapply(out_folds, function(x) neff_ratio(x))
  
  if("low" %in% unlist(lapply(n_eff_ratios, function(x) mcmc_neff_data(x)$rating))){
    print("effective sample size too small")
  }
  
  # effective number of samples in bulk and tail:
  out_sum <- lapply(out_folds, summary)
  n_eff_ratio_bulk <- lapply(out_sum, function(x) x$fixed$Bulk_ESS/x$total_ndraws)
  n_eff_ratio_tail <- lapply(out_sum, function(x) x$fixed$Tail_ESS/x$total_ndraws)
  
  if(min(unlist(n_eff_ratio_bulk)) < 0.1 |  min(unlist(n_eff_ratio_tail)) < 0.1 ){
    print("effective sample size in bulk or tail too small")
  } else {
    print("effective sample size in bulk and tail fine")
  }
  
  # are R-hat values fine:
  rhats <- lapply(out_folds, function(x) bayesplot::rhat(x$fit))
  if(max(unlist(rhats)) > 1.01){
    print(paste("R-hat > 1.01 for", length(which(unlist(rhats) > 1.01)), "parameters. R-hat max.", max(unlist(rhats))))
  } else {
    print("R-hat < 1.01 for all parameters")
  }
  
  # which folds are problematic:
  n_eff_rat <- which(lapply(lapply(n_eff_ratios, function(x) mcmc_neff_data(x)$rating), function(x) "low" %in% x) == TRUE)
  n_eff_bulk <- which(lapply(n_eff_ratio_bulk, function(x) min(unlist(x))) < 0.1)
  n_eff_tail <- which(lapply(n_eff_ratio_tail, function(x) min(unlist(x))) < 0.1)
  rhat_large <- which(lapply(rhats, function(x) max(unlist(x))) >= 1.02)
  folds_MCMCfail <- sort(unique(c(n_eff_rat, n_eff_bulk, n_eff_tail, rhat_large)))
  print(paste("Folds for which MCMC fails:", paste(folds_MCMCfail, collapse = ", ")))

  spec_folds_MCMC_fail[[i]] <- folds_MCMCfail
  
  }
sink(file = NULL)

save(spec_folds_MCMC_fail, 
     file = file.path(results_dir, "check_output", paste0("MCMC_check_CV_buffer_",  buffer_km, "_MCMC_fail.RData")))



## second fitting round (single folds, single species): ----

results_dir <- file.path("M:", "Documents", "DEBTs", "analysis", "Schifferle_BBS_occupancy_models_2023",
                         "results", "CV_cluster", "round2_2000_2000")

buffer_km <- 750

MCMC_check_file <- file(file.path(results_dir, "check_output", paste0("MCMC_check_buffer_",  buffer_km, "_final.txt")), open = "wt") # write console output here
sink(MCMC_check_file, type = "output")

print(paste("buffer distance:", buffer_km))
print(paste("model outputs:", results_dir))


species_set <- gsub(pattern = "(out_)|(_CV_fold..RData)", x = list.files(results_dir, pattern = "out_"), replacement = "")

files <- list.files(results_dir, pattern = "out_")

for(i in 1:length(files)){
  
  print(i)
  print(files[i])
  
  # load fitted model:
  skip_to_next <- FALSE
  #tryCatch(print(load(file = file.path(results_dir, paste0("out_", spec, "_temp_val_5yrs.RData")))),
  #         error = function(e) { skip_to_next <<- TRUE})
  tryCatch(print(load(file = file.path(results_dir, files[i]))),
           error = function(e) { skip_to_next <<- TRUE})
  if(skip_to_next) { next }
  
  
  # check MCMC:
  
  # check for divergent transitions:
  hmc_diagnostics <- nuts_params(out)
  div_trans <- sum(subset(hmc_diagnostics, Parameter == "divergent__")$Value)
  print(paste("divergent transitions:", div_trans))
  
  # is effective number of MCMC samples large enough: 
  n_eff_ratios <- neff_ratio(out)
  if("low" %in% mcmc_neff_data(n_eff_ratios)$rating){
    print("effective sample size too small")
  }
  
  # effective number of samples in bulk and tail:
  out_sum <- summary(out)
  n_eff_ratio_bulk <- out_sum$fixed$Bulk_ESS/out_sum$total_ndraws
  n_eff_ratio_tail <- out_sum$fixed$Tail_ESS/out_sum$total_ndraws
  
  if(min(n_eff_ratio_bulk) < 0.1 |  min(n_eff_ratio_tail) < 0.1 ){
    print("effective sample size in bulk or tail too small")
  } else {
    print("effective sample size in bulk and tail fine")
  }
  
  # are R-hat values fine:
  rhats <- bayesplot::rhat(out$fit)
  if(max(rhats) > 1.01){
    print(paste("R-hat > 1.01 for", length(which(rhats > 1.01)), "parameters.  R-hat max.", max(rhats)))
  } else {
    print("R-hat < 1.01 for all parameters")
  }
}
sink(file = NULL)
