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

# results_dir <- file.path("//NAS-2-P-SN-01.ibb.uni-potsdam.de", "users$", "schifferle1", "Documents", "DEBTs", "analysis",
#                          "Schifferle_BBS_occupancy_models_2023", "results", "CV_buffer750km")
# results_dir <- file.path("M:", "Documents", "DEBTs", "analysis", "Schifferle_BBS_occupancy_models_2023",
#                          "results", "CV_buffer750km")
results_dir <- file.path("M:", "Documents", "DEBTs", "analysis", "Schifferle_BBS_occupancy_models_2023",
                         "results", "CV_buffer750km", "refit_2000_2000")
buffer_km <- 750

# log output / store initial information as text file:

if(!dir.exists(file.path(results_dir, "check_output"))){
  dir.create(file.path(results_dir, "check_output"))
}

MCMC_check_file <- file(file.path(results_dir, "check_output", "MCMC_check_res.txt"), open = "wt") # write console output here
sink(MCMC_check_file, type = "output")

print(paste("buffer distance:", buffer_km))
print(paste("model outputs:", results_dir))

# iterate over species: ----

# selected species, sorted by ecoregion:
#load(file = file.path("data", "final_species_selection_eco_sorted.RData")) # final_species_eco_sorted; output of 1_2_species_selection.R

species_set <- unique(gsub(pattern = "(out_)|(_CV_fold..RData)", x = list.files(results_dir, pattern = "out_"), replacement = ""))

# species and folds for which MCMC failed:
spec_folds_MCMC_fail <- vector(mode = "list", length = length(species_set))
names(spec_folds_MCMC_fail) <- species_set

for(i in 1:length(species_set)){
  
  spec <- species_set[i] # 1
  
  print(paste(i, spec))
  
  # load model fits of each fold:
  out_folds <- vector(mode = "list", length = 5)
  
  # # test whether species data are there:
  # skip_to_next <- FALSE
  # tryCatch(print(load(file = file.path(results_dir, paste0("test_preds_", spec, "_CV_fold5.RData")))),
  #          error = function(e) { skip_to_next <<- TRUE})
  # if(skip_to_next) { next }
  
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
  
  if(any(div_trans != 0)) {
    spec_folds_MCMC_fail[[i]] <- which(div_trans != 0)
  }
  
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
  if(max(unlist(rhats)) > 1.02){
    print(paste("R-hat > 1.02 for", length(which(unlist(rhats) > 1.02)), "parameters. R-hat max.", max(unlist(rhats))))
  } else {
    print("R-hat <= 1.02 for all parameters")
  }
  
  # which folds are problematic:
  n_eff_rat <- which(lapply(lapply(n_eff_ratios, function(x) mcmc_neff_data(x)$rating), function(x) "low" %in% x) == TRUE)
  n_eff_bulk <- which(lapply(n_eff_ratio_bulk, function(x) min(unlist(x))) < 0.1)
  n_eff_tail <- which(lapply(n_eff_ratio_tail, function(x) min(unlist(x))) < 0.1)
  rhat_large <- which(lapply(rhats, function(x) max(unlist(x))) >= 1.02)
  folds_MCMCfail <- sort(unique(c(n_eff_rat, n_eff_bulk, n_eff_tail, rhat_large)))
  print(paste("Folds for which MCMC fails:", paste(folds_MCMCfail, collapse = ", ")))

  spec_folds_MCMC_fail[[i]] <- unique(c(spec_folds_MCMC_fail[[i]], folds_MCMCfail))
  
  }
sink(file = NULL)

save(spec_folds_MCMC_fail, 
     file = file.path(results_dir, "check_output", "specs_folds_MCMC_failed.RData"))

names(which(lengths(spec_folds_MCMC_fail) != 0))


# have species been excluded in previous steps (fitting full model, fitting model for temporal validation):
# temporal validation:

# refit for species not excluded in the previous step:

if(results_dir == "M:/Documents/DEBTs/analysis/Schifferle_BBS_occupancy_models_2023/results/CV_buffer750km"){
  
  load(file.path("M:", "Documents", "DEBTs", "analysis", "Schifferle_BBS_occupancy_models_2023",
                 "results", "fm_buffer750km", "refit_2000_2000", "check_output", "specs_MCMC_failed.RData"))
  specs_discard_fm <- specs_MCMC_failed
  
  load(file.path("M:", "Documents", "DEBTs", "analysis", "Schifferle_BBS_occupancy_models_2023",
                 "results", "temp_val_buffer_750_10yrs", "refit_2000_2000", "check_output", "specs_MCMC_failed.RData"))
  specs_discard_tv <- specs_MCMC_failed
  
  load(file.path(results_dir, "check_output", "specs_folds_MCMC_failed.RData"))
  specs_MCMC_failed_cv <- names(which(lengths(spec_folds_MCMC_fail) != 0))
    
  species_refit <- subset(specs_MCMC_failed_cv, !specs_MCMC_failed_cv %in% c(specs_discard_fm, specs_discard_tv))
  species_refit
  save(species_refit, file = file.path(results_dir, "check_output", "species_refit.RData"))
}
