# check for multiple models whether MCMC worked:
# 1) save text file with minimal information for each model
# 2) call quarto script to generate pdfs with more details for each model


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
# results_dir <- file.path("M:", "Documents", "DEBTs", "analysis", "Schifferle_BBS_occupancy_models_2023",
#                          "results", "temp_val", "round2_2000_2000")
results_dir <- file.path("M:", "Documents", "DEBTs", "analysis", "Schifferle_BBS_occupancy_models_2023",
                         "results", "full_model")
# results_dir <- file.path("//NAS-2-P-SN-01.ibb.uni-potsdam.de", "users$", "schifferle1", "Documents", "DEBTs", "analysis", 
#                          "Schifferle_BBS_occupancy_models_2023", "results", "full_model", "round2_warmup1000")
# results_dir <- file.path("M:", "Documents", "DEBTs", "analysis", "Schifferle_BBS_occupancy_models_2023",
#                          "results", "CV_cluster", "round2_2000_2000")
# results_dir <- file.path("M:", "Documents", "DEBTs", "analysis", "Schifferle_BBS_occupancy_models_2023",
#                          "results", "attribution")
# results_dir <- file.path("C:", "Users", "schifferle1", "Documents", "DEBTs", "analysis", "Schifferle_BBS_occupancy_models_2023",
#                          "results", "full_model")
buffer_km <- 750


# log output / store initial information as text file:
if(!dir.exists(file.path(results_dir, "check_output"))){
  dir.create(file.path(results_dir, "check_output"))
}

specs_MCMC_failed <- vector(mode = "character")

MCMC_check_file <- file(file.path(results_dir, "check_output", "MCMC_check3.txt"), open = "wt") # write console output here
sink(MCMC_check_file, type = "output")

print(paste("buffer distance:", buffer_km))
print(paste("model outputs:", results_dir))

# iterate over species / files: ----

# selected species, sorted by ecoregion:
#load(file = file.path("data", "final_species_selection_eco_sorted.RData")) # final_species_eco_sorted; output of 1_2_species_selection.R

# species for which models are fit:
#species_set <- final_species_eco_sorted
#species_set <- gsub(pattern = "(out_)|(_fm_buffer250.RData)", x = list.files(results_dir, pattern = "out_"), replacement = "")
#species_set <- gsub(pattern = "(out_)|(_temp_val_5yrs.RData)", x = list.files(results_dir, pattern = "out_"), replacement = "")

# models in folder:
model_file <- list.files(results_dir, pattern = "out_")

for(i in 1:length(model_file)){ 

  spec <- unlist(strsplit(model_file[i], split = "_"))[2]

  # check whether species ran already:
  if(file.exists(file.path(results_dir, "check_output", paste0("check_out_", spec, "_fm_buffer", buffer_km, ".pdf")))){
    next
  }

  print(paste(i, model_file[i]))
  
  load(file = file.path(results_dir, model_file[i]))
  
  # check MCMC:
  
  # check for divergent transitions:
  hmc_diagnostics <- nuts_params(out)
  div_trans <- sum(subset(hmc_diagnostics, Parameter == "divergent__")$Value)
  print(paste("divergent transitions:", div_trans))
  
  if(div_trans != 0) {specs_MCMC_failed <- c(specs_MCMC_failed, spec)}
  
  # is effective number of MCMC samples large enough: 
  n_eff_ratios <- neff_ratio(out)
  if("low" %in% mcmc_neff_data(n_eff_ratios)$rating){
    print("effective sample size too small")
    specs_MCMC_failed <- c(specs_MCMC_failed, spec)
  }
  
  # effective number of samples in bulk and tail:
  out_sum <- summary(out)
  n_eff_ratio_bulk <- out_sum$fixed$Bulk_ESS/out_sum$total_ndraws
  n_eff_ratio_tail <- out_sum$fixed$Tail_ESS/out_sum$total_ndraws
  
  if(min(n_eff_ratio_bulk) < 0.1 |  min(n_eff_ratio_tail) < 0.1 ){
    print("effective sample size in bulk or tail too small")
    specs_MCMC_failed <- c(specs_MCMC_failed, spec)
  } else {
    print("effective sample size in bulk and tail fine")
  }
  
  # are R-hat values fine:
  rhats <- bayesplot::rhat(out$fit)
  if(max(rhats) > 1.01){
    print(paste("R-hat > 1.01 for", length(which(rhats > 1.01)), "parameters.  R-hat max.", max(rhats)))
    specs_MCMC_failed <- c(specs_MCMC_failed, spec)
  } else {
    print("R-hat < 1.01 for all parameters")
  }
  
  # run quarto for more detailed report:
  
  # dir_name <- file.path("C:", "Users", "schifferle1", "Documents", "DEBTs", "analysis", 
  #                       "Schifferle_BBS_occupancy_models_2023", "results", "full_model", "check_output")
  # dir_name <- file.path("M:", "Documents", "DEBTs", "analysis", 
  #                       "Schifferle_BBS_occupancy_models_2023", "results", "full_model", "check_output")
  # dir_name <- file.path("results", "full_model", "check_output")
  
  #file_name <-  paste0("check_", spec, "_fm_buffer", buffer_km, ".pdf")
  file_name <-  paste0("check_", gsub(pattern = ".RData", replacement = "", x = model_file[i]), ".pdf")
  
  quarto::quarto_render("2_2_DOM_check_flocker_fit_details.qmd",
                        output_file = file_name,
                        output_format = "pdf",
                        execute_params = list(spec = spec,
                                              buffer_km = buffer_km,
                                              results_dir = results_dir,
                                              model_file = model_file[i]),
                        #quarto_args = c("output-dir" = dir_name) # didn't work (permissions denied)
                        quiet = TRUE)
}
sink(file = NULL)
specs_MCMC_failed <- unique(specs_MCMC_failed)

save(specs_MCMC_failed, file = file.path(results_dir, "check_output", "specs_MCMC_failed.RData"))
specs_MCMC_failed # 16 species
sort(specs_MCMC_failed)
# which part failed:
# Acadian Flycatcher: effective sample size too small
# American Kestrel fine (R-hat < 1.02)
# Bald Eagle: effective sample size too small
# Black-billed Magpie: divergent transitions
# Broad-winged Hawk: effective sample size too small, too large R-hat
# Brown-headed Nuthatch: divergent transitions
# Chipping Sparrow: effective sample size too small, too large R-hat
# Cooper's Hawk: effective sample size too small, too large R-hat
# Evening Grosbeak: effective sample size too small
# Golden Eagle: effective sample size too small, too large R-hat
# Hooded Warbler: effective sample size too small 
# Lincoln's Sparrow: fine (Rhat < 1.02)
# Louisiana Waterthrush: effective sample size too small, too large R-hat
# Prairie Warbler: effective sample size too small
# Sharp-shinned Hawk: effective sample size too small, too large R-hat
# Winter Wren: effective sample size too small, too large R-hat

specs_refit <- sort(specs_MCMC_failed[-which(specs_MCMC_failed %in% c("American Kestrel", "Lincoln's Sparrow"))])
save(specs_refit, file = file.path(results_dir, "check_output", "specs_refit.RData"))
