# Script:   2_3a_fit_DOMs_check_fit.R
# Purpose:  Check whether dynamic occupancy model fitting was successful based on MCMC diagnostics
# Inputs:   results/fm_buffer750km/out_<species>_fm_buffer750.RData
#           results/fm_buffer750km/refit_2000_2000/out_<species>_fm_buffer750.RData
#           results/temp_val_buffer_750_10yrs/out_<species>_temp_val_10yrs_buffer_750.RData
#           results/temp_val_buffer_750_10yrs/refit_2000_2000/out_<species>_temp_val_10yrs_buffer_750.RData
# Outputs:  results/fm_buffer750km/check_output/specs_MCMC_failed.RData
#           results/fm_buffer750km/refit_2000_2000/check_output/specs_MCMC_failed.RData
#           results/fm_buffer750km/check_output/check_<input-file-name>.pdf (one file per species)
#           results/fm_buffer750km/refit_2000_2000/check_output/check_<input-file-name>.pdf
#           results/temp_val_buffer_750_10yrs/check_output/specs_MCMC_failed.RData
#           results/temp_val_buffer_750_10yrs/refit_2000_2000/check_output/specs_MCMC_failed.RData
#           results/temp_val_buffer_750_10yrs/check_output/check_<input-file-name>.pdf (one file per species)
#           results/temp_val_buffer_750_10yrs/refit_2000_2000/check_output/check_<input-file-name>.pdf
# Runs on:  Local
# Notes:    this script calls 2_3b_fit_DOMs_check_fit_details.qmd, which generates pdfs with details on MCMC diagnostics for model for each species
#           this script is run four times: 
#           1) check fit of dynamic occupancy models fitted with all data and 1000 iterations
#           2) check fit of dynamic occupancy models fitted with all data and 2000 iterations
#           3) check fit of dynamic occupancy models fitted with subset of years for temporal validation and 1000 iterations
#           3) check fit of dynamic occupancy models fitted with subset of years for temporal validation and 2000 iterations
#           -> set results directory accordingly:


source(file.path("scripts", "0_paths.R"))


# set results directory:

# models fitted with all data:
results_dir <- file.path(dir, "results", "fm_buffer750km") # first fitting round
# results_dir <- file.path(dir, "results", "fm_buffer750km", "refit_2000_2000") # second fitting round:

# models fitted for temporal validation:
# results_dir <- file.path(dir, "results", "temp_val_buffer_750_10yrs") # first fitting round
# results_dir <- file.path(dir, "results", "temp_val_buffer_750_10yrs", "refit_2000_2000") # second fitting round


# packages: --------------------------------------------------------------------

library(dplyr)
library(flocker)
library(sf)
library(ggplot2)
library(cmdstanr)
library(bayesplot)
library(brms)


# directories: -----------------------------------------------------------------

set_cmdstan_path(path = NULL) # for HPC; local: set_cmdstan_path("C:/Users/schifferle1/Documents/cmdstan-2.34.1")


# MCMC diagnostics: ------------------------------------------------------------

# save results of MCMC diagostic checks as text file:
if(!dir.exists(file.path(results_dir, "check_output"))){dir.create(file.path(results_dir, "check_output"))}
MCMC_check_file <- file(file.path(results_dir, "check_output", "MCMC_check_res.txt"), open = "wt")
sink(MCMC_check_file, type = "output")

# save species names with fitting issues:
specs_MCMC_failed <- vector(mode = "character")

print(results_dir)

# fitted models:
model_file <- list.files(results_dir, pattern = "out_")

# buffer size:
buffer_km <- 750


# iterate over species / files: 

for(i in 1:length(model_file)){ 

  spec <- unlist(strsplit(model_file[i], split = "_"))[2]

  # check whether check ran already:
  if(file.exists(file.path(results_dir, "check_output", paste0("check_out_", spec, "_fm_buffer", buffer_km, ".pdf")))){
    next
  }

  print(paste(i, spec))
  
  # load model:
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
  if(max(rhats) > 1.02){
    print(paste("R-hat > 1.02 for", length(which(rhats > 1.02)), "parameters.  R-hat max.", max(rhats)))
    specs_MCMC_failed <- c(specs_MCMC_failed, spec)
  } else {
    print("R-hat < 1.02 for all parameters")
  }
  
  # run quarto script that generated pdf with more detailed report:

  file_name <-  paste0("check_", gsub(pattern = ".RData", replacement = "", x = model_file[i]), ".pdf")

  quarto::quarto_render("scripts/2_3b_fit_DOMs_check_fit_details.qmd",
                        output_file = file_name,
                        output_format = "pdf",
                        execute_params = list(spec = spec,
                                              buffer_km = buffer_km,
                                              results_dir = results_dir,
                                              model_file = model_file[i]),
                        #quarto_args = c("output-dir" = results_dir), # didn't work (permissions denied) -> file saved in working directory
                        quiet = FALSE)
}
sink(file = NULL)

# save species names with fitting issues:
specs_MCMC_failed <- unique(specs_MCMC_failed)
save(specs_MCMC_failed, file = file.path(results_dir, "check_output", "specs_MCMC_failed.RData"))

# session info:
writeLines(capture.output(sessionInfo()), file.path(dir, "results", "sessionInfo", "2_3a_fit_DOMs_check_fit.txt"))
