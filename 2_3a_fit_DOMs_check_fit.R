# check whether model fitting was successful based on MCMC diagnostics

# 1) saves text file with check results regarding 
# - divergent transitions
# - effective number of MCMC samples, also in bulk and tail, 
# - R-hat values < 1.02 
# 2) calls 2_3b_fit_DOMs_check_fit_details.qmd which generates pdfs with more details on MCMC diagnostics


# packages: --------------------------------------------------------------------

library(dplyr)
library(flocker)
library(sf)
library(ggplot2)
library(cmdstanr)
library(bayesplot)
library(brms)


# directories: -----------------------------------------------------------------

set_cmdstan_path("C:/Users/schifferle1/Documents/cmdstan-2.34.1")

# project directory:
dir <- file.path("//NAS-2-P-SN-01.ibb.uni-potsdam.de", "daten$", "AG26", "Transfer", "Schifferle_BBS_occupancy_models_2023")

# full models:
# first fitting round:
# results_dir <- file.path(dir, "results", "fm_buffer750km")
# second fitting round:
# results_dir <- file.path(dir, "results", "fm_buffer750km", "refit_2000_2000")

# models fitted for temporal validation:
# first fitting round:
# results_dir <- file.path(dir, "results", "temp_val_buffer_750_10yrs")
# second fitting round:
results_dir <- file.path(dir, "results", "temp_val_buffer_750_10yrs", "refit_2000_2000")


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

  quarto::quarto_render("2_3b_fit_DOMs_check_fit_details.qmd",
                        output_file = file_name,
                        output_format = "pdf",
                        execute_params = list(spec = spec,
                                              buffer_km = buffer_km,
                                              results_dir = results_dir,
                                              model_file = model_file[i]),
                        #quarto_args = c("output-dir" = results_dir), # didn't work (permissions denied) -> file saved in working directory
                        quiet = TRUE)
}
sink(file = NULL)

# save species names with fitting issues:
specs_MCMC_failed <- unique(specs_MCMC_failed)
save(specs_MCMC_failed, file = file.path(results_dir, "check_output", "specs_MCMC_failed.RData"))