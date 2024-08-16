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
results_dir <- file.path("M:", "Documents", "DEBTs", "analysis", "Schifferle_BBS_occupancy_models_2023", 
                         "results", "full_model")

buffer_km <- 750
#fold <- 1

# log output / store initial information as text file:
MCMC_check_file <- file(file.path(results_dir, "check_output", paste0("MCMC_check_buffer_",  buffer_km, ".txt")), open = "wt") # write console output here
sink(MCMC_check_file, type = "output")

print(paste("buffer distance:", buffer_km))
print(paste("model outputs:", results_dir))

# iterate over species: ----

# selected species, sorted by ecoregion:
load(file = file.path("data", "final_species_selection_eco_sorted.RData")) # final_species_eco_sorted; output of 1_2_species_selection.R

for(i in 1:140){#length(final_species_eco_sorted)){
  
  spec <- final_species_eco_sorted[i] # 1
  
  print(paste(i, spec))
  
  # load fitted model:
  skip_to_next <- FALSE
  tryCatch(print(load(file = file.path(results_dir, paste0("out_", spec, "_fm_buffer", buffer_km, ".RData")))),
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
    print(paste("R-hat > 1.01 for", length(which(rhats > 1.01)), "parameters."))
  } else {
    print("R-hat < 1.01 for all parameters")
  }
  
  # run quarto for more detailed report:
  
  # dir_name <- file.path("C:", "Users", "schifferle1", "Documents", "DEBTs", "analysis", 
  #                       "Schifferle_BBS_occupancy_models_2023", "results", "full_model", "check_output")
  # dir_name <- file.path("M:", "Documents", "DEBTs", "analysis", 
  #                       "Schifferle_BBS_occupancy_models_2023", "results", "full_model", "check_output")
  # dir_name <- file.path("results", "full_model", "check_output")
  
  file_name <-  paste0("check_", spec, "_fm_buffer", buffer_km, ".pdf")

  quarto::quarto_render("2_2_DOM_check_flocker_fit_details.qmd", 
                        output_file = file_name, 
                        output_format = "pdf",
                        execute_params = list(spec = spec,
                                              buffer_km = buffer_km,
                                              results_dir = results_dir), 
                        #quarto_args = c("output-dir" = dir_name) # didn't work (permissions denied)
                        quiet = TRUE)
}
sink(file = NULL)