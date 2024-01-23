# occupancy null model for 5 test species and 7 versions of routes included

# packages: ----

library(dplyr)
library(tidyr)
library(collapse)
library(jagsUI)
library(doParallel)

# load data: ----

## BBS data:
load(file = file.path("data", "BBS_data_merged.RData")) # output of DEBTs\analysis\Schifferle_BBS_explorations_2023\BBS_data_prep.R

# all identified species:
species_common <- sort(unique(bbs_dt$English_Common_Name)) # 571

# reduce dataset to necessary columns:
bbs_dt_occ <- bbs_dt %>% 
  select(English_Common_Name, AOU, RTENO, Latitude, Longitude, BCR, Year, paste0("Count", seq(10, 50, 10)), Month, Day, ObsN) %>% 
  # convert month and date to day of year:
  mutate(date = lubridate::ymd(paste(Year, Month, Day, sep = "/"))) %>% 
  mutate(doy = lubridate::yday(date)) %>% 
  select(-c(Month, Day, date)) %>% 
  # add column on whether route was surveyed (needed later):
  mutate(Surveyed = 1) %>% 
  # site needs to be numeric:
  mutate(RTENO = as.numeric(RTENO))

# reformat BBS data for occupancy modelling:

# expand data to have one row per route and year:
route_dt <- tidyr::expand_grid(RTENO = unique(bbs_dt_occ$RTENO),
                                  Year = min(bbs_dt_occ$Year):max(bbs_dt_occ$Year)) %>% # 224'124
  # join route data:
  collapse::join(bbs_dt_occ[, c("RTENO", "Latitude", "Longitude", "BCR")], on = c("RTENO"), how = "left") %>% 
  
  # add observer and date when route was surveyed:
  collapse::join(bbs_dt_occ[, c("RTENO", "Year", "ObsN", "doy")], on = c("RTENO", "Year"), how = "left") %>%
  # all route-year combinations without date / observer haven't been surveyed:
  mutate(Surveyed = if_else(is.na(doy), 0, 1))  %>% 
  # use only later years to check whether model fitting works then:
  filter(Year >= 1995 & Year <= 2019) # no data in 2020 due to covid

nyears <- length(unique(route_dt$Year)) 
nsurveys <- 5

# only use data of routes that provide enough temporal and approx. even spatial coverage: in first section xx
load(file.path("data", "route_selection_25ys_surv_beg_end_max_5y_miss_max_30_r_per_BCR_v2.RData")) # output of 1_route_selection.R
sel_routes_final 

route_dt_ss <- route_dt %>% 
  filter(RTENO %in% sel_routes_final)

nsites <- length(unique(route_dt_ss$RTENO))

# model definition: -----

# simplest dynamic occupancy model: I assume that all parameters are constant across time and space, no covariates
# constant colonisation, extinction, detection

cat(file = "dynoccmod1.txt"," #xxDEBTs/analysis/Schifferle_BBS_occupancy_models_2023/
model {
  
  # Specify priors ---------------------------------
  
  # Initial occupancy priors:
  
	psi1 ~ dunif(0, 1)
  
  # Colonisation and extinction priors:
  
  eps ~ dunif(0, 1)
  gamma ~ dunif(0, 1)

  # Detection priors:
 	p ~ dunif(0, 1)

  # Likelihood: ---------
  
  # Ecological submodel:
  
  # Define occupancy state conditional on parameters
  # state transitions:
  
  for (s in 1:nsites){
    z[s,1] ~ dbern(psi1)
    
    for (t in 2:nyears){
      z[s,t] ~ dbern(z[s,t-1]*(1-eps) + (1-z[s,t-1])*gamma)
    }
  }

  # Observation model:
  
  for (s in 1:nsites){
      for (j in 1:nsurveys){
        for (t in 1:nyears){
              y[s,j,t] ~ dbern(p * z[s,t])
        }
      }
  }

  # Derived parameters:
  
  # Population occupancy in each year:
  psi[1] <- psi1                
  for (t in 2:nyears){
    psi[t] <- psi[t-1]*(1-eps) + (1-psi[t-1])*gamma
  }
}
")


# iterate over species: ----

testspecs <- c("American Goldfinch",
               "White-breasted Nuthatch",
               "Winter Wren",
               "House Finch",
               "Eastern Meadowlark")

# register cores for parallel computation:
ncores <- 15 # 15: 5 species * 3 chains? # xx
cl <- makeCluster(ncores, setup_timeout = 0.5)
registerDoParallel(cl)

# new log file (clean existing one):
logfile <- file.path("job_log.txt")
writeLines(c(""), logfile)

res_metrics_all <- foreach(spec = testspecs,
                           .combine = rbind,
                           .packages = c("dplyr", "collapse", "jagsUI"),
                           .errorhandling = "remove",
                           .verbose = TRUE) %dopar% {

  ## assemble data: ----
  
  occ_dt_spec <- route_dt_ss %>% 
    # add species information:
    mutate(English_Common_Name = spec, 
           AOU = unique(bbs_dt_occ$AOU[which(bbs_dt_occ$English_Common_Name == spec)])) %>%
    # add observations:
    collapse::join(bbs_dt_occ[, c("RTENO", "Year", "AOU", paste0("Count", seq(10, 50, 10)))], 
                   on = c("RTENO", "Year", "AOU"), how = "left") %>% 
    # if route was surveyed but species not observed, replace NA with 0:
    mutate(across(Count10:Count50, ~ 
                    case_when(Surveyed == 1 & is.na(.) ~ 0,
                              .default = .))) %>% # TRUE ~ . for older dplyr version
    # convert bird counts to presence / absence:
    mutate(across(Count10:Count50, ~ 
                    case_when(. > 1 ~ 1,
                              .default = .)))

  ###
  # to fasten model testing subsample data by using random 30 % (following Kery and Royle): xx
  # occ_dt_spec_route_ss_full <- occ_dt_spec_route_ss                      # Make a copy of full data set
  # prop.data <- 0.3                     # Proportion of data to be used
  # ncase <- nrow(occ_dt_spec_route_ss)                   # 116204
  # set.seed(1)                          # Ensures you get the same subset
  # sel.cases <- sort(sample(1:ncase, ncase * prop.data))
  # occ_dt_spec_route_ss <- occ_dt_spec_route_ss[sel.cases,] # Smaller data set
  ### 

  
  # data frame to store results:
  res_metrics <- data.frame("routes_subset" = "25ys_surv_beg_end_max_5y_miss_max_30_r_per_BCR_v2",
                            "species" = spec,
                            "Rhat_fine" = NA,
                            "n_routes" = NA,
                            "psiobs" = NA,
                            "DIC" = NA, "pD" = NA, 
                            "eps_mean" = NA, "eps_sd" = NA, "eps_2.5" = NA, "eps_97.5" = NA,
                            "gamma_mean" = NA, "gamma_sd" = NA, "gamma_2.5" = NA, "gamma_97.5" = NA,
                            "p_mean" = NA, "p_sd" = NA, "p_2.5" = NA, "p_97.5" = NA,
                            "psi1_mean" = NA, "psi1_sd" = NA, "psi1_2.5" = NA, "psi1_97.5" = NA)
  
  # compute observed occupancy:
  obs_occupancy <- occ_dt_spec %>% 
    mutate(presence = rowSums(pick(Count10:Count50))) %>%
    mutate(presence = if_else(presence > 1, 1, 0)) %>% 
    summarise(psiobs = mean(presence, na.rm = TRUE)) %>%  # number of route-year comb. with presence / number of all surveyed route-year comb.
    pull(psiobs)
  
  res_metrics$psiobs <- round(obs_occupancy,2)
  
  # reformat observations to array for model:
  years <- seq(min(occ_dt_spec$Year), max(occ_dt_spec$Year))
  y_array <- array(NA, dim = c(nsites, nsurveys, nyears))
  for (t in 1:length(years)){
    y_array[1:nsites, 1:nsurveys, t] <- as.matrix(occ_dt_spec[which(occ_dt_spec$Year == years[t]), 
                                                                       c(paste0("Count", seq(10, 50, 10)))])
  }
  
  # fit model: --------------------------------------------------------------

  # bundle data:
  str(bdata <- list(y = y_array, 
                    nsurveys = nsurveys,
                    nsites = nsites, 
                    nyears = nyears))
  
  # Initial values
  # for z = occupancy status at site s at year t (see model definition)
  # use observed occupancy status per route and year as initial values for z[s, t]:
  zst <- apply(bdata$y, c(1, 3), max) # max. over all surveys
  inits <- function(){ list(z = zst)}
  
  # Parameters monitored
  params <- c("psi1", "psi", "eps", "gamma", "p")
  
  out <- jags(data = bdata, inits = inits, parameters.to.save = params, model.file = "dynoccmod1.txt", #xx
              n.adapt = 1000, n.chains = 3, n.thin = 2, n.iter = 3000, n.burnin = 500,
              parallel = TRUE)
  
  if(!all(unlist(out$Rhat) < 1.1)){
    out <- update(out, n.iter = 3000)
    } # consider adding while loop
  
  # store results:
  res_all <- out
  
  # monitor progress by logging completed iterations:
  cat(paste("completed fitting model for", spec,"\n"), file = logfile, append = T)
  
  save(res_all, 
       file = file.path("data",
                        paste0("res_nullmodel_", gsub(" ", "_", spec), "_25yr_be5_BCR_ss.RData")))
  
  res_metrics$Rhat_fine <- ifelse(all(unlist(out$Rhat) < 1.1), 1, 0)
  res_metrics$DIC <- out$DIC
  res_metrics$pD <- out$pD
  res_metrics$eps_mean <- round(out$summary["eps", "mean"],2)
  res_metrics$eps_sd <- round(out$summary["eps", "sd"],2)
  res_metrics$eps_2.5 <- round(out$summary["eps", "2.5%"],2)
  res_metrics$eps_97.5 <- round(out$summary["eps", "97.5%"],2)
  res_metrics$gamma_mean <- round(out$summary["gamma", "mean"],2)
  res_metrics$gamma_sd <- round(out$summary["gamma", "sd"],2)
  res_metrics$gamma_2.5 <- round(out$summary["gamma", "2.5%"],2)
  res_metrics$gamma_97.5 <- round(out$summary["gamma", "97.5%"],2)
  res_metrics$p_mean <- round(out$summary["p", "mean"],2)
  res_metrics$p_sd <- round(out$summary["p", "sd"],2)
  res_metrics$p_2.5 <- round(out$summary["p", "2.5%"],2)
  res_metrics$p_97.5 <- round(out$summary["p", "97.5%"],2)
  res_metrics$psi1_mean <- round(out$summary["psi1", "mean"],2)
  res_metrics$psi1_sd <- round(out$summary["psi1", "sd"],2)
  res_metrics$psi1_2.5 <- round(out$summary["psi1", "2.5%"],2)
  res_metrics$psi1_97.5 <- round(out$summary["psi1", "97.5%"],2)
  res_metrics$n_routes <- nsites
  
  res_metrics
  
}

stopCluster(cl)

save(res_metrics_all, 
     file = file.path("data",
                      "res_nullmodel_testspecs_25yr_be5_BCR_ss.RData"))