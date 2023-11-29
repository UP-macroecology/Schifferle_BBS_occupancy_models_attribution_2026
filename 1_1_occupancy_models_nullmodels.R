# occupancy null model for 5 test species and 7 versions of routes included

# packages: ----

library(dplyr)
library(tidyr)
library(collapse)
library(jagsUI)
library(doParallel)

# load data: ----

# results of DEBTs\analysis\Schifferle_BBS_explorations_2023\BBS_data_prep.R

## BBS data:
load(file = file.path("data", "BBS_data_merged.RData"))
#load(file = file.path("DEBTs", "analysis", "Schifferle_BBS_occupancy_models_2023", "data", "BBS_data_merged.RData"))
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


# reformat BBS data for occupancy modelling: ----

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
  filter(Year >= 1980 & Year <= 2019) # no data in 2020 due to covid

nyears <- length(unique(route_dt$Year)) 
nsurveys <- 5


# df to create route subsets: --------------------------------------------------------

route_subsets <- c("each year", 
                   "min. every 2nd y.", 
                   "min. every 3rd y.", 
                   "1 y. missing", 
                   "2 y. missing", 
                   "5 y. missing", 
                   "10 y. missing")

# number of years routes are surveyed and max. number of years missed between consecutive surveys:
routes_survey_inf <- route_dt %>% 
  select(RTENO, Year, Surveyed) %>%
  arrange(RTENO, Year) %>%
  # number of missing years between consecutive surveys:
  group_by(RTENO, grp = with(rle(Surveyed), rep(seq_along(lengths), times = lengths))) %>% 
  mutate(counter = seq_along(grp)) %>% 
  mutate(gap = if_else(Surveyed == 1, 0, counter)) %>%
  # maximum number of missing years between consecutive surveys and number of years route was surveyed in total:
  ungroup() %>% 
  group_by(RTENO) %>% 
  summarise(max_gap = max(gap), years_surveyed = sum(Surveyed))


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
#registerDoParallel(cores = 6) # 15: 5 species * 3 chains? # xx
#getDoParWorkers() # check registered number of cores

res_metrics_all <- foreach(spec = testspecs,
                           .combine = rbind,
                           .packages = c("dplyr", "collapse", "jagsUI"),
                           .errorhandling = "remove",
                           .verbose = TRUE) %dopar% {

  # assemble data:
  
  occ_dt_spec <- route_dt %>% 
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
  
  # store model results:
  res_all <- vector(mode = "list", length = 7)
  
  # data frame to store results:
  res_metrics <- data.frame("routes_subset" = rep(NA, 7),
                            "species" = spec,
                            "Rhat_fine" = NA,
                            "n_routes" = NA,
                            "psiobs" = NA,
                            "DIC" = NA, "pD" = NA, 
                            "eps_mean" = NA, "eps_sd" = NA, "eps_2.5" = NA, "eps_97.5" = NA,
                            "gamma_mean" = NA, "gamma_sd" = NA, "gamma_2.5" = NA, "gamma_97.5" = NA,
                            "p_mean" = NA, "p_sd" = NA, "p_2.5" = NA, "p_97.5" = NA,
                            "psi1_mean" = NA, "psi1_sd" = NA, "psi1_2.5" = NA, "psi1_97.5" = NA)
  
  ## iterate over route subsets: ----
  
  for(i in 2:3){ #1:length(route_subsets)){
    
    print(paste(i, "of", length(route_subsets)))
    
    if(i == 1){
      # 1) routes surveyed in each year of the considered time period:
      routes_subset <- routes_survey_inf %>% 
        filter(years_surveyed == nyears) %>% 
        pull(RTENO)
    }
    
    if(i == 2){
      # 2) routes surveyed at least every second year (= max. gap 1 year between consecutive surveys):
      routes_subset <- routes_survey_inf %>% 
        filter(max_gap <= 1) %>% 
        pull(RTENO)
    }
    
    if(i == 3){
      # 3) routes surveyed at least every third year (= max. gap 2 year between consecutive surveys):
      routes_subset <- routes_survey_inf %>% 
        filter(max_gap <= 2) %>% 
        pull(RTENO)
    }
    
    if(i == 4){
      # 4) routes surveyed in minimum all except for 1 year:
      routes_subset <- routes_survey_inf %>% 
        filter(years_surveyed >= nyears-1) %>% 
        pull(RTENO)
    }
    
    if(i == 5){
      # 5) routes surveyed in minimum all except for 2 year:
      routes_subset <- routes_survey_inf %>% 
        filter(years_surveyed >= nyears-2) %>% 
        pull(RTENO)
    }
    
    if(i == 6){
      # 6) routes surveyed in minimum all except for 5 year:
      routes_subset <- routes_survey_inf %>% 
        filter(years_surveyed >= nyears-5) %>% 
        pull(RTENO)
    }
    
    if(i == 7){
      # 7) routes surveyed in minimum all except for 10 year:
      routes_subset <- routes_survey_inf %>% 
        filter(years_surveyed >= nyears-10) %>% 
        pull(RTENO)
    }
    
    # subsample / thin data:xx
    # should I do this after route selection regarding temporal coverage? (yes - first delete routes
    # too rarely surveyed, then delete routes not surveyed if BCR already contains too many routes compared to other BCRs)
    # how many BCRs have more than 30 routes after accounting for temporal coverage:
    
    # think about this tomorrow with app!

    # assemble data:
    
    # use only route subset:
    occ_dt_spec_route_ss <- occ_dt_spec %>% 
      filter(RTENO %in% routes_subset)
    
    ###
    
    # # for testing: Randomly thin out the data set by subsampling 30%: xx
    # occ_dt_spec_route_ss_full <- occ_dt_spec_route_ss                      # Make a copy of full data set
    # prop.data <- 0.3                     # Proportion of data to be used
    # ncase <- nrow(occ_dt_spec_route_ss)                   # 116204
    # set.seed(1)                          # Ensures you get the same subset
    # sel.cases <- sort(sample(1:ncase, ncase * prop.data))
    # occ_dt_spec_route_ss <- occ_dt_spec_route_ss[sel.cases,] # Smaller data set
    #   
    ###
    
    nsites <- length(unique(occ_dt_spec_route_ss$RTENO))
    
    # compute observed occupancy:
    obs_occupancy <- occ_dt_spec_route_ss %>% 
      mutate(presence = rowSums(pick(Count10:Count50))) %>%
      mutate(presence = if_else(presence > 1, 1, 0)) %>% 
      summarise(psiobs = mean(presence, na.rm = TRUE)) %>%  # number of route-year comb. with presence / number of all surveyed route-year comb.
      pull(psiobs)
    
    res_metrics$psiobs[i] <- round(obs_occupancy,2)
    
    # reformat observations to array for model:
    years <- seq(min(occ_dt_spec_route_ss$Year), max(occ_dt_spec_route_ss$Year))
    y_array <- array(NA, dim = c(nsites, nsurveys, nyears))
    for (t in 1:length(years)){
      y_array[1:nsites, 1:nsurveys, t] <- as.matrix(occ_dt_spec_route_ss[which(occ_dt_spec_route_ss$Year == years[t]), 
                                                                         c(paste0("Count", seq(10, 50, 10)))])
    }
    

    # scale date later
    
    ## estimate model: --------------------------------------------------------------

    # bundle data:
    str(bdata <- list(y = y_array, 
                      nsurveys = nsurveys,
                      #site = occ_dt_spec$RTENO, 
                      #year = occ_dt_spec$Year-(min(occ_dt_spec$Year)-1), 
                      #date = occ_dt_spec$doy,
                      nsites = nsites, 
                      nyears = nyears))
    
    # Initial values
    #inits <- function(){list(z = zst)} # z same z as in model definition = occupancy status at site s at year t
    inits <- function(){ list(z = matrix(1, nrow = nsites, ncol = nyears))} # all sites and years have occupancy status 1 xx change
    
    # Parameters monitored
    params <- c("psi1", "psi", "eps", "gamma", "p")
    
    out <- jags(data = bdata, inits = inits, parameters.to.save = params, model.file = "dynoccmod1.txt", #xx
                n.adapt = 1000, n.chains = 3, n.thin = 2, n.iter = 3000, n.burnin = 500,
                parallel = TRUE)
    
    if(!all(unlist(out$Rhat) < 1.1)){
      out <- update(out, n.iter = 3000)
      } # consider adding while loop
    
    # store results:
    res_all[[i]] <- out
    
    res_metrics$routes_subset[i] <- route_subsets[i]
    res_metrics$Rhat_fine[i] <- ifelse(all(unlist(out$Rhat) < 1.1), 1, 0)
    res_metrics$DIC[i] <- out$DIC
    res_metrics$pD[i] <- out$pD
    res_metrics$eps_mean[i] <- round(out$summary["eps", "mean"],2)
    res_metrics$eps_sd[i] <- round(out$summary["eps", "sd"],2)
    res_metrics$eps_2.5[i] <- round(out$summary["eps", "2.5%"],2)
    res_metrics$eps_97.5[i] <- round(out$summary["eps", "97.5%"],2)
    res_metrics$gamma_mean[i] <- round(out$summary["gamma", "mean"],2)
    res_metrics$gamma_sd[i] <- round(out$summary["gamma", "sd"],2)
    res_metrics$gamma_2.5[i] <- round(out$summary["gamma", "2.5%"],2)
    res_metrics$gamma_97.5[i] <- round(out$summary["gamma", "97.5%"],2)
    res_metrics$p_mean[i] <- round(out$summary["p", "mean"],2)
    res_metrics$p_sd[i] <- round(out$summary["p", "sd"],2)
    res_metrics$p_2.5[i] <- round(out$summary["p", "2.5%"],2)
    res_metrics$p_97.5[i] <- round(out$summary["p", "97.5%"],2)
    res_metrics$psi1_mean[i] <- round(out$summary["psi1", "mean"],2)
    res_metrics$psi1_sd[i] <- round(out$summary["psi1", "sd"],2)
    res_metrics$psi1_2.5[i] <- round(out$summary["psi1", "2.5%"],2)
    res_metrics$psi1_97.5[i] <- round(out$summary["psi1", "97.5%"],2)
    res_metrics$n_routes[i] <- nsites
    
  }
  
  save(res_all, res_metrics, 
       file = file.path("data",
                        paste0("res_nullmodel_", gsub(" ", "_", spec), "_route_subsets.RData")))
  
  res_metrics
}

save(res_metrics_all, 
     file = file.path("data",
                      paste0("res_nullmodel_route_subsets_BCR_thinned.RData")))
