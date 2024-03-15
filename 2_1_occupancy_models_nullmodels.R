# occupancy null model for test species:

# packages: ----

library(dplyr)
library(jagsUI)

sink("output_nullmodel_House_Finch_GoF.txt") # write console output here

# directories: ---- 

cluster_import <- file.path("/import", "ecoc9z", "data-zurell", "schifferle", "BBS_occupancy_models_2023")

# load data: ----

# species for which to fit model:
spec <- "House Finch"

# selected routes and focal years matched to environmental data:
# merged route, year, environment data:
load(file = file.path("data", "route_year_env_data.RData"))

## scale covariates:
route_sel_env_dt_scaled <- route_sel_env_dt_final %>% 
  mutate(across(bio14:primn_3yrs, ~ (scale(.)) %>% as.vector()))

# route-year-species information (only surveyed)
load(file = file.path("data", "BBS_for_occ_spec_records.RData")) # output of 1_0_reformat_BBS_data.R

# selected species:
load(file = file.path("data", "final_species_selection.RData")) # output of 1_2_species_selection.R


# model definition: -----

# simplest dynamic occupancy model: I assume that all parameters are constant across time and space, no covariates
# constant colonisation, extinction, detection

model_file <- "dynoccmod1.txt"

cat(file = model_file," #xxDEBTs/analysis/Schifferle_BBS_occupancy_models_2023/
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
  
  
  # GoF computation part of code: ----

  # Draw a replicate data set under the fitted model
  for (s in 1:nsites){
    for (t in 1:nyears){
      for (j in 1:nsurveys){
        yrep[s,j,t] ~ dbern(p * z[s,t])
      }
    }
  }

  # (a) Computations for the GoF of the open part of the model
  # (based on number of state transitions)

  # Compute observed z matrix for observed and replicated data
  for (s in 1:nsites){
    for (t in 1:nyears){
      zobs[s,t] <- max(y[s,,t])       # For observed data
      zobsrep[s,t] <- max(yrep[s,,t]) # For replicated data
    }

  # Identify extinctions, persistence, colonization and non-colonizations
  for (t in 2:nyears){

    # for observed data:
    ext[s,(t-1)] <- equals(zobs[s,t],0) * equals(zobs[s,t-1],1)
    nonext[s,(t-1)] <- equals(zobs[s,t],1) * equals(zobs[s,t-1],1)
    colo[s,(t-1)] <- equals(zobs[s,t],1) * equals(zobs[s,t-1],0)
    noncolo[s,(t-1)] <- equals(zobs[s,t],0) * equals(zobs[s,t-1],0)

    # for replicated data:
    extrep[s,(t-1)] <- equals(zobsrep[s,t],0) * equals(zobsrep[s,t-1],1)
    nonextrep[s,(t-1)] <- equals(zobsrep[s,t],1) * equals(zobsrep[s,t-1],1)
    colorep[s,(t-1)] <- equals(zobsrep[s,t],1) * equals(zobsrep[s,t-1],0)
    noncolorep[s,(t-1)] <- equals(zobsrep[s,t],0)*equals(zobsrep[s,t-1],0)
    }
  }

  # Tally up number of transitions and put into a matrix for each year
  for(t in 1:(nyears-1)){

    # for observed data:
    tm[1,1,t] <- sum(noncolo[,t]) # transition mat for obs. data
    tm[1,2,t] <- sum(colo[,t])
    tm[2,1,t] <- sum(ext[,t])
    tm[2,2,t] <- sum(nonext[,t])

    # for replicated data:
    tmrep[1,1,t] <- sum(noncolorep[,t]) # transition mat for rep. data
    tmrep[1,2,t] <- sum(colorep[,t])
    tmrep[2,1,t] <- sum(extrep[,t])
    tmrep[2,2,t] <- sum(nonextrep[,t])
  }

  # Compute expected numbers of transitions under the model
  # Probability of each individual transition
  for(s in 1:nsites){
   for(t in 1:(nyears-1)){
    noncolo.exp[s,t] <- (1-psi[t]) * (1-gamma)
    colo.exp[s,t] <- (1-psi[t]) * gamma
    ext.exp[s,t] <- psi[t] * eps
    nonext.exp[s,t] <- psi[t] * (1-eps)
   }
  }

  # Sum up over sites to obtain the expected number of those transitions
  for(t in 1:(nyears-1)){
    Etm[1,1,t] <- sum(noncolo.exp[,t])
    Etm[1,2,t] <- sum(colo.exp[,t])
    Etm[2,1,t] <- sum(ext.exp[,t])
    Etm[2,2,t] <- sum(nonext.exp[,t])
  }

  # Compute Chi-square discrepancy ~~~ see Errata 2021-10-09 (KS: haven't found them)
  for(t in 1:(nyears-1)){

    # for observed data:
    x2Open[1,1,t] <- pow((tm[1,1,t] - Etm[1,1,t]), 2) / (Etm[1,1,t]+e)
    x2Open[1,2,t] <- pow((tm[1,2,t] - Etm[1,2,t]), 2) / (Etm[1,2,t]+e)
    x2Open[2,1,t] <- pow((tm[2,1,t] - Etm[2,1,t]), 2) / (Etm[2,1,t]+e)
    x2Open[2,2,t] <- pow((tm[2,2,t] - Etm[2,2,t]), 2) / (Etm[2,2,t]+e)

    # for replicated data:
    x2repOpen[1,1,t] <- pow((tmrep[1,1,t]-Etm[1,1,t]),2)/(Etm[1,1,t]+e)
    x2repOpen[1,2,t] <- pow((tmrep[1,2,t]-Etm[1,2,t]),2)/(Etm[1,2,t]+e)
    x2repOpen[2,1,t] <- pow((tmrep[2,1,t]-Etm[2,1,t]),2)/(Etm[2,1,t]+e)
    x2repOpen[2,2,t] <- pow((tmrep[2,2,t]-Etm[2,2,t]),2)/(Etm[2,2,t]+e)
  }

  # Add up overall test statistic and compute fit stat ratio (open part)
  Chi2Open <- sum(x2Open[,,])       # Chisq. statistic for observed data
  Chi2repOpen <- sum(x2repOpen[,,]) # Chisq. statistic for replicated data
  Chi2ratioOpen <- Chi2Open / Chi2repOpen


  # (b) Computations for the GoF of the closed part of the model
  # (based on the number of times detected, i.e., detection freqiencies)

  # Compute detection frequencies for observed and replicated data
  for (s in 1:nsites){
    for (t in 1:nyears){

      # Det. frequencies for observed and replicated data
      detfreq[s,t] <- sum(y[s,,t])
      detfreqrep[s,t] <- sum(yrep[s,,t])

      # Expected detection frequencies under the model
      for (j in 1:nsurveys){
        tmp[s,j,t] <- p * z[s, t]
      }
      E[s,t] <- sum(tmp[s,,t])     # Expected number of detections

      # Chi-square and Freeman-Tukey discrepancy measures
      # for actual data set:
      x2Closed[s,t] <- pow((detfreq[s,t] - E[s,t]),2) / (E[s,t]+e)
      ftClosed[s,t] <- pow((sqrt(detfreq[s,t]) - sqrt(E[s,t])),2)

      # for replicated data set:
      x2repClosed[s,t] <- pow((detfreqrep[s,t] - E[s,t]),2) / (E[s,t]+e)
      ftrepClosed[s,t] <- pow((sqrt(detfreqrep[s,t]) - sqrt(E[s,t])),2)
    }
  }

  # Add up Chi-square and FT discrepancies and compute fit stat ratio (closed part)
  Chi2Closed <- sum(x2Closed[,])
  FTClosed <- sum(ftClosed[,])
  Chi2repClosed <- sum(x2repClosed[,])
  FTrepClosed <- sum(ftrepClosed[,])
  Chi2ratioClosed <- Chi2Closed / Chi2repClosed
  FTratioClosed <- FTClosed / FTrepClosed
  
}
")


# assemble data: ----

nyears <- length(unique(route_sel_env_dt_final$Year)) # 25
nsurveys <- 5
nsites <- length(unique(route_sel_env_dt_final$RTENO)) # 476

presences_spec <- bbs_dt_occ %>% 
  select(c(English_Common_Name, RTENO, Year, paste0("Count", seq(10, 50, 10)))) %>% 
  filter(English_Common_Name == spec)

# match to routes-year-env:
occ_dt_spec <- route_sel_env_dt_scaled %>% 
  # add observations:
  collapse::join(presences_spec, on = c("RTENO", "Year"), how = "left") %>% 
  # if route was surveyed but species not observed, replace NA with 0:
  mutate(across(Count10:Count50, ~ 
                  case_when(Surveyed == 1 & is.na(.) ~ 0,
                            .default = .))) %>%
  # convert bird counts to presence / absence:
  mutate(across(Count10:Count50, ~ 
                  case_when(. > 1 ~ 1,
                            .default = .)))

# reformat observations to array for model:
years <- seq(min(occ_dt_spec$Year), max(occ_dt_spec$Year))
y_array <- array(NA, dim = c(nsites, nsurveys, nyears))
for (t in 1:nyears){
  y_array[1:nsites, 1:nsurveys, t] <- as.matrix(occ_dt_spec[which(occ_dt_spec$Year == years[t]), 
                                                            c(paste0("Count", seq(10, 50, 10)))])
}

# data bundle:
bdata <- list(y = y_array, 
              nsites = nsites, 
              nyears = nyears,
              nsurveys = nsurveys,
              e = 0.0001 # to avoid division by zero in GoF part
)
str(bdata)


# fit model: ----

# parameters to monitor:
params <- c("psi1", "psi", "eps", "gamma","p",
            "Chi2Open", "Chi2repOpen",
            "Chi2ratioOpen", "Chi2Closed", "Chi2repClosed", "Chi2ratioClosed",
            "FTClosed", "FTrepClosed", "FTratioClosed", "tm", "tmrep", "Etm"
           )

# initial values:
# for z = occupancy status at site s at year t (see model definition)
# use observed occupancy status per route and year as initial values for z[s, t]:
zst <- apply(bdata$y, c(1, 3), max) # max. over all surveys
inits <- function(){list(z = zst)}

# fit model:
# out <- jags(data = bdata, inits = inits, parameters.to.save = params, 
#             model.file = "dynoccmod_cov_quadr_det_sect_lasso.txt",
#             n.adapt = 1000, n.chains = 3, n.thin = 2, n.iter = 2000, n.burnin = 500,
#             parallel = TRUE,
#             n.cores = 3)

out <- autojags(data = bdata, inits = inits, parameters.to.save = params, 
                model.file = model_file,
                #n.adapt = 1000, # default is NULL, which will result in the function running groups of 100 adaptation iterations (to a max of 10,000) until JAGS reports adaptation is sufficient
                n.chains = 3, 
                n.thin = 5, 
                iter.increment = 3000, 
                n.burnin = 500,
                parallel = TRUE,
                n.cores = 3,
                #save.all.iter = TRUE,
                Rhat.limit = 1.1,
                max.iter = 150000) # Rushing: 50000, Mingjian: 20000

sink(type = "message") # write console output to file

save(out, file = file.path(cluster_import, "data", paste0(gsub(" ", "_", spec), "_nullmodel_GoF.RData")))

