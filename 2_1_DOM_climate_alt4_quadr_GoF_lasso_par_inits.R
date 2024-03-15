# Fit occupancy models including goodness-of-fit measure as following Kéry and Royle 2021 AHM:

# quadratic effects of bioclim covariates: bio1, bio2, bio3, bio12, bio15,  lasso prior, p intercept only

# packages: ----

library(dplyr)
library(jagsUI)
library(doParallel)

print(tempdir())

# directories: ----
dir <- file.path("/import", "ecoc9z", "data-zurell", "schifferle", "BBS_occupancy_models_2023")
#dir <- getwd()

# load data: ----

# selected routes and focal years matched to environmental data:
# merged route, year, environment data:
load(file = file.path("data", "route_year_env_data.RData"))

# scale covariates:
route_sel_env_dt_scaled <- route_sel_env_dt_final %>% 
  mutate(across(bio2:pr_winter_3yrs, ~ (scale(.)) %>% as.vector()))

# route-year-species information (only surveyed)
load(file = file.path("data", "BBS_for_occ_spec_records.RData")) # output of 1_0_reformat_BBS_data.R

# selected species:
load(file = file.path("data", "final_species_selection.RData")) # output of 1_2_species_selection.R


# model definition: -----

model_file <- "dynoccmod_cl_alt4_quadr_lasso_GoF.txt"

cat(file = model_file,"
model {
  
## Declare priors ------------------------------------------------------------------------------
  
  # Initial occupancy Priors
	
	# Simple uniform priors for intercepts
	alpha_psi <- logit(psi_intercept)
	psi_intercept ~ dunif(0.0001, 0.9999)
	
  for (i in 1:16) { 
	  beta_psi[i] ~ ddexp(0, tau_psi)
	}
	tau_psi <- 1/lambda_psi 
	lambda_psi ~ dt(0,1,1)T(0,)

  # Colonisation and Extinction Priors
	
	# Sets a uniform prior for intercepts on epsilon
	alpha_eps <- logit(eps_intercept)
 	eps_intercept ~ dunif(0.0001,0.9999)

	# Sets a uniform prior for intercepts on gamma
 	alpha_gamma <- logit(gamma_intercept)
  gamma_intercept ~ dunif(0.0001,0.9999)
    	
	# For each covariate on Extinction, set an Lasso prior for the coefficients only
	for (i in 1:16) {
		beta_eps[i] ~ ddexp(0, tau_eps) 
	}
  tau_eps <- 1/lambda_eps 
	lambda_eps ~ dt(0,1,1)T(0,)

	# For each covariate on Colonisation, set an Lasso prior for the coefficients only
	for (i in 1:16) { 
		beta_gamma[i] ~ ddexp(0, tau_gamma)
	}
  tau_gamma <- 1/lambda_gamma 
	lambda_gamma ~ dt(0,1,1)T(0,)


  # Detection Priors
	p ~ dunif(0.0001,0.9999)

  # Likelihood: ---------
  
  ## Ecological submodel: ----
  
  # Iterate over sites:
	for (s in 1:nsites){
		
		# Sets out the formula for Psi (initial occupancy) at the site level
		logit(psi1[s]) <- alpha_psi + beta_psi[1]*bio1_3yrs[s] + beta_psi[2]*bio2_3yrs[s] + 
		                  beta_psi[3]*bio3_3yrs[s] + beta_psi[4]*pr_spring_3yrs[s] + 
		                  beta_psi[5]*pr_summer_3yrs[s] + beta_psi[6]*pr_autumn_3yrs[s] + 
		                  beta_psi[7]*pr_winter_3yrs[s] + beta_psi[8]*bio15_3yrs[s] +
		                  beta_psi[9]*pow(bio1_3yrs[s],2) + beta_psi[10]*pow(bio2_3yrs[s],2) +
		                  beta_psi[11]*pow(bio3_3yrs[s],2) + beta_psi[12]*pow(pr_spring_3yrs[s],2) +
		                  beta_psi[13]*pow(pr_summer_3yrs[s],2) + beta_psi[14]*pow(pr_autumn_3yrs[s],2) +
		                  beta_psi[15]*pow(pr_winter_3yrs[s],2) + beta_psi[16]*pow(bio15_3yrs[s],2)

		# Based on Psi, determine Year1 occupancy
		z[s,1] ~ dbern(psi1[s]) 
	
    # State transitions
		# Iterated over years:
		for (t in 2:nyears){
      
      # Set formula for Epsilon (Extinction)
			logit(eps[s, t-1]) <- alpha_eps + beta_eps[1]*bio1[s, t] + beta_eps[2]*bio2[s, t] + 
		                        beta_eps[3]*bio3[s, t] + beta_eps[4]*pr_spring[s, t] + 
		                        beta_eps[5]*pr_summer[s, t] + beta_eps[6]*pr_autumn[s, t] + 
		                        beta_eps[7]*pr_winter[s, t] + beta_eps[8]*bio15[s, t] +
		                        beta_eps[9]*pow(bio1[s, t],2) + beta_eps[10]*pow(bio2[s, t],2) +
		                        beta_eps[11]*pow(bio3[s, t],2) + beta_eps[12]*pow(pr_spring[s, t],2) +
		                        beta_eps[13]*pow(pr_summer[s, t],2) + beta_eps[14]*pow(pr_autumn[s, t],2) +
		                        beta_eps[15]*pow(pr_winter[s, t],2) + beta_eps[16]*pow(bio15[s, t],2)
			
			# Set formula for Gamma (Colonisation)
			logit(gamma[s, t-1]) <- alpha_gamma + beta_gamma[1]*bio1[s, t] + beta_gamma[2]*bio2[s, t] + 
		                          beta_gamma[3]*bio3[s, t] + beta_gamma[4]*pr_spring[s, t] + 
		                          beta_gamma[5]*pr_summer[s, t] + beta_gamma[6]*pr_autumn[s, t] + 
		                          beta_gamma[7]*pr_winter[s, t] + beta_gamma[8]*bio15[s, t] +
		                          beta_gamma[9]*pow(bio1[s, t],2) + beta_gamma[10]*pow(bio2[s, t],2) +
		                          beta_gamma[11]*pow(bio3[s, t],2) + beta_gamma[12]*pow(pr_spring[s, t],2) +
		                          beta_gamma[13]*pow(pr_summer[s, t],2) + beta_gamma[14]*pow(pr_autumn[s, t],2) +
		                          beta_gamma[15]*pow(pr_winter[s, t],2) + beta_gamma[16]*pow(bio15[s, t],2)
			
			# Based on these and prior season state, determine occupancy state
			z[s,t] ~ dbern(z[s,t-1]*(1-eps[s,t-1]) + (1-z[s,t-1])*gamma[s,t-1])
		}
	}

  ## Observation model: ----
  
	# Iterate over Year-Site-Survey combinations
	for (s in 1:nsites){
		for (j in 1:nsurveys){
			for (t in 1:nyears){
				# Declare detection formula
				y[s,j,t] ~ dbern(z[s,t]*p)
			}
		}
	}

  ## Derived parameters: ----
  
  # occupancy at each route and in each year:
  for(s in 1:nsites){
    
    psi[s,1] <- psi1[s]    

    for (t in 2:nyears){
      psi[s, t] <- psi[s, t-1]*(1-eps[s, t-1]) + (1-psi[s, t-1])*gamma[s, t-1]
  }
  }
  
  
  # GoF computation part of code: ----
  
  # Draw a replicate data set under the fitted model
  for (s in 1:nsites){
    for (t in 1:nyears){
      for (j in 1:nsurveys){
        yrep[s,j,t] ~ dbern(z[s,t] * p)
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
    noncolo.exp[s,t] <- (1-psi[s,t]) * (1-gamma[s,t])
    colo.exp[s,t] <- (1-psi[s,t]) * gamma[s,t]
    ext.exp[s,t] <- psi[s,t] * eps[s,t]
    nonext.exp[s,t] <- psi[s,t] * (1-eps[s,t])
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
        tmp[s,j,t] <- z[s, t] * p
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

# check whether species has already run: xx

testspecs <- c("American Goldfinch",
               "Black-capped Chickadee", 
               "Eastern Phoebe",
               "Yellow Warbler",
               "Eastern Kingbird",
               "American Robin",
               "Black-billed Cuckoo",
               "House Finch",
               "Warbling Vireo")

# register cores for parallel computation:
ncores <- 27 # n species * 3 chains? 
cl <- makeCluster(ncores, setup_timeout = 0.5)
registerDoParallel(cl)


foreach(spec = testspecs,
        .packages = c("dplyr", "collapse", "jagsUI"),
        .errorhandling = "remove",
        .verbose = TRUE) %dopar% {
          
          
          # assemble data: ----
          
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
                        bio1_3yrs = occ_dt_spec$bio1_3yrs[seq(1, nrow(occ_dt_spec), nyears)], # one value per site
                        bio2_3yrs = occ_dt_spec$bio2_3yrs[seq(1, nrow(occ_dt_spec), nyears)],
                        bio3_3yrs = occ_dt_spec$bio3_3yrs[seq(1, nrow(occ_dt_spec), nyears)],
                        pr_spring_3yrs = occ_dt_spec$pr_spring_3yrs[seq(1, nrow(occ_dt_spec), nyears)],
                        pr_summer_3yrs = occ_dt_spec$pr_summer_3yrs[seq(1, nrow(occ_dt_spec), nyears)],
                        pr_autumn_3yrs = occ_dt_spec$pr_autumn_3yrs[seq(1, nrow(occ_dt_spec), nyears)],
                        pr_winter_3yrs = occ_dt_spec$pr_winter_3yrs[seq(1, nrow(occ_dt_spec), nyears)],
                        bio15_3yrs = occ_dt_spec$bio15_3yrs[seq(1, nrow(occ_dt_spec), nyears)],
                        bio1 = matrix(occ_dt_spec$bio1, nrow = nsites, ncol = nyears, byrow = TRUE), # matrix [site, year]
                        bio2 = matrix(occ_dt_spec$bio2, nrow = nsites, ncol = nyears, byrow = TRUE),
                        bio3 = matrix(occ_dt_spec$bio3, nrow = nsites, ncol = nyears, byrow = TRUE),
                        pr_spring = matrix(occ_dt_spec$pr_spring, nrow = nsites, ncol = nyears, byrow = TRUE),
                        pr_summer = matrix(occ_dt_spec$pr_summer, nrow = nsites, ncol = nyears, byrow = TRUE),
                        pr_autumn = matrix(occ_dt_spec$pr_autumn, nrow = nsites, ncol = nyears, byrow = TRUE),
                        pr_winter = matrix(occ_dt_spec$pr_winter, nrow = nsites, ncol = nyears, byrow = TRUE),
                        bio15 = matrix(occ_dt_spec$bio15, nrow = nsites, ncol = nyears, byrow = TRUE),
                        nsites = nsites, 
                        nyears = nyears,
                        nsurveys = nsurveys,
                        e = 0.0001 # to avoid division by zero in GoF part
          )
          str(bdata)
          
          # fit model: ----
          
          # parameters to monitor:
          params <- c("psi1", "psi", "alpha_psi", "beta_psi",
                      "eps", "alpha_eps", "beta_eps", 
                      "gamma", "alpha_gamma", "beta_gamma",
                      "p",
                      "Chi2Open", "Chi2repOpen",
                      "Chi2ratioOpen", "Chi2Closed", "Chi2repClosed", "Chi2ratioClosed",
                      "FTClosed", "FTrepClosed", "FTratioClosed", "tm", "tmrep", "Etm")
          
          # initial values:
          # for z = occupancy status at site s at year t (see model definition)
          # use observed occupancy status per route and year as initial values for z[s, t]:
          zst <- apply(bdata$y, c(1, 3), max) # max. over all surveys
          
          # model parameters: psi_intercept, eps_intercept, gamma_intercept, beta_psi[1:10], 
          # lambda_psi, beta_eps[1:10], lambda_eps, beta_gamma[1:10], lambda_gamma, p
          
          # initial values for psi_intercept, eps_intercept, gamma_intercept, p?
          
          # raw mean extinction and colonisation probability across sites for each year:
          raw_ext <- numeric(length=nyears - 1)
          raw_col <- numeric(length=nyears - 1)
          
          for(t in 1:nyears - 1){
            transitions <- zst[, t+1] - zst[, t]
            raw_ext[t] <- length(which(transitions == -1))/nsites
            raw_col[t] <- length(which(transitions == 1))/nsites
          }
          
          # raw detection probability:
          # mean per site and year:
          raw_p <- numeric(length=nyears)
          for(t in 1:nyears){
            raw_p[t] <- mean(rowSums(bdata$y[,,t], na.rm = TRUE)/nsurveys) # year t, mean across sites
          }
          
          inits <- function(){  
            list(z = zst,
                 psi_intercept = sum(zst[,1])/nsites, # raw psi first year across sites
                 eps_intercept = runif(n = 1, min = max(0, mean(raw_ext)-2*sd(raw_ext)), 
                                       max = min(1, mean(raw_ext)+2*sd(raw_ext))),
                 gamma_intercept = runif(n = 1,min = max(0, mean(raw_col)-2*sd(raw_col)), 
                                         max = min(1, mean(raw_col)+2*sd(raw_col))),
                 p_intercept = runif(n = 1, 
                                     min = max(0, mean(raw_p)-2*sd(raw_p)), 
                                     max = min(1, mean(raw_p)+2*sd(raw_p)))
            )  
          }
          
          
          # fit model:
          # out <- jags(data = bdata, inits = inits, parameters.to.save = params, 
          #             model.file = model_file,
          #             n.adapt = 1000, n.chains = 3, n.thin = 2, n.iter = 2000, n.burnin = 500,
          #             parallel = TRUE,
          #             n.cores = 3)
          
          sink(paste0("out_cl_alt4_quadr_lasso_", gsub(" ", "_", spec),"_GoF_inits.txt")) # write console output here
          sink(type = "message")
          
          print(spec)
          
          print(inits()[2:5])
          
          start.time <- Sys.time()
          
          out <- autojags(data = bdata, inits = inits, parameters.to.save = params, 
                          model.file = model_file,
                          #n.adapt = 1000, # default is NULL, which will result in the function running groups of 100 adaptation iterations (to a max of 10,000) until JAGS reports adaptation is sufficient
                          n.chains = 3, 
                          n.thin = 2, 
                          iter.increment = 1000, 
                          n.burnin = 500,
                          parallel = TRUE,
                          n.cores = 3,
                          Rhat.limit = 1.1,
                          max.iter = 40000) # Rushing: 50000, Mingjian: 20000
          
          save(out, file = file.path(dir, "data", paste0(gsub(" ", "_", spec), "_cl_alt4_quadr_lasso_GoF_inits.RData")))
          
          end.time <- Sys.time()
          time.taken <- round(end.time - start.time,2)
          print(time.taken)
          
          sink(file = NULL)
        }

stopCluster(cl)

rm(list=ls())
gc()
