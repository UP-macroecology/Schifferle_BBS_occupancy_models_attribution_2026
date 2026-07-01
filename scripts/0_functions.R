# Script: 0_functions.R
# Purpose: Defines functions used throughout the analyses
# Inputs:  -
# Outputs: -
# Runs on: HPC (NAS Potsdam)


# function for spatial thinning of routes:

thin <- function(sf, thin_dist = 3000, runs = 10, ncores = 10){
  
  require(sf, quietly = TRUE)
  require(purrr, quietly = TRUE)
  require(furrr, quietly = TRUE)
  
  sample.vec <- function(x, ...) x[sample(length(x), ...)]
  
  sf_buffer <- st_buffer(sf, thin_dist)
  buff_int <- st_intersects(sf, sf_buffer) 
  buff_int <- setNames(buff_int, 1:length(buff_int))
  
  n_int <- map_dbl(buff_int, length)
  
  plan(multisession, workers = ncores)
  
  seeds <- sample.int(n = runs)
  results_runs <- future_map(seeds, function(i){
    
    set.seed(i)
    while (max(n_int) > 1) {
      max_neighbors <- names(which(n_int == max(n_int)))
      
      # remove point with max neighbors
      sampled_id <- sample.vec(max_neighbors, 1)
      
      pluck(buff_int, sampled_id) <- NULL
      buff_int <- map(buff_int, function(x) setdiff(x, as.numeric(sampled_id)))
      n_int <- map_dbl(buff_int, length)
    }
    
    unlist(buff_int) %>% unique()
    
  })
  
  lengths <- map_dbl(results_runs, length)
  
  selected_run <- results_runs[[sample.vec(which(lengths == max(lengths)), 1)]]
  
  out <- sf[selected_run,]
  
  out
}


# function to choose most important variable of those that have spearman correlation > 0.7: 

select07 <- function(imp, X, threshold = 0.7, method="spearman") { # imp = importance
  
  cm <- cor(X, method=method) # correlation matrix
  sort.imp <- colnames(X)[order(imp,decreasing=T)] # sort variables according to importance based on PCA
  
  pairs <- which(abs(cm)>= threshold, arr.ind=TRUE) # identifies correlated variable pairs
  index <- which(pairs[,1]==pairs[,2]) # removes entry on diagonal
  pairs <- pairs[-index,]                       
  
  exclude <- NULL
  
  # iterate over variables, start with most important:
  for (i in seq_len(length(sort.imp))){
    
    if ((sort.imp[i] %in% row.names(pairs))&
        ((sort.imp[i] %in% exclude)==F)) { # if variable is in correlated pairs and not yet excluded
      
      cv <- cm[setdiff(row.names(cm),exclude), sort.imp[i]] # correlations between focal variable and all unexcluded variables
      cv <- cv[setdiff(names(cv), sort.imp[1:i])] # of these only those that has not been iterated over so far
      
      exclude <- c(exclude, names(which((abs(cv)>=threshold)))) # exclude those above threshold -> these are too highly correlated with focal variable and less important
    }
  }
  sort.imp[!(sort.imp %in% unique(exclude)), drop = FALSE]
}


# function to subsample routes so that there is a maximum of 30 routes in each Bird Conservation Region:

# from: https://gist.github.com/danlwarren/271288d5bab45d2da549:

# Input is:
# x, a data frame containing the columns to be used to calculate distances along with whatever other data you need
# cols, a vector of column names or indices to use for calculating distances
# npoints, the number of rarefied points to spit out
#
# e.g., thin.max(my.data, c("latitude", "longitude"), 200)


thin.max <- function(x, cols, npoints){
  #Create empty vector for output
  inds <- vector(mode="numeric")
  
  #Create distance matrix
  this.dist <- as.matrix(dist(x[,cols], upper=TRUE))
  
  #Draw first index at random
  inds <- c(inds, as.integer(runif(1, 1, length(this.dist[,1]))))
  
  #Get second index from maximally distant point from first one
  #Necessary because apply needs at least two columns or it'll barf
  #in the next bit
  inds <- c(inds, which.max(this.dist[,inds]))
  
  while(length(inds) < npoints){
    #For each point, find its distance to the closest point that's already been selected
    min.dists <- apply(this.dist[,inds], 1, min)
    
    #Select the point that is furthest from everything we've already selected
    this.ind <- which.max(min.dists)
    
    #Get rid of ties, if they exist
    if(length(this.ind) > 1){
      print("Breaking tie...")
      this.ind <- this.ind[1]
    }
    inds <- c(inds, this.ind)
  }
  
  return(x[inds,])
}


# function get presence-absence data for a species from BBS data, as input for occupancy models:

BBS_pres_abs_spec <- function(species, 
                              BBS_spec_data = bbs_dt_occ, 
                              BBS_route_data = route_sel_dt){
  

  # presences only:
  presences_spec <- BBS_spec_data %>% 
    select(c(English_Common_Name, RTENO, Year, paste0("Count", seq(10, 50, 10)))) %>% 
    filter(English_Common_Name == species)

  if(nrow(presences_spec) == 0) stop("Species name not found.")
  
  # presences-absences:
  occ_dt_spec <- BBS_route_data %>%  
    # add observations:
    collapse::join(presences_spec, on = c("RTENO", "Year"), how = "left") %>% 
    mutate(English_Common_Name = species) %>% 
    # if route was surveyed but species not observed, replace NA with 0:
    mutate(across(Count10:Count50, ~ 
                    case_when(Surveyed == 1 & is.na(.) ~ 0,
                              .default = .))) %>%
    # convert bird counts to presence / absence:
    mutate(across(Count10:Count50, ~ 
                    case_when(. > 1 ~ 1,
                              .default = .))) %>% 
    # presence on route across all sections:
    mutate(presence = rowSums(across(paste0("Count", seq(10, 50, 10))))) %>%
    mutate(presence = ifelse(presence >= 1, 1, 0)) %>% 
    #mutate(presence = factor(presence)) %>% 
    arrange(RTENO)
  
  # for the following species there are reported presences on routes far outside the breeding range
  # these are considered absences (not used for models): # see 1_3_dataprep_BBS_outlier_check_selected_species.R.R
  if(species %in% c("Yellow-throated Warbler", "Red Crossbill")){
    if(species == "Yellow-throated Warbler"){
      occ_dt_spec$presence[which(occ_dt_spec$RTENO == 84069052)] <- 0
    }
    if(species == "Red Crossbill"){
      occ_dt_spec$presence[which(occ_dt_spec$RTENO == 84007015)] <- 0
    }
  }
  return(occ_dt_spec)
}


# function to get BBS routes that are within a certain buffer around the routes where
# a species was recorded, as input for occupancy models:

training_routes <- function(species, buffer_km, BBS_spec_data = bbs_dt_occ,
                            routes_sf = routes_sel_sf,
                            output = c("RTENOs", "buffer")){
  
  # presences only:
  presences_spec <- BBS_spec_data %>% 
    select(c(English_Common_Name, RTENO, Year, paste0("Count", seq(10, 50, 10)))) %>% 
    filter(English_Common_Name == species)
  
  if(nrow(presences_spec) == 0) stop("Species name not found.")
  
  # presence - absences:
  occ_dt_spec <- BBS_pres_abs_spec(species = species)
  
  # routes with presences (sf):
  occ_spec_sf <- routes_sf %>%
    left_join(occ_dt_spec, by = c("RTENO_BBS" = "RTENO")) %>%
    # summarise presence on route across all years:
    group_by(RTENO_BBS) %>%
    summarise(presence_summarised = max(presence, na.rm=TRUE)) %>%
    mutate(presence_summarised = factor(presence_summarised, levels = c(1,0)))
  
  # buffer presences:
  pres_buffer <- occ_spec_sf %>% 
    filter(presence_summarised == 1) %>%
    st_buffer(dist = buffer_km*1000) %>% # 750000
    st_union
  
  if(output == "buffer"){
    return(pres_buffer)
  } else {
    # routes within buffer:
    routes_within <- occ_spec_sf %>% 
      st_filter(., y = pres_buffer, join = st_within) %>% 
      pull(RTENO_BBS) %>% 
      sort()

    return(routes_within)
  }
  
}


# modified function of R package flocker to get Z matrix of occupancy probabilities from flocker models faster:

# adapted version to use data frame as new data:
 get_Z_mod <- function(flocker_fit, 
                       draw_ids = NULL, 
                       new_data = NULL){ # as data frame
   
   new_data <- as.data.frame(new_data)
   
   if (is.null(draw_ids)) {
     n_iter <- brms::ndraws(flocker_fit)
   }
   else {
     n_iter <- length(draw_ids)
   }
   
   # initial occ. prob., colonisation and extinction prob.:
   
   lps2 <- fitted_flocker(flocker_fit, 
                          components = c("occ", "colo", "ex"), 
                          draw_ids = draw_ids, 
                          new_data = new_data)
     
   # reformat from long df to array:
   
   nsites <- length(unique(new_data$cellID))
   nyears <- length(unique(new_data$year))
   
   init <- lps2$linpred_occ[1:nsites, ]
   
   colo <- array(NA, dim = c(nsites, nyears, n_iter))
   for(t in 1:nyears){
     colo[,t,] <- lps2$linpred_col[(t-1) * nsites + (1:nsites),]
   }

   ex <- array(NA, dim = c(nsites, nyears, n_iter))
   for(t in 1:nyears){
     ex[,t,] <- lps2$linpred_ex[(t-1) * nsites + (1:nsites),]
   }

   # calculate occ. probs for each site and year:
   
   out <- array(NA, dim = dim(colo))

   for (i in 1:nsites) {
     
     print(i)
     
     for (j in 1:n_iter) {

       length_out <- length(colo[i, , j])

       if (!identical(colo[i, , j], NA)) { # if there is at least one NA in colo 
         colo[i, , j] <- colo[i, , j][1:max(which(!is.na(colo[i, , j])))] # removes NAs at the end only?
         ex[i, , j] <- ex[i, , j][1:max(which(!is.na(colo[i, , j])))]
       }
       
       time_series <- rep(NA, length_out)
       
       time_series[1] <- init[i, j]
         
       if (length(colo[i, , j]) > 1) {
         
         for (y in 2:length_out) {
           time_series[y] <- (1 - time_series[y - 1]) * colo[i, y, j] + (time_series[y - 1]) * (1 - ex[i, y, j])
         }
       }
       
       out[i, , j] <- time_series
     }
   }
  out # Z matrix
  
 }
 
 
 # get sum of routes across the conterminous USA with observations of a species for each year between 1995 and 2019:
 
 obs_time_series <- function(spec){
   
   rel_routes <- training_routes(species = spec, buffer_km = 750, output = "RTENOs")
   occ_dt_spec <- BBS_pres_abs_spec(species = spec) %>%
     filter(RTENO %in% rel_routes)
   
   # route-level presence:
   # sum all routes for each year (temporal trend)
   obs_temp_trend <- occ_dt_spec %>%
     group_by(Year) %>%
     summarise(pres_sum = sum(presence, na.rm = TRUE))
   
   return(obs_temp_trend) 
 }
