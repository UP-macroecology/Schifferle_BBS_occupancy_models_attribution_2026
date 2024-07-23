# function for spatial thinning of routes (fom Anna):

# function for "fast" spatial thinning"
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
  
  # out <- sf_to_df(out)[,3:4] %>%
  #   rename("lon" = "x", "lat" = "y") # xx
  
  out
}


# Function select07 (from Damaris, comments from me):
select07 <- function(imp, X, threshold=0.7, method="spearman") # imp = importance
{
  
  # choose most important variable of those that have spearman correlation > 0.7
  
  cm <- cor(X, method=method) # correlation matrix
  sort.imp <- colnames(X)[order(imp,decreasing=T)] # sort variables according to importance based on PCA
  
  pairs <- which(abs(cm)>= threshold, arr.ind=TRUE) # identifies correlated variable pairs
  index <- which(pairs[,1]==pairs[,2])           # removes entry on diagonal
  pairs <- pairs[-index,]                        # -"-
  
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

# Function to subsample routes so that there is a maximum of 30 routes in each Bird Conservation Region:

# 2 versions:

## 1) first remove one of the two closest points, 
## then calculate distances again and again remove one of the two closest points...

subsample1 <- function(data, npoints_keep = 30){
  
  # data = sf object
  
  # number of points to remove:
  remove_n <- nrow(data) - npoints_keep
  
  output_data <- data
  
  for(r in 1:remove_n){
    
    #print(r)
    
    # distance matrix:
    dist_mat <- output_data %>% 
      st_distance()
    
    # which pair of points are closest:
    dist_mat[dist_mat==units::as_units(0, "m")] <- NA  # exclude 0 (distance to itself)
    ind <- arrayInd(which.min(dist_mat), dim(dist_mat))
    
    # corresponding RTENO:
    remove_this <- output_data %>% 
      slice(ind[sample(c(1,2), 1)]) %>% # arbitrarily remove one of two closest points
      pull(RTENO)
    
    # remove that route:
    output_data <- output_data %>% 
      filter(RTENO != remove_this)
  }
  return(output_data)
}

## 2) from: https://gist.github.com/danlwarren/271288d5bab45d2da549:

# randomly choose one point, calculate distance to all others, choose the most distant one,
# then find distance of all not chosen points to the chosen points and choose again the most distant point to all chosen points, repeat
# input:
# x, a data frame containing the columns to be used to calculate distances along with whatever other data you need
# cols, a vector of column names or indices to use for calculating distances
# npoints, the number of rarefied points to spit out
# e.g. subsample2(my.data, c("latitude", "longitude"), 200)

# advantage: random drawing -> repeat to get different sets -> what to use? xx

subsample2 <- function(x, cols, npoints){
  
  # create empty vector for output:
  inds <- vector(mode = "numeric")
  
  # create distance matrix:
  this.dist <- as.matrix(dist(x[, cols], upper=TRUE))
  
  # draw first index at random:
  inds <- c(inds, as.integer(runif(1, 1, length(this.dist[,1]))))
  
  # get second index from maximally distant point from first one (necessary for apply ftc)
  inds <- c(inds, which.max(this.dist[,inds]))
  
  while(length(inds) < npoints){
    # for each point, find its distance to the closest point that's already been selected:
    min.dists <- apply(this.dist[,inds], 1, min)
    
    # select the point that is furthest from everything already selected:
    this.ind <- which.max(min.dists)
    
    # get rid of ties, if they exist:
    if(length(this.ind) > 1){
      print("Breaking tie...")
      this.ind <- this.ind[1]
    }
    inds <- c(inds, this.ind)
  }
  return(x[inds,])
}



# species specific presence-absence data as input for occupancy models, based on BBS data:

BBS_pres_abs_spec <- function(species, 
                              BBS_spec_data = bbs_dt_occ, 
                              BBS_route_data = route_sel_dt){
  

  # presences only:
  presences_spec <- BBS_spec_data %>% 
    select(c(English_Common_Name, RTENO, Year, paste0("Count", seq(10, 50, 10)))) %>% 
    filter(English_Common_Name == species)

  if(nrow(presences_spec) == 0) stop("Species name not found.")
  
  # presences-absences:
  occ_dt_spec <- BBS_route_data %>%  #route_sel_env_dt_scaled %>% 
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
    mutate(presence = ifelse(presence >= 1, 1, 0))
  
  return(occ_dt_spec)
}


# returns BBS routes that are within a certain buffer distance around the routes where
# a species was recorded; these are considered to fit dynamic occupancy models:


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
      st_filter(., y = pres_buffer, join = st_within)
    return(routes_within$RTENO_BBS)
  }
  
}

