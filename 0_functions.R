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