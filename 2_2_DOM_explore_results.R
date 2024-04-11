# explore fitted dynamic occupancy models:

# packages: ----

library(jagsUI)
library(dplyr)
library(ggplot2)
library(stringr)
library(flocker)
library(brms)

# load fitted model(s):

load(file.path("results", "American_Goldfinch_cl_alt4_lu_quadr_p_sect_lasso_GoF_inits_less_pars.RData"))
out

# general explorations: ------

## JAGS: ----

print(out)
summary(out)
out$summary
str(out)
str(out$samples) # samples per chain; original output object from the rjags package, as class mcmc.list
out$sims.list$p # actual draws from posterior for p
densityplot(out, parameters = "beta_psi")
jagshelper::plotdens(out, p="beta_psi[2]")
traceplot(x = out, parameters = "psi1")
traceplot(out, parameters = "gamma", Rhat_min =  1.1)
whiskerplot(out, parameters = "beta_psi")

#library(jagshelper)
jagshelper::check_Rhat(jags_out) # proportion of Rhats below a threshold of 1.1
jagshelper::traceworstRhat(jags_out, parmfrow=c(3,2))  # trace plots for least-converged nodes
out_df <- jagshelper::jags_df(jags_out)
str(out_df)

# possibly relevant covariates:
which(!out$overlap0$beta_psi)
which(!out$overlap0$beta_gamma)
which(!out$overlap0$beta_eps)

# function to return covariates with mixing issues:

pars_conv_issue <- function(jags_output){
  
  pars_no_conv <- names(which(jags_output$summary[,"Rhat"] > 1.1))
  # parameters with mixing issues:
  par_names_no_conv <- unique(str_extract(pars_no_conv, ".*(?=\\[)"))
  
  return(par_names_no_conv)
}


# routes with JAGS mixing issues:

# function to return sites (RTENO number) with mixing issues (in any parameter):
sites_conv_issue <- function(jags_output){
  
  # parameters with mixing issues:
  pars_no_conv <- names(which(jags_output$summary[,"Rhat"] > 1.1))

  # selected routes and years:
  load(file = file.path("data", "route_year_env_data.RData"))
  nsites = 476
  nyears = 25
  route_nrs <- matrix(route_sel_env_dt_final$RTENO, nrow = nsites, ncol = nyears, byrow = TRUE)[,1]
  # site numbers with mixing issues:
  routes_no_conv <- unique(as.numeric(str_extract(pars_no_conv, "((?<=\\[).*(?=,))|((?<=\\[).*(?=\\]))"))) # (?<=\\[) looks for what's behind [, (?=\\])) looks for what's before ]
  # corresponding RTENOs:
  RTENOs_no_conv <- route_nrs[routes_no_conv]
  return(RTENOs_no_conv)
}

RTENOs_no_conv <- sites_conv_issue(jags_output)
  
# load routes:
routes_sel_sf <- sf::st_read(file.path("data", "route_selection_1991_2015_surv_beg_end_max_5y_miss_v2_spat_thin_100km_max_30_r_per_BCR_centroids.shp")) # output of 1_1_route_selection.R
# plots:
bio_sf <- routes_sel_sf %>% 
  left_join(route_sel_env_dt_final[, c("RTENO", "bio15_3yrs")], by = c("RTENO_BBS" = "RTENO")) %>% 
  distinct()
plot(bio_sf["bio15_3yrs"],reset = FALSE)
routes_sel_sf %>% 
  filter(RTENO_BBS %in% RTENOs_no_conv) %>% 
  plot(add = TRUE, col = "black", reset = FALSE, pch = 4)


## flocker: ----

load(file.path("results", "American_Goldfinch_cl_lu_p_Stan_prior.RData"))
flocker_out <- multi_colex_cl_lu_p_prior # 33 div.

flocker_out$fit
flocker_out$prior
str(flocker_out)
#plot(flocker_out)
plot(flocker_out, variable = "b_occ_bio1_3yrs")
bayesplot::mcmc_combo(flocker_out, pars = "b_occ_bio1_3yrs")
posterior <- as.array(flocker_out)
bayesplot::mcmc_intervals(posterior, pars = c("b_occ_bio1_3yrs"))
bayesplot::mcmc_dens(posterior, pars = c("b_occ_bio1_3yrs"))


# JAGS goodness-of-fit following Kéry and Royle AHM 2021: ----

load(file.path("results", "American_Goldfinch_cl_alt4_lu_quadr_p_sect_uninf_GoF_inits_less_pars_butGoF.RData"))
jags_out <- out

# Plots of expected versus observed value of fit stats:

# Open part pf model = occupancy submodel:
pl <- range(c(jags_out$sims.list$Chi2Open, jags_out$sims.list$Chi2repOpen))
plot(jags_out$sims.list$Chi2Open, jags_out$sims.list$Chi2repOpen,
     xlab = "Chi2 observed data", ylab = "Chi2 expected data",
     main = "Open part of model ", xlim = pl, ylim = pl, frame.plot = FALSE)
abline(0, 1, lwd = 2)
text(200, 500, paste('Bpv = ', round(mean(jags_out$sims.list$Chi2repOpen >
                                            jags_out$sims.list$Chi2Open), 2)), cex = 1) # Bayesian p-value 0 proportion of points above the line

# Closed part of model = detection submodel:
pl <- range(c(jags_out$sims.list$Chi2Closed, jags_out$sims.list$Chi2repClosed))
plot(jags_out$sims.list$Chi2Closed, jags_out$sims.list$Chi2repClosed,
     xlab = "Chi2 observed data", ylab = "Chi2 expected data",
     main = "Closed part of model (Chi-squared)", xlim = pl, ylim = pl,
     frame.plot = FALSE)
abline(0, 1, lwd = 2)
text(3500, 5000, paste('Bpv = ', round(mean(jags_out$sims.list$Chi2repClosed >
                                              jags_out$sims.list$Chi2Closed), 2)), cex = 1)

# GoF for open part of the model looks better



# compare JAGS and brms/flocker model fits: ----

# load fitted models:
load(file.path("results", "American_Goldfinch_cl_alt4_lu_quadr_p_sect_lasso_GoF_inits_less_pars.RData"))
jags_out <- out
load(file.path("results", "American_Goldfinch_cl_lu_p_Stan_prior.RData"))
flocker_out <- multi_colex_cl_lu_p_prior # 33 div.


# df matching JAGS parameter names to flocker parameter names:
jags_pars <- rownames(jags_out$summary)[which(grepl("(alpha)|(beta)|(p\\[)", rownames(jags_out$summary)))]
stan_pars <- rownames(fixef(flocker_out, summary = TRUE))
par_names_jags_stan <- data.frame("JAGS_pars" = jags_pars,
                                  "Stan_pars" = stan_pars[c(2, 9:25, 29, 26:28, 30, 34, 31:33, # psi
                                                            4, 61:77, 81, 78:80, 82, 86, 83:85, # eps
                                                            3, 35:51, 55, 52:54, 56, 60, 57:59, # gamma
                                                            1, 5:8)]) # p1-5 stan intercept = p1??

# function to compare output for full model:

compare_Jags_Stan_fm <- function(jags_output, stan_output){
  
  stan_estimates <- as.data.frame(fixef(stan_output, summary = TRUE, robust = FALSE, 
                                        probs = c(0.025, 0.975)))
  stan_estimates$par <- row.names(stan_estimates)
  
  jags_estimates <- as.data.frame(jags_output$summary[which(grepl("(alpha)|(beta)|(p)", rownames(jags_output$summary))),
                                                      c("mean", "sd", "2.5%", "97.5%")])
  jags_estimates$par <- row.names(jags_estimates)
  
  comparison <- par_names_jags_stan %>% 
    left_join(jags_estimates, by = c(JAGS_pars = "par")) %>% 
    left_join(stan_estimates, by = c(Stan_pars = "par")) %>% 
    rename("JAGS_mean" = mean, "JAGS_sd" = sd, "JAGS_Q2.5" = "2.5%", "JAGS_Q97.5" = "97.5%",
           "Stan_mean" = Estimate , "Stan_error" = Est.Error, "Stan_Q2.5" = "Q2.5", "Stan_Q97.5" = "Q97.5")
  
  comparison <- comparison[ , c(1, 2, 3, 7, 4, 8,5,9, 6, 10)]
  
  comparison$Stan_sign <- ifelse(comparison$Stan_Q2.5 < 0 & comparison$Stan_Q97.5 > 0, 0, 1)
  comparison$Jags_sign <- ifelse(comparison$JAGS_Q2.5 < 0 & comparison$JAGS_Q97.5 > 0, 0, 1)
  
  comparison$sign_match <- ifelse((comparison$Stan_sign == 1 & comparison$Jags_sign == 1) |
                                    (comparison$Stan_sign == 0 & comparison$Jags_sign == 0), 1, 0)
  
  return(comparison)
  
}

AG_fm_comp <- compare_Jags_Stan_fm(jags_output = jags_out, stan_output = flocker_out)


# compare posterior distributions of JAGS and flocker: ----

# function to compare JAGS and Stan posterior density for single parameters:
compare_post_dens <- function(jags_output, stan_output, jags_par, stan_par) {
  
  jags_draws <- unlist(jags_output$samples[ , jags_par])
  brms_draws <- unlist(lapply(stan_output$fit@sim$samples, "[", stan_par))
  
  colors <- c("JAGS" = "blue", "flocker" = "red")
  
  ggplot() + 
    geom_density(aes(x = brms_draws, colour = "flocker"))  +
    geom_density(aes(x = jags_draws, colour = "JAGS")) +
    theme_bw() +
    labs(x = "estimate",
         y = "density",
         color = "Legend") +
    ggtitle(stan_par) +
    scale_color_manual(values = colors) +
    theme(legend.position="bottom") +
    geom_vline(aes(xintercept = mean(brms_draws), colour = "flocker"), linetype = "dashed") +
    geom_vline(aes(xintercept = mean(jags_draws), colour = "JAGS"), linetype = "dashed") +
    geom_vline(xintercept = 0)
}

compare_post_dens(jags_output = jags_out, stan_output = flocker_out,
                  jags_par = "beta_psi[2]", stan_par = "b_occ_bio2_3yrs")


# compare posterior distributions using different priors: ----

## compare JAGS lasso, JAGS normal and flocker horseshoe: ----

load(file.path("results", "American_Goldfinch_cl_alt4_lu_quadr_p_sect_uninf_GoF_inits_less_pars_butGoF.RData"))
norm_prior <- out
load(file.path("results", "American_Goldfinch_cl_alt4_lu_quadr_p_sect_lasso_GoF_inits_less_pars_butGoF.RData"))
lasso_prior <- out
load(file.path("results", "American_Goldfinch_cl_lu_p_Stan_prior_1000its.RData"))
hs_prior <- multi_colex_cl_lu_p_prior

for(p in 1:nrow(par_names_jags_stan)){
  
  print(p)
  
  lasso_draws <- unlist(lasso_prior$samples[ , par_names_jags_stan$JAGS_pars[p]]) # 3000
  norm_draws <- unlist(norm_prior$samples[ , par_names_jags_stan$JAGS_pars[p]]) # 3000
  hs_draws <- unlist(lapply(hs_prior$fit@sim$samples, "[", paste0("b_", par_names_jags_stan$Stan_pars[p]))) # 4000
  
  
  p_df <- tibble(par = par_names_jags_stan$Stan_pars[p], 
                 lasso_draws, 
                 norm_draws, "hs_draws" = sample(hs_draws, size = 3000, replace = FALSE)) #%>% # use only 3000 of the 4000 samples
  #tidyr::pivot_longer(cols = 2:3, values_to = "draw", names_to = "prior") %>% 
  #mutate(prior = ifelse(prior == "lasso_draws", "lasso", "norm"))
  
  if(p == 1){
    all_p_df <- p_df
  } else{
    all_p_df <- rbind(all_p_df, p_df)
  }
}

colors <- c("norm(0, 0.5)" = "blue", "lasso" = "red", "flocker_hs" = "#4daf4a")
ggplot(data = all_p_df) + 
  geom_density(aes(x = lasso_draws, colour = "lasso"))  +
  geom_density(aes(x = norm_draws, colour = "norm(0, 0.5)")) +
  geom_density(aes(x = hs_draws, colour = "flocker_hs")) +
  ggforce::facet_wrap_paginate(facets = ~par, scales = "free", nrow = 3, ncol = 4, page = 8) +
  geom_vline(xintercept = 0) +
  theme_bw() +
  labs(x = "estimate",
       y = "density",
       color = "Legend") +
  scale_color_manual(values = colors) +
  theme(legend.position="bottom")


## compare flocker normal prior vs. horsehoe prior: ----

load(file.path("results", "American_Goldfinch_cl_lu_p_Stan_no_hs.RData"))
norm_prior <- multi_colex_cl_lu_p_prior
load(file.path("results", "American_Goldfinch_cl_lu_p_Stan_prior_1000its.RData"))
hs_prior <- multi_colex_cl_lu_p_prior


stan_pars <- rownames(fixef(norm_prior, summary = TRUE))

for(p in 1:length(stan_pars)){
  
  print(p)
  
  fl_norm_draws <- unlist(lapply(norm_prior$fit@sim$samples, "[", paste0("b_", stan_pars[p]))) # 4000
  fl_hs_draws <- unlist(lapply(hs_prior$fit@sim$samples, "[", paste0("b_", stan_pars[p]))) # 4000
  
  p_df <- tibble(par = stan_pars[p], 
                 fl_norm_draws, 
                 fl_hs_draws) 
  
  if(p == 1){
    all_stan_p_df <- p_df
  } else{
    all_stan_p_df <- rbind(all_stan_p_df, p_df)
  }
}

colors <- c("norm(0, 2)" = "blue", "flocker_hs" = "#4daf4a")
ggplot(data = all_stan_p_df) + 
  geom_density(aes(x = fl_norm_draws, colour = "norm(0, 2)")) +
  geom_density(aes(x = fl_hs_draws, colour = "flocker_hs")) +
  ggforce::facet_wrap_paginate(facets = ~par, scales = "free", nrow = 3, ncol = 4, page = 1) +
  geom_vline(xintercept = 0) +
  theme_bw() +
  labs(x = "estimate",
       y = "density",
       color = "Legend") +
  scale_color_manual(values = colors) +
  theme(legend.position="bottom")
