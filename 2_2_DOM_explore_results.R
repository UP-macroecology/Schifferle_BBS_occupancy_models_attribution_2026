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
                                            jags_out$sims.list$Chi2Open), 2)), cex = 1) # Bayesian p-value = proportion of points above the line

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





# explore flocker horseshoe prior settings: ----

# parameter par_ratio:

# load fitted model(s):

## normal priors for comparison:
load(file.path("results", "American_Goldfinch_cl_lu_p_Stan_no_hs.RData"))
norm_prior <- multi_colex_cl_lu_p_prior

load(file.path("results", "American_Goldfinch_cl_lu_p_Stan_prior_1000its.RData"))
parrat_03 <- multi_colex_cl_lu_p_prior

load(file.path("results", "AG_cl_lu_p_flock_horseshoe(df = 1, par_ratio = 0.01, df_global = 1, scale_slab = 2, df_slab = 4).RData"))
parrat_001 <- multi_colex_cl_lu_p_prior

load(file.path("results", "AG_cl_lu_p_flock_horseshoe(df = 1, par_ratio = 0.1, df_global = 1, scale_slab = 2, df_slab = 4).RData"))
parrat_01 <- multi_colex_cl_lu_p_prior

load(file.path("results", "AG_cl_lu_p_flock_horseshoe(df = 1, par_ratio = 0.2, df_global = 1, scale_slab = 2, df_slab = 4).RData"))
parrat_02 <- multi_colex_cl_lu_p_prior

load(file.path("results", "AG_cl_lu_p_flock_horseshoe(df = 1, par_ratio = 1, df_global = 1, scale_slab = 2, df_slab = 4).RData"))
parrat_1 <- multi_colex_cl_lu_p_prior

load(file.path("results", "AG_cl_lu_p_flock_horseshoe(df = 1, par_ratio = 3, df_global = 1, scale_slab = 2, df_slab = 4).RData"))
parrat_3 <- multi_colex_cl_lu_p_prior

load(file.path("results", "AG_cl_lu_p_flock_horseshoe(df = 1, par_ratio = 99, df_global = 1, scale_slab = 2, df_slab = 4).RData"))
parrat_99 <- multi_colex_cl_lu_p_prior


# plot posteriors for coefficients:

stan_pars <- rownames(fixef(parrat_01, summary = TRUE))

for(p in 1:length(stan_pars)){
  
  print(p)
  
  fl_norm_prior_draws <- unlist(lapply(norm_prior$fit@sim$samples, "[", paste0("b_", stan_pars[p]))) # 4000
  fl_parrat_001_draws <- unlist(lapply(parrat_001$fit@sim$samples, "[", paste0("b_", stan_pars[p]))) # 4000
  fl_parrat_01_draws <- unlist(lapply(parrat_01$fit@sim$samples, "[", paste0("b_", stan_pars[p]))) # 4000
  fl_parrat_02_draws <- unlist(lapply(parrat_02$fit@sim$samples, "[", paste0("b_", stan_pars[p]))) # 4000
  fl_parrat_03_draws <- unlist(lapply(parrat_03$fit@sim$samples, "[", paste0("b_", stan_pars[p]))) # 4000
  fl_parrat_1_draws <- unlist(lapply(parrat_1$fit@sim$samples, "[", paste0("b_", stan_pars[p]))) # 4000
  fl_parrat_3_draws <- unlist(lapply(parrat_3$fit@sim$samples, "[", paste0("b_", stan_pars[p]))) # 4000
  fl_parrat_99_draws <- unlist(lapply(parrat_99$fit@sim$samples, "[", paste0("b_", stan_pars[p]))) # 4000
  
  p_df <- tibble(par = stan_pars[p], 
                 fl_norm_prior_draws,
                 fl_parrat_001_draws,
                 fl_parrat_01_draws, 
                 fl_parrat_02_draws,
                 fl_parrat_03_draws,
                 fl_parrat_1_draws,
                 fl_parrat_3_draws,
                 fl_parrat_99_draws) 
  
  if(p == 1){
    all_stan_p_df <- p_df
  } else{
    all_stan_p_df <- rbind(all_stan_p_df, p_df)
  }
}

colors <- c("norm(0, 2)" = "grey", "par.rat 0.01" = "yellow2", "par.rat 0.1" = "orange", "par.rat 0.2" = "red", 
            "par.rat 0.3" = "green", "par.rat 1" = "darkgreen", "par.rat 3" = "blue", 
            "par.rat 99" = "darkviolet")

# means:
means <- all_stan_p_df %>% 
  group_by(par) %>% 
  summarise(mean = across(fl_norm_prior_draws:fl_parrat_99_draws, mean))

for(i in 1:8){
  
  png(file=file.path("plots", "hs_explorations", paste0("par_ratio_", i, ".png")),
      width=1200, height=900)
  
  p <- ggplot(data = all_stan_p_df) + 
    geom_density(aes(x = fl_norm_prior_draws, colour = "norm(0, 2)")) +
    geom_density(aes(x = fl_parrat_001_draws, colour = "par.rat 0.01")) +
    geom_density(aes(x = fl_parrat_01_draws, colour = "par.rat 0.1")) +
    geom_density(aes(x = fl_parrat_02_draws, colour = "par.rat 0.2")) +
    geom_density(aes(x = fl_parrat_03_draws, colour = "par.rat 0.3")) +
    geom_density(aes(x = fl_parrat_1_draws, colour = "par.rat 1")) +
    geom_density(aes(x = fl_parrat_3_draws, colour = "par.rat 3")) +
    geom_density(aes(x = fl_parrat_99_draws, colour = "par.rat 99")) +
    ggforce::facet_wrap_paginate(facets = ~par, scales = "free", nrow = 3, ncol = 4, page = i) +
    geom_vline(xintercept = 0) +
    geom_vline(data = means, aes(xintercept = mean$fl_norm_prior_draws, colour = "norm(0, 2)"), linetype = "dashed") +
    geom_vline(data = means, aes(xintercept = mean$fl_parrat_001_draws, colour = "par.rat 0.01"), linetype = "dashed") +
    geom_vline(data = means, aes(xintercept = mean$fl_parrat_01_draws, colour = "par.rat 0.1"), linetype = "dashed") +
    geom_vline(data = means, aes(xintercept = mean$fl_parrat_02_draws, colour = "par.rat 0.2"), linetype = "dashed") +
    geom_vline(data = means, aes(xintercept = mean$fl_parrat_03_draws, colour = "par.rat 0.3"), linetype = "dashed") +
    geom_vline(data = means, aes(xintercept = mean$fl_parrat_1_draws, colour = "par.rat 1"), linetype = "dashed") +
    geom_vline(data = means, aes(xintercept = mean$fl_parrat_3_draws, colour = "par.rat 3"), linetype = "dashed") +
    geom_vline(data = means, aes(xintercept = mean$fl_parrat_99_draws, colour = "par.rat 99"), linetype = "dashed") +
    theme_bw() +
    labs(x = "estimate",
         y = "density",
         color = "Legend") +
    scale_color_manual(values = colors) +
    theme(legend.position="bottom", text=element_text(size=21))
  
  print(p)
  
  dev.off()
  
}



# parameter df_global:

# load fitted models:
load(file.path("results", "American_Goldfinch_cl_lu_p_Stan_prior_1000its.RData"))
dfglob_1 <- multi_colex_cl_lu_p_prior

load(file.path("results", "AG_cl_lu_p_flock_horseshoe(df = 1, par_ratio = 0.3, df_global = 0.1, scale_slab = 2, df_slab = 4).RData"))
dfglob_01 <- multi_colex_cl_lu_p_prior

load(file.path("results", "AG_cl_lu_p_flock_horseshoe(df = 1, par_ratio = 0.3, df_global = 0.5, scale_slab = 2, df_slab = 4).RData"))
dfglob_05 <- multi_colex_cl_lu_p_prior

load(file.path("results", "AG_cl_lu_p_flock_horseshoe(df = 1, par_ratio = 0.3, df_global = 2, scale_slab = 2, df_slab = 4).RData"))
dfglob_2 <- multi_colex_cl_lu_p_prior

load(file.path("results", "AG_cl_lu_p_flock_horseshoe(df = 1, par_ratio = 0.3, df_global = 5, scale_slab = 2, df_slab = 4).RData"))
dfglob_5 <- multi_colex_cl_lu_p_prior

# plot posteriors for coefficients:

stan_pars <- rownames(fixef(dfglob_1, summary = TRUE))

for(p in 1:length(stan_pars)){
  
  print(p)
  
  fl_norm_prior_draws <- unlist(lapply(norm_prior$fit@sim$samples, "[", paste0("b_", stan_pars[p]))) # 4000
  fl_dfglob_01_draws <- unlist(lapply(dfglob_01$fit@sim$samples, "[", paste0("b_", stan_pars[p]))) # 4000
  fl_dfglob_05_draws <- unlist(lapply(dfglob_05$fit@sim$samples, "[", paste0("b_", stan_pars[p]))) # 4000
  fl_dfglob_1_draws <- unlist(lapply(dfglob_1$fit@sim$samples, "[", paste0("b_", stan_pars[p]))) # 4000
  fl_dfglob_2_draws <- unlist(lapply(dfglob_2$fit@sim$samples, "[", paste0("b_", stan_pars[p]))) # 4000
  fl_dfglob_5_draws <- unlist(lapply(dfglob_5$fit@sim$samples, "[", paste0("b_", stan_pars[p]))) # 4000

  p_df <- tibble(par = stan_pars[p], 
                 fl_norm_prior_draws,
                 fl_dfglob_01_draws,
                 fl_dfglob_05_draws,
                 fl_dfglob_1_draws, 
                 fl_dfglob_2_draws,
                 fl_dfglob_5_draws) 
  
  if(p == 1){
    all_stan_p_df <- p_df
  } else{
    all_stan_p_df <- rbind(all_stan_p_df, p_df)
  }
}

# means:
means <- all_stan_p_df %>% 
  group_by(par) %>% 
  summarise(mean = across(fl_norm_prior_draws:fl_dfglob_5_draws, mean))

colors <- c("norm(0, 2)" = "grey", "df_glob 0.1" = "yellow2","df_glob 0.5" = "orange",
            "df_glob 1" = "red", "df_glob 2" = "blue", "df_glob 5" = "darkgreen")

for(i in 1:8){
  
  png(file=file.path("plots", "hs_explorations", paste0("df_global_", i, ".png")),
      width=1200, height=900)
  
  p <- ggplot(data = all_stan_p_df) + 
    geom_density(aes(x = fl_norm_prior_draws, colour = "norm(0, 2)")) +
    geom_density(aes(x = fl_dfglob_01_draws, colour = "df_glob 0.1")) +
    geom_density(aes(x = fl_dfglob_05_draws, colour = "df_glob 0.5")) +
    geom_density(aes(x = fl_dfglob_1_draws, colour = "df_glob 1")) +
    geom_density(aes(x = fl_dfglob_2_draws, colour = "df_glob 2")) +
    geom_density(aes(x = fl_dfglob_5_draws, colour = "df_glob 5")) +
    ggforce::facet_wrap_paginate(facets = ~par, scales = "free", nrow = 3, ncol = 4, page = i) +
    geom_vline(xintercept = 0) +
    geom_vline(data = means, aes(xintercept = mean$fl_norm_prior_draws, colour = "norm(0, 2)"), linetype = "dashed") +
    geom_vline(data = means, aes(xintercept = mean$fl_dfglob_01_draws, colour = "df_glob 0.1"), linetype = "dashed") +
    geom_vline(data = means, aes(xintercept = mean$fl_dfglob_05_draws, colour = "df_glob 0.5"), linetype = "dashed") +
    geom_vline(data = means, aes(xintercept = mean$fl_dfglob_1_draws, colour = "df_glob 1"), linetype = "dashed") +
    geom_vline(data = means, aes(xintercept = mean$fl_dfglob_2_draws, colour = "df_glob 2"), linetype = "dashed") +
    geom_vline(data = means, aes(xintercept = mean$fl_dfglob_5_draws, colour = "df_glob 5"), linetype = "dashed") +
    theme_bw() +
    labs(x = "estimate",
         y = "density",
         color = "Legend") +
    scale_color_manual(values = colors) +
    theme(legend.position="bottom", text=element_text(size=21))
  
  print(p)
  
  dev.off()
  
}

# parameter scale_slab:

# load fitted models:
load(file.path("results", "American_Goldfinch_cl_lu_p_Stan_prior_1000its.RData"))
scaleslab_2 <- multi_colex_cl_lu_p_prior

load(file.path("results", "AG_cl_lu_p_flock_horseshoe(df = 1, par_ratio = 0.3, df_global = 1, scale_slab = 0.5, df_slab = 4).RData"))
scaleslab_05 <- multi_colex_cl_lu_p_prior

load(file.path("results", "AG_cl_lu_p_flock_horseshoe(df = 1, par_ratio = 0.3, df_global = 1, scale_slab = 1, df_slab = 4).RData"))
scaleslab_1 <- multi_colex_cl_lu_p_prior

load(file.path("results", "AG_cl_lu_p_flock_horseshoe(df = 1, par_ratio = 0.3, df_global = 1, scale_slab = 3, df_slab = 4).RData"))
scaleslab_3 <- multi_colex_cl_lu_p_prior


# plot posteriors for coefficients:

stan_pars <- rownames(fixef(dfglob_1, summary = TRUE))

for(p in 1:length(stan_pars)){
  
  print(p)
  
  fl_norm_prior_draws <- unlist(lapply(norm_prior$fit@sim$samples, "[", paste0("b_", stan_pars[p]))) # 4000
  fl_scaleslab_2_draws <- unlist(lapply(scaleslab_2$fit@sim$samples, "[", paste0("b_", stan_pars[p]))) # 4000
  fl_scaleslab_05_draws <- unlist(lapply(scaleslab_05$fit@sim$samples, "[", paste0("b_", stan_pars[p]))) # 4000
  fl_scaleslab_1_draws <- unlist(lapply(scaleslab_1$fit@sim$samples, "[", paste0("b_", stan_pars[p]))) # 4000
  fl_scaleslab_3_draws <- unlist(lapply(scaleslab_3$fit@sim$samples, "[", paste0("b_", stan_pars[p]))) # 4000
  
  p_df <- tibble(par = stan_pars[p], 
                 fl_norm_prior_draws,
                 fl_scaleslab_2_draws, 
                 fl_scaleslab_05_draws,
                 fl_scaleslab_1_draws,
                 fl_scaleslab_3_draws) 
  
  if(p == 1){
    all_stan_p_df <- p_df
  } else{
    all_stan_p_df <- rbind(all_stan_p_df, p_df)
  }
}

# means:
means <- all_stan_p_df %>% 
  group_by(par) %>% 
  summarise(mean = across(fl_norm_prior_draws:fl_scaleslab_3_draws, mean))

colors <- c("norm(0, 2)" = "grey", "scale_slab 0.5" = "orange", "scale_slab 1" = "red", 
            "scale_slab 2" = "blue", "scale_slab 3" = "darkgreen")

for(i in 1:8){
  
  png(file=file.path("plots", "hs_explorations", paste0("scale_slab_", i, ".png")),
      width=1200, height=900)
  
  p <- ggplot(data = all_stan_p_df) + 
    geom_density(aes(x = fl_norm_prior_draws, colour = "norm(0, 2)")) +
    geom_density(aes(x = fl_scaleslab_2_draws, colour = "scale_slab 2")) +
    geom_density(aes(x = fl_scaleslab_05_draws, colour = "scale_slab 0.5")) +
    geom_density(aes(x = fl_scaleslab_1_draws, colour = "scale_slab 1")) +
    geom_density(aes(x = fl_scaleslab_3_draws, colour = "scale_slab 3")) +
    ggforce::facet_wrap_paginate(facets = ~par, scales = "free", nrow = 3, ncol = 4, page = i) +
    geom_vline(xintercept = 0) +
    geom_vline(data = means, aes(xintercept = mean$fl_norm_prior_draws, colour = "norm(0, 2)"), linetype = "dashed") +
    geom_vline(data = means, aes(xintercept = mean$fl_scaleslab_2_draws, colour = "scale_slab 2"), linetype = "dashed") +
    geom_vline(data = means, aes(xintercept = mean$fl_scaleslab_05_draws, colour = "scale_slab 0.5"), linetype = "dashed") +
    geom_vline(data = means, aes(xintercept = mean$fl_scaleslab_1_draws, colour = "scale_slab 1"), linetype = "dashed") +
    geom_vline(data = means, aes(xintercept = mean$fl_scaleslab_3_draws, colour = "scale_slab 3"), linetype = "dashed") +
    theme_bw() +
    labs(x = "estimate",
         y = "density",
         color = "Legend") +
    scale_color_manual(values = colors) +
    theme(legend.position="bottom", text=element_text(size=21))
  
  print(p)
  
  dev.off()

}

# parameter df_slab:

# load fitted models:
load(file.path("results", "American_Goldfinch_cl_lu_p_Stan_prior_1000its.RData"))
dfslab_4 <- multi_colex_cl_lu_p_prior

load(file.path("results", "AG_cl_lu_p_flock_horseshoe(df = 1, par_ratio = 0.3, df_global = 1, scale_slab = 2, df_slab = 2).RData"))
dfslab_2 <- multi_colex_cl_lu_p_prior

# plot posteriors for coefficients:

stan_pars <- rownames(fixef(dfglob_1, summary = TRUE))

for(p in 1:length(stan_pars)){
  
  print(p)
  
  fl_norm_prior_draws <- unlist(lapply(norm_prior$fit@sim$samples, "[", paste0("b_", stan_pars[p]))) # 4000
  fl_dfslab_2_draws <- unlist(lapply(dfslab_2$fit@sim$samples, "[", paste0("b_", stan_pars[p]))) # 4000
  fl_dfslab_4_draws <- unlist(lapply(dfslab_4$fit@sim$samples, "[", paste0("b_", stan_pars[p]))) # 4000

  p_df <- tibble(par = stan_pars[p], 
                 fl_norm_prior_draws,
                 fl_dfslab_2_draws, 
                 fl_dfslab_4_draws) 
  
  if(p == 1){
    all_stan_p_df <- p_df
  } else{
    all_stan_p_df <- rbind(all_stan_p_df, p_df)
  }
}

# means:
means <- all_stan_p_df %>% 
  group_by(par) %>% 
  summarise(mean = across(fl_norm_prior_draws:fl_dfslab_4_draws, mean))

colors <- c("norm(0, 2)" = "grey", "df_slab 2" = "blue", "df_slab 4" = "red")

for(i in 1:8){
  
  png(file=file.path("plots", "hs_explorations", paste0("df_slab_", i, ".png")),
      width=1200, height=900)

  p <- ggplot(data = all_stan_p_df) + 
    geom_density(aes(x = fl_norm_prior_draws, colour = "norm(0, 2)")) +
    geom_density(aes(x = fl_dfslab_2_draws, colour = "df_slab 2")) +
    geom_density(aes(x = fl_dfslab_4_draws, colour = "df_slab 4")) +
    ggforce::facet_wrap_paginate(facets = ~par, scales = "free", nrow = 3, ncol = 4, page = i) +
    geom_vline(xintercept = 0) +
    geom_vline(data = means, aes(xintercept = mean$fl_norm_prior_draws, colour = "norm(0, 2)"), linetype = "dashed") +
    geom_vline(data = means, aes(xintercept = mean$fl_dfslab_2_draws, colour = "df_slab 2"), linetype = "dashed") +
    geom_vline(data = means, aes(xintercept = mean$fl_dfslab_4_draws, colour = "df_slab 4"), linetype = "dashed") +
    theme_bw() +
    labs(x = "estimate",
         y = "density",
         color = "Legend") +
    scale_color_manual(values = colors) +
    theme(legend.position="bottom", text=element_text(size=21))
  
  print(p)
  
  dev.off()

}


# parameter df:

# load fitted models:
load(file.path("results", "American_Goldfinch_cl_lu_p_Stan_prior_1000its.RData"))
df_1 <- multi_colex_cl_lu_p_prior

load(file.path("results", "AG_cl_lu_p_flock_horseshoe(df = 0.5, par_ratio = 0.3, df_global = 1, scale_slab = 2, df_slab = 4).RData"))
df_05 <- multi_colex_cl_lu_p_prior

load(file.path("results", "AG_cl_lu_p_flock_horseshoe(df = 2, par_ratio = 0.3, df_global = 1, scale_slab = 2, df_slab = 4).RData"))
df_2 <- multi_colex_cl_lu_p_prior

load(file.path("results", "AG_cl_lu_p_flock_horseshoe(df = 3, par_ratio = 0.3, df_global = 1, scale_slab = 2, df_slab = 4).RData"))
df_3 <- multi_colex_cl_lu_p_prior

# plot posteriors for coefficients:

stan_pars <- rownames(fixef(df_1, summary = TRUE))

for(p in 1:length(stan_pars)){
  
  print(p)
  
  fl_norm_prior_draws <- unlist(lapply(norm_prior$fit@sim$samples, "[", paste0("b_", stan_pars[p]))) # 4000
  fl_df_1_draws <- unlist(lapply(df_1$fit@sim$samples, "[", paste0("b_", stan_pars[p]))) # 4000
  fl_df_05_draws <- unlist(lapply(df_05$fit@sim$samples, "[", paste0("b_", stan_pars[p]))) # 4000
  fl_df_2_draws <- unlist(lapply(df_2$fit@sim$samples, "[", paste0("b_", stan_pars[p]))) # 4000
  fl_df_3_draws <- unlist(lapply(df_3$fit@sim$samples, "[", paste0("b_", stan_pars[p]))) # 4000
  
  p_df <- tibble(par = stan_pars[p], 
                 fl_norm_prior_draws,
                 fl_df_1_draws, 
                 fl_df_05_draws,
                 fl_df_2_draws,
                 fl_df_3_draws) 
  
  if(p == 1){
    all_stan_p_df <- p_df
  } else{
    all_stan_p_df <- rbind(all_stan_p_df, p_df)
  }
}

# means:
means <- all_stan_p_df %>% 
  group_by(par) %>% 
  summarise(mean = across(fl_norm_prior_draws:fl_df_3_draws, mean))

colors <- c("norm(0, 2)" = "grey", "df 0.5" = "yellow2", "df 1" = "red",
            "df 2" = "blue", "df 3" = "darkgreen")

for(i in 1:8){
  
  png(file=file.path("plots", "hs_explorations", paste0("df_", i, ".png")),
      width=1200, height=900)
  
  p <- ggplot(data = all_stan_p_df) + 
    geom_density(aes(x = fl_norm_prior_draws, colour = "norm(0, 2)")) +
    geom_density(aes(x = fl_df_05_draws, colour = "df 0.5")) +
    geom_density(aes(x = fl_df_1_draws, colour = "df 1")) +
    geom_density(aes(x = fl_df_2_draws, colour = "df 2")) +
    geom_density(aes(x = fl_df_3_draws, colour = "df 3")) +
    ggforce::facet_wrap_paginate(facets = ~par, scales = "free", nrow = 3, ncol = 4, page = i) +
    geom_vline(xintercept = 0) +
    geom_vline(data = means, aes(xintercept = mean$fl_norm_prior_draws, colour = "norm(0, 2)"), linetype = "dashed") +
    geom_vline(data = means, aes(xintercept = mean$fl_df_05_draws, colour = "df 0.5"), linetype = "dashed") +
    geom_vline(data = means, aes(xintercept = mean$fl_df_1_draws, colour = "df 1"), linetype = "dashed") +
    geom_vline(data = means, aes(xintercept = mean$fl_df_2_draws, colour = "df 2"), linetype = "dashed") +
    geom_vline(data = means, aes(xintercept = mean$fl_df_3_draws, colour = "df 3"), linetype = "dashed") +
    theme_bw() +
    labs(x = "estimate",
         y = "density",
         color = "Legend") +
    scale_color_manual(values = colors) +
    theme(legend.position="bottom", text=element_text(size=21))
  
  print(p)
  
  dev.off()
  
}

# check whether which variable is significant changes depending on parameter values:

# models to compare:
model_fits <- list(norm_prior, parrat_03, parrat_001, parrat_01, parrat_02, parrat_1, parrat_3, parrat_99, 
                dfglob_1, dfglob_01, dfglob_05, dfglob_2, dfglob_5, 
                scaleslab_2, scaleslab_05, scaleslab_1, scaleslab_3, 
                dfslab_4, dfslab_2, 
                df_1, df_05, df_2, df_3)
names(model_fits) <- c("norm_prior", "parrat_03", "parrat_001", "parrat_01", "parrat_02", "parrat_1", "parrat_3", "parrat_99", 
                      "dfglob_1", "dfglob_01", "dfglob_05", "dfglob_2", "dfglob_5", 
                      "scaleslab_2", "scaleslab_05", "scaleslab_1", "scaleslab_3", 
                      "dfslab_4", "dfslab_2", 
                      "df_1", "df_05", "df_2", "df_3")

get_sign_pars <- function(x){
  
  sign_pars_df <- as_tibble(fixef(model_fits[[1]], summary = TRUE), rownames = "par") %>% 
    mutate(sign = if_else(Q2.5 * Q97.5 > 0, 1, 0)) %>% 
    select(par, sign)
}

as.data.frame(lapply(model_fits, FUN = get_sign_pars)) %>% 
  select(c(1, seq(2, ncol(.), by = 2))) %>% 
  mutate(total = rowSums(across(where(is.numeric)))) %>% 
  distinct(total)
# all 23 models agree

# explore divergent transitions: ----

# following Gabry et al. 2019:
library(bayesplot)
post <- brms::as_draws(dfglob_01) # df_1
hmc_diagnostics <- nuts_params(dfglob_01) # df_1

color_scheme_set("darkgray")
div_style <- parcoord_style_np(div_color = "green", div_size = 0.15, div_alpha = 0.4)

mcmc_parcoord(post, 
              regex_pars = c("b_occ"),
              np = hmc_diagnostics,
              np_style = div_style)

mcmc_parcoord(post, 
              regex_pars = c("b_col"),
              np = hmc_diagnostics,
              np_style = div_style)

mcmc_parcoord(post, 
              regex_pars = c("b_ex"),
              np = hmc_diagnostics,
              np_style = div_style)

mcmc_parcoord(post, 
              regex_pars = c("b_route"),
              np = hmc_diagnostics,
              np_style = div_style)

# according to Gabry 2019 this pattern of divergences may indicate that they are false positives

div_style <- scatter_style_np(div_color = "green", div_size = 2.5, div_alpha = 0.75)
mcmc_scatter(
  post,
  size = 1.5,
  alpha = 2/3,
  pars = c("b_occ_Ipr_spring_3yrsE2", "b_ex_Intercept"), 
  np = hmc_diagnostics,
  np_style = div_style
  )


## compare flocker normal priors vs. horsehoe priors: ----

load(file.path("results", "AG_cl_lu_p_flock_logistic(0,1)_normal(0,1).RData"))
log01_norm01 <- out
load(file.path("results", "AG_cl_lu_p_flock_logistic(0,1)_normal(0,3).RData"))
log01_norm03 <- out
load(file.path("results", "AG_cl_lu_p_flock_logistic(0,2)_normal(0,2).RData"))
log02_norm02 <- out
load(file.path("results", "AG_cl_lu_p_flock_logistic(0,3)_normal(0,2).RData"))
log03_norm02 <- out
load(file.path("results", "AG_cl_lu_p_flock_logistic(0,2)_horseshoe(df = 3, par_ratio = 0.3, df_global = 1, scale_slab = 2, df_slab = 4).RData"))
log02_hs3 <- out # 3 divergent transitions

stan_pars <- rownames(fixef(log01_norm01, summary = TRUE))

for(p in 1:length(stan_pars)){
  
  print(p)
  
  l01n01_draws <- unlist(lapply(log01_norm01$fit@sim$samples, "[", paste0("b_", stan_pars[p]))) # 4000
  l01n03_draws <- unlist(lapply(log01_norm03$fit@sim$samples, "[", paste0("b_", stan_pars[p]))) # 4000
  l02n02_draws <- unlist(lapply(log02_norm02$fit@sim$samples, "[", paste0("b_", stan_pars[p]))) # 4000
  l03n02_draws <- unlist(lapply(log03_norm02$fit@sim$samples, "[", paste0("b_", stan_pars[p]))) # 4000
  l02hs3_draws <- unlist(lapply(log02_hs3$fit@sim$samples, "[", paste0("b_", stan_pars[p]))) # 4000
  
  p_df <- tibble(par = stan_pars[p], 
                 l01n01_draws, 
                 l01n03_draws,
                 l02n02_draws,
                 l03n02_draws,
                 l02hs3_draws) 
  
  if(p == 1){
    all_stan_p_df <- p_df
  } else{
    all_stan_p_df <- rbind(all_stan_p_df, p_df)
  }
}

colors <- c("logistic(0,1)_normal(0,1)" = "yellow2", 
            "logistic(0,1)_normal(0,3)" = "orange",
            "logistic(0,2)_normal(0,2)" = "red",
            "logistic(0,3)_normal(0,2)" = "purple",
            "logistic(0,2)_horseshoe(3)" = "blue"
            )

ggplot(data = all_stan_p_df) + 
  geom_density(aes(x = l01n01_draws, colour = "logistic(0,1)_normal(0,1)")) +
  geom_density(aes(x = l01n03_draws, colour = "logistic(0,1)_normal(0,3)")) +
  geom_density(aes(x = l02n02_draws, colour = "logistic(0,2)_normal(0,2)")) +
  geom_density(aes(x = l03n02_draws, colour = "logistic(0,3)_normal(0,2)")) +
  geom_density(aes(x = l02hs3_draws, colour = "logistic(0,2)_horseshoe(3)")) +
  ggforce::facet_wrap_paginate(facets = ~par, scales = "free", nrow = 3, ncol = 4, page = 8) +
  geom_vline(xintercept = 0) +
  theme_bw() +
  labs(x = "estimate",
       y = "density",
       color = "Legend") +
  scale_color_manual(values = colors) +
  theme(legend.position="bottom")


# compare prior and posterior distribution: ------------------------------------

# "fit" model again using only prior (not data):


# load data:
load(file = file.path("data", "route_year_env_data.RData"))
# scale covariates:
route_sel_env_dt_scaled <- route_sel_env_dt_final %>% 
  mutate(across(bio2:pr_winter_3yrs, ~ (scale(.)) %>% as.vector()))
# route-year-species information (only surveyed)
load(file = file.path("data", "BBS_for_occ_spec_records.RData")) # output of 1_0_reformat_BBS_data.R

# assemble data:
spec <- "American Goldfinch"
nyears <- length(unique(route_sel_env_dt_final$Year)) # 25
nsurveys <- 5
nsites <- length(unique(route_sel_env_dt_final$RTENO)) # 476

# species presences:
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

# reformat obs. as array sites x surveys x years (same as for JAGS):
years <- seq(min(occ_dt_spec$Year), max(occ_dt_spec$Year))
y_array <- array(NA, dim = c(nsites, nsurveys, nyears))
for (t in 1:nyears){
  y_array[1:nsites, 1:nsurveys, t] <- as.matrix(occ_dt_spec[which(occ_dt_spec$Year == years[t]), c(paste0("Count", seq(10, 50, 10)))])
}

# reformat environmental covariates:
route_sel_env_dt_scaled
env_cov <- vector("list", length = nyears)
for (t in 1:nyears){
  env_cov[[t]] <- route_sel_env_dt_scaled[which(route_sel_env_dt_scaled$Year == years[t]), 
                                          c("bio1", "bio2", "bio3", "pr_spring", "pr_summer","pr_autumn", 
                                            "pr_winter", "bio15","bio1_3yrs", "bio2_3yrs", "bio3_3yrs", 
                                            "pr_spring_3yrs", "pr_summer_3yrs", "pr_autumn_3yrs", "pr_winter_3yrs", 
                                            "bio15_3yrs",
                                            "sum_annual_crops", "secdn","pastr", "urban", "primn",
                                            "sum_annual_crops_3yrs", "secdn_3yrs", "pastr_3yrs", "urban_3yrs", "primn_3yrs")]
}

# covariate for detection probability:
det_cov <- vector("list", length = 1)
names(det_cov) <- "route_section"
det_cov$route_section <- array(NA, dim = c(nsites, nsurveys, nyears))
det_cov$route_section[ , , 1:nyears] <- matrix(rep(c("Sect1", "Sect2", "Sect3", "Sect4", "Sect5"), nsites), nsites, byrow = TRUE)

# make flocker data:
fd <- make_flocker_data_dynamic(
  obs = y_array,
  unit_covs = env_cov, 
  event_covs = det_cov, 
  quiet = FALSE
)

# fit model with prior only:
out_prior <- flock(
  f_occ = ~ bio1_3yrs + bio2_3yrs + bio3_3yrs + pr_spring_3yrs + pr_summer_3yrs + pr_autumn_3yrs + 
    pr_winter_3yrs + bio15_3yrs + I(bio1_3yrs^2) + I(bio2_3yrs^2) + I(bio3_3yrs^2) + 
    I(pr_spring_3yrs^2) + I(pr_summer_3yrs^2) + I(pr_autumn_3yrs^2) + I(pr_winter_3yrs^2) + I(bio15_3yrs^2) +
    sum_annual_crops_3yrs + secdn_3yrs + pastr_3yrs + urban_3yrs + primn_3yrs +
    I(sum_annual_crops_3yrs^2) + I(secdn_3yrs^2) + I(pastr_3yrs^2) + I(urban_3yrs^2) + I(primn_3yrs^2),
  f_det = ~ route_section,
  f_col = ~ bio1 + bio2 + bio3 + pr_spring + pr_summer + pr_autumn + pr_winter + bio15 +
    I(bio1^2) + I(bio2^2) + I(bio3^2) + I(pr_spring^2) + I(pr_summer^2) + I(pr_autumn^2) + 
    I(pr_winter^2) + I(bio15^2) + sum_annual_crops + secdn + pastr + urban + primn +
    I(sum_annual_crops^2) + I(secdn^2) + I(pastr^2) + I(urban^2) + I(primn^2),
  f_ex = ~ bio1 + bio2 + bio3 + pr_spring + pr_summer + pr_autumn + pr_winter + bio15 +
    I(bio1^2) + I(bio2^2) + I(bio3^2) + I(pr_spring^2) + I(pr_summer^2) + I(pr_autumn^2) + 
    I(pr_winter^2) + I(bio15^2) + sum_annual_crops + secdn + pastr + urban + primn +
    I(sum_annual_crops^2) + I(secdn^2) + I(pastr^2) + I(urban^2) + I(primn^2),
  flocker_data = fd,
  prior = c(brms::set_prior("logistic(0,3)", class = "Intercept") + 
              brms::set_prior("logistic(0,3)", class = "Intercept", dpar = "occ"),
            brms::set_prior("logistic(0,3)", class = "Intercept", dpar = "colo"),
            brms::set_prior("logistic(0,3)", class = "Intercept", dpar = "ex"),
            brms::set_prior("normal(0,2)", class = "b"),
            brms::set_prior("normal(0,2)", dpar = "occ", class = "b"),
            brms::set_prior("normal(0,2)", dpar = "colo", class = "b"),
            brms::set_prior("normal(0,2)", dpar = "ex", class = "b")),
  multiseason = "colex",
  multi_init = "explicit",
  backend = "cmdstanr",
  cores = 4,
  chains = 4,
  warmup = 250,
  iter = 250 + 1000,
  sample_prior = "only"
)
out_prior

# load model fitted with data:

load(file.path("results", "AG_cl_lu_p_flock_logistic(0,3)_normal(0,2).RData"))
log03_norm02 <- out

# compare prior and posterior distributions:

stan_pars <- rownames(fixef(out_prior, summary = TRUE))

for(p in 1:length(stan_pars)){
  
  print(p)
  
  posterior_draws <- unlist(lapply(log03_norm02$fit@sim$samples, "[", paste0("b_", stan_pars[p]))) # 4000
  prior_draws <- unlist(lapply(out_prior$fit@sim$samples, "[", paste0("b_", stan_pars[p]))) # 4000

  p_df <- tibble(par = stan_pars[p], 
                 posterior_draws, 
                 prior_draws) 
  
  if(p == 1){
    all_stan_p_df <- p_df
  } else{
    all_stan_p_df <- rbind(all_stan_p_df, p_df)
  }
}
all_stan_p_df

colors <- c("prior" = "cornflowerblue", 
            "posterior" = "red3")




for(i in 1:8){
  
  png(file=file.path("plots", "prior_post_comparison", paste0("logistic(0,3)_normal(0,2)_", i, ".png")),
      width=1200, height=900)
  
  p <- ggplot(data = all_stan_p_df) + 
    geom_density(aes(x = prior_draws, fill = "prior"), alpha = 0.6, colour = "transparent") +
    geom_density(aes(x = posterior_draws, fill = "posterior"), alpha = 0.6, colour = "transparent") +
    ggforce::facet_wrap_paginate(facets = ~par, scales = "free", nrow = 3, ncol = 4, page = i) +
    geom_vline(xintercept = 0) +
    theme_bw() +
    labs(x = "estimate",
         y = "density",
         color = "Legend") +
    scale_fill_manual(values = colors) +
    theme(legend.position="bottom") +
    theme(legend.position="bottom", text=element_text(size=21))
  
  print(p)
  
  dev.off()
  
}

## same for horseshoe prior: ----

# fit model with prior only:
out_prior_hs <- flock(
  f_occ = ~ bio1_3yrs + bio2_3yrs + bio3_3yrs + pr_spring_3yrs + pr_summer_3yrs + pr_autumn_3yrs + 
    pr_winter_3yrs + bio15_3yrs + I(bio1_3yrs^2) + I(bio2_3yrs^2) + I(bio3_3yrs^2) + 
    I(pr_spring_3yrs^2) + I(pr_summer_3yrs^2) + I(pr_autumn_3yrs^2) + I(pr_winter_3yrs^2) + I(bio15_3yrs^2) +
    sum_annual_crops_3yrs + secdn_3yrs + pastr_3yrs + urban_3yrs + primn_3yrs +
    I(sum_annual_crops_3yrs^2) + I(secdn_3yrs^2) + I(pastr_3yrs^2) + I(urban_3yrs^2) + I(primn_3yrs^2),
  f_det = ~ route_section,
  f_col = ~ bio1 + bio2 + bio3 + pr_spring + pr_summer + pr_autumn + pr_winter + bio15 +
    I(bio1^2) + I(bio2^2) + I(bio3^2) + I(pr_spring^2) + I(pr_summer^2) + I(pr_autumn^2) + 
    I(pr_winter^2) + I(bio15^2) + sum_annual_crops + secdn + pastr + urban + primn +
    I(sum_annual_crops^2) + I(secdn^2) + I(pastr^2) + I(urban^2) + I(primn^2),
  f_ex = ~ bio1 + bio2 + bio3 + pr_spring + pr_summer + pr_autumn + pr_winter + bio15 +
    I(bio1^2) + I(bio2^2) + I(bio3^2) + I(pr_spring^2) + I(pr_summer^2) + I(pr_autumn^2) + 
    I(pr_winter^2) + I(bio15^2) + sum_annual_crops + secdn + pastr + urban + primn +
    I(sum_annual_crops^2) + I(secdn^2) + I(pastr^2) + I(urban^2) + I(primn^2),
  flocker_data = fd,
  prior = c(brms::set_prior("logistic(0,2)", class = "Intercept") + 
              brms::set_prior("logistic(0,2)", class = "Intercept", dpar = "occ"),
            brms::set_prior("logistic(0,2)", class = "Intercept", dpar = "colo"),
            brms::set_prior("logistic(0,2)", class = "Intercept", dpar = "ex"),
            brms::set_prior("horseshoe(df = 3, par_ratio = 0.3, df_global = 1, scale_slab = 2, df_slab = 4)", class = "b"),
            brms::set_prior("horseshoe(df = 3, par_ratio = 0.3, df_global = 1, scale_slab = 2, df_slab = 4)", dpar = "occ", class = "b"),
            brms::set_prior("horseshoe(df = 3, par_ratio = 0.3, df_global = 1, scale_slab = 2, df_slab = 4)", dpar = "colo", class = "b"),
            brms::set_prior("horseshoe(df = 3, par_ratio = 0.3, df_global = 1, scale_slab = 2, df_slab = 4)", dpar = "ex", class = "b")),
  multiseason = "colex",
  multi_init = "explicit",
  backend = "cmdstanr",
  cores = 4,
  chains = 4,
  warmup = 250,
  iter = 250 + 1000,
  sample_prior = "only"
)
out_prior_hs

# load model fitted with data:

load(file.path("results", "AG_cl_lu_p_flock_logistic(0,2)_horseshoe(df = 3, par_ratio = 0.3, df_global = 1, scale_slab = 2, df_slab = 4).RData"))
log02_hs3 <- out

# compare prior and posterior distributions:

stan_pars <- rownames(fixef(out_prior_hs, summary = TRUE))

for(p in 1:length(stan_pars)){
  
  print(p)
  
  posterior_draws <- unlist(lapply(log02_hs3$fit@sim$samples, "[", paste0("b_", stan_pars[p]))) # 4000
  prior_draws <- unlist(lapply(out_prior_hs$fit@sim$samples, "[", paste0("b_", stan_pars[p]))) # 4000
  
  p_df <- tibble(par = stan_pars[p], 
                 posterior_draws, 
                 prior_draws) 
  
  if(p == 1){
    all_stan_p_df <- p_df
  } else{
    all_stan_p_df <- rbind(all_stan_p_df, p_df)
  }
}
all_stan_p_df

colors <- c("prior" = "cornflowerblue", 
            "posterior" = "red3")




for(i in 1:8){
  
  png(file=file.path("plots", "prior_post_comparison", paste0("logistic(0,2)_hs(3)", i, ".png")),
      width=1200, height=900)
  
  p <- ggplot(data = all_stan_p_df) + 
    geom_density(aes(x = prior_draws, fill = "prior"), alpha = 0.6, colour = "transparent") +
    geom_density(aes(x = posterior_draws, fill = "posterior"), alpha = 0.6, colour = "transparent") +
    ggforce::facet_wrap_paginate(facets = ~par, scales = "free", nrow = 3, ncol = 4, page = i) +
    geom_vline(xintercept = 0) +
    theme_bw() +
    labs(x = "estimate",
         y = "density",
         color = "Legend") +
    scale_fill_manual(values = colors) +
    theme(legend.position="bottom") +
    theme(legend.position="bottom", text=element_text(size=21)) +
    ylim(c(0, 10))
  
  print(p)
  
  dev.off()
  
}
