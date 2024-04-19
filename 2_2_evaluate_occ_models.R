# explore and evaluate occupancy models:

# packages: ----

library(jagsUI)
library(dplyr)
library(precrec)
library(ggplot2)
library(stringr)

# load data: ----

# route-year-species information (only surveyed)
load(file = file.path("data", "BBS_for_occ_spec_records.RData")) # output of 1_0_reformat_BBS_data.R
bbs_dt_occ
nrow(bbs_dt_occ)

# selected species:
load(file = file.path("data", "final_species_selection.RData")) # output of 1_2_species_selection.R
species_selection_final

# # selected routes and years:
# load(file = file.path("data", "route_year_env_data.RData"))
# route_sel_env_dt_final



# output of model fitting: ----

## quick checks regarding convergence: ----

#load(file.path("results", "Grasshopper_Sparrow_cl_alt3_quadr_lasso_GoF_inits.RData"))
#load(file.path("results", "Black-billed_Cuckoo_cl_alt3_quadr_lasso_GoF_inits.RData"))
#load(file.path("results", "American_Robin_cl_alt3_quadr_lasso_GoF_inits.RData"))
#load(file.path("results", "Black-capped_Chickadee_cl_alt3_quadr_lasso_GoF_inits.RData"))
#load(file.path("results", "Eastern_Phoebe_cl_alt3_quadr_lasso_GoF_inits.RData"))
#load(file.path("results", "Black-capped_Chickadee_cl_alt4_quadr_lasso_GoF_inits.RData"))
#load(file.path("results", "Black-billed_Cuckoo_cl_alt4_quadr_lasso_GoF_inits.RData"))
#load(file.path("results", "Eastern_Phoebe_cl_alt4_quadr_lasso_GoF_inits.RData"))
#load(file.path("results", "American_Goldfinch_cl_alt4_quadr_lasso_GoF_inits.RData"))
#load(file.path("results", "American_Goldfinch_cl_alt3_quadr_lasso_GoF_inits.RData"))
#load(file.path("results", "Black-capped_Chickadee_cl_alt4_lu_quadr_p_sect_lasso_GoF_inits_less_pars_butGoF.RData"))
#load(file.path("results", "American_Goldfinch_cl_alt4_lu_quadr_p_sect_lasso_GoF_inits_less_pars_butGoF.RData"))
#load(file.path("results", "Eastern_Kingbird_cl_alt4_lu_quadr_p_sect_lasso_GoF_inits_less_pars_butGoF.RData"))
#load(file.path("results", "Grasshopper_Sparrow_cl_alt4_lu_quadr_p_sect_lasso_GoF_inits_less_pars_butGoF.RData"))
load(file.path("results", "Black-billed_Cuckoo_cl_alt4_lu_quadr_p_sect_lasso_GoF_inits_less_pars_butGoF.RData"))

# psi, gamma, eps:
chain_range <- apply(X = out$sims.list$psi, MARGIN = 2:3, FUN = range)
#dim(chain_range) # 2 476 24
min_df <- as.data.frame(chain_range[1,,])
max_df <- as.data.frame(chain_range[2,,])
diff_df <- max_df-min_df
#dim(diff_df)
diff_nomix_df <- diff_df[which(out$Rhat$psi > 1.1, arr.ind = TRUE)] # difference only where chains haven't mixed
sort(diff_nomix_df, decreasing = TRUE)

# psi1:
chain_range <- apply(X = out$sims.list$psi1, MARGIN = 2, FUN = range)
min_df <- as.data.frame(chain_range[1,])
max_df <- as.data.frame(chain_range[2,])
diff_df <- max_df-min_df
diff_nomix_df <- diff_df[which(out$Rhat$psi1 > 1.1, arr.ind = TRUE)] # difference only where chains haven't mixed
sort(diff_nomix_df, decreasing = TRUE)


# evaluate spatial predictive ability: -----------------------------------------

# model output:
out <- HF_quad_res_6000

# AUC over all years or per year and track change over time?

# assemble data: true status (xx) -- probability of occupancy -- probability of non-occupancy:

AUC_dt <- data.frame("site" = route_sel_env_dt_final$RTENO, 
                     "year" = route_sel_env_dt_final$Year, 
                     "obs_occ" = NA, 
                     "prob_occ" = NA, 
                     "prob_non_occ" = NA)

# observed status:

# species:
spec <- "House Finch"

# data:
presences_spec <- bbs_dt_occ %>% 
  select(c(English_Common_Name, RTENO, Year, paste0("Count", seq(10, 50, 10)))) %>% 
  filter(English_Common_Name == spec)

# match to routes-year-env:
occ_dt_spec <- route_sel_env_dt_final %>% 
  select(c(RTENO, Year, Latitude, Longitude, Surveyed)) %>% 
  # add observations:
  collapse::join(presences_spec, on = c("RTENO", "Year"), how = "left") %>% 
  # if route was surveyed but species not observed, replace NA with 0:
  mutate(across(Count10:Count50, ~ 
                  case_when(Surveyed == 1 & is.na(.) ~ 0,
                            .default = .))) %>%
  # convert bird counts to presence / absence:
  mutate(across(Count10:Count50, ~ 
                  case_when(. > 1 ~ 1,
                            .default = .))) %>% 
  # sum presence/absence on route:
  mutate(occ_route = rowSums(across(Count10:Count50))) %>%
  # convert presence absence:
  mutate(occ_route = ifelse(occ_route > 0, 1, 0))

AUC_dt$obs_occ <- as.factor(occ_dt_spec$occ_route)
nrow(AUC_dt)

dim(out$sims.list$psi) # first dimension = values sampled from the posterior distributions of each monitored parameter = draws posterior * 3 (chains)?

AUC_dt2 <- AUC_dt[order(AUC_dt$year, decreasing = FALSE), ] # do site number match? looks good

# predicted mean probability of occupancy: xx

AUC_dt2$prob_occ <- out$summary[grepl(pattern = "^psi\\[", row.names(out$summary)),"mean"]
AUC_dt2$prob_non_occ <- 1 - AUC_dt2$prob_occ
AUC_dt3 <- AUC_dt2[complete.cases(AUC_dt2),]# only complete rows:

# calculate overall AUC:
sscurves <- evalmod(scores = AUC_dt3$prob_occ, labels = AUC_dt3$obs_occ, posclass = 1)

# Show ROC and Precision-Recall plots
plot(sscurves)

overall_auc <- auc(sscurves)
overall_auc




# compare_Jags_Stan_cl_only <- function(jags_output, stan_output){
#   
#   stan_estimates <- fixef(stan_output, summary = TRUE, robust = FALSE, probs = c(0.025, 0.975))
#   jags_estimates <- jags_output$summary[which(grepl("(alpha)|(beta)", rownames(jags_output$summary))), c("mean", "2.5%", "97.5%")]
#   
#   jags_estimates
#   
#   
#   jags_reformatted <- rbind(c(NA,NA, NA),
#                             jags_estimates["alpha_psi", c("mean", "2.5%", "97.5%")],
#                             jags_estimates["alpha_gamma", c("mean", "2.5%", "97.5%")],
#                             jags_estimates["alpha_eps", c("mean", "2.5%", "97.5%")],
#                             jags_estimates["beta_psi[1]", c("mean", "2.5%", "97.5%")],
#                             jags_estimates["beta_psi[2]", c("mean", "2.5%", "97.5%")],
#                             jags_estimates["beta_psi[3]", c("mean", "2.5%", "97.5%")],
#                             jags_estimates["beta_psi[4]", c("mean", "2.5%", "97.5%")],
#                             jags_estimates["beta_psi[5]", c("mean", "2.5%", "97.5%")],
#                             jags_estimates["beta_psi[6]", c("mean", "2.5%", "97.5%")],
#                             jags_estimates["beta_psi[7]", c("mean", "2.5%", "97.5%")],
#                             jags_estimates["beta_psi[8]", c("mean", "2.5%", "97.5%")],
#                             jags_estimates["beta_psi[9]", c("mean", "2.5%", "97.5%")],
#                             jags_estimates["beta_psi[10]", c("mean", "2.5%", "97.5%")],
#                             jags_estimates["beta_psi[11]", c("mean", "2.5%", "97.5%")],
#                             jags_estimates["beta_psi[12]", c("mean", "2.5%", "97.5%")],
#                             jags_estimates["beta_psi[13]", c("mean", "2.5%", "97.5%")],
#                             jags_estimates["beta_psi[14]", c("mean", "2.5%", "97.5%")],
#                             jags_estimates["beta_psi[15]", c("mean", "2.5%", "97.5%")],
#                             jags_estimates["beta_psi[16]", c("mean", "2.5%", "97.5%")],
#                             jags_estimates["beta_gamma[1]", c("mean", "2.5%", "97.5%")],
#                             jags_estimates["beta_gamma[2]", c("mean", "2.5%", "97.5%")],
#                             jags_estimates["beta_gamma[3]", c("mean", "2.5%", "97.5%")],
#                             jags_estimates["beta_gamma[4]", c("mean", "2.5%", "97.5%")],
#                             jags_estimates["beta_gamma[5]", c("mean", "2.5%", "97.5%")],
#                             jags_estimates["beta_gamma[6]", c("mean", "2.5%", "97.5%")],
#                             jags_estimates["beta_gamma[7]", c("mean", "2.5%", "97.5%")],
#                             jags_estimates["beta_gamma[8]", c("mean", "2.5%", "97.5%")],
#                             jags_estimates["beta_gamma[9]", c("mean", "2.5%", "97.5%")],
#                             jags_estimates["beta_gamma[10]", c("mean", "2.5%", "97.5%")],
#                             jags_estimates["beta_gamma[11]", c("mean", "2.5%", "97.5%")],
#                             jags_estimates["beta_gamma[12]", c("mean", "2.5%", "97.5%")],
#                             jags_estimates["beta_gamma[13]", c("mean", "2.5%", "97.5%")],
#                             jags_estimates["beta_gamma[14]", c("mean", "2.5%", "97.5%")],
#                             jags_estimates["beta_gamma[15]", c("mean", "2.5%", "97.5%")],
#                             jags_estimates["beta_gamma[16]", c("mean", "2.5%", "97.5%")],
#                             jags_estimates["beta_eps[1]", c("mean", "2.5%", "97.5%")],
#                             jags_estimates["beta_eps[2]", c("mean", "2.5%", "97.5%")],
#                             jags_estimates["beta_eps[3]", c("mean", "2.5%", "97.5%")],
#                             jags_estimates["beta_eps[4]", c("mean", "2.5%", "97.5%")],
#                             jags_estimates["beta_eps[5]", c("mean", "2.5%", "97.5%")],
#                             jags_estimates["beta_eps[6]", c("mean", "2.5%", "97.5%")],
#                             jags_estimates["beta_eps[7]", c("mean", "2.5%", "97.5%")],
#                             jags_estimates["beta_eps[8]", c("mean", "2.5%", "97.5%")],
#                             jags_estimates["beta_eps[9]", c("mean", "2.5%", "97.5%")],
#                             jags_estimates["beta_eps[10]", c("mean", "2.5%", "97.5%")],
#                             jags_estimates["beta_eps[11]", c("mean", "2.5%", "97.5%")],
#                             jags_estimates["beta_eps[12]", c("mean", "2.5%", "97.5%")],
#                             jags_estimates["beta_eps[13]", c("mean", "2.5%", "97.5%")],
#                             jags_estimates["beta_eps[14]", c("mean", "2.5%", "97.5%")],
#                             jags_estimates["beta_eps[15]", c("mean", "2.5%", "97.5%")],
#                             jags_estimates["beta_eps[16]", c("mean", "2.5%", "97.5%")])
#   
#   comparison <- cbind(as.data.frame(stan_estimates)[,c(1,3,4)], jags_reformatted)
#   colnames(comparison) <- c("Stan_mean", "Stan_2.5%", "Stan_97.5%",
#                             "JAGS_mean", "JAGS_2.5%", "JAGS_97.5%")
#   comparison <- comparison[,c(1,4,2,5,3,6)]
#   
#   comparison$Stan_sign <- ifelse(comparison$`Stan_2.5%` < 0 & comparison$`Stan_97.5%` > 0, 0, 1)
#   comparison$Jags_sign <- ifelse(comparison$`JAGS_2.5%` < 0 & comparison$`JAGS_97.5%` > 0, 0, 1)
#   
#   return(comparison)
#   
# }


# # function to compare output for full model:
# 
# compare_Jags_Stan_fm <- function(jags_output, stan_output){
#   
#   stan_estimates <- fixef(stan_output, summary = TRUE, robust = FALSE, probs = c(0.025, 0.975))
#   jags_estimates <- jags_output$summary[which(grepl("(alpha)|(beta)|(p)", rownames(jags_output$summary))), c("mean", "2.5%", "97.5%")]
#   
#   jags_reformatted <- rbind(c(NA,NA, NA),
#                             jags_estimates["alpha_psi", c("mean", "2.5%", "97.5%")],
#                             jags_estimates["alpha_gamma", c("mean", "2.5%", "97.5%")],
#                             jags_estimates["alpha_eps", c("mean", "2.5%", "97.5%")],
#                             jags_estimates["p[1]", c("mean", "2.5%", "97.5%")],
#                             jags_estimates["p[2]", c("mean", "2.5%", "97.5%")],
#                             jags_estimates["p[3]", c("mean", "2.5%", "97.5%")],
#                             jags_estimates["p[4]", c("mean", "2.5%", "97.5%")],
#                             jags_estimates["p[5]", c("mean", "2.5%", "97.5%")],
#                             jags_estimates["beta_psi[1]", c("mean", "2.5%", "97.5%")],
#                             jags_estimates["beta_psi[2]", c("mean", "2.5%", "97.5%")],
#                             jags_estimates["beta_psi[3]", c("mean", "2.5%", "97.5%")],
#                             jags_estimates["beta_psi[4]", c("mean", "2.5%", "97.5%")],
#                             jags_estimates["beta_psi[5]", c("mean", "2.5%", "97.5%")],
#                             jags_estimates["beta_psi[6]", c("mean", "2.5%", "97.5%")],
#                             jags_estimates["beta_psi[7]", c("mean", "2.5%", "97.5%")],
#                             jags_estimates["beta_psi[8]", c("mean", "2.5%", "97.5%")],
#                             jags_estimates["beta_psi[9]", c("mean", "2.5%", "97.5%")],
#                             jags_estimates["beta_psi[10]", c("mean", "2.5%", "97.5%")],
#                             jags_estimates["beta_psi[11]", c("mean", "2.5%", "97.5%")],
#                             jags_estimates["beta_psi[12]", c("mean", "2.5%", "97.5%")],
#                             jags_estimates["beta_psi[13]", c("mean", "2.5%", "97.5%")],
#                             jags_estimates["beta_psi[14]", c("mean", "2.5%", "97.5%")],
#                             jags_estimates["beta_psi[15]", c("mean", "2.5%", "97.5%")],
#                             jags_estimates["beta_psi[16]", c("mean", "2.5%", "97.5%")],
#                             jags_estimates["beta_psi[17]", c("mean", "2.5%", "97.5%")],
#                             jags_estimates["beta_psi[19]", c("mean", "2.5%", "97.5%")],
#                             jags_estimates["beta_psi[20]", c("mean", "2.5%", "97.5%")],
#                             jags_estimates["beta_psi[21]", c("mean", "2.5%", "97.5%")],
#                             jags_estimates["beta_psi[18]", c("mean", "2.5%", "97.5%")],
#                             jags_estimates["beta_psi[22]", c("mean", "2.5%", "97.5%")],
#                             jags_estimates["beta_psi[24]", c("mean", "2.5%", "97.5%")],
#                             jags_estimates["beta_psi[25]", c("mean", "2.5%", "97.5%")],
#                             jags_estimates["beta_psi[26]", c("mean", "2.5%", "97.5%")],
#                             jags_estimates["beta_psi[23]", c("mean", "2.5%", "97.5%")],
#                             jags_estimates["beta_gamma[1]", c("mean", "2.5%", "97.5%")],
#                             jags_estimates["beta_gamma[2]", c("mean", "2.5%", "97.5%")],
#                             jags_estimates["beta_gamma[3]", c("mean", "2.5%", "97.5%")],
#                             jags_estimates["beta_gamma[4]", c("mean", "2.5%", "97.5%")],
#                             jags_estimates["beta_gamma[5]", c("mean", "2.5%", "97.5%")],
#                             jags_estimates["beta_gamma[6]", c("mean", "2.5%", "97.5%")],
#                             jags_estimates["beta_gamma[7]", c("mean", "2.5%", "97.5%")],
#                             jags_estimates["beta_gamma[8]", c("mean", "2.5%", "97.5%")],
#                             jags_estimates["beta_gamma[9]", c("mean", "2.5%", "97.5%")],
#                             jags_estimates["beta_gamma[10]", c("mean", "2.5%", "97.5%")],
#                             jags_estimates["beta_gamma[11]", c("mean", "2.5%", "97.5%")],
#                             jags_estimates["beta_gamma[12]", c("mean", "2.5%", "97.5%")],
#                             jags_estimates["beta_gamma[13]", c("mean", "2.5%", "97.5%")],
#                             jags_estimates["beta_gamma[14]", c("mean", "2.5%", "97.5%")],
#                             jags_estimates["beta_gamma[15]", c("mean", "2.5%", "97.5%")],
#                             jags_estimates["beta_gamma[16]", c("mean", "2.5%", "97.5%")],
#                             jags_estimates["beta_gamma[17]", c("mean", "2.5%", "97.5%")],
#                             jags_estimates["beta_gamma[19]", c("mean", "2.5%", "97.5%")],
#                             jags_estimates["beta_gamma[20]", c("mean", "2.5%", "97.5%")],
#                             jags_estimates["beta_gamma[21]", c("mean", "2.5%", "97.5%")],
#                             jags_estimates["beta_gamma[18]", c("mean", "2.5%", "97.5%")],
#                             jags_estimates["beta_gamma[22]", c("mean", "2.5%", "97.5%")],
#                             jags_estimates["beta_gamma[24]", c("mean", "2.5%", "97.5%")],
#                             jags_estimates["beta_gamma[25]", c("mean", "2.5%", "97.5%")],
#                             jags_estimates["beta_gamma[26]", c("mean", "2.5%", "97.5%")],
#                             jags_estimates["beta_gamma[23]", c("mean", "2.5%", "97.5%")],
#                             jags_estimates["beta_eps[1]", c("mean", "2.5%", "97.5%")],
#                             jags_estimates["beta_eps[2]", c("mean", "2.5%", "97.5%")],
#                             jags_estimates["beta_eps[3]", c("mean", "2.5%", "97.5%")],
#                             jags_estimates["beta_eps[4]", c("mean", "2.5%", "97.5%")],
#                             jags_estimates["beta_eps[5]", c("mean", "2.5%", "97.5%")],
#                             jags_estimates["beta_eps[6]", c("mean", "2.5%", "97.5%")],
#                             jags_estimates["beta_eps[7]", c("mean", "2.5%", "97.5%")],
#                             jags_estimates["beta_eps[8]", c("mean", "2.5%", "97.5%")],
#                             jags_estimates["beta_eps[9]", c("mean", "2.5%", "97.5%")],
#                             jags_estimates["beta_eps[10]", c("mean", "2.5%", "97.5%")],
#                             jags_estimates["beta_eps[11]", c("mean", "2.5%", "97.5%")],
#                             jags_estimates["beta_eps[12]", c("mean", "2.5%", "97.5%")],
#                             jags_estimates["beta_eps[13]", c("mean", "2.5%", "97.5%")],
#                             jags_estimates["beta_eps[14]", c("mean", "2.5%", "97.5%")],
#                             jags_estimates["beta_eps[15]", c("mean", "2.5%", "97.5%")],
#                             jags_estimates["beta_eps[16]", c("mean", "2.5%", "97.5%")],
#                             jags_estimates["beta_eps[17]", c("mean", "2.5%", "97.5%")],
#                             jags_estimates["beta_eps[19]", c("mean", "2.5%", "97.5%")],
#                             jags_estimates["beta_eps[20]", c("mean", "2.5%", "97.5%")],
#                             jags_estimates["beta_eps[21]", c("mean", "2.5%", "97.5%")],
#                             jags_estimates["beta_eps[18]", c("mean", "2.5%", "97.5%")],
#                             jags_estimates["beta_eps[22]", c("mean", "2.5%", "97.5%")],
#                             jags_estimates["beta_eps[24]", c("mean", "2.5%", "97.5%")],
#                             jags_estimates["beta_eps[25]", c("mean", "2.5%", "97.5%")],
#                             jags_estimates["beta_eps[26]", c("mean", "2.5%", "97.5%")],
#                             jags_estimates["beta_eps[23]", c("mean", "2.5%", "97.5%")])
#   
#   stan_estimates2 <- rbind(stan_estimates[1:4,c(1,3,4)], c(NA, NA, NA), stan_estimates[5:nrow(stan_estimates),c(1,3,4)])
#   row.names(stan_estimates2)[5] <- "routeSectionSect1"
#   comparison <- as.data.frame(cbind(stan_estimates2, jags_reformatted))
#   colnames(comparison) <- c("Stan_mean", "Stan_Q2.5", "Stan_Q97.5",
#                             "JAGS_mean", "JAGS_Q2.5", "JAGS_Q97.5")
#   comparison <- comparison[,c(1,4,2,5,3,6)]
#   
#   comparison$Stan_sign <- ifelse(comparison$Stan_Q2.5 < 0 & comparison$Stan_Q97.5 > 0, 0, 1)
#   comparison$Jags_sign <- ifelse(comparison$JAGS_Q2.5 < 0 & comparison$JAGS_Q97.5 > 0, 0, 1)
#   
#   return(comparison)
#   
# }


# AUC per year: ----

# same as later in Briscoe modified function!

years <- unique(AUC_dt3$year)
yearly_AUC <- vector(mode = "numeric", length = length(years))

for(t in 1:length(years)){
  
  year_dt <- AUC_dt3 %>% 
    filter(year == years[t])
  sscurves <- evalmod(scores = year_dt$prob_occ, labels = year_dt$obs_occ, posclass = 1)
  yearly_AUC[t] <- auc(sscurves)$aucs[1]
  
}
yearly_AUC

# plot AUC over time:
ggplot(data = data.frame("year" = years, "yearly_AUC" = yearly_AUC),
       aes(x = year, y = yearly_AUC)) +
  geom_line() +
  geom_point() +
  ylim(c(0.5, 1)) +
  theme_bw() +
  ggtitle(spec) +
  xlab("Year") + ylab("AUC") +
  theme(text = element_text(size=15))


# Briscoe: ----
# Function to calculate deviance (use as is!)
calc.deviance<- function (obs, pred, weights = rep(1, length(obs)), family = "binomial", 
                          calc.mean = FALSE) 
{
  if (length(obs) != length(pred)) {
    stop("observations and predictions must be of equal length")
  }
  y_i <- obs
  u_i <- pred
  family = tolower(family)
  if (family == "binomial" | family == "bernoulli") {
    deviance.contribs <- (y_i * log(u_i)) + ((1 - y_i) * 
                                               log(1 - u_i))
    
    deviance <- -2 * sum(deviance.contribs * weights)
  }
  else if (family == "poisson") {
    deviance.contribs <- ifelse(y_i == 0, 0, (y_i * log(y_i/u_i))) - 
      (y_i - u_i)
    deviance <- 2 * sum(deviance.contribs * weights)
  }
  else if (family == "laplace") {
    deviance <- sum(abs(y_i - u_i))
  }
  else if (family == "gaussian") {
    deviance <- sum((y_i - u_i) * (y_i - u_i))
  }
  else {
    stop("unknown family, should be one of: \"binomial\", \"bernoulli\", \"poisson\", \"laplace\", \"gaussian\"")
  }
  if (calc.mean) 
    deviance <- deviance/length(obs)
  return(deviance)
}


# Function to get evaluation statistics **individually for all years**.

# eval.f() of briscoe changed to my data:

x = 2000

AUC_dt3

eval.f.2 <-function(x) { # x = year, 1991:2015
  
  # eval.yr<-pred.st[[x-1999]] # pred.st - raster stack of predictions ([space, time]?)
  # spObs.yr<-spObs_sdf[spObs_sdf$year==x & !is.na(spObs_sdf$collDet),] #spObs - SpatialPointsDataFrame with species observations; colldet = data aggregated across all survey visits;whether the species was ever detected at the site in each year)
  # eval.pred<-raster::extract(eval.yr,spObs.yr) # occ. prob. prediction at species occurrence point
  
  dt_yr <- AUC_dt3[which(AUC_dt3$year == x),] # only complete cases
  eval.pred <- dt_yr$prob_occ
  eval.pred[eval.pred==0] <- 0.0000000000000001 #Replace 0s with very small number to avoid NaN errors
  eval.pred[eval.pred==1] <- (1-0.0000000000000001) #Replace 0s with 1-very small number to avoid NaN errors

  # observations for one year:
  # spObs.yr<-spObs.yr$collDet #colldet = data aggregated across all survey visits (1 if species was observed in this site and year at least once)
  spObs.yr <- as.numeric(as.character(dt_yr$obs_occ))
  
  eval_obj <- precrec::evalmod(scores = eval.pred, labels = spObs.yr) #ROC and Precision-Recall curves
  
  prev <- sum(spObs.yr)/length(spObs.yr) # prevalence, sites at which species was observed / number of sites?
  
  dev_mod <- calc.deviance(obs=spObs.yr, 
                           pred=eval.pred, 
                           weights = rep(1,length(spObs.yr)), 
                           family="binomial", 
                           calc.mean = FALSE) # deviance explained by the model
  dev_null <- calc.deviance(obs=spObs.yr,
                            pred=rep(prev,length(spObs.yr)), # same prevalence for each site
                            weights = rep(1,length(spObs.yr)), 
                            family="binomial", 
                            calc.mean = FALSE) # null model that estimates species prevalence in each year (i.e. the proportion of sites where the species was detected)
  # % explained deviance:
  dev_exp<-((dev_null - dev_mod)/dev_null) *100
  # constrain to 0 to 100:
  if(dev_exp < 0) dev_exp <- 0
  
  eval.out<-c(auc(eval_obj)[1,4], # area under ROC curve
              auc(eval_obj)[2,4], # area under precision-recall curve
              dev_exp) # explained deviance
  
  return(eval.out)
}

## yeary ROC-AUC, PR-AUC, explained deviance: ----

years <- sort(unique(AUC_dt3$year))

res_yearly_df <- data.frame("year" = years,
                            "AUC_ROC" = rep(NA, length(years)), 
                            "AUC_PR" = rep(NA, length(years)), 
                            "dev_expl" = rep(NA, length(years)))
for(y in 1:length(years)){
  print(years[y])
  res_yearly_df[y, c(2:4)] <- eval.f.2(x = years[y])
}
res_yearly_df

ggplot(data = res_yearly_df,
       aes(x = year, y = AUC_ROC)) +
  geom_line() +
  geom_point() +
  ylim(c(0.5, 1)) +
  theme_bw() +
  ggtitle(spec) +
  xlab("Year") + ylab("AUC ROC") +
  theme(text = element_text(size=15))
ggplot(data = res_yearly_df,
       aes(x = year, y = AUC_PR)) +
  geom_line() +
  geom_point() +
  ylim(c(0.5, 1)) +
  theme_bw() +
  ggtitle(spec) +
  xlab("Year") + ylab("AUC PR") +
  theme(text = element_text(size=15))
ggplot(data = res_yearly_df,
       aes(x = year, y = dev_expl)) +
  geom_line() +
  geom_point() +
  ylim(c(0, max(res_yearly_df$dev_expl))) +
  theme_bw() +
  ggtitle(spec) +
  xlab("Year") + ylab("% explained deviance") +
  theme(text = element_text(size=15))

# but this was also for single years?

# continue here: xx 
# Function to get evaluation statistics pooled across all years#####
#Assumes have the following loaded into the workspace:
# pred.st - raster stack of predictions
# spObs - SpatialPointsDataFrame with species observations 

eval_all_years.f<-function(y){
  
  years<-y
  eval.pred.all<-NULL
  spObs.all<-NULL
  prev.all<-NULL
  
  for(x in years){
    # print(x)
    eval.yr<-pred.st[[x-1999]]
    spObs.yr<-spObs[spObs$year==x & !is.na(spObs$collDet),]
    
    #add check so skipped if no obs data for that species in that year
    if(length(spObs.yr@data[,1])>0){
      eval.pred<-raster::extract(eval.yr,spObs.yr)
      
      #Add check to get rid of any NA predictions (due to missing NDVI?)
      ex<-which(is.na(eval.pred==TRUE))
      if(length(ex)!=0){
        eval.pred<-eval.pred[-ex]
      }
      
      eval.pred[eval.pred==0]<-0.0000000000000001 #Replace 0s with very small number to avoid NaN errors
      eval.pred[eval.pred==1]<-1-0.0000000000000001 #Replace 0s with 1-very small number to avoid NaN errors
      spObs.yr<-spObs.yr$collDet
      if(length(ex)!=0){
        spObs.yr<-spObs.yr[-ex]
      }
      prev<-rep(sum(spObs.yr)/length(spObs.yr),length(spObs.yr))
      
      eval.pred.all<-c(eval.pred.all,eval.pred)
      prev.all<-c(prev.all,prev)
      spObs.all<-c(spObs.all, spObs.yr)
    }
  }
  if(length(unique(spObs.all))==1){ #All 1s or all 0s, can't assess
    eval.out<-c(NA,NA,NA)
  }else{
    eval_obj <- precrec::evalmod(scores = eval.pred.all, labels = spObs.all)
    #  prev<-sum(spObs.all)/length(spObs.all)
    dev_mod<-calc.deviance(obs=spObs.all, pred=eval.pred.all, weights = rep(1,length(spObs.all)),
                           family="binomial", calc.mean = FALSE)
    dev_null<-calc.deviance(obs=spObs.all, pred=prev.all, weights = rep(1,length(spObs.all)),
                            family="binomial", calc.mean = FALSE)
    
    dev_exp<-((dev_null - dev_mod)/dev_null) *100
    #roc <- rocprc(eval_obj)
    eval.out<-c(auc(eval_obj)[1,4],auc(eval_obj)[2,4],dev_exp)
  }
  return(eval.out)
}






# temporal predictive ability: ----

# Briscoe et al.
# compare expected area of occupancy to number of sites observed to be occupied (observed AOO)
# calculate C-index values (Hmisc::rcorr.cens) [[C-index]]

##Calculate temporal AUC


temp.cor.out<-as.data.frame(matrix(ncol=13,nrow=5))
temp.cor.out.all<-temp.cor.out[0,]

data.yrs<-c(3,5,10)

for(s in 1:length(species)){
  speciesName<-species[s]
  
  for(dat.yrs in data.yrs){
    
    for(m in 1:length(methods)){
      method<-methods1[m]
      
      pred.dat<-SumPrOcc_long[SumPrOcc_long$species==speciesName 
                              & SumPrOcc_long$data_yrs==dat.yrs & 
                                SumPrOcc_long$method==method,]
      obs.dat<-mhb.trends[mhb.trends$CommonName==speciesName & 
                            mhb.trends$Time %in% pred.dat$year,]
      obs.dat<-obs.dat[order(obs.dat$Time),]
      pred.dat<-pred.dat[order(pred.dat$year),]
      
      pred.dat$SumPrOcc_st<-pred.dat$SumPrOcc/pred.dat$SumPrOcc[1] #
      
      
      out<-Hmisc::rcorr.cens(pred.dat$SumPrOcc_st, obs.dat$Imp.) # c-index 
      
      
      out.cor<-cor(pred.dat$SumPrOcc_st, obs.dat$Imp.)
      
      temp.cor.out[m,1]<-speciesName
      temp.cor.out[m,2]<-as.character(method)
      temp.cor.out[m,3]<-as.character(dat.yrs)
      temp.cor.out[m,4:12]<-out
      temp.cor.out[m,13]<-out.cor
    }
    
    temp.cor.out.all<-rbind(temp.cor.out.all,temp.cor.out)
  }
}

####
# KS:

# expected area of occupancy:
# =  yearly sums of estimated probability of occupancy across the landscape
exp_AOO <- AUC_dt2 %>% 
  group_by(year) %>% 
  summarise(sum_prob_occ = sum(prob_occ)) # rather with predictions to all US cells xx!?

# number of sites observed to be occupied:
obs_AOO <- AUC_dt2 %>% 
  group_by(year) %>% 
  summarise(sum_obs_occ = sum(as.numeric(as.character(obs_occ)), na.rm = TRUE))

Hmisc::rcorr.cens(exp_AOO$sum_prob_occ, obs_AOO$sum_obs_occ) # ??
x <- round(rnorm(200))
y <- rnorm(200)
c_index <- Hmisc::rcorr.cens(x, y, outx=F)[1]
c_index

ggplot(data = cbind(obs_AOO, exp_AOO[,2]),
       aes(x = year)) +
  geom_line(aes(y = sum_obs_occ, colour = "blue")) +
  geom_point(aes(y = sum_obs_occ), colour = "blue") +
  geom_line(aes(y = sum_prob_occ, colour = "red")) +
  geom_point(aes(y = sum_prob_occ), colour = "red") + 
  scale_color_identity(guide = "legend", labels = c("observed AOO", "expected AOO"), name = "") +
  theme_bw() +
  ggtitle(spec, paste("C-index =", round(c_index, 2))) +
  xlab("Year") + ylab("AOO") +
  theme(text = element_text(size=15), legend.position="bottom")


# same for independent test routes