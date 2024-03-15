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

# selected routes and years:
load(file = file.path("data", "route_year_env_data.RData"))
route_sel_env_dt_final



# output of model fitting: ----

load(file.path("results", "Grasshopper_Sparrow_cl_alt3_quadr_lasso_GoF_inits.RData"))
BC <- out
summary(BC)
pars_no_conv <- names(which(BC$summary[,"Rhat"] > 1.1))
traceplot(BC, parameters = "beta_eps", Rhat_min =  1.1)
traceplot(BC, parameters = "beta_gamma", Rhat_min =  1.1)
sort(unique(as.numeric(str_extract(pars_no_conv, "((?<=\\[).*(?=,))|((?<=\\[).*(?=\\]))"))))
unique(str_extract(pars_no_conv, ".*(?=\\[)"))
traceplot(BC, parameters = "beta_gamma", Rhat_min =  1.1)


load(file.path("results", "Eastern_Phoebe_cl_alt3_quadr_lasso_GoF_inits.RData"))  # model output
EP_inits <- out
summary(EP_inits)
EP_inits$summary %>% View
traceplot(EP_inits, parameters = "p")
which(unlist(EP_inits$Rhat) > 1.1)
traceplot(EP_inits, parameters = "psi")
traceplot(EP_inits, Rhat_min =  1.1)
traceplot(EP_inits, parameters = "gamma")
densityplot(EP_inits, parameters = "gamma")
traceplot(EP_inits, parameters = "alpha_gamma")
densityplot(EP_inits, parameters = "alpha_gamma")
traceplot(EP_inits, parameters = "beta_gamma")
densityplot(EP_inits, parameters = "beta_gamma")
densityplot(EP_inits, parameters = "psi1")
# get RTENOs of routes with convergence problems:

pars_no_conv <- names(which(EP_inits$summary[,"Rhat"] > 1.1))
nsites = 476
nyears = 25
route_nrs <- matrix(route_sel_env_dt_final$RTENO, nrow = nsites, ncol = nyears, byrow = TRUE)[,1]
route_nrs


# load(file.path("results", "House_Finch_cl_lu_noC3C4_quadr_p_sect_lasso.RData"))  # model output
# HF_cl_lu_quadr_noC3C4_p_sect <- out
routes_no_conv <- unique(as.numeric(str_extract(pars_no_conv, "((?<=\\[).*(?=,))|((?<=\\[).*(?=\\]))"))) # (?<=\\[) looks for what's behind [, (?=\\])) looks for what's before ]
routes_no_conv
RTENOs_no_conv <- route_nrs[routes_no_conv]
# load routes:
library(sf)
routes_sel_sf <- st_read(file.path("data", "route_selection_1991_2015_surv_beg_end_max_5y_miss_v2_spat_thin_100km_max_30_r_per_BCR_centroids.shp")) # output of 1_1_route_selection.R
bio_sf <- routes_sel_sf %>% 
  left_join(route_sel_env_dt_final[, c("RTENO", "bio15_3yrs")], by = c("RTENO_BBS" = "RTENO")) %>% 
  distinct()
plot(bio_sf["bio15_3yrs"],reset = FALSE)

routes_sel_sf %>% 
  filter(RTENO_BBS %in% RTENOs_no_conv) %>% 
  plot(add = TRUE, col = "black", reset = FALSE, pch = 4)



summary(out)
length(which(unlist(HF_quad_res$Rhat) > 1.1)) # 137
str(HF_quad_res)
HF_quad_res$mcmc.info$n.samples
HF_quad_res$mcmc.info$n.iter

load(file.path("results", "House_Finch_quadr_GoF_update.RData")) 
HF_quad_res_update <- out2
summary(HF_quad_res_update) # convergence failure -> 3000 iterations not enough
length(which(unlist(HF_quad_res_update$Rhat) > 1.1)) # 157

load(file.path("results", "House_Finch_quadr_GoF_update5000.RData")) 
HF_quad_res_update2 <- out3
summary(HF_quad_res_update2) # convergence failure -> 5000 iterations not enough
length(which(unlist(HF_quad_res_update2$Rhat) > 1.1)) # 73

load(file.path("results", "House_Finch_quadr_GoF_p_site_re_lasso.RData")) 
HF_quad_res_6000<- out
summary(HF_quad_res_6000) 
length(which(unlist(HF_quad_res_6000$Rhat) > 1.1)) # 866

sort(names(unlist(HF_quad_res_6000$Rhat)[which(unlist(HF_quad_res_6000$Rhat) > 1.1)]), decreasing = TRUE)

# other species:
load(file.path("results", "Eastern_Meadowlark_quadr_GoF_p_site_re_lasso.RData")) 
EM_quad_res_6000 <- out
summary(EM_quad_res_6000) 
length(which(unlist(EM_quad_res_6000$Rhat) > 1.1)) # 1707

load(file.path("results", "Winter_Wren_quadr_GoF_p_site_re_lasso.RData")) 
WW_quad_res_6000 <- out
summary(WW_quad_res_6000) 
length(which(unlist(WW_quad_res_6000$Rhat) > 1.1)) # 15839

#load(file.path("results", "Brown_headed_Cowbird_test1.RData"))

HF_quad_res_update2$sims.list$Chi2ratioOpen
HF_quad_res_update2$sims.list$Chi2ratioClosed
mean(HF_quad_res_update2$sims.list$Chi2ratioOpen)
mean(HF_quad_res_update2$sims.list$Chi2ratioClosed)

load(file.path("results", "Yellow_Warbler_cl_alt2_quadr_lasso_GoF.RData")) 
YW_bioalt2 <- out



# explore output files: ----
summary(HF_cl_alt_quadr_lasso) # all iterations before the final incremental step are considered burnin
# 2500 iteration in total (1500+1000)
# 1000 iteration not discarded
# thin rate = 5 -> 200 draws from the posterior
summary(HF_cl_quadr)

print(HF_cl_alt_quadr_lasso)
HF_cl_alt_quadr_lasso$parameters
HF_cl_alt_quadr_lasso$mcmc.info
HF_cl_alt_quadr_lasso$summary %>% View

HF_cl_alt_quadr_lasso$summary[grep("beta_psi", rownames(HF_cl_alt_quadr_lasso$summary)),]

traceplot(x = HF_cl_alt_quadr_lasso, parameters = "psi1")
traceplot(x = HF_cl_alt_quadr_lasso, parameters = "alpha_psi")
densityplot(HF_cl_alt_quadr_lasso, parameters = "alpha_psi")
HF_cl_alt_quadr_lasso$summary[grep("alpha_psi", rownames(HF_cl_alt_quadr_lasso$summary)),]
densityplot(HF_cl_alt_quadr_lasso, parameters = "beta_psi")
HF_cl_alt_quadr_lasso$model
HF_cl_alt_quadr_lasso$sims.list$p # actual draws from posterior for p
hist(HF_cl_alt_quadr_lasso$sims.list$p, xlab="Value", main = "p posterior")
densityplot(HF_cl_alt_quadr_lasso, parameters = "p")
str(HF_cl_alt_quadr_lasso)

whiskerplot(HF_cl_alt_quadr_lasso, parameters = "beta_psi")
whiskerplot(HF_cl_alt_quadr_lasso, parameters = "beta_eps")
whiskerplot(HF_cl_alt_quadr_lasso, parameters = "beta_gamma")

var(HF_cl_alt_quadr_lasso$sims.list$deviance)/2 # roughly the same as pD (6340 instead of 6356)
str(HF_cl_alt_quadr_lasso$samples) # samples per chain; original output object from the rjags package, as class mcmc.list
str(HF_cl_alt_quadr_lasso$samples[[1]])
HF_cl_alt_quadr_lasso$samples[[1]][1,1:5]
HF_cl_alt_quadr_lasso$sims.list$psi1[1, 1:5] # same

# explore goodness-of-fit following Kéry and Royle AHM 2021: -------------------

# Plots of expected versus observed value of fit stats
# Open part
res <- YW_bioalt2#HF_cl_alt_quadr_lasso#HF_cl_quadr_unif#HF_cl_quadr_p_sect#HF_cl_quadr#HF_lu_quadr_noC3C4_p_sect#HF_cl_quadr_p_sect#HF_cl_lu_quadr_noC3C4_p_sect#HF_cl_lu_quadr_noC3C4#HF_nm#EM_quad_res_6000# HF_quad_res_6000 #HF_quad_res # HF_lin_res

pl <- range(c(res$sims.list$Chi2Open, res$sims.list$Chi2repOpen))
plot(res$sims.list$Chi2Open, res$sims.list$Chi2repOpen,
     xlab = "Chi2 observed data", ylab = "Chi2 expected data",
     main = "Open part of model ", xlim = pl, ylim = pl, frame.plot = FALSE)
abline(0, 1, lwd = 2)
text(200, 500, paste('Bpv = ', round(mean(res$sims.list$Chi2repOpen >
                                            res$sims.list$Chi2Open), 2)), cex = 1) # Bayesian p-value 0 proportion of points above the line

# Closed part of model: Chi-squared
pl <- range(c(res$sims.list$Chi2Closed, res$sims.list$Chi2repClosed))
plot(res$sims.list$Chi2Closed, res$sims.list$Chi2repClosed,
     xlab = "Chi2 observed data", ylab = "Chi2 expected data",
     main = "Closed part of model (Chi-squared)", xlim = pl, ylim = pl,
     frame.plot = FALSE)
abline(0, 1, lwd = 2)
text(3500, 5000, paste('Bpv = ', round(mean(res$sims.list$Chi2repClosed >
                                              res$sims.list$Chi2Closed), 2)), cex = 1)





#library(jagshelper)
#check_Rhat(HF_quad_res_update2) # proportion of Rhats below a threshold of 1.1
#traceworstRhat(HF_quad_res_update2, parmfrow=c(3,2))  # trace plots for least-converged nodes
#traceworstRhat(HF_quad_res_update2)[1]  # trace plots for least-converged nodes
#out_df <- jags_df(HF_quad_res_update2)
#str(out_df)



# fit for open part of the model looks better when including quadratic terms (but not converged yet!)


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