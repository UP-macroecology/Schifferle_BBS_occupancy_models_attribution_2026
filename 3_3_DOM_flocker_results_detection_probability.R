# gather estimated detection probability for each species
# plot detection probability posterior distributions for each species
# explore relationship of detection probability and traits and model performance

# packages:

library(dplyr)
library(flocker)
library(cmdstanr)
#set_cmdstan_path(path = NULL)
set_cmdstan_path("C:/Users/schifferle1/Documents/cmdstan-2.34.1")
library(ggplot2)

# directory of results:
#results_dir <- file.path("/mnt", "ibb_share", "zurell_transfer", "Schifferle_BBS_occupancy_models_2023", "results", "fm_buffer750km")
results_dir <- file.path("T:", "Schifferle_BBS_occupancy_models_2023", "results", "fm_buffer750km")

# directory to save plots:
plot_dir <- file.path("plots", "detection_probability")
if(!dir.exists(plot_dir)){dir.create(plot_dir, recursive = TRUE)}

# selected species:
load(file = file.path("data", "species_set_analysis.RData"))
final_species

# save detection probability:
det_prob_df <- data.frame(species = final_species, 
                          det_prob_median_mean = NA, 
                          CI90_low_mean = NA, 
                          CI90_high_mean = NA,
                          det_prob_median_sect1 = NA,
                          det_prob_CI90_low_sect1 = NA,
                          det_prob_CI90_high_sect1 = NA,
                          det_prob_median_sect2 = NA,
                          det_prob_CI90_low_sect2 = NA,
                          det_prob_CI90_high_sect2 = NA,
                          det_prob_median_sect3 = NA,
                          det_prob_CI90_low_sect3 = NA,
                          det_prob_CI90_high_sect3 = NA,
                          det_prob_median_sect4 = NA,
                          det_prob_CI90_low_sect4 = NA,
                          det_prob_CI90_high_sect4 = NA,
                          det_prob_median_sect5 = NA,
                          det_prob_CI90_low_sect5 = NA,
                          det_prob_CI90_high_sect5 = NA)

# extract and store detection probability: ----

for(i in 1:length(final_species)){
  
  spec <- final_species[i]
  
  print(paste(i, spec))
      
  # check where to look for model output (did MCMC fitting work with less or only with more iterations?)
  if(file.exists(file.path(results_dir, "refit_2000_2000", paste0("out_", spec, "_fm_buffer750.RData")))){
    output_dir <- file.path(results_dir, "refit_2000_2000")
  } else {
    output_dir <- results_dir
  }
  
  # model predictions:
  load(file.path(output_dir, paste0("postproc_", spec, "_fm_buffer750.RData")))

  # plot posterior distribution:
  
  # reformat for plotting:
  det_prob_samples <- res_list$fitted$linpred_det[1,,1,]
  rm(res_list)
  det_prob_samples2 <- t(det_prob_samples)

  plot_title <- ggtitle(paste(spec, "\nDetection probability (median + 90% intervals)"))
  p <- bayesplot::mcmc_areas(det_prob_samples2,prob = 0.9) +
    plot_title +
    xlim(0, 1) +
    geom_vline(xintercept = 0.5, linetype = "dashed", colour = "grey60") +
    geom_vline(xintercept = 1/3, linetype = "dashed", colour = "grey60") +
    geom_vline(xintercept = 2/3, linetype = "dashed", colour = "grey60") +
    theme(text = element_text(size = 22))
  
  jpeg(file = file.path(plot_dir, paste0(spec,"_detection_prob_fm.jpg")), 
       width = 1000, height = 700, quality = 100)
  print(p)
  dev.off()
  
  # median:
  det_prob_route <- apply(det_prob_samples, MAR = 1, FUN = median)

  det_prob_df$det_prob_median_sect1[which(det_prob_df$species == spec)] <- det_prob_route[1]
  det_prob_df$det_prob_median_sect2[which(det_prob_df$species == spec)] <- det_prob_route[2]
  det_prob_df$det_prob_median_sect3[which(det_prob_df$species == spec)] <- det_prob_route[3]
  det_prob_df$det_prob_median_sect4[which(det_prob_df$species == spec)] <- det_prob_route[4]
  det_prob_df$det_prob_median_sect5[which(det_prob_df$species == spec)] <- det_prob_route[5]
  
  # 90 % CI (percentile interval):
  ci_visit1 <- bayestestR::ci(det_prob_samples[1,], ci = 0.9, method = "ETI")
  ci_visit2 <- bayestestR::ci(det_prob_samples[2,], ci = 0.9, method = "ETI")
  ci_visit3 <- bayestestR::ci(det_prob_samples[3,], ci = 0.9, method = "ETI")
  ci_visit4 <- bayestestR::ci(det_prob_samples[4,], ci = 0.9, method = "ETI")
  ci_visit5 <- bayestestR::ci(det_prob_samples[5,], ci = 0.9, method = "ETI")

  det_prob_df$det_prob_CI90_low_sect1[which(det_prob_df$species == spec)] <- ci_visit1$CI_low
  det_prob_df$det_prob_CI90_high_sect1[which(det_prob_df$species == spec)] <- ci_visit1$CI_high
  det_prob_df$det_prob_CI90_low_sect2[which(det_prob_df$species == spec)] <- ci_visit2$CI_low
  det_prob_df$det_prob_CI90_high_sect2[which(det_prob_df$species == spec)] <- ci_visit2$CI_high
  det_prob_df$det_prob_CI90_low_sect3[which(det_prob_df$species == spec)] <- ci_visit3$CI_low
  det_prob_df$det_prob_CI90_high_sect3[which(det_prob_df$species == spec)] <- ci_visit3$CI_high
  det_prob_df$det_prob_CI90_low_sect4[which(det_prob_df$species == spec)] <- ci_visit4$CI_low
  det_prob_df$det_prob_CI90_high_sect4[which(det_prob_df$species == spec)] <- ci_visit4$CI_high
  det_prob_df$det_prob_CI90_low_sect5[which(det_prob_df$species == spec)] <- ci_visit5$CI_low
  det_prob_df$det_prob_CI90_high_sect5[which(det_prob_df$species == spec)] <- ci_visit5$CI_high
  
  # mean across route sections for median:
  det_prob_df$det_prob_median_mean[which(det_prob_df$species == spec)] <- mean(det_prob_route)
  
  # mean across route sections for lower bound:
  mean_CI_low <- mean(ci_visit1$CI_low, ci_visit2$CI_low, ci_visit3$CI_low, ci_visit4$CI_low, ci_visit5$CI_low)
  mean_CI_high <- mean(ci_visit1$CI_high, ci_visit2$CI_high, ci_visit3$CI_high, ci_visit4$CI_high, ci_visit5$CI_high)
  
  det_prob_df$CI90_low_mean[which(det_prob_df$species == spec)] <- mean_CI_low
  det_prob_df$CI90_high_mean[which(det_prob_df$species == spec)] <- mean_CI_high
  }

#write.csv(det_prob_df, file.path("data", "detection_probability_fm.csv"))
summary(det_prob_df)


# explore detection probability and traits: ------------------------------------

load(file.path("data", "BBS_data_merged.RData"))

det_prob_traits <- bbs_dt %>% 
  select(English_Common_Name:Range.Size) %>% 
  distinct() %>% 
  right_join(det_prob_df, by = c(English_Common_Name = "species")) %>% 
  rename("species" = English_Common_Name)

det_prob_traits %>% 
  select(species, det_prob_median_mean) %>%  View

plot(det_prob_traits$det_prob_median_mean ~ det_prob_traits$Mass)
summary(lm(det_prob_median_mean ~ Trophic.Level + Trophic.Niche + Mass + Primary.Lifestyle + Habitat, data = det_prob_traits))

# boxplots for categorial traits:

plot_dir2 <- file.path(plot_dir, "explorative_plots")
if(!dir.exists(plot_dir2)){dir.create(plot_dir2)}

category <- c("ORDER", "Habitat", "Migration", "Trophic.Level",
              "Trophic.Niche", "Primary.Lifestyle")

# to add number of cases above boxplots:
get_box_stats <- function(y, upper_limit = 1.15) {
  return(data.frame(
    y = 0.95 * upper_limit,
    label = length(y)
  ))
}

for(c in 1:length(category)){
  
  print(category[c])
  
  jpeg(file = file.path(plot_dir2, paste0("det_prob_", category[c], ".jpg")), 
       width = 1200, height = 900, quality = 100)
  print(
    
    ggplot(det_prob_traits, 
           aes(x = forcats::fct_reorder(.data[[category[c]]], det_prob_median_mean, .fun = median, .desc =TRUE), y = det_prob_median_mean)) +
      geom_boxplot() +
      ggtitle(paste("Detection probability and", category[c])) +
      ylab("mean detection probability across route sections") +
      xlab("") +
      stat_summary(fun.data = get_box_stats, geom = "text", hjust = 0.5, vjust = 0.9, size = 8) +
      theme_bw() +
      theme(#axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1), 
            text = element_text(size = 25))
  )
  dev.off()
  
}


# is model performance related to detection probability: -----------------------

# okay in time:
load(file.path("M:", "Documents", "DEBTs", "analysis", "Schifferle_BBS_occupancy_models_2023",
               "results", "temp_val_buffer_750_10yrs", "temp_eval", "10_years", "spec_set_temp_val_ok1.RData")) # 3_2_quant_temp_performance.R
spec_temp_okay <- specs_thresh

# okay in space:
CV_eval_summary <- read.csv(file = file.path("M:", "Documents", "DEBTs", "analysis", "Schifferle_BBS_occupancy_models_2023",
                                             "results", "CV_buffer750km", "CV_eval", "CV_eval_summary.csv")) # 3_1_DOM_CV_evaluation_metrics.R
spec_spat_okay <- CV_eval_summary %>% 
  filter(occ_spat_auc_mean >= 0.7) %>% 
  pull(species)

# okay in both:
spec_okay <- intersect(spec_temp_okay, spec_spat_okay) # 80

det_prob_traits <- det_prob_traits %>% 
  mutate(perf_fine = factor(ifelse(species %in% spec_okay, 1, 0))) %>% 
  mutate(temp_perf_fine = factor(ifelse(species %in% spec_temp_okay, 1, 0))) %>% 
  mutate(spat_perf_fine = factor(ifelse(species %in% spec_spat_okay, 1, 0))) %>% 
  left_join(CV_eval_summary %>%  select(species, occ_spat_auc_mean))
  
ggplot(det_prob_traits, 
       aes(x = perf_fine, y = det_prob_median_mean)) +
  geom_boxplot() +
  stat_summary(fun.data = get_box_stats, geom = "text", hjust = 0.5, vjust = 0.9, size = 8) +
  xlab("overall performance") +
  ylab("mean detection probability across route sections") +
  theme_bw()
# no difference
ggplot(det_prob_traits, 
       aes(x = temp_perf_fine, y = det_prob_median_mean)) +
  geom_boxplot() +
  stat_summary(fun.data = get_box_stats, geom = "text", hjust = 0.5, vjust = 0.9, size = 8) +
  xlab("temp. performance") +
  ylab("mean detection probability across route sections") +
  theme_bw()
# no difference
ggplot(det_prob_traits, 
       aes(x = spat_perf_fine, y = det_prob_median_mean)) +
  geom_boxplot() +
  stat_summary(fun.data = get_box_stats, geom = "text", hjust = 0.5, vjust = 0.9, size = 8) +
  xlab("spat. performance") +
  ylab("mean detection probability across route sections") +
  theme_bw()

# spatial AUC and detection probability:
ggplot(det_prob_traits) +
  geom_point(aes(x = det_prob_median_mean, y = occ_spat_auc_mean)) +
  xlab("mean detection probability across route sections") +
  ylab("mean spatial AUC") +
  geom_smooth(aes(x = det_prob_median_mean, y = occ_spat_auc_mean), method = "lm") +
  theme_bw() +
  theme(text = element_text(size = 16))

# temporal performance:
load(file.path("M:", "Documents", "DEBTs", "analysis", "Schifferle_BBS_occupancy_models_2023", 
               "results", "temp_val_buffer_750_10yrs", "temp_eval", "10_years", "temp_val_metrics_final.RData")) # 3_1_DOM_CV_evaluation_metrics.R
temp_val_metrics

det_prob_traits %>% 
  left_join(temp_val_metrics) %>% 
  ggplot() +
  geom_point(aes(x = det_prob_median_mean, y = rmse)) +
  xlab("mean detection probability across route sections") +
  #ylab("mean absolute percentage error time series") +
  #ylab("temporal C index") +
  ylab("root mean squared error") +
  geom_smooth(aes(x = det_prob_median_mean, y = rmse), method = "lm") +
  theme_bw() +
  theme(text = element_text(size = 16))
