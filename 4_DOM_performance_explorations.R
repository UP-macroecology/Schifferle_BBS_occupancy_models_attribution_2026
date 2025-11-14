# explore results of spatially blocked cross-validation of DOMs:
# for what kind of species is predictive performance in space/time/for future data better or worse:

# packages: ----

library(dplyr)
library(ggplot2)
library(ggrepel)

# directories: ----

results_dir <- file.path("M:", "Documents", "DEBTs", "analysis", "Schifferle_BBS_occupancy_models_2023", 
                         "results")
dir <- getwd()


# functions: -----

source("0_functions.R")


# load data: ----

# trait data, merged with BBS data:
load(file.path("data", "BBS_data_merged.RData")) # bbs_dt; output of Schifferle_BBS_explorations_2023/BBS_data_prep.R

# selected species:
load(file = file.path("data", "species_set_analysis.RData"))
final_species

# ecoregions:
load(file = file.path("data", "species_ecoregions.RData")) # spec_eco_df; output of 1_2_species_selection.R

# performance metrics:

# cross validation:
CV_eval_summary <- read.csv(file = file.path(results_dir,  "CV_buffer750km", "CV_eval", "CV_eval_summary.csv")) # output of 3_1_DOM_CV_evaluation_metrics.R
# temporal validation:
load(file = file.path(results_dir, "temp_val_buffer_750_10yrs", "temp_eval", "10_years", "temp_val_metrics_final.RData")) # output of 3_1_DOM_temp_evaluation_metrics.R
temp_val_metrics


# merge data:
spec_traits_df <- spec_eco_df %>% 
  filter(species %in% final_species) %>% 
  left_join(bbs_dt[, c("English_Common_Name", "ORDER", "Family", "Genus", "HWI", "Mass", "Habitat", "Migration", "Trophic.Level",
                       "Trophic.Niche", "Primary.Lifestyle", "Min.Latitude", "Max.Latitude", "Centroid.Latitude", "Centroid.Longitude",
                       "Range.Size")], by = c("species" = "English_Common_Name"), multiple = "first")

spec_traits_perf_df <- spec_traits_df %>% 
  left_join(CV_eval_summary) %>% 
  left_join(temp_val_metrics)


# explorations: ----

range(spec_traits_perf_df$y_spattemp_C - spec_traits_perf_df$y_spattemp_auc, na.rm = TRUE)
range(spec_traits_perf_df$occ_spattemp_C - spec_traits_perf_df$occ_spattemp_auc, na.rm = TRUE)
# -> spatio-temporal: AUC and C-Index identical
range(spec_traits_perf_df$y_spat_C_mean - spec_traits_perf_df$y_spat_auc_mean, na.rm = TRUE)
range(spec_traits_perf_df$occ_spat_C_mean - spec_traits_perf_df$occ_spat_auc_mean, na.rm = TRUE)
# -> mean spatial AUC and mean spatial C-index identical


## spatial vs. temporal performance (CV based): ----

# boxplots validation metrics:
jpeg(file = file.path("plots", "performance_explorations", "CV_metrics_boxplot.jpg"), 
     width = 1000, height = 1200, quality = 100)
CV_eval_summary %>% 
  select(grep(pattern = "y_", x = colnames(CV_eval_summary), value = TRUE)) %>% 
  tidyr::pivot_longer(cols = y_spattemp_C:y_temp_C_c, names_to = "metric", values_to = "value") %>% 
  ggplot() + 
  geom_boxplot(aes(group = metric, y = value, x = metric)) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1), 
        text = element_text(size = 30)) +
  geom_hline(yintercept = 0.5, linetype = "dashed") +
  scale_x_discrete(labels = c("y spat. AUC", "y spat. AUC change", "y spat. C-index", "y spat. C-index change",
                              "y spat.-temp. AUC", "y spat.-temp. AUC change", "y spat.-temp. C-index", "y spat.-temp. C-index change",
                              "y temp. C-index", "y temp. C-index change"))
dev.off()

# scatterplot spatial. vs. temporal:
ggplot(spec_traits_perf_df) +
  geom_point(aes(x = y_spat_auc_mean, y = y_temp_C)) +
  theme_bw() +
  ggtitle("CV") +
  ylim(c(0.25, 1)) +
  xlim(c(0.45, 1)) +
  geom_hline(yintercept = 0.5, linetype = "dashed") +
  geom_vline(xintercept = 0.5, linetype = "dashed")
# -> spatial predictive performance for all species better than random guessing (between 0.66 and 0.98)
summary(spec_traits_perf_df$y_spat_auc_mean)
summary(spec_traits_perf_df$occ_spat_auc_mean)
length(which(spec_traits_perf_df$y_temp_C > 0.5))
# -> temporal predictive performance for 123 of 159 species better than random guessing (between 0.26 and 0.99)
summary(spec_traits_perf_df$y_temp_C)
# -> temporal predictive performance varies more across species than spatial predictive performance
cor.test(spec_traits_perf_df$occ_spat_auc_mean, spec_traits_perf_df$y_temp_C, use = "complete.obs")
# -> if model performance better in space it tends to perform worse in time and vice versa


## all routes vs. only routes with change in occupancy: ----

# is model only right at sites that are always / never occupied?:

# scatterplot spatialtemp. performance all routes vs. only routes with change:
jpeg(file = file.path("plots", "performance_explorations", "CV_all_vs_change_routes.jpg"), 
     width = 1400, height = 1000, quality = 100)
ggplot(data = CV_eval_summary) +
  geom_point(aes(x = y_spat_auc_mean, y = y_spat_auc_mean_c), size = 2) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed") +
  ylim(c(0.5, 1)) +
  xlim(0.5, 1) +
  geom_text_repel(data = CV_eval_summary,
                  aes(x = y_spat_auc_mean, y = y_spat_auc_mean_c, label = species),
                  size = 6,
                  force = 0.5,
                  max.overlaps = 12) +
  ylab("mean yearly spat.temp. AUC, routes with change") +
  xlab("mean yearly spat.temp. AUC, all routes") +
  theme_bw() +
  theme(text = element_text(size = 30))
dev.off()

# same, but spatial only and add information which species passed temporal validation criteria:
# species okay in time:
load(file.path("M:", "Documents", "DEBTs", "analysis", "Schifferle_BBS_occupancy_models_2023",
               "results", "temp_val_buffer_750_10yrs", "temp_eval", "10_years", "spec_set_temp_val_ok1.RData"))
spec_temp_okay <- specs_thresh

jpeg(file = file.path("plots", "performance_explorations", "spat_CV_all_vs_change_routes_plus_temp.jpg"), 
     width = 1400, height = 1000, quality = 100)
CV_eval_summary %>% 
  mutate(temp_val_fine = ifelse(species %in% spec_temp_okay, "yes", "no")) %>% 
  ggplot() +
  geom_point(aes(x = y_spat_auc_mean, y = y_spat_auc_mean_c, colour = as.factor(temp_val_fine)), size = 3) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed") +
  ylim(c(0.5, 0.85)) +
  xlim(0.6, 1) +
  scale_colour_manual(values = c("grey70", "dodgerblue3"), name = "temporal performance fine:") +
  geom_text_repel(data = CV_eval_summary,
                  aes(x = y_spat_auc_mean, y = y_spat_auc_mean_c, label = species),
                  size = 6,
                  force = 0.5,
                  max.overlaps = 10) +
  ylab("mean yearly AUC, routes with change") +
  xlab("mean yearly AUC, all routes") +
  theme_bw() +
  theme(text = element_text(size = 30), legend.position = "top")
dev.off()


## test some thresholds / criteria for validation: ----


# scatterplot spat. vs temp.; selection citeria incl. temporal validation
jpeg(file = file.path("plots", "performance_explorations", "CV_plus_temp_val_thresholds.jpg"), 
     width = 1200, height = 900, quality = 100)
ggplot(data = spec_traits_perf_df) +
  geom_point(aes(x = y_spat_C_mean, y = y_temp_C, 
                 colour = as.factor((y_temp_C >= 0.7 | trend_corr_CV == 1) & 
                                      (C_ind_10yrs_preds >= 0.7 | trend_corr == 1))),
             size = 3) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed") +
  labs(colour = "CV temp. C > 0.7 or trend correct and\ntemp. val. C > 0.7 or trend correct") +
  ylim(c(0.5, 1)) +
  xlim(c(0.5, 1)) +
  ylab("y temp. C-index") +
  xlab("y spat. C-index") +
  geom_hline(yintercept = 0.7, linetype = "dashed", col = "blue") +
  geom_vline(xintercept = 0.7, linetype = "dashed", col = "blue") +
  geom_text_repel(data = spec_traits_perf_df,
                  aes(x = y_spat_C_mean, y = y_temp_C, label = species),
                  size = 3, force = 0.5, max.overlaps = 10) +
  theme_bw() +
  theme(text = element_text(size = 20))
dev.off()

# spec_traits_perf_df %>%
#   filter(trend_corr == 1 | C_ind_10yrs_preds >= 0.7) %>%
#   filter(y_spattemp_C >= 0.7) %>%
#   filter(y_temp_C >= 0.7 | trend_corr_CV == 1) %>%
#   filter(y_spattemp_C_c >= 0.6) %>%
#   filter(y_temp_C_c >= 0.6) %>%
#   pull(species)
# # 
# species_fine <- spec_traits_perf_df %>% 
#   filter(y_spat_auc_mean >= 0.7) %>% 
#   filter((y_temp_C >= 0.7 | trend_corr_CV > 0) & 
#            (C_ind_10yrs_preds >= 0.7 | trend_corr > 0)) %>% 
#   pull(species) # 76
# species_fine
# 
# # also requirements for sites with change:
# species_fine2 <- spec_traits_perf_df %>% 
#   filter(y_spat_auc_mean >= 0.7 & y_spat_auc_mean_c >= 0.6) %>% 
#   filter((y_temp_C >= 0.7 | trend_corr_CV > 0) & 
#            (C_ind_10yrs_preds >= 0.7 | trend_corr > 0)) %>% 
#   pull(species) # 65
# species_fine2
# species_fine[!species_fine %in% species_fine2]
# 
# species_discard <- final_species[!final_species %in% species_fine]
# 
# spec_traits_perf_df %>% 
#   mutate(keep = ifelse(species %in% species_discard, 0, 1)) %>% 
#   #filter(species %in% species_discard) %>% 
#   select(species, y_spattemp_auc, y_spat_auc_mean, y_temp_C, trend_corr_CV, C_ind_10yrs_preds, 
#          trend_corr, trend_diff, keep) %>% 
#   View
# 
# species_fine2 <- spec_traits_perf_df %>%
#   filter(y_spattemp_C >= 0.7) %>%
#   filter(y_temp_C >= 0.8 | (y_temp_C >= 0.7 | trend_corr_CV > 0) |
#            (C_ind_10yrs_preds >= 0.8 | (C_ind_10yrs_preds >= 0.7 | trend_corr > 0))) %>%
#   pull(species)
# # 136
# 
# species_fine3 <- spec_traits_perf_df %>%
#   filter(y_spattemp_C >= 0.7) %>%
#   filter((y_temp_C >= 0.8 | (y_temp_C >= 0.7 | trend_corr_CV > 0)) &
#            (C_ind_10yrs_preds >= 0.8 | (C_ind_10yrs_preds >= 0.7 | trend_corr > 0))) %>%
#   pull(species)
# species_fine3 # 77
# 
# # # one temp val good, or the other: 136
# # # both: 77
# species_fine4 <- spec_traits_perf_df %>%
#   filter(y_spattemp_C >= 0.7) %>%
#   filter(y_temp_C >= 0.8 | C_ind_10yrs_preds >= 0.8 |
#            ((y_temp_C >= 0.6 | trend_corr_CV > 0) & (C_ind_10yrs_preds >= 0.6 | trend_corr > 0))) %>%
#   pull(species)
# species_fine4 # 86
# # either one C very good
# # or the other C very good
# # or either one C medium or trend captured and either other C medium or trend captured
# species_fine4[!species_fine4 %in% species_fine]
# species_fine[!species_fine %in% species_fine4]



## number of raw colonisation and extinction events and performance: ----

# raw colonisation and extinction events per species between 1995 and 2019:
  
raw_col_ext_df <- tibble("species" = final_species,
                         "raw_col_events" = NA,
                         "raw_ext_events" = NA)

for(i in 1:nrow(raw_col_ext_df)){
  
  spec <- raw_col_ext_df$species[i]
  
  # observations:
  occ_dt_spec <- BBS_pres_abs_spec(species = spec) %>% 
    select(RTENO, Year, presence) %>% 
    group_by(RTENO) %>% 
    mutate(switch = presence - lag(presence)) %>% 
    mutate(col_event = ifelse(switch == 1, 1, 0),
           ext_event = ifelse(switch == -1, 1, 0))
  
  raw_col_ext_df$raw_col_events[i] <- sum(occ_dt_spec$col_event, na.rm = TRUE)
  raw_col_ext_df$raw_ext_events[i] <- sum(occ_dt_spec$ext_event, na.rm = TRUE)
}

raw_col_ext_df
#save(raw_col_ext_df, file = file.path("data", "raw_col_ext_event_final_species_1995_2019.RData"))
#write.csv(raw_col_ext_df, file = file.path("data", "raw_col_ext_event_final_species_1995_2019.csv"), row.names = FALSE)

# is number of raw colonisation and extinction events correlated with DOM performance?
eval_raw_dt <- all_eval_metrics %>% 
  left_join(raw_col_ext_df) %>% 
  select(species, y_spattemp_auc, y_spattemp_auc_c, y_spat_auc_mean, y_spat_auc_mean_c, 
         y_temp_C, y_temp_C_c, C_ind_10yrs_preds, trend_corr, raw_col_events, raw_ext_events)

range(eval_raw_dt$raw_col_events)
range(eval_raw_dt$raw_ext_events)
# between ~ 60 and ~ 1300

plot(eval_raw_dt$raw_col_events, eval_raw_dt$y_spat_auc_mean)
plot(eval_raw_dt$raw_col_events, eval_raw_dt$y_spat_auc_mean_c)
plot(eval_raw_dt$raw_col_events, eval_raw_dt$y_temp_C)
plot(eval_raw_dt$raw_col_events, eval_raw_dt$C_ind_10yrs_preds)

cor <- cor(eval_raw_dt[, -1], method = "s")
corrplot::corrplot(cor, method = "square", type = "upper", order = "original", 
                   tl.col = "black", tl.srt = 45)

cor[, c("raw_col_events", "raw_ext_events")]
# not much correlation between number of raw dynamic events and model performance
# max. abs. correlation 0.39 -> more col. / ext. events -> lower spatial AUC
# but: considering only routes with change: 0.24 -> more col. / ext. events, slightly higher spatial AUC
# less correlation regarding temporal performance




## CV performance and traits: ----

if(!dir.exists(file.path("plots", "performance_explorations"))){
  dir.create(file.path("plots", "performance_explorations"))
}

category <- c("eco", "ORDER", "Habitat", "Migration", "Trophic.Level",
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
  
  # scatterplot spatial vs. temporal performance:
  jpeg(file = file.path("plots", "performance_explorations", paste0("CV_", category[c], ".jpg")), 
       width = 1200, height = 900, quality = 100)

  print(ggplot(spec_traits_perf_df) +
          geom_hline(yintercept = 0.5, linetype = "dashed", color = "grey40") +
          geom_vline(xintercept = 0.5, linetype = "dashed", color = "grey40") +
          geom_point(aes(x = y_spat_auc_mean, y = y_temp_C, colour = .data[[category[c]]]), size = 4) +
          scale_color_brewer(palette = "Paired", na.value = "grey50") +
          theme_bw() +
          ggtitle("CV") +
          ylim(c(0.25, 1)) +
          xlim(c(0.45, 1)) +
          theme(text = element_text(size = 35))
        )
  dev.off()
  
  # boxplots spatial vs. temporal performance:
  jpeg(file = file.path("plots", "performance_explorations", paste0("CV_spat_temp_", category[c], ".jpg")), 
       width = 1200, height = 900, quality = 100)
  print(spec_traits_perf_df %>% 
          tidyr::pivot_longer(cols = y_spattemp_C:occ_temp_C_c, 
                              names_to = "metric", values_to = "value") %>% 
          filter(metric %in% c("y_spat_auc_mean", "y_temp_C")) %>% 
          
          ggplot(aes(x = .data[[category[c]]], y = value, 
                   fill = as.factor(metric))) +
          geom_boxplot() +
          #scale_fill_brewer(palette = "Paired", na.value = "grey50") +
          scale_fill_discrete(name = "C index", labels = c("spatial", "temporal")) +
          ggtitle("CV") +
          stat_summary(fun.data = get_box_stats, geom = "text", hjust = 0.5, vjust = 0.9, size = 10) +
          theme_bw() +
          theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1), text = element_text(size = 35))
  )
  dev.off()
  
}




## CV performance and families: ----

Family <- names(which(table(spec_traits_perf_df$Family) > 1))

for(f in 1:length(Family)){
  
  print(Family[f])
  
  jpeg(file = file.path("plots", "performance_explorations", paste0("CV_", Family[f], ".jpg")), 
       width = 1200, height = 900, quality = 100)
  
  print(ggplot(spec_traits_perf_df) +
          geom_hline(yintercept = 0.5, linetype = "dashed", color = "grey40") +
          geom_vline(xintercept = 0.5, linetype = "dashed", color = "grey40") +
          geom_point(aes(x = y_spat_auc_mean, y = y_temp_C), size = 5) +
          geom_point(data = spec_traits_perf_df[which(spec_traits_perf_df$Family == Family[f]), ], 
                     aes(x = y_spat_auc_mean, y = y_temp_C), colour="goldenrod1", size=6) +
          theme_bw() +
          ggtitle(paste("CV", Family[f])) +
          ylim(c(0.25, 1)) +
          xlim(c(0.45, 1)) +
          ylab("y temp. C-index") +
          xlab("y spat. AUC") +
          theme(text = element_text(size = 30))
  )
  dev.off()
  
}


## temp. predictive performance based on CV vs. performance of predicting future data: ----


ggplot(spec_traits_perf_df) +
  geom_point(aes(x = y_temp_C, y = C_ind_10yrs_preds)) +
  theme_bw() 
cor.test(spec_traits_perf_df$y_temp_C, 
         spec_traits_perf_df$C_ind_10yrs_preds, use = "complete.obs")
# -> correlated, but rather weakly

# boxplots temp. predictive performance CV vs. future data:

category <- c("eco", "ORDER", "Habitat", "Migration", "Trophic.Level",
              "Trophic.Niche", "Primary.Lifestyle")

for(c in 1:length(category)){
  
  print(category[c])
  
  jpeg(file = file.path("plots", "performance_explorations", paste0("temp_val_vs_CV_temp_", category[c], ".jpg")), 
       width = 1200, height = 900, quality = 100)
  
  print(spec_traits_perf_df %>% 
          tidyr::pivot_longer(cols = y_spattemp_C:trend_diff, 
                              names_to = "metric", values_to = "value") %>% 
          filter(metric %in% c("y_temp_C", "C_ind_10yrs_preds")) %>%
          ggplot(aes(x = .data[[category[c]]], y = value, 
                     fill = as.factor(metric))) +
          geom_boxplot() +
          scale_fill_discrete(name = "C index", labels = c("temp. val.", "CV temp.")) +
          theme_bw() +
          ggtitle("temporal predictive performance") +
          stat_summary(fun.data = get_box_stats, geom = "text", hjust = 0.5, vjust = 0.9, size = 10) +
          theme_bw() +
          theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1), text = element_text(size = 35))
  )
  dev.off()
  
}

## some species groups: ----

# Wrens:

jpeg(file = file.path("plots", "performance_explorations", "CV_Wrens.jpg"), 
     width = 1200, height = 900, quality = 100)
ggplot(spec_traits_perf_df) +
  geom_point(aes(x = y_spat_auc_mean, y = y_temp_C), colour = "grey50", size = 4) +
  geom_point(data = spec_traits_perf_df[which(spec_traits_perf_df$species == "Carolina Wren"), ],
             aes(x = y_spat_auc_mean, y = y_temp_C), colour="red", size=5) +
  geom_point(data = spec_traits_perf_df[which(spec_traits_perf_df$species == "Bewick's Wren"), ],
             aes(x = y_spat_auc_mean, y = y_temp_C), colour="blue", size=5) +
  geom_point(data = spec_traits_perf_df[which(spec_traits_perf_df$species == "House Wren"), ],
             aes(x = y_spat_auc_mean, y = y_temp_C), colour="coral4", size=5) +
  geom_point(data = spec_traits_perf_df[which(spec_traits_perf_df$species == "Marsh Wren"), ],
             aes(x = y_spat_auc_mean, y = y_temp_C), colour="olivedrab", size=5)+
  geom_point(data = spec_traits_perf_df[which(spec_traits_perf_df$species == "Rock Wren"), ],
             aes(x = y_spat_auc_mean, y = y_temp_C), colour="grey20", size=5) +
  geom_point(data = spec_traits_perf_df[which(spec_traits_perf_df$species == "Sedge Wren"), ],
             aes(x = y_spat_auc_mean, y = y_temp_C), colour="darkgoldenrod2", size=5) +
  geom_point(data = spec_traits_perf_df[which(spec_traits_perf_df$species == "Winter Wren"), ],
             aes(x = y_spat_auc_mean, y = y_temp_C), colour="cadetblue1", size=5) +
  geom_point(data = spec_traits_perf_df[which(spec_traits_perf_df$species == "Canyon Wren"), ],
             aes(x = y_spat_auc_mean, y = y_temp_C), colour="orange", size=5) +
  geom_label_repel(data = subset(spec_traits_perf_df, Family == "Troglodytidae"),
                aes(x = y_spat_auc_mean, y = y_temp_C, label = species),
                nudge_y = seq(-0.2, 0.2, length = 7),
                size = 6,
                box.padding = 1.2,
                point.padding = 0.5,
                force = 100,
                segment.size  = 0.2,
                segment.color = "grey50",
                direction = "x") +
  ylab("y temp. C-index") + xlab("y spat. AUC") +
  theme_bw() +
  ggtitle("Wrens") +
  ylim(c(0.25, 1)) + xlim(c(0.55, 1)) +
  geom_hline(yintercept = 0.5, linetype = "dashed") +
  geom_vline(xintercept = 0.5, linetype = "dashed") +
  theme(text = element_text(size = 25))
dev.off()


# birds of prey (Accipitridae, Falconidae, Pandionidae):

jpeg(file = file.path("plots", "performance_explorations", "CV_Birds_of_prey.jpg"), 
     width = 1200, height = 900, quality = 100)
ggplot(spec_traits_perf_df) +
  geom_point(aes(x = y_spat_auc_mean, y = y_temp_C), colour = "grey50", size = 4) +
  geom_point(data = spec_traits_perf_df[which(spec_traits_perf_df$Family == "Accipitridae"), ],
             aes(x = y_spat_auc_mean, y = y_temp_C), colour="red", size=5) +
  geom_point(data = spec_traits_perf_df[which(spec_traits_perf_df$Family == "Falconidae"), ],
             aes(x = y_spat_auc_mean, y = y_temp_C), colour="blue", size=5) +
  geom_point(data = spec_traits_perf_df[which(spec_traits_perf_df$Family == "Pandionidae"), ],
             aes(x = y_spat_auc_mean, y = y_temp_C), colour="coral4", size=5) +
  geom_label_repel(data = subset(spec_traits_perf_df, Family == "Accipitridae" | Family == "Falconidae"| Family == "Pandionidae"),
                   aes(x = y_spat_auc_mean, y = y_temp_C, label = species),
                   nudge_y = seq(-0.1, 0.1, length = 6),
                   size = 8,
                   box.padding = 1.2,
                   point.padding = 0.5,
                   force = 100,
                   segment.size  = 0.2,
                   segment.color = "grey50",
                   direction = "x") +
  ylab("y temp. C-index") + xlab("y spat. AUC") +
  theme_bw() +
  ggtitle("Birds of prey") +
  ylim(c(0.25, 1)) + xlim(c(0.55, 1)) +
  geom_hline(yintercept = 0.5, linetype = "dashed") +
  geom_vline(xintercept = 0.5, linetype = "dashed") +
  theme(text = element_text(size = 25))
dev.off()


# Woodpeckers:

jpeg(file = file.path("plots", "performance_explorations", "CV_Woodpeckers.jpg"), 
     width = 1200, height = 900, quality = 100)
ggplot(spec_traits_perf_df) +
  geom_point(aes(x = y_spat_auc_mean, y = y_temp_C), colour = "grey50", size = 4) +
  geom_point(data = spec_traits_perf_df[which(spec_traits_perf_df$Genus == "Dryobates"), ],
             aes(x = y_spat_auc_mean, y = y_temp_C), colour="lightblue", size=5) +
  geom_point(data = spec_traits_perf_df[which(spec_traits_perf_df$Genus == "Dryocopus"), ],
             aes(x = y_spat_auc_mean, y = y_temp_C), colour="black", size=5) +
  geom_point(data = spec_traits_perf_df[which(spec_traits_perf_df$Genus == "Melanerpes"), ],
             aes(x = y_spat_auc_mean, y = y_temp_C), colour="red", size=5) +
  geom_point(data = spec_traits_perf_df[which(spec_traits_perf_df$Genus == "Sphyrapicus"), ],
             aes(x = y_spat_auc_mean, y = y_temp_C), colour="yellow2", size=5) +
  geom_label_repel(data = subset(spec_traits_perf_df, Genus == "Dryobates" | Genus == "Melanerpes" | Genus == "Sphyrapicus" | Genus == "Dryocopus"),
                   aes(x = y_spat_auc_mean, y = y_temp_C, label = species),
                   nudge_y = seq(-0.1, 0.1, length = 6),
                   size = 6,
                   box.padding = 1.2,
                   point.padding = 0.5,
                   force = 100,
                   segment.size  = 0.2,
                   segment.color = "grey50",
                   direction = "x") +
  ylab("y temp. C-index") + xlab("y spat. AUC") +
  theme_bw() +
  ggtitle("Woodpeckers") +
  ylim(c(0.25, 1)) + xlim(c(0.55, 1)) +
  geom_hline(yintercept = 0.5, linetype = "dashed") +
  geom_vline(xintercept = 0.5, linetype = "dashed") +
  theme(text = element_text(size = 25))
dev.off()


# Vireos:

jpeg(file = file.path("plots", "performance_explorations", "CV_Vireos.jpg"), 
     width = 1200, height = 900, quality = 100)
ggplot(spec_traits_perf_df) +
  geom_point(aes(x = y_spat_auc_mean, y = y_temp_C), colour = "grey50", size = 4) +
  geom_point(data = spec_traits_perf_df[which(spec_traits_perf_df$Family == "Vireonidae"), ],
             aes(x = y_spat_auc_mean, y = y_temp_C), colour="green", size=5) +
  geom_label_repel(data = subset(spec_traits_perf_df, Family == "Vireonidae"),
                   aes(x = y_spat_auc_mean, y = y_temp_C, label = species),
                   nudge_y = seq(-0.1, 0.1, length = 7),
                   size = 6,
                   box.padding = 1.2,
                   point.padding = 0.5,
                   force = 100,
                   segment.size  = 0.2,
                   segment.color = "grey50",
                   direction = "x") +
  ylab("y temp. C-index") + xlab("y spat. AUC") +
  theme_bw() +
  ggtitle("Vireos") +
  ylim(c(0.25, 1)) + xlim(c(0.55, 1)) +
  geom_hline(yintercept = 0.5, linetype = "dashed") +
  geom_vline(xintercept = 0.5, linetype = "dashed") +
  theme(text = element_text(size = 25))
dev.off()


# cardinals:

jpeg(file = file.path("plots", "performance_explorations", "CV_Cardinals.jpg"), 
     width = 1200, height = 900, quality = 100)
ggplot(spec_traits_perf_df) +
  geom_point(aes(x = y_spat_auc_mean, y = y_temp_C), colour = "grey50", size = 4) +
  geom_point(data = spec_traits_perf_df[which(spec_traits_perf_df$Family == "Cardinalidae"), ],
             aes(x = y_spat_auc_mean, y = y_temp_C, colour = Genus), size=5) +
  geom_label_repel(data = subset(spec_traits_perf_df, Family == "Cardinalidae"),
                   aes(x = y_spat_auc_mean, y = y_temp_C, label = species),
                   size = 6,
                   box.padding = 1.2,
                   point.padding = 0.5,
                   force = 100,
                   segment.size  = 0.2,
                   segment.color = "grey50",
                   ) +
  ylab("y temp. C-index") + xlab("y spat. AUC") +
  theme_bw() +
  ggtitle("Cardinals") +
  ylim(c(0.25, 1)) +
  xlim(c(0.55, 1)) +
  geom_hline(yintercept = 0.5, linetype = "dashed") +
  geom_vline(xintercept = 0.5, linetype = "dashed") +
  theme(text = element_text(size = 25))
dev.off()


# swallows etc.:

jpeg(file = file.path("plots", "performance_explorations", "CV_Swallows.jpg"), 
     width = 1200, height = 900, quality = 100)
ggplot(spec_traits_perf_df) +
  geom_point(aes(x = y_spat_auc_mean, y = y_temp_C), colour = "grey50", size = 4) +
  geom_point(data = spec_traits_perf_df[which(spec_traits_perf_df$Family == "Hirundinidae"), ],
             aes(x = y_spat_auc_mean, y = y_temp_C, colour = Habitat), size=5) +
  geom_label_repel(data = subset(spec_traits_perf_df, Family == "Hirundinidae"),
                   aes(x = y_spat_auc_mean, y = y_temp_C, label = species),
                   size = 6,
                   box.padding = 1.2,
                   point.padding = 0.5,
                   force = 100,
                   max.overlaps = 20, 
                   segment.size  = 0.2,
                   segment.color = "grey50") +
  ylab("y temp. C-index") + xlab("y spat. AUC") +
  theme_bw() +
  ggtitle("Swallows") +
  ylim(c(0.25, 1)) + xlim(c(0.55, 1)) +
  geom_hline(yintercept = 0.5, linetype = "dashed") +
  geom_vline(xintercept = 0.5, linetype = "dashed") +
  theme(text = element_text(size = 20))
dev.off()





## some species: ----

ggplot(spec_traits_perf_df) +
  geom_hline(yintercept = 0.5, linetype = "dashed") +
  geom_vline(xintercept = 0.5, linetype = "dashed") +
  geom_point(aes(x = y_spat_C_mean, y = y_temp_C)) +
  geom_point(data = spec_traits_perf_df[which(spec_traits_perf_df$species == "Carolina Wren"), ], 
             aes(x = y_spat_C_mean, y = y_temp_C), colour="blue", size=5) +
  geom_point(data = spec_traits_perf_df[which(spec_traits_perf_df$species == "Northern Cardinal"), ], 
             aes(x = y_spat_C_mean, y = y_temp_C), colour="limegreen", size=5) +
  geom_point(data = spec_traits_perf_df[which(spec_traits_perf_df$species == "Carolina Chickadee"), ], 
             aes(x = y_spat_C_mean, y = y_temp_C), colour="yellow2", size=5) +
  geom_point(data = spec_traits_perf_df[which(spec_traits_perf_df$species == "Red-bellied Woodpecker"), ], 
             aes(x = y_spat_C_mean, y = y_temp_C), colour="orange", size=5) +
  geom_point(data = spec_traits_perf_df[which(spec_traits_perf_df$species == "Tufted Titmouse"), ], 
             aes(x = y_spat_C_mean, y = y_temp_C), colour="pink", size=5) +
  geom_point(data = spec_traits_perf_df[which(spec_traits_perf_df$species == "American Redstart"), ], 
             aes(x = y_spat_C_mean, y = y_temp_C), colour="red", size=5) +
  geom_point(data = spec_traits_perf_df[which(spec_traits_perf_df$species == "Wild Turkey"), ], 
             aes(x = y_spat_C_mean, y = y_temp_C), colour="purple", size=5) +
  geom_point(data = spec_traits_perf_df[which(spec_traits_perf_df$species == "Scarlet Tanager"), ], 
             aes(x = y_spat_C_mean, y = y_temp_C), colour="tomato", size=5) +
  geom_point(data = spec_traits_perf_df[which(spec_traits_perf_df$species == "Lark Bunting"), ], 
             aes(x = y_spat_C_mean, y = y_temp_C), colour="lightblue", size=5) +
  geom_point(data = spec_traits_perf_df[which(spec_traits_perf_df$species == "Chestnut-sided Warbler"), ], 
             aes(x = y_spat_C_mean, y = y_temp_C), colour="brown4", size=5) +
  geom_point(data = spec_traits_perf_df[which(spec_traits_perf_df$species == "Eastern Meadowlark"), ], 
             aes(x = y_spat_C_mean, y = y_temp_C), colour="lightgreen", size=5) +
  geom_point(data = spec_traits_perf_df[which(spec_traits_perf_df$species == "Least Flycatcher"), ], 
             aes(x = y_spat_C_mean, y = y_temp_C), colour="darkgreen", size=5) +
  geom_point(data = spec_traits_perf_df[which(spec_traits_perf_df$species == "Northern Mockingbird"), ], 
             aes(x = y_spat_C_mean, y = y_temp_C), colour="grey", size=5) +
  geom_point(data = spec_traits_perf_df[which(spec_traits_perf_df$species == "Pine Warbler"), ], 
             aes(x = y_spat_C_mean, y = y_temp_C), colour="violet", size=5) +
  geom_point(data = spec_traits_perf_df[which(spec_traits_perf_df$species == "Prairie Warbler"), ], 
             aes(x = y_spat_C_mean, y = y_temp_C), colour="khaki2", size=5) +
  geom_point(data = spec_traits_perf_df[which(spec_traits_perf_df$species == "Prothonotary Warbler"), ], 
             aes(x = y_spat_C_mean, y = y_temp_C), colour="aquamarine", size=5) +
  geom_point(data = spec_traits_perf_df[which(spec_traits_perf_df$species == "Purple Finch"), ], 
             aes(x = y_spat_C_mean, y = y_temp_C), colour="deeppink4", size=5) +
  theme_bw()


## some lm explorations on performance ~ trait/range characteristics etc.: ----

summary(lm(y_spat_auc_mean  ~ eco, data = spec_traits_perf_df))
summary(lm(y_temp_C  ~ eco, data = spec_traits_perf_df))

summary(lm(y_spat_auc_mean ~ Range.Size, data = spec_traits_perf_df))
# slightly higher spatial predictive performance for smaller range size
summary(lm(y_temp_C ~ Range.Size, data = spec_traits_perf_df))

# latitudinal range centroid:
summary(lm(y_spat_auc_mean ~ Centroid.Latitude, data = spec_traits_perf_df))
summary(lm(y_temp_C ~ Centroid.Latitude, data = spec_traits_perf_df))

summary(lm(y_spat_auc_mean ~ Max.Latitude, data = spec_traits_perf_df))
summary(lm(y_temp_C ~ Max.Latitude, data = spec_traits_perf_df))
summary(lm(y_spat_auc_mean ~ Min.Latitude, data = spec_traits_perf_df))
# slightly higher spatial performance for species with ranges extending less to the south
summary(lm(y_temp_C ~ Min.Latitude, data = spec_traits_perf_df))
# slightly higher temporal performance for species with ranges extending further south


## (change of env. data over time regarding single species):----

# load data:

# env data:
load(file = file.path("data", "route_year_env_data.RData")) # route_sel_env_dt_final
# routes-years:
load(file = file.path("data", "BBS_for_occ_selection.RData")) # route_sel_dt; output of 1_3_match_BBS_to_env_data.R 
# selected routes spatial data (to buffer presences):
routes_sel_sf <- st_read(file.path("data", "route_selection_1995_2019_surv_beg_end_max_5y_miss_v2_spat_thin_100km_max_30_r_per_BCR_centroids.shp")) # output of 1_1_route_selection.R
# route-year-species information (only surveyed)
load(file = file.path("data", "BBS_for_occ_spec_records.RData")) # bbs_dt_occ; output of 1_0_reformat_BBS_data.R

# species:
spec <- "American Kestrel"

# species presences-absences:
occ_dt_spec <- BBS_pres_abs_spec(species = spec)

# relevant routes, within distance of 750 km of species records:
rel_routes <- training_routes(species = spec, buffer_km = 750, output = "RTENOs")

dt <- occ_dt_spec %>% 
  filter(RTENO %in% rel_routes) %>%   
  group_by(RTENO) %>% 
  # add whether colonization or extinction event happened on a route in a given year:
  mutate(colo = as.factor(presence - lag(presence))) %>% 
  left_join(route_sel_env_dt_final, by = c("RTENO", "Year"))

# plots:

dt %>% 
  filter(colo == 1) %>% 
  ggplot() +
  geom_point(aes(x = Year, y = pr_mean_winter)) +
  geom_smooth(aes(x = Year, y = pr_mean_winter), method = "lm")
# winter precipitation decreased over time on routes that got colonized
dt %>% 
  filter(colo == -1) %>% 
  ggplot() +
  geom_point(aes(x = Year, y = pr_mean_winter)) +
  geom_smooth(aes(x = Year, y = pr_mean_winter), method = "lm")
# winter precipitation decreased over time on routes where species went (temporarily) extinct
dt %>% 
  #filter(colo %in% c(0, 1)) %>% 
  ggplot() +
  geom_point(aes(x = Year, y = pr_mean_winter, color = colo)) +
  geom_smooth(aes(x = Year, y = pr_mean_winter, group = colo, colour = colo), method = "lm")

dt %>% 
  #filter(colo == 1) %>% 
  ggplot() +
  geom_point(aes(x = Year, y = bio1, color = colo)) +
  geom_smooth(aes(x = Year, y = bio1, group = colo, colour = colo), method = "lm")


## ((Elton traits xx)): ----

## temporal performance according to criteria set and species: ----

load(file = file.path("data", "species_DOM_val_okay.RData"))
spec_okay


final_val_expl_df <- tibble(species = final_species) %>% 
  mutate(val_okay = factor(ifelse(species %in% spec_okay, "yes", "no"))) %>% 
  # add traits and phylogenetic information:
  left_join(bbs_dt %>% 
              select(English_Common_Name, Family, ORDER, HWI, Mass, Habitat, Migration, Trophic.Level, Trophic.Niche, Primary.Lifestyle, Range.Size, Centroid.Latitude, Centroid.Longitude) %>% distinct,
            by = c(species = "English_Common_Name"))
final_val_expl_df

ggplot(final_val_expl_df) +
  geom_boxplot(aes(x = val_okay, y = Range.Size)) +
  theme_bw()
# range size. centroid latitude and longitude not important
# HWI and Mass not important

final_val_expl_df %>% 
  group_by(Family, val_okay) %>% 
  summarise(n = n()) %>% 
  ggplot() +
  geom_linerange(aes(x = Family, 
                     ymin = 0, ymax = n)) +
  geom_point(aes(y = n, x = Family, colour = val_okay)) +
  coord_flip() +
  theme_bw()

final_val_expl_df %>% 
  group_by(Primary.Lifestyle, val_okay) %>% 
  summarise(n = n()) %>% 
  ggplot() +
  geom_linerange(aes(x = Primary.Lifestyle, 
                     ymin = 0, ymax = n)) +
  geom_point(aes(y = n, x = Primary.Lifestyle, colour = val_okay)) +
  coord_flip() +
  theme_bw()
# Habitat, Migration, Trophic.Level, Trophic.Niche, Primary.Lifestyle
# no clear difference!