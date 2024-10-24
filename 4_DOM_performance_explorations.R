# explore results of spatially blocked cross-validation of DOMs:
# for what kind of species is predictive performance in space/time/for future data better or worse:

# packages: ----

library(dplyr)
library(ggplot2)
library(ggrepel)

# functions: -----

source("0_functions.R")

# dir <- file.path("//NAS-2-P-SN-01.ibb.uni-potsdam.de", "users$", "schifferle1", "Documents", "DEBTs", "analysis", 
#                                 "Schifferle_BBS_occupancy_models_2023") 
dir <- getwd()


# load data:

# trait data, merged with BBS data:
load(file.path("data", "BBS_data_merged.RData")) # bbs_dt; output of Schifferle_BBS_explorations_2023/BBS_data_prep.R

# selected species:
load(file = file.path("data", "final_species_selection_eco_sorted.RData")) # final_species_eco_sorted; output of 1_2_species_selection.R
final_species_eco_sorted

# ecoregions:
load(file = file.path("data", "species_ecoregions.RData")) # spec_eco_df; output of 1_2_species_selection.R

# performance metrics:
CV_eval_summary <- read.csv(file = file.path(dir, "results", "CV_cluster", "CV_eval2", "CV_eval_summary.csv")) # output of 3_1_DOM_CV_evaluation_metrics.R
load(file = file.path("data", "C_temp_val.RData")) # C_temp_val_df; output of 3_1_DOM_temp_val_eval.R


# merge data:

spec_traits_perf_df <- spec_traits_df %>% 
  left_join(CV_eval_summary) %>% 
  left_join(C_temp_val_df)

spec_traits_df <- spec_eco_df %>% 
  left_join(bbs_dt[, c("English_Common_Name", "ORDER", "Family", "Genus", "HWI", "Mass", "Habitat", "Migration", "Trophic.Level",
                       "Trophic.Niche", "Primary.Lifestyle", "Min.Latitude", "Max.Latitude", "Centroid.Latitude", "Centroid.Longitude",
                       "Range.Size")], by = c("species" = "English_Common_Name"), multiple = "first")

spec_traits_df



range(spec_traits_perf_df$y_spattemp_C - spec_traits_perf_df$y_spattemp_auc, na.rm = TRUE)
range(spec_traits_perf_df$occ_spattemp_C - spec_traits_perf_df$occ_spattemp_auc, na.rm = TRUE)
# -> spatio-temporal: AUC and C-Index identical
range(spec_traits_perf_df$y_spatial_C_mean - spec_traits_perf_df$y_spatial_auc_mean, na.rm = TRUE)
range(spec_traits_perf_df$occ_spatial_C_mean - spec_traits_perf_df$occ_spatial_auc_mean, na.rm = TRUE)
# -> mean spatial AUC and mean spatial C-index almost identical

ggplot(spec_traits_perf_df) +
  geom_point(aes(x = y_spatial_auc_mean, y = y_temp_C)) +
  theme_bw() +
  ggtitle("CV") +
  ylim(c(0.25, 1)) +
  xlim(c(0.45, 1)) +
  geom_hline(yintercept = 0.5, linetype = "dashed") +
  geom_vline(xintercept = 0.5, linetype = "dashed")
# -> spatial predictive performance for all species better than random guessing (between 0.6 and 0.97)
summary(spec_traits_perf_df$y_spatial_auc_mean)
summary(spec_traits_perf_df$occ_spatial_auc_mean)
length(which(spec_traits_perf_df$y_temp_C > 0.5))
# -> temporal predictive performance for 129 of 174 species better than random guessing (between 0.26 and 0.97)
summary(spec_traits_perf_df$y_temp_C)
summary(spec_traits_perf_df$occ_temp_C)
# -> temporal predictive performance varies more across species than spatial predictive performance
cor.test(spec_traits_perf_df$occ_spatial_auc_mean, spec_traits_perf_df$y_temp_C, use = "complete.obs")
# -> if model performance better in space it tends to perform worse in time and wise versa


# CV performance and traits: ----

dir.create(file.path("plots", "CV_performance_explorations"))

category <- c("eco", "ORDER", "Habitat", "Migration", "Trophic.Level",
              "Trophic.Niche", "Primary.Lifestyle")

for(c in 1:length(category)){
  
  print(category[c])
  
  jpeg(file = file.path("plots", "CV_traits", paste0("CV_", category[c], ".jpg")), 
       width = 1200, height = 900, quality = 100)

  print(ggplot(spec_traits_perf_df) +
          geom_hline(yintercept = 0.5, linetype = "dashed", color = "grey40") +
          geom_vline(xintercept = 0.5, linetype = "dashed", color = "grey40") +
          geom_point(aes(x = y_spatial_auc_mean, y = y_temp_C, colour = .data[[category[c]]]), size = 4) +
          scale_color_brewer(palette = "Paired", na.value = "grey50") +
          theme_bw() +
          ggtitle("CV") +
          ylim(c(0.25, 1)) +
          xlim(c(0.45, 1)) +
          theme(text = element_text(size = 20))
        )
  dev.off()
  
}


# CV performance and families: ----

Family <- names(which(table(spec_traits_perf_df$Family) > 1))

for(f in 1:length(Family)){
  
  print(Family[f])
  
  jpeg(file = file.path("plots", "CV_traits", paste0("CV_", Family[f], ".jpg")), 
       width = 1200, height = 900, quality = 100)
  
  print(ggplot(spec_traits_perf_df) +
          geom_hline(yintercept = 0.5, linetype = "dashed", color = "grey40") +
          geom_vline(xintercept = 0.5, linetype = "dashed", color = "grey40") +
          geom_point(aes(x = y_spatial_auc_mean, y = y_temp_C), size = 4) +
          geom_point(data = spec_traits_perf_df[which(spec_traits_perf_df$Family == Family[f]), ], 
                     aes(x = y_spatial_auc_mean, y = y_temp_C), colour="goldenrod1", size=5) +
          theme_bw() +
          ggtitle(paste("CV", Family[f])) +
          ylim(c(0.25, 1)) +
          xlim(c(0.45, 1)) +
          theme(text = element_text(size = 20))
  )
  dev.off()
  
}

# some interesting groups: ----


# Wrens:

jpeg(file = file.path("plots", "CV_traits", "CV_Wrens.jpg"), 
     width = 1200, height = 900, quality = 100)
ggplot(spec_traits_perf_df) +
  geom_point(aes(x = y_spatial_auc_mean, y = y_temp_C), colour = "grey50", size = 4) +
  geom_point(data = spec_traits_perf_df[which(spec_traits_perf_df$species == "Carolina Wren"), ],
             aes(x = y_spatial_auc_mean, y = y_temp_C), colour="red", size=5) +
  geom_point(data = spec_traits_perf_df[which(spec_traits_perf_df$species == "Bewick's Wren"), ],
             aes(x = y_spatial_auc_mean, y = y_temp_C), colour="blue", size=5) +
  geom_point(data = spec_traits_perf_df[which(spec_traits_perf_df$species == "House Wren"), ],
             aes(x = y_spatial_auc_mean, y = y_temp_C), colour="coral4", size=5) +
  geom_point(data = spec_traits_perf_df[which(spec_traits_perf_df$species == "Marsh Wren"), ],
             aes(x = y_spatial_auc_mean, y = y_temp_C), colour="olivedrab", size=5)+
  geom_point(data = spec_traits_perf_df[which(spec_traits_perf_df$species == "Rock Wren"), ],
             aes(x = y_spatial_auc_mean, y = y_temp_C), colour="grey20", size=5) +
  geom_point(data = spec_traits_perf_df[which(spec_traits_perf_df$species == "Sedge Wren"), ],
             aes(x = y_spatial_auc_mean, y = y_temp_C), colour="darkgoldenrod2", size=5) +
  geom_point(data = spec_traits_perf_df[which(spec_traits_perf_df$species == "Winter Wren"), ],
             aes(x = y_spatial_auc_mean, y = y_temp_C), colour="cadetblue1", size=5) +
  geom_label_repel(data = subset(spec_traits_perf_df, Family == "Troglodytidae"),
                aes(x = y_spatial_auc_mean, y = y_temp_C, label = species),
                nudge_y = seq(-0.2, 0.2, length = 7),
                size = 8,
                box.padding = 1.2,
                point.padding = 0.5,
                force = 100,
                segment.size  = 0.2,
                segment.color = "grey50",
                direction = "x") +

  theme_bw() +
  ggtitle("Wrens") +
  ylim(c(0.25, 1)) +
  xlim(c(0.45, 1)) +
  geom_hline(yintercept = 0.5, linetype = "dashed") +
  geom_vline(xintercept = 0.5, linetype = "dashed") +
  theme(text = element_text(size = 20))
dev.off()


# birds of prey (Accipitridae, Falconidae, Pandionidae):

jpeg(file = file.path("plots", "CV_traits", "CV_Birds_of_prey.jpg"), 
     width = 1200, height = 900, quality = 100)
ggplot(spec_traits_perf_df) +
  geom_point(aes(x = y_spatial_auc_mean, y = y_temp_C), colour = "grey50", size = 4) +
  geom_point(data = spec_traits_perf_df[which(spec_traits_perf_df$Family == "Accipitridae"), ],
             aes(x = y_spatial_auc_mean, y = y_temp_C), colour="red", size=5) +
  geom_point(data = spec_traits_perf_df[which(spec_traits_perf_df$Family == "Falconidae"), ],
             aes(x = y_spatial_auc_mean, y = y_temp_C), colour="blue", size=5) +
  geom_point(data = spec_traits_perf_df[which(spec_traits_perf_df$Family == "Pandionidae"), ],
             aes(x = y_spatial_auc_mean, y = y_temp_C), colour="coral4", size=5) +
  geom_label_repel(data = subset(spec_traits_perf_df, Family == "Accipitridae" | Family == "Falconidae"| Family == "Pandionidae"),
                   aes(x = y_spatial_auc_mean, y = y_temp_C, label = species),
                   nudge_y = seq(-0.1, 0.1, length = 13),
                   #nudge_y       = 0.9 - subset(spec_traits_perf_df, Family == "Accipitridae" | Family == "Falconidae"| Family == "Pandionidae")$y_temp_C,
                   size = 8,
                   box.padding = 1.2,
                   point.padding = 0.5,
                   force = 100,
                   segment.size  = 0.2,
                   segment.color = "grey50",
                   direction = "x") +
  
  theme_bw() +
  ggtitle("Birds of prey") +
  ylim(c(0.25, 1)) +
  xlim(c(0.45, 1)) +
  geom_hline(yintercept = 0.5, linetype = "dashed") +
  geom_vline(xintercept = 0.5, linetype = "dashed") +
  theme(text = element_text(size = 20))
dev.off()

# Woodpeckers:

jpeg(file = file.path("plots", "CV_traits", "CV_Woodpeckers.jpg"), 
     width = 1200, height = 900, quality = 100)
ggplot(spec_traits_perf_df) +
  geom_point(aes(x = y_spatial_auc_mean, y = y_temp_C), colour = "grey50", size = 4) +
  geom_point(data = spec_traits_perf_df[which(spec_traits_perf_df$Genus == "Dryobates"), ],
             aes(x = y_spatial_auc_mean, y = y_temp_C), colour="lightblue", size=5) +
  geom_point(data = spec_traits_perf_df[which(spec_traits_perf_df$Genus == "Dryocopus"), ],
             aes(x = y_spatial_auc_mean, y = y_temp_C), colour="black", size=5) +
  geom_point(data = spec_traits_perf_df[which(spec_traits_perf_df$Genus == "Melanerpes"), ],
             aes(x = y_spatial_auc_mean, y = y_temp_C), colour="red", size=5) +
  geom_point(data = spec_traits_perf_df[which(spec_traits_perf_df$Genus == "Sphyrapicus"), ],
             aes(x = y_spatial_auc_mean, y = y_temp_C), colour="yellow2", size=5) +
  geom_label_repel(data = subset(spec_traits_perf_df, Genus == "Dryobates" | Genus == "Melanerpes" | Genus == "Sphyrapicus" | Genus == "Dryocopus"),
                   aes(x = y_spatial_auc_mean, y = y_temp_C, label = species),
                   nudge_y = seq(-0.1, 0.1, length = 7),
                   size = 8,
                   box.padding = 1.2,
                   point.padding = 0.5,
                   force = 100,
                   segment.size  = 0.2,
                   segment.color = "grey50",
                   direction = "x") +
  
  theme_bw() +
  ggtitle("Woodpeckers") +
  ylim(c(0.25, 1)) +
  xlim(c(0.45, 1)) +
  geom_hline(yintercept = 0.5, linetype = "dashed") +
  geom_vline(xintercept = 0.5, linetype = "dashed") +
  theme(text = element_text(size = 20))
dev.off()

# Vireos:

jpeg(file = file.path("plots", "CV_traits", "CV_Vireos.jpg"), 
     width = 1200, height = 900, quality = 100)
ggplot(spec_traits_perf_df) +
  geom_point(aes(x = y_spatial_auc_mean, y = y_temp_C), colour = "grey50", size = 4) +
  geom_point(data = spec_traits_perf_df[which(spec_traits_perf_df$Family == "Vireonidae"), ],
             aes(x = y_spatial_auc_mean, y = y_temp_C), colour="green", size=5) +
  geom_label_repel(data = subset(spec_traits_perf_df, Family == "Vireonidae"),
                   aes(x = y_spatial_auc_mean, y = y_temp_C, label = species),
                   nudge_y = seq(-0.1, 0.1, length = 6),
                   size = 8,
                   box.padding = 1.2,
                   point.padding = 0.5,
                   force = 100,
                   segment.size  = 0.2,
                   segment.color = "grey50",
                   direction = "x") +
  
  theme_bw() +
  ggtitle("Vireos") +
  ylim(c(0.25, 1)) +
  xlim(c(0.45, 1)) +
  geom_hline(yintercept = 0.5, linetype = "dashed") +
  geom_vline(xintercept = 0.5, linetype = "dashed") +
  theme(text = element_text(size = 20))
dev.off()


# Cardinals:

jpeg(file = file.path("plots", "CV_traits", "CV_Cardinals.jpg"), 
     width = 1200, height = 900, quality = 100)
ggplot(spec_traits_perf_df) +
  geom_point(aes(x = y_spatial_auc_mean, y = y_temp_C), colour = "grey50", size = 4) +
  geom_point(data = spec_traits_perf_df[which(spec_traits_perf_df$Family == "Cardinalidae"), ],
             aes(x = y_spatial_auc_mean, y = y_temp_C, colour = Genus), size=5) +
  geom_label_repel(data = subset(spec_traits_perf_df, Family == "Cardinalidae"),
                   aes(x = y_spatial_auc_mean, y = y_temp_C, label = species),
                   #nudge_y = seq(-0.1, 0.1, length = 11),
                   size = 8,
                   box.padding = 1.2,
                   point.padding = 0.5,
                   force = 100,
                   segment.size  = 0.2,
                   segment.color = "grey50",
                   #direction = "x"
                   ) +
  
  theme_bw() +
  ggtitle("Cardinals") +
  ylim(c(0.25, 1)) +
  xlim(c(0.45, 1)) +
  geom_hline(yintercept = 0.5, linetype = "dashed") +
  geom_vline(xintercept = 0.5, linetype = "dashed") +
  theme(text = element_text(size = 20))
dev.off()


# swallows etc.:

jpeg(file = file.path("plots", "CV_traits", "CV_Swallows.jpg"), 
     width = 1200, height = 900, quality = 100)
ggplot(spec_traits_perf_df) +
  geom_point(aes(x = y_spatial_auc_mean, y = y_temp_C), colour = "grey50", size = 4) +
  geom_point(data = spec_traits_perf_df[which(spec_traits_perf_df$Family == "Hirundinidae"), ],
             aes(x = y_spatial_auc_mean, y = y_temp_C, colour = Habitat), size=5) +
  geom_label_repel(data = subset(spec_traits_perf_df, Family == "Hirundinidae"),
                   aes(x = y_spatial_auc_mean, y = y_temp_C, label = species),
                   size = 8,
                   box.padding = 1.2,
                   point.padding = 0.5,
                   force = 100,
                   max.overlaps = 20, 
                   segment.size  = 0.2,
                   segment.color = "grey50") +
  theme_bw() +
  ggtitle("Swallows") +
  ylim(c(0.25, 1)) +
  xlim(c(0.45, 1)) +
  geom_hline(yintercept = 0.5, linetype = "dashed") +
  geom_vline(xintercept = 0.5, linetype = "dashed") +
  theme(text = element_text(size = 20))
dev.off()



# temporal predictive performance based on CV vs. performance of predicting future data: ----

ggplot(spec_traits_perf_df) +
  geom_point(aes(x = y_temp_C, y = C_ind_5yrs_preds)) +
  theme_bw() 
cor.test(spec_traits_perf_df$y_temp_C, 
         spec_traits_perf_df$C_ind_5yrs_preds, use = "complete.obs")
# -> correlated, but rather weakly

ggplot(spec_traits_perf_df) +
  geom_point(aes(x = y_spatial_auc_mean, y = C_ind_5yrs_preds)) +
  geom_point(aes(x = y_spatial_auc_mean, y = y_temp_C), color = "green") +
  theme_bw() 



# boxplots: ----

category <- c("eco", "ORDER", "Habitat", "Migration", "Trophic.Level",
              "Trophic.Niche", "Primary.Lifestyle")

# spatial predictive performance:

for(c in 1:length(category)){
  
  print(category[c])
  
  jpeg(file = file.path("plots", "CV_traits", paste0("CV_", category[c], "_spat_boxplot.jpg")), 
       width = 1200, height = 900, quality = 100)
  
  print(ggplot(spec_traits_perf_df) +
          geom_boxplot(aes(y = y_spatial_auc_mean, x = .data[[category[c]]])) +
          theme_bw() +
          ggtitle("CV") +
          theme(text = element_text(size = 20)) +
          annotate("text",
                   x = 1:length(table(spec_traits_perf_df[[category[c]]])[which(table(spec_traits_perf_df[[category[c]]]) > 0)]),
                   y = 1.05,
                   label = table(spec_traits_perf_df[[category[c]]])[which(table(spec_traits_perf_df[[category[c]]]) > 0)],
                   size = 10)
  )
  dev.off()
  
}

# temporal predictive performance:

for(c in 1:length(category)){
  
  print(category[c])
  
  jpeg(file = file.path("plots", "CV_traits", paste0("CV_", category[c], "_temp_boxplot.jpg")), 
       width = 1200, height = 900, quality = 100)
  
  print(ggplot(spec_traits_perf_df) +
          geom_boxplot(aes(y = y_temp_C, x = .data[[category[c]]])) +
          theme_bw() +
          ggtitle("CV") +
          theme(text = element_text(size = 20)) +
          annotate("text",
                   x = 1:length(table(spec_traits_perf_df[[category[c]]])[which(table(spec_traits_perf_df[[category[c]]]) > 0)]),
                   y = 1.05,
                   label = table(spec_traits_perf_df[[category[c]]])[which(table(spec_traits_perf_df[[category[c]]]) > 0)],
                   size = 10)
  )
  dev.off()
  
}


# predictive performance on future data:

for(c in 1:length(category)){
  
  print(category[c])
  
  jpeg(file = file.path("plots", "CV_traits", paste0("temp_val_", category[c], "_boxplot.jpg")), 
       width = 1200, height = 900, quality = 100)
  
  print(ggplot(spec_traits_perf_df) +
          geom_boxplot(aes(y = C_ind_5yrs_preds, x = .data[[category[c]]])) +
          theme_bw() +
          ggtitle("CV") +
          theme(text = element_text(size = 20)) +
          annotate("text",
                   x = 1:length(table(spec_traits_perf_df[[category[c]]])[which(table(spec_traits_perf_df[[category[c]]]) > 0)]),
                   y = 1.05,
                   label = table(spec_traits_perf_df[[category[c]]])[which(table(spec_traits_perf_df[[category[c]]]) > 0)],
                   size = 10)
  )
  dev.off()
  
}


# some species: ----

ggplot(spec_traits_perf_df) +
  geom_hline(yintercept = 0.5, linetype = "dashed") +
  geom_vline(xintercept = 0.5, linetype = "dashed") +
  geom_point(aes(x = y_spatial_C_mean, y = y_temp_C)) +
  geom_point(data = spec_traits_perf_df[which(spec_traits_perf_df$species == "Carolina Wren"), ], 
             aes(x = y_spatial_C_mean, y = y_temp_C), colour="blue", size=5) +
  geom_point(data = spec_traits_perf_df[which(spec_traits_perf_df$species == "Northern Cardinal"), ], 
             aes(x = y_spatial_C_mean, y = y_temp_C), colour="limegreen", size=5) +
  geom_point(data = spec_traits_perf_df[which(spec_traits_perf_df$species == "Carolina Chickadee"), ], 
             aes(x = y_spatial_C_mean, y = y_temp_C), colour="yellow2", size=5) +
  geom_point(data = spec_traits_perf_df[which(spec_traits_perf_df$species == "Red-bellied Woodpecker"), ], 
             aes(x = y_spatial_C_mean, y = y_temp_C), colour="orange", size=5) +
  geom_point(data = spec_traits_perf_df[which(spec_traits_perf_df$species == "Tufted Titmouse"), ], 
             aes(x = y_spatial_C_mean, y = y_temp_C), colour="pink", size=5) +
  geom_point(data = spec_traits_perf_df[which(spec_traits_perf_df$species == "American Redstart"), ], 
             aes(x = y_spatial_C_mean, y = y_temp_C), colour="red", size=5) +
  geom_point(data = spec_traits_perf_df[which(spec_traits_perf_df$species == "Wild Turkey"), ], 
             aes(x = y_spatial_C_mean, y = y_temp_C), colour="purple", size=5) +
  geom_point(data = spec_traits_perf_df[which(spec_traits_perf_df$species == "Scarlet Tanager"), ], 
             aes(x = y_spatial_C_mean, y = y_temp_C), colour="tomato", size=5) +
  geom_point(data = spec_traits_perf_df[which(spec_traits_perf_df$species == "Lark Bunting"), ], 
             aes(x = y_spatial_C_mean, y = y_temp_C), colour="lightblue", size=5) +
  geom_point(data = spec_traits_perf_df[which(spec_traits_perf_df$species == "Chestnut-sided Warbler"), ], 
             aes(x = y_spatial_C_mean, y = y_temp_C), colour="brown4", size=5) +
  geom_point(data = spec_traits_perf_df[which(spec_traits_perf_df$species == "Eastern Meadowlark"), ], 
             aes(x = y_spatial_C_mean, y = y_temp_C), colour="lightgreen", size=5) +
  geom_point(data = spec_traits_perf_df[which(spec_traits_perf_df$species == "Least Flycatcher"), ], 
             aes(x = y_spatial_C_mean, y = y_temp_C), colour="darkgreen", size=5) +
  geom_point(data = spec_traits_perf_df[which(spec_traits_perf_df$species == "Northern Mockingbird"), ], 
             aes(x = y_spatial_C_mean, y = y_temp_C), colour="grey", size=5) +
  geom_point(data = spec_traits_perf_df[which(spec_traits_perf_df$species == "Pine Warbler"), ], 
             aes(x = y_spatial_C_mean, y = y_temp_C), colour="violet", size=5) +
  geom_point(data = spec_traits_perf_df[which(spec_traits_perf_df$species == "Prairie Warbler"), ], 
             aes(x = y_spatial_C_mean, y = y_temp_C), colour="khaki2", size=5) +
  geom_point(data = spec_traits_perf_df[which(spec_traits_perf_df$species == "Prothonotary Warbler"), ], 
             aes(x = y_spatial_C_mean, y = y_temp_C), colour="aquamarine", size=5) +
  geom_point(data = spec_traits_perf_df[which(spec_traits_perf_df$species == "Purple Finch"), ], 
             aes(x = y_spatial_C_mean, y = y_temp_C), colour="deeppink4", size=5) +
  theme_bw()


# some lms: ----

summary(lm(y_spatial_auc_mean  ~ eco, data = spec_traits_perf_df))
summary(lm(y_temp_C  ~ eco, data = spec_traits_perf_df))

summary(lm(y_spatial_auc_mean ~ Range.Size, data = spec_traits_perf_df))
# slightly higher spatial predictive performance for smaller range size
summary(lm(y_temp_C ~ Range.Size, data = spec_traits_perf_df))
# slightly higher temporal predictive performance for smaller range size

# latitudinal range centroid:
summary(lm(y_spatial_auc_mean ~ Centroid.Latitude, data = spec_traits_perf_df))
summary(lm(y_temp_C ~ Centroid.Latitude, data = spec_traits_perf_df))

summary(lm(y_spatial_auc_mean ~ Max.Latitude, data = spec_traits_perf_df))
summary(lm(y_temp_C ~ Max.Latitude, data = spec_traits_perf_df))
summary(lm(y_spatial_auc_mean ~ Min.Latitude, data = spec_traits_perf_df))
# slightly higher spatial performance for species with ranges extending less to the south
summary(lm(y_temp_C ~ Min.Latitude, data = spec_traits_perf_df))
# slightly higher temporal performance for species with ranges extending further south


# explore for single species how env. conditions changed over time: ----

# load data:

# env data:
load(file = file.path("data", "route_year_env_data.RData")) # route_sel_env_dt_final
route_sel_env_dt_final
# routes-years:
load(file = file.path("data", "BBS_for_occ_selection.RData")) # route_sel_dt; output of 1_3_match_BBS_to_env_data.R 
# selected routes spatial data (to buffer presences):
routes_sel_sf <- st_read(file.path("data", "route_selection_1991_2015_surv_beg_end_max_5y_miss_v2_spat_thin_100km_max_30_r_per_BCR_centroids.shp")) # output of 1_1_route_selection.R
# route-year-species information (only surveyed)
load(file = file.path("data", "BBS_for_occ_spec_records.RData")) # bbs_dt_occ; output of 1_0_reformat_BBS_data.R

# species:
spec <- "Bald Eagle"

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
  geom_point(aes(x = Year, y = pr_winter)) +
  geom_smooth(aes(x = Year, y = pr_winter), method = "lm")
# winter precipitation decreased over time on routes that got colonized
dt %>% 
  filter(colo == -1) %>% 
  ggplot() +
  geom_point(aes(x = Year, y = pr_winter)) +
  geom_smooth(aes(x = Year, y = pr_winter), method = "lm")
# winter precipitation decreased over time on routes where species went (temporarily) extinct
dt %>% 
  #filter(colo %in% c(0, 1)) %>% 
  ggplot() +
  geom_point(aes(x = Year, y = pr_winter, color = colo)) +
  geom_smooth(aes(x = Year, y = pr_winter, group = colo, colour = colo), method = "lm")

dt %>% 
  #filter(colo == 1) %>% 
  ggplot() +
  geom_point(aes(x = Year, y = secdf, color = colo)) +
  geom_smooth(aes(x = Year, y = secdf, group = colo, colour = colo), method = "lm")

dt %>% 
  #filter(colo == 1) %>% 
  ggplot() +
  geom_point(aes(x = Year, y = bio1, color = colo)) +
  geom_smooth(aes(x = Year, y = bio1, group = colo, colour = colo), method = "lm")


# Elton traits...