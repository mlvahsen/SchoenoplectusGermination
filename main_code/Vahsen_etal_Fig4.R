# Create Figure 4 in Vahsen et al. germination paper
# Also calculates predicted means and 95% CI for in-text
# 19 Feb 2021

## Preliminaries ####
# Load libraries
library(tidyverse);library(patchwork);library(ggmcmc)
library(EnvStats);library(cowplot); library(here)

# Set ggplot theme
theme_set(theme_classic())

# Read in data that was used to fit the models (same as original data file that
# takes out grouping by Assay_Num)
germ_all_final <- read_csv(here("outputs","germ_all_data_forModel.csv"))

# Read in coda objects from all four models
Model3_coda <- read_rds(here("outputs","Model3_coda_forPlotting.rds"))

## Calculate predicted means for each treatment ####  

# This averages across all other treatments and predicts means at depth = 0.

# Extract regression slopes from best model
ggs(Model3_coda) %>% 
  filter(Parameter %in% c("b0", "b1", "b4", "b5", "b6", "b7", "b8", "b9", "b10", "b11")) %>% 
  spread(key = Parameter, value = value) -> pred_betas

# Get predicted standardized depths
ggs(Model3_coda, family = "xa_std") %>% 
  group_by(Parameter) %>% 
  summarize(mean = mean(value)) %>% 
  mutate(depth = germ_all_final$Depth_Top) -> pred_xa_std

# Get the minimum depth from the standardized depths vector. This is to predict
# means across treatments at depth = 0
xa_std_min <- min(pred_xa_std$mean)

# Create empty vectors to store the predicted values for each level of each
# environmental factor

# Temperature
temp1_pred <- NULL
temp2_pred <- NULL
temp3_pred <- NULL
temp4_pred <- NULL
# Media
media1_pred <- NULL
media2_pred <- NULL
media3_pred <- NULL
# Treatment
treatment0_pred <- NULL
treatment1_pred <- NULL
# Photoperiod
photoperiod1_pred <- NULL
photoperiod2_pred <- NULL
photoperiod3_pred <- NULL

# Calculate average values across other treatments except for temperature --
# note this is not a weighted mean, but just showing differences across
# treatments

for(i in 1:nrow(pred_betas)){
  temp3_pred[i] <- mean(# photoperiod = 3
                        plogis(pred_betas$b0[i]+pred_betas$b1[i]*xa_std_min+pred_betas$b4[i]), # Media 2
                        plogis(pred_betas$b0[i]+pred_betas$b1[i]*xa_std_min+pred_betas$b5[i]), # Media 3
                        plogis(pred_betas$b0[i]+pred_betas$b1[i]*xa_std_min), 
                        plogis(pred_betas$b0[i]+pred_betas$b1[i]*xa_std_min+pred_betas$b4[i]+pred_betas$b9[i]),# Media 2 + Treatment 1
                        plogis(pred_betas$b0[i]+pred_betas$b1[i]*xa_std_min+pred_betas$b5[i]+pred_betas$b9[i]),# Media 3 + Treatment 1
                        plogis(pred_betas$b0[i]+pred_betas$b1[i]*xa_std_min+pred_betas$b9[i]),# Treatment 1
                        # photoperiod = 1
                        plogis(pred_betas$b0[i]+pred_betas$b1[i]*xa_std_min+pred_betas$b4[i]+pred_betas$b10[i]),# Media 2 + Photo 1
                        plogis(pred_betas$b0[i]+pred_betas$b1[i]*xa_std_min+pred_betas$b5[i]+pred_betas$b10[i]),# Media 3 + Photo 1
                        plogis(pred_betas$b0[i]+pred_betas$b1[i]*xa_std_min+pred_betas$b10[i]), # Photo 1
                        plogis(pred_betas$b0[i]+pred_betas$b1[i]*xa_std_min+pred_betas$b4[i]+pred_betas$b9[i]+pred_betas$b10[i]), # Media 2 + Treatment 1 + Photo 1
                        plogis(pred_betas$b0[i]+pred_betas$b1[i]*xa_std_min+pred_betas$b5[i]+pred_betas$b9[i]+pred_betas$b10[i]), # Media 3 + Treatment 1 + Photo 1
                        plogis(pred_betas$b0[i]+pred_betas$b1[i]*xa_std_min+pred_betas$b9[i]+pred_betas$b10[i]), # Treatment 1 + Photo 1
                        # photoperiod = 2
                        plogis(pred_betas$b0[i]+pred_betas$b1[i]*xa_std_min+pred_betas$b4[i]+pred_betas$b11[i]), # Media 2 + Photo 2
                        plogis(pred_betas$b0[i]+pred_betas$b1[i]*xa_std_min+pred_betas$b5[i]+pred_betas$b11[i]),# Media 3 + Photo 2
                        plogis(pred_betas$b0[i]+pred_betas$b1[i]*xa_std_min+pred_betas$b11[i]),# Photo 2
                        plogis(pred_betas$b0[i]+pred_betas$b1[i]*xa_std_min+pred_betas$b4[i]+pred_betas$b9[i]+pred_betas$b11[i]),# Media 2 + Treatment 1 + Photo 2
                        plogis(pred_betas$b0[i]+pred_betas$b1[i]*xa_std_min+pred_betas$b5[i]+pred_betas$b9[i]+pred_betas$b11[i]),# Media 2 + Treatment 1 + Photo 2
                        plogis(pred_betas$b0[i]+pred_betas$b1[i]*xa_std_min+pred_betas$b9[i]+pred_betas$b11[i]))# Treatment 1 + Photo 2
  
  temp1_pred[i] <- mean(# photoperiod = 3
    plogis(pred_betas$b0[i]+pred_betas$b1[i]*xa_std_min+pred_betas$b4[i]+pred_betas$b6[i]),# Media 2 + Temp 1
    plogis(pred_betas$b0[i]+pred_betas$b1[i]*xa_std_min+pred_betas$b5[i]+pred_betas$b6[i]),# Media 3 + Temp 1
    plogis(pred_betas$b0[i]+pred_betas$b1[i]*xa_std_min+pred_betas$b6[i]),# Temp 1
    plogis(pred_betas$b0[i]+pred_betas$b1[i]*xa_std_min+pred_betas$b4[i]+pred_betas$b9[i]+pred_betas$b6[i]),# Media 2 + Treatment 1 + Temp 1
    plogis(pred_betas$b0[i]+pred_betas$b1[i]*xa_std_min+pred_betas$b5[i]+pred_betas$b9[i]+pred_betas$b6[i]),# Media 3 + Treatment 1 + Temp 1
    plogis(pred_betas$b0[i]+pred_betas$b1[i]*xa_std_min+pred_betas$b9[i]+pred_betas$b6[i]),# Treatment 1 + Temp 1
    # photoperiod = 1
    plogis(pred_betas$b0[i]+pred_betas$b1[i]*xa_std_min+pred_betas$b4[i]+pred_betas$b10[i]+pred_betas$b6[i]),# Media 2 + Photo 1 + Temp 1
    plogis(pred_betas$b0[i]+pred_betas$b1[i]*xa_std_min+pred_betas$b5[i]+pred_betas$b10[i]+pred_betas$b6[i]),# Media 3 + Photo 1 + Temp 1
    plogis(pred_betas$b0[i]+pred_betas$b1[i]*xa_std_min+pred_betas$b10[i]+pred_betas$b6[i]),# Photo 1 + Temp 1
    plogis(pred_betas$b0[i]+pred_betas$b1[i]*xa_std_min+pred_betas$b4[i]+pred_betas$b9[i]+pred_betas$b10[i]+pred_betas$b6[i]), # Media 2 + Treatment 1 + Photo 1 + Temp 1
    plogis(pred_betas$b0[i]+pred_betas$b1[i]*xa_std_min+pred_betas$b5[i]+pred_betas$b9[i]+pred_betas$b10[i]+pred_betas$b6[i]),# Media 3 + Treatment 1 + Photo 1 + Temp 1
    plogis(pred_betas$b0[i]+pred_betas$b1[i]*xa_std_min+pred_betas$b9[i]+pred_betas$b10[i]+pred_betas$b6[i]),# Treatment 1 + Photo 1 + Temp 1
    # photoperiod = 2
    plogis(pred_betas$b0[i]+pred_betas$b1[i]*xa_std_min+pred_betas$b4[i]+pred_betas$b11[i]+pred_betas$b6[i]),# Media 2 + Photo 2 + Temp 1
    plogis(pred_betas$b0[i]+pred_betas$b1[i]*xa_std_min+pred_betas$b5[i]+pred_betas$b11[i]+pred_betas$b6[i]),# Media 3 + Photo 2 + Temp 1
    plogis(pred_betas$b0[i]+pred_betas$b1[i]*xa_std_min+pred_betas$b11[i]+pred_betas$b6[i]),# Photo 1 + Temp 1
    plogis(pred_betas$b0[i]+pred_betas$b1[i]*xa_std_min+pred_betas$b4[i]+pred_betas$b9[i]+pred_betas$b11[i]+pred_betas$b6[i]),# Media 2 + Treatment 1 + Photo 2 + Temp 1
    plogis(pred_betas$b0[i]+pred_betas$b1[i]*xa_std_min+pred_betas$b5[i]+pred_betas$b9[i]+pred_betas$b11[i]+pred_betas$b6[i]),# Media 3 + Treatment 1 + Photo 2 + Temp 1
    plogis(pred_betas$b0[i]+pred_betas$b1[i]*xa_std_min+pred_betas$b9[i]+pred_betas$b11[i]+pred_betas$b6[i]))# Treatment 1 + Photo 2 + Temp 1
  
  temp2_pred[i] <- mean(# photoperiod = 3
    plogis(pred_betas$b0[i]+pred_betas$b1[i]*xa_std_min+pred_betas$b4[i]+pred_betas$b7[i]), # Media 2 + Temp 2
    plogis(pred_betas$b0[i]+pred_betas$b1[i]*xa_std_min+pred_betas$b5[i]+pred_betas$b7[i]),# Media 3 + Temp 2
    plogis(pred_betas$b0[i]+pred_betas$b1[i]*xa_std_min+pred_betas$b7[i]),# Temp 2
    plogis(pred_betas$b0[i]+pred_betas$b1[i]*xa_std_min+pred_betas$b4[i]+pred_betas$b9[i]+pred_betas$b7[i]),# Media 2 +Treatment 1+ Temp 2
    plogis(pred_betas$b0[i]+pred_betas$b1[i]*xa_std_min+pred_betas$b5[i]+pred_betas$b9[i]+pred_betas$b7[i]),# Media 3 +Treatment 1+ Temp 2
    plogis(pred_betas$b0[i]+pred_betas$b1[i]*xa_std_min+pred_betas$b9[i]+pred_betas$b7[i]),# Treatment 1+ Temp 2
    # photoperiod = 1
    plogis(pred_betas$b0[i]+pred_betas$b1[i]*xa_std_min+pred_betas$b4[i]+pred_betas$b10[i]+pred_betas$b7[i]),# Media 2 + Photo 1 + Temp 2
    plogis(pred_betas$b0[i]+pred_betas$b1[i]*xa_std_min+pred_betas$b5[i]+pred_betas$b10[i]+pred_betas$b7[i]),# Media 3 + Photo 1 + Temp 2
    plogis(pred_betas$b0[i]+pred_betas$b1[i]*xa_std_min+pred_betas$b10[i]+pred_betas$b7[i]),# Photo 1 + Temp 2
    plogis(pred_betas$b0[i]+pred_betas$b1[i]*xa_std_min+pred_betas$b4[i]+pred_betas$b9[i]+pred_betas$b10[i]+pred_betas$b7[i]),# Media 2 + Treatment 1 + Photo 1 + Temp 2
    plogis(pred_betas$b0[i]+pred_betas$b1[i]*xa_std_min+pred_betas$b5[i]+pred_betas$b9[i]+pred_betas$b10[i]+pred_betas$b7[i]),# Media 3 + Treatment 1 + Photo 1 + Temp 2
    plogis(pred_betas$b0[i]+pred_betas$b1[i]*xa_std_min+pred_betas$b9[i]+pred_betas$b10[i]+pred_betas$b7[i]),# Media 2 + Treatment 1 + Photo 1 + Temp 2
    # photoperiod = 2
    plogis(pred_betas$b0[i]+pred_betas$b1[i]*xa_std_min+pred_betas$b4[i]+pred_betas$b11[i]+pred_betas$b7[i]),# Media 2 + Photo 2 + Temp 2
    plogis(pred_betas$b0[i]+pred_betas$b1[i]*xa_std_min+pred_betas$b5[i]+pred_betas$b11[i]+pred_betas$b7[i]),# Media 3 + Photo 2 + Temp 2
    plogis(pred_betas$b0[i]+pred_betas$b1[i]*xa_std_min+pred_betas$b11[i]+pred_betas$b7[i]),# Photo 2 + Temp 2
    plogis(pred_betas$b0[i]+pred_betas$b1[i]*xa_std_min+pred_betas$b4[i]+pred_betas$b9[i]+pred_betas$b11[i]+pred_betas$b7[i]),# Media 2 + Treatment 1 + Photo 2 + Temp 2
    plogis(pred_betas$b0[i]+pred_betas$b1[i]*xa_std_min+pred_betas$b5[i]+pred_betas$b9[i]+pred_betas$b11[i]+pred_betas$b7[i]),# Media 3 + Treatment 1 + Photo 2 + Temp 2
    plogis(pred_betas$b0[i]+pred_betas$b1[i]*xa_std_min+pred_betas$b9[i]+pred_betas$b11[i]+pred_betas$b7[i])) # Treatment 1 + Photo 2 + Temp 2
  
  temp4_pred[i] <- mean(# photoperiod = 3
    plogis(pred_betas$b0[i]+pred_betas$b1[i]*xa_std_min+pred_betas$b4[i]+pred_betas$b8[i]),# Media 2 + Temp 4
    plogis(pred_betas$b0[i]+pred_betas$b1[i]*xa_std_min+pred_betas$b5[i]+pred_betas$b8[i]),# Media 3 + Temp 4
    plogis(pred_betas$b0[i]+pred_betas$b1[i]*xa_std_min+pred_betas$b8[i]), # Temp 4
    plogis(pred_betas$b0[i]+pred_betas$b1[i]*xa_std_min+pred_betas$b4[i]+pred_betas$b9[i]+pred_betas$b8[i]),# Media 2 +Treatment 1+ Temp 4
    plogis(pred_betas$b0[i]+pred_betas$b1[i]*xa_std_min+pred_betas$b5[i]+pred_betas$b9[i]+pred_betas$b8[i]),# Media 3 +Treatment 1+ Temp 4
    plogis(pred_betas$b0[i]+pred_betas$b1[i]*xa_std_min+pred_betas$b9[i]+pred_betas$b8[i]),# Treatment 1+ Temp 2
    # photoperiod = 1
    plogis(pred_betas$b0[i]+pred_betas$b1[i]*xa_std_min+pred_betas$b4[i]+pred_betas$b10[i]+pred_betas$b8[i]),# Media 2 +Photo 1+ Temp 4
    plogis(pred_betas$b0[i]+pred_betas$b1[i]*xa_std_min+pred_betas$b5[i]+pred_betas$b10[i]+pred_betas$b8[i]),# Media 3 +Photo 1+ Temp 4
    plogis(pred_betas$b0[i]+pred_betas$b1[i]*xa_std_min+pred_betas$b10[i]+pred_betas$b8[i]),# Photo 1+ Temp 4
    plogis(pred_betas$b0[i]+pred_betas$b1[i]*xa_std_min+pred_betas$b4[i]+pred_betas$b9[i]+pred_betas$b10[i]+pred_betas$b8[i]),# Media 2 + Treatment 1 + Photo 1+ Temp 4
    plogis(pred_betas$b0[i]+pred_betas$b1[i]*xa_std_min+pred_betas$b5[i]+pred_betas$b9[i]+pred_betas$b10[i]+pred_betas$b8[i]),# Media 3 + Treatment 1 + Photo 1+ Temp 4
    plogis(pred_betas$b0[i]+pred_betas$b1[i]*xa_std_min+pred_betas$b9[i]+pred_betas$b10[i]+pred_betas$b8[i]),# Treatment 1 + Photo 1+ Temp 4
    # photoperiod = 2
    plogis(pred_betas$b0[i]+pred_betas$b1[i]*xa_std_min+pred_betas$b4[i]+pred_betas$b11[i]+pred_betas$b8[i]),# Media 2 +Photo 2+ Temp 4
    plogis(pred_betas$b0[i]+pred_betas$b1[i]*xa_std_min+pred_betas$b5[i]+pred_betas$b11[i]+pred_betas$b8[i]),# Media 3 +Photo 2+ Temp 4
    plogis(pred_betas$b0[i]+pred_betas$b1[i]*xa_std_min+pred_betas$b11[i]+pred_betas$b8[i]),# Photo 2+ Temp 4
    plogis(pred_betas$b0[i]+pred_betas$b1[i]*xa_std_min+pred_betas$b4[i]+pred_betas$b9[i]+pred_betas$b11[i]+pred_betas$b8[i]),# Media 2 + Treatment 1 + Photo 2+ Temp 4
    plogis(pred_betas$b0[i]+pred_betas$b1[i]*xa_std_min+pred_betas$b5[i]+pred_betas$b9[i]+pred_betas$b11[i]+pred_betas$b8[i]),# Media 3 + Treatment 1 + Photo 2+ Temp 4
    plogis(pred_betas$b0[i]+pred_betas$b1[i]*xa_std_min+pred_betas$b9[i]+pred_betas$b11[i]+pred_betas$b8[i]))# Treatment 1 + Photo 2+ Temp 4

}

# Calculate average values across other treatments except for media
for(i in 1:nrow(pred_betas)){
  media3_pred[i] <- mean(plogis(pred_betas$b0[i]+pred_betas$b1[i]*xa_std_min), 
                        plogis(pred_betas$b0[i]+pred_betas$b1[i]*xa_std_min+pred_betas$b6[i]), # Temp = 1
                        plogis(pred_betas$b0[i]+pred_betas$b1[i]*xa_std_min+pred_betas$b7[i]), # Temp = 2
                        plogis(pred_betas$b0[i]+pred_betas$b1[i]*xa_std_min+pred_betas$b8[i]),# Temp = 4
                        plogis(pred_betas$b0[i]+pred_betas$b1[i]*xa_std_min+pred_betas$b9[i]),# Treatment = 1
                        plogis(pred_betas$b0[i]+pred_betas$b1[i]*xa_std_min+pred_betas$b9[i]+pred_betas$b6[i]),# Treatment = 1, Temp = 1
                        plogis(pred_betas$b0[i]+pred_betas$b1[i]*xa_std_min+pred_betas$b9[i]+pred_betas$b7[i]),# Treatment = 1, Temp = 2
                        plogis(pred_betas$b0[i]+pred_betas$b1[i]*xa_std_min+pred_betas$b9[i]+pred_betas$b8[i]),# Treatment = 1, Temp = 4
                        # photoperiod = 1
                        plogis(pred_betas$b0[i]+pred_betas$b1[i]*xa_std_min+pred_betas$b10[i]),# Photo = 1
                        plogis(pred_betas$b0[i]+pred_betas$b1[i]*xa_std_min+pred_betas$b6[i]+pred_betas$b10[i]),# Temp = 1, Photo = 1
                        plogis(pred_betas$b0[i]+pred_betas$b1[i]*xa_std_min+pred_betas$b7[i]+pred_betas$b10[i]),# Temp = 2, Photo = 1
                        plogis(pred_betas$b0[i]+pred_betas$b1[i]*xa_std_min+pred_betas$b8[i]+pred_betas$b10[i]),# Temp = 4, Photo = 1
                        plogis(pred_betas$b0[i]+pred_betas$b1[i]*xa_std_min+pred_betas$b9[i]+pred_betas$b10[i]),# Treatment = 1, Photo = 1
                        plogis(pred_betas$b0[i]+pred_betas$b1[i]*xa_std_min+pred_betas$b9[i]+pred_betas$b6[i]+pred_betas$b10[i]),# Treatment = 1, Temp = 1, Photo = 1
                        plogis(pred_betas$b0[i]+pred_betas$b1[i]*xa_std_min+pred_betas$b9[i]+pred_betas$b7[i]+pred_betas$b10[i]),# Treatment = 1, Temp = 2, Photo = 1
                        plogis(pred_betas$b0[i]+pred_betas$b1[i]*xa_std_min+pred_betas$b9[i]+pred_betas$b8[i]+pred_betas$b10[i]),# Treatment = 1, Temp = 4, Photo = 1
                        # photoperiod = 2
                        plogis(pred_betas$b0[i]+pred_betas$b1[i]*xa_std_min+pred_betas$b11[i]),# Photo = 2
                        plogis(pred_betas$b0[i]+pred_betas$b1[i]*xa_std_min+pred_betas$b6[i]+pred_betas$b11[i]),# Temp = 1, Photo = 2
                        plogis(pred_betas$b0[i]+pred_betas$b1[i]*xa_std_min+pred_betas$b7[i]+pred_betas$b11[i]),# Temp = 2, Photo = 2
                        plogis(pred_betas$b0[i]+pred_betas$b1[i]*xa_std_min+pred_betas$b8[i]+pred_betas$b11[i]),# Temp = 4, Photo = 2
                        plogis(pred_betas$b0[i]+pred_betas$b1[i]*xa_std_min+pred_betas$b9[i]+pred_betas$b11[i]),# Treatment = 1, Photo = 2
                        plogis(pred_betas$b0[i]+pred_betas$b1[i]*xa_std_min+pred_betas$b9[i]+pred_betas$b6[i]+pred_betas$b11[i]),# Treatment = 1, Temp = 1, Photo = 2
                        plogis(pred_betas$b0[i]+pred_betas$b1[i]*xa_std_min+pred_betas$b9[i]+pred_betas$b7[i]+pred_betas$b11[i]),# Treatment = 1, Temp = 2, Photo = 2
                        plogis(pred_betas$b0[i]+pred_betas$b1[i]*xa_std_min+pred_betas$b9[i]+pred_betas$b8[i]+pred_betas$b11[i]))# Treatment = 1, Temp = 4, Photo = 2
  
  media1_pred[i] <- mean(plogis(pred_betas$b0[i]+pred_betas$b1[i]*xa_std_min+pred_betas$b4[i]),# Media = 2
                         plogis(pred_betas$b0[i]+pred_betas$b1[i]*xa_std_min+pred_betas$b6[i]+pred_betas$b4[i]),# Temp = 1, Media = 2
                         plogis(pred_betas$b0[i]+pred_betas$b1[i]*xa_std_min+pred_betas$b7[i]+pred_betas$b4[i]),# Temp = 2, Media = 2
                         plogis(pred_betas$b0[i]+pred_betas$b1[i]*xa_std_min+pred_betas$b8[i]+pred_betas$b4[i]),# Temp = 4, Media = 2
                         plogis(pred_betas$b0[i]+pred_betas$b1[i]*xa_std_min+pred_betas$b9[i]+pred_betas$b4[i]),# Treatment = 1, Media = 2
                         plogis(pred_betas$b0[i]+pred_betas$b1[i]*xa_std_min+pred_betas$b9[i]+pred_betas$b6[i]+pred_betas$b4[i]),# Treatment = 1, Temp = 1,Media = 2
                         plogis(pred_betas$b0[i]+pred_betas$b1[i]*xa_std_min+pred_betas$b9[i]+pred_betas$b7[i]+pred_betas$b4[i]),# Treatment = 1, Temp = 2,Media = 2
                         plogis(pred_betas$b0[i]+pred_betas$b1[i]*xa_std_min+pred_betas$b9[i]+pred_betas$b8[i]+pred_betas$b4[i]),# Treatment = 1, Temp = 4,Media = 2
                         # photoperiod = 1
                         plogis(pred_betas$b0[i]+pred_betas$b1[i]*xa_std_min+pred_betas$b10[i]+pred_betas$b4[i]),# Photo = 1, Media = 2
                         plogis(pred_betas$b0[i]+pred_betas$b1[i]*xa_std_min+pred_betas$b6[i]+pred_betas$b10[i]+pred_betas$b4[i]),# Temp = 1, Photo = 1, Media = 2
                         plogis(pred_betas$b0[i]+pred_betas$b1[i]*xa_std_min+pred_betas$b7[i]+pred_betas$b10[i]+pred_betas$b4[i]),# Temp = 2, Photo = 1, Media = 2
                         plogis(pred_betas$b0[i]+pred_betas$b1[i]*xa_std_min+pred_betas$b8[i]+pred_betas$b10[i]+pred_betas$b4[i]),# Temp = 4, Photo = 1, Media = 2
                         plogis(pred_betas$b0[i]+pred_betas$b1[i]*xa_std_min+pred_betas$b9[i]+pred_betas$b10[i]+pred_betas$b4[i]),# Treatment = 1, Photo = 1, Media = 2
                         plogis(pred_betas$b0[i]+pred_betas$b1[i]*xa_std_min+pred_betas$b9[i]+pred_betas$b6[i]+pred_betas$b10[i]+pred_betas$b4[i]),# Treatment = 1, Temp = 1, Photo = 1, Media = 2
                         plogis(pred_betas$b0[i]+pred_betas$b1[i]*xa_std_min+pred_betas$b9[i]+pred_betas$b7[i]+pred_betas$b10[i]+pred_betas$b4[i]),# Treatment = 1, Temp = 2, Photo = 1, Media = 2
                         plogis(pred_betas$b0[i]+pred_betas$b1[i]*xa_std_min+pred_betas$b9[i]+pred_betas$b8[i]+pred_betas$b10[i]+pred_betas$b4[i]),# Treatment = 1, Temp = 4, Photo = 1, Media = 2
                         # photoperiod = 2
                         plogis(pred_betas$b0[i]+pred_betas$b1[i]*xa_std_min+pred_betas$b11[i]+pred_betas$b4[i]),# Photo = 2, Media = 2
                         plogis(pred_betas$b0[i]+pred_betas$b1[i]*xa_std_min+pred_betas$b6[i]+pred_betas$b11[i]+pred_betas$b4[i]),# Temp = 1, Photo = 2, Media = 2
                         plogis(pred_betas$b0[i]+pred_betas$b1[i]*xa_std_min+pred_betas$b7[i]+pred_betas$b11[i]+pred_betas$b4[i]),# Temp = 2, Photo = 2, Media = 2
                         plogis(pred_betas$b0[i]+pred_betas$b1[i]*xa_std_min+pred_betas$b8[i]+pred_betas$b11[i]+pred_betas$b4[i]),# Temp = 4, Photo = 2, Media = 2
                         plogis(pred_betas$b0[i]+pred_betas$b1[i]*xa_std_min+pred_betas$b9[i]+pred_betas$b11[i]+pred_betas$b4[i]),# Treatment = 1, Photo = 2, Media = 2
                         plogis(pred_betas$b0[i]+pred_betas$b1[i]*xa_std_min+pred_betas$b9[i]+pred_betas$b6[i]+pred_betas$b11[i]+pred_betas$b4[i]),# Treatment = 1, Temp = 1, Photo = 2, Media = 2
                         plogis(pred_betas$b0[i]+pred_betas$b1[i]*xa_std_min+pred_betas$b9[i]+pred_betas$b7[i]+pred_betas$b11[i]+pred_betas$b4[i]),# Treatment = 1, Temp = 2, Photo = 2, Media = 2
                         plogis(pred_betas$b0[i]+pred_betas$b1[i]*xa_std_min+pred_betas$b9[i]+pred_betas$b8[i]+pred_betas$b11[i]+pred_betas$b4[i]))# Treatment = 1, Temp = 4, Photo = 2, Media = 2
  
  media2_pred[i] <- mean(plogis(pred_betas$b0[i]+pred_betas$b1[i]*xa_std_min+pred_betas$b5[i]), # Media = 3
                         plogis(pred_betas$b0[i]+pred_betas$b1[i]*xa_std_min+pred_betas$b6[i]+pred_betas$b5[i]), # Temp = 1, Media = 3
                         plogis(pred_betas$b0[i]+pred_betas$b1[i]*xa_std_min+pred_betas$b7[i]+pred_betas$b5[i]),# Temp = 2, Media = 3
                         plogis(pred_betas$b0[i]+pred_betas$b1[i]*xa_std_min+pred_betas$b8[i]+pred_betas$b5[i]),# Temp = 4, Media = 3
                         plogis(pred_betas$b0[i]+pred_betas$b1[i]*xa_std_min+pred_betas$b9[i]+pred_betas$b5[i]),# Treatment = 1, Media = 3
                         plogis(pred_betas$b0[i]+pred_betas$b1[i]*xa_std_min+pred_betas$b9[i]+pred_betas$b6[i]+pred_betas$b5[i]),# Treatment = 1, Temp = 1, Media = 3
                         plogis(pred_betas$b0[i]+pred_betas$b1[i]*xa_std_min+pred_betas$b9[i]+pred_betas$b7[i]+pred_betas$b5[i]),# Treatment = 1, Temp = 2, Media = 3
                         plogis(pred_betas$b0[i]+pred_betas$b1[i]*xa_std_min+pred_betas$b9[i]+pred_betas$b8[i]+pred_betas$b5[i]),# Treatment = 1, Temp = 4, Media = 3
                         # photoperiod = 1
                         plogis(pred_betas$b0[i]+pred_betas$b1[i]*xa_std_min+pred_betas$b10[i]+pred_betas$b5[i]),# Photo = 1, Media = 3
                         plogis(pred_betas$b0[i]+pred_betas$b1[i]*xa_std_min+pred_betas$b6[i]+pred_betas$b10[i]+pred_betas$b5[i]),# Temp = 1, Photo = 1, Media = 3
                         plogis(pred_betas$b0[i]+pred_betas$b1[i]*xa_std_min+pred_betas$b7[i]+pred_betas$b10[i]+pred_betas$b5[i]),# Temp = 2, Photo = 1, Media = 3
                         plogis(pred_betas$b0[i]+pred_betas$b1[i]*xa_std_min+pred_betas$b8[i]+pred_betas$b10[i]+pred_betas$b5[i]),# Temp = 4, Photo = 1, Media = 3
                         plogis(pred_betas$b0[i]+pred_betas$b1[i]*xa_std_min+pred_betas$b9[i]+pred_betas$b10[i]+pred_betas$b5[i]),# Treatment = 1, Photo = 1, Media = 3
                         plogis(pred_betas$b0[i]+pred_betas$b1[i]*xa_std_min+pred_betas$b9[i]+pred_betas$b6[i]+pred_betas$b10[i]+pred_betas$b5[i]),# Treatment = 1,Temp = 1, Photo = 1, Media = 3
                         plogis(pred_betas$b0[i]+pred_betas$b1[i]*xa_std_min+pred_betas$b9[i]+pred_betas$b7[i]+pred_betas$b10[i]+pred_betas$b5[i]),# Treatment = 1,Temp = 2, Photo = 1, Media = 3
                         plogis(pred_betas$b0[i]+pred_betas$b1[i]*xa_std_min+pred_betas$b9[i]+pred_betas$b8[i]+pred_betas$b10[i]+pred_betas$b5[i]),# Treatment = 1,Temp = 4, Photo = 1, Media = 3
                         # photoperiod = 2
                         plogis(pred_betas$b0[i]+pred_betas$b1[i]*xa_std_min+pred_betas$b11[i]+pred_betas$b5[i]),# Photo = 2, Media = 3
                         plogis(pred_betas$b0[i]+pred_betas$b1[i]*xa_std_min+pred_betas$b6[i]+pred_betas$b11[i]+pred_betas$b5[i]),# Temp = 1, Photo = 2, Media = 3
                         plogis(pred_betas$b0[i]+pred_betas$b1[i]*xa_std_min+pred_betas$b7[i]+pred_betas$b11[i]+pred_betas$b5[i]),# Temp = 2, Photo = 2, Media = 3
                         plogis(pred_betas$b0[i]+pred_betas$b1[i]*xa_std_min+pred_betas$b8[i]+pred_betas$b11[i]+pred_betas$b5[i]),# Temp = 4, Photo = 2, Media = 3
                         plogis(pred_betas$b0[i]+pred_betas$b1[i]*xa_std_min+pred_betas$b9[i]+pred_betas$b11[i]+pred_betas$b5[i]),# Treatment = 1, Photo = 2, Media = 3
                         plogis(pred_betas$b0[i]+pred_betas$b1[i]*xa_std_min+pred_betas$b9[i]+pred_betas$b6[i]+pred_betas$b11[i]+pred_betas$b5[i]),# Treatment = 1, Temp = 1, Photo = 2, Media = 3
                         plogis(pred_betas$b0[i]+pred_betas$b1[i]*xa_std_min+pred_betas$b9[i]+pred_betas$b7[i]+pred_betas$b11[i]+pred_betas$b5[i]),# Treatment = 1, Temp = 2, Photo = 2, Media = 3
                         plogis(pred_betas$b0[i]+pred_betas$b1[i]*xa_std_min+pred_betas$b9[i]+pred_betas$b8[i]+pred_betas$b11[i]+pred_betas$b5[i]))# Treatment = 1, Temp = 4, Photo = 2, Media = 3
  

}

# Calculate average values across other treatments except for pre-treatment
for(i in 1:nrow(pred_betas)){
  treatment0_pred[i] <- mean(plogis(pred_betas$b0[i]+pred_betas$b1[i]*xa_std_min),
                         plogis(pred_betas$b0[i]+pred_betas$b1[i]*xa_std_min+pred_betas$b6[i]), # Temp = 1
                         plogis(pred_betas$b0[i]+pred_betas$b1[i]*xa_std_min+pred_betas$b7[i]), # Temp = 2
                         plogis(pred_betas$b0[i]+pred_betas$b1[i]*xa_std_min+pred_betas$b8[i]),# Temp = 4
                         plogis(pred_betas$b0[i]+pred_betas$b1[i]*xa_std_min+pred_betas$b4[i]), # Media = 2
                         plogis(pred_betas$b0[i]+pred_betas$b1[i]*xa_std_min+pred_betas$b4[i]+pred_betas$b6[i]), # Media = 2, Temp = 1
                         plogis(pred_betas$b0[i]+pred_betas$b1[i]*xa_std_min+pred_betas$b4[i]+pred_betas$b7[i]),# Media = 2, Temp = 2
                         plogis(pred_betas$b0[i]+pred_betas$b1[i]*xa_std_min+pred_betas$b4[i]+pred_betas$b8[i]),# Media = 2, Temp = 4
                         plogis(pred_betas$b0[i]+pred_betas$b1[i]*xa_std_min+pred_betas$b5[i]),# Media = 3
                         plogis(pred_betas$b0[i]+pred_betas$b1[i]*xa_std_min+pred_betas$b5[i]+pred_betas$b6[i]),# Media = 3, Temp = 1
                         plogis(pred_betas$b0[i]+pred_betas$b1[i]*xa_std_min+pred_betas$b5[i]+pred_betas$b7[i]),# Media = 3, Temp = 2
                         plogis(pred_betas$b0[i]+pred_betas$b1[i]*xa_std_min+pred_betas$b5[i]+pred_betas$b8[i]),# Media = 3, Temp = 4
                         # photoperiod = 1
                         plogis(pred_betas$b0[i]+pred_betas$b1[i]*xa_std_min+pred_betas$b10[i]),# Photo = 1
                         plogis(pred_betas$b0[i]+pred_betas$b1[i]*xa_std_min+pred_betas$b6[i]+pred_betas$b10[i]),# Temp = 1, Photo = 1
                         plogis(pred_betas$b0[i]+pred_betas$b1[i]*xa_std_min+pred_betas$b7[i]+pred_betas$b10[i]),# Temp = 2, Photo = 1
                         plogis(pred_betas$b0[i]+pred_betas$b1[i]*xa_std_min+pred_betas$b8[i]+pred_betas$b10[i]),# Temp = 4, Photo = 1
                         plogis(pred_betas$b0[i]+pred_betas$b1[i]*xa_std_min+pred_betas$b4[i]+pred_betas$b10[i]),# Media = 2, Photo = 1
                         plogis(pred_betas$b0[i]+pred_betas$b1[i]*xa_std_min+pred_betas$b4[i]+pred_betas$b6[i]+pred_betas$b10[i]),# Media = 2,Temp = 1, Photo = 1
                         plogis(pred_betas$b0[i]+pred_betas$b1[i]*xa_std_min+pred_betas$b4[i]+pred_betas$b7[i]+pred_betas$b10[i]),# Media = 2,Temp = 2, Photo = 1
                         plogis(pred_betas$b0[i]+pred_betas$b1[i]*xa_std_min+pred_betas$b4[i]+pred_betas$b8[i]+pred_betas$b10[i]),# Media = 2,Temp = 4, Photo = 1
                         plogis(pred_betas$b0[i]+pred_betas$b1[i]*xa_std_min+pred_betas$b5[i]+pred_betas$b10[i]),# Media = 3, Photo = 1
                         plogis(pred_betas$b0[i]+pred_betas$b1[i]*xa_std_min+pred_betas$b5[i]+pred_betas$b6[i]+pred_betas$b10[i]),# Media = 3, Temp = 1, Photo = 1
                         plogis(pred_betas$b0[i]+pred_betas$b1[i]*xa_std_min+pred_betas$b5[i]+pred_betas$b7[i]+pred_betas$b10[i]),# Media = 3, Temp = 2, Photo = 1
                         plogis(pred_betas$b0[i]+pred_betas$b1[i]*xa_std_min+pred_betas$b5[i]+pred_betas$b8[i]+pred_betas$b10[i]),# Media = 3, Temp = 3, Photo = 1
                         # photoperiod = 2
                         plogis(pred_betas$b0[i]+pred_betas$b1[i]*xa_std_min+pred_betas$b11[i]), # Photo = 2
                         plogis(pred_betas$b0[i]+pred_betas$b1[i]*xa_std_min+pred_betas$b6[i]+pred_betas$b11[i]), # Temp = 1, Photo = 2
                         plogis(pred_betas$b0[i]+pred_betas$b1[i]*xa_std_min+pred_betas$b7[i]+pred_betas$b11[i]),# Temp = 2, Photo = 2
                         plogis(pred_betas$b0[i]+pred_betas$b1[i]*xa_std_min+pred_betas$b8[i]+pred_betas$b11[i]),# Temp = 3, Photo = 2
                         plogis(pred_betas$b0[i]+pred_betas$b1[i]*xa_std_min+pred_betas$b4[i]+pred_betas$b11[i]),# Media = 2, Photo = 2
                         plogis(pred_betas$b0[i]+pred_betas$b1[i]*xa_std_min+pred_betas$b4[i]+pred_betas$b6[i]+pred_betas$b11[i]), # Media = 2, Temp = 1, Photo = 2
                         plogis(pred_betas$b0[i]+pred_betas$b1[i]*xa_std_min+pred_betas$b4[i]+pred_betas$b7[i]+pred_betas$b11[i]),# Media = 2, Temp = 2, Photo = 2
                         plogis(pred_betas$b0[i]+pred_betas$b1[i]*xa_std_min+pred_betas$b4[i]+pred_betas$b8[i]+pred_betas$b11[i]),# Media = 2, Temp = 4, Photo = 2
                         plogis(pred_betas$b0[i]+pred_betas$b1[i]*xa_std_min+pred_betas$b5[i]+pred_betas$b11[i]),# Media = 3, Photo = 2
                         plogis(pred_betas$b0[i]+pred_betas$b1[i]*xa_std_min+pred_betas$b5[i]+pred_betas$b6[i]+pred_betas$b11[i]),# Media = 3, Temp = 1, Photo = 2
                         plogis(pred_betas$b0[i]+pred_betas$b1[i]*xa_std_min+pred_betas$b5[i]+pred_betas$b7[i]+pred_betas$b11[i]),# Media = 3, Temp = 2, Photo = 2
                         plogis(pred_betas$b0[i]+pred_betas$b1[i]*xa_std_min+pred_betas$b5[i]+pred_betas$b8[i]+pred_betas$b11[i]))# Media = 3, Temp = 4, Photo = 2
  
  treatment1_pred[i] <- mean(plogis(pred_betas$b0[i]+pred_betas$b1[i]*xa_std_min+pred_betas$b9[i]), # Treatment = 1
                             plogis(pred_betas$b0[i]+pred_betas$b1[i]*xa_std_min+pred_betas$b6[i]+pred_betas$b9[i]),# Temp = 1, Treatment = 1
                             plogis(pred_betas$b0[i]+pred_betas$b1[i]*xa_std_min+pred_betas$b7[i]+pred_betas$b9[i]),# Temp = 2, Treatment = 1
                             plogis(pred_betas$b0[i]+pred_betas$b1[i]*xa_std_min+pred_betas$b8[i]+pred_betas$b9[i]),# Temp = 4, Treatment = 1
                             plogis(pred_betas$b0[i]+pred_betas$b1[i]*xa_std_min+pred_betas$b4[i]+pred_betas$b9[i]),# Media = 2, Treatment = 1
                             plogis(pred_betas$b0[i]+pred_betas$b1[i]*xa_std_min+pred_betas$b4[i]+pred_betas$b6[i]+pred_betas$b9[i]),# Media = 2, Temp = 1, Treatment = 1
                             plogis(pred_betas$b0[i]+pred_betas$b1[i]*xa_std_min+pred_betas$b4[i]+pred_betas$b7[i]+pred_betas$b9[i]),# Media = 2, Temp = 2, Treatment = 1
                             plogis(pred_betas$b0[i]+pred_betas$b1[i]*xa_std_min+pred_betas$b4[i]+pred_betas$b8[i]+pred_betas$b9[i]),# Media = 2, Temp = 4, Treatment = 1
                             plogis(pred_betas$b0[i]+pred_betas$b1[i]*xa_std_min+pred_betas$b5[i]+pred_betas$b9[i]),# Media = 3, Treatment = 1
                             plogis(pred_betas$b0[i]+pred_betas$b1[i]*xa_std_min+pred_betas$b5[i]+pred_betas$b6[i]+pred_betas$b9[i]),# Media = 3, Temp = 1, Treatment = 1
                             plogis(pred_betas$b0[i]+pred_betas$b1[i]*xa_std_min+pred_betas$b5[i]+pred_betas$b7[i]+pred_betas$b9[i]),# Media = 3, Temp = 2, Treatment = 1
                             plogis(pred_betas$b0[i]+pred_betas$b1[i]*xa_std_min+pred_betas$b5[i]+pred_betas$b8[i]+pred_betas$b9[i]),# Media = 3, Temp = 4, Treatment = 1
                             # photoperiod = 1
                             plogis(pred_betas$b0[i]+pred_betas$b1[i]*xa_std_min+pred_betas$b10[i]+pred_betas$b9[i]),# Photo = 1, Treatment = 1
                             plogis(pred_betas$b0[i]+pred_betas$b1[i]*xa_std_min+pred_betas$b6[i]+pred_betas$b10[i]+pred_betas$b9[i]),# Temp = 1, Photo = 1, Treatment = 1
                             plogis(pred_betas$b0[i]+pred_betas$b1[i]*xa_std_min+pred_betas$b7[i]+pred_betas$b10[i]+pred_betas$b9[i]),# Temp = 2, Photo = 1, Treatment = 1
                             plogis(pred_betas$b0[i]+pred_betas$b1[i]*xa_std_min+pred_betas$b8[i]+pred_betas$b10[i]+pred_betas$b9[i]),# Temp = 4, Photo = 1, Treatment = 1
                             plogis(pred_betas$b0[i]+pred_betas$b1[i]*xa_std_min+pred_betas$b4[i]+pred_betas$b10[i]+pred_betas$b9[i]),# Media = 2, Photo = 1, Treatment = 1
                             plogis(pred_betas$b0[i]+pred_betas$b1[i]*xa_std_min+pred_betas$b4[i]+pred_betas$b6[i]+pred_betas$b10[i]+pred_betas$b9[i]),# Media = 2, Temp = 1, Photo = 1, Treatment = 1
                             plogis(pred_betas$b0[i]+pred_betas$b1[i]*xa_std_min+pred_betas$b4[i]+pred_betas$b7[i]+pred_betas$b10[i]+pred_betas$b9[i]),# Media = 2, Temp = 2, Photo = 1, Treatment = 1
                             plogis(pred_betas$b0[i]+pred_betas$b1[i]*xa_std_min+pred_betas$b4[i]+pred_betas$b8[i]+pred_betas$b10[i]+pred_betas$b9[i]),# Media = 2, Temp = 4, Photo = 1, Treatment = 1
                             plogis(pred_betas$b0[i]+pred_betas$b1[i]*xa_std_min+pred_betas$b5[i]+pred_betas$b10[i]+pred_betas$b9[i]),# Media = 3, Photo = 1, Treatment = 1
                             plogis(pred_betas$b0[i]+pred_betas$b1[i]*xa_std_min+pred_betas$b5[i]+pred_betas$b6[i]+pred_betas$b10[i]+pred_betas$b9[i]),# Media = 3, Temp = 1, Photo = 1, Treatment = 1
                             plogis(pred_betas$b0[i]+pred_betas$b1[i]*xa_std_min+pred_betas$b5[i]+pred_betas$b7[i]+pred_betas$b10[i]+pred_betas$b9[i]),# Media = 3, Temp = 2, Photo = 1, Treatment = 1
                             plogis(pred_betas$b0[i]+pred_betas$b1[i]*xa_std_min+pred_betas$b5[i]+pred_betas$b8[i]+pred_betas$b10[i]+pred_betas$b9[i]),# Media = 3, Temp = 4, Photo = 1, Treatment = 1
                             # photoperiod = 2
                             plogis(pred_betas$b0[i]+pred_betas$b1[i]*xa_std_min+pred_betas$b11[i]+pred_betas$b9[i]),# Photo = 2, Treatment = 2
                             plogis(pred_betas$b0[i]+pred_betas$b1[i]*xa_std_min+pred_betas$b6[i]+pred_betas$b11[i]+pred_betas$b9[i]), # Temp = 1, Photo = 2, Treatment = 1
                             plogis(pred_betas$b0[i]+pred_betas$b1[i]*xa_std_min+pred_betas$b7[i]+pred_betas$b11[i]+pred_betas$b9[i]),# Temp = 2, Photo = 2, Treatment = 1
                             plogis(pred_betas$b0[i]+pred_betas$b1[i]*xa_std_min+pred_betas$b8[i]+pred_betas$b11[i]+pred_betas$b9[i]),# Temp = 4, Photo = 2, Treatment = 1
                             plogis(pred_betas$b0[i]+pred_betas$b1[i]*xa_std_min+pred_betas$b4[i]+pred_betas$b11[i]+pred_betas$b9[i]),# Media = 2, Photo = 2, Treatment = 1
                             plogis(pred_betas$b0[i]+pred_betas$b1[i]*xa_std_min+pred_betas$b4[i]+pred_betas$b6[i]+pred_betas$b11[i]+pred_betas$b9[i]),# Media = 2, Temp = 1, Photo = 2, Treatment = 1
                             plogis(pred_betas$b0[i]+pred_betas$b1[i]*xa_std_min+pred_betas$b4[i]+pred_betas$b7[i]+pred_betas$b11[i]+pred_betas$b9[i]),# Media = 2, Temp = 2, Photo = 2, Treatment = 1
                             plogis(pred_betas$b0[i]+pred_betas$b1[i]*xa_std_min+pred_betas$b4[i]+pred_betas$b8[i]+pred_betas$b11[i]+pred_betas$b9[i]),# Media = 2, Temp = 4, Photo = 2, Treatment = 1
                             plogis(pred_betas$b0[i]+pred_betas$b1[i]*xa_std_min+pred_betas$b5[i]+pred_betas$b11[i]+pred_betas$b9[i]),# Media = 3, Photo = 2, Treatment = 1
                             plogis(pred_betas$b0[i]+pred_betas$b1[i]*xa_std_min+pred_betas$b5[i]+pred_betas$b6[i]+pred_betas$b11[i]+pred_betas$b9[i]),# Media = 3, Temp = 1, Photo = 2, Treatment = 1
                             plogis(pred_betas$b0[i]+pred_betas$b1[i]*xa_std_min+pred_betas$b5[i]+pred_betas$b7[i]+pred_betas$b11[i]+pred_betas$b9[i]),# Media = 3, Temp = 2, Photo = 2, Treatment = 1
                             plogis(pred_betas$b0[i]+pred_betas$b1[i]*xa_std_min+pred_betas$b5[i]+pred_betas$b8[i]+pred_betas$b11[i]+pred_betas$b9[i]))# Media = 3, Temp = 4, Photo = 2, Treatment = 1
  
  
}

# Calculate average values across other treatments except for photoperiod
for(i in 1:nrow(pred_betas)){
  photoperiod3_pred[i] <- mean(plogis(pred_betas$b0[i]+pred_betas$b1[i]*xa_std_min),
                         plogis(pred_betas$b0[i]+pred_betas$b1[i]*xa_std_min+pred_betas$b6[i]),# Temp = 1
                         plogis(pred_betas$b0[i]+pred_betas$b1[i]*xa_std_min+pred_betas$b7[i]),# Temp = 2
                         plogis(pred_betas$b0[i]+pred_betas$b1[i]*xa_std_min+pred_betas$b8[i]),# Temp = 4
                         plogis(pred_betas$b0[i]+pred_betas$b1[i]*xa_std_min+pred_betas$b9[i]),# Treatment = 1
                         plogis(pred_betas$b0[i]+pred_betas$b1[i]*xa_std_min+pred_betas$b9[i]+pred_betas$b6[i]),# Treatment = 1, Temp = 1
                         plogis(pred_betas$b0[i]+pred_betas$b1[i]*xa_std_min+pred_betas$b9[i]+pred_betas$b7[i]),# Treatment = 1, Temp = 2
                         plogis(pred_betas$b0[i]+pred_betas$b1[i]*xa_std_min+pred_betas$b9[i]+pred_betas$b8[i]),# Treatment = 1, Temp = 4
                         # media = 2
                         plogis(pred_betas$b0[i]+pred_betas$b1[i]*xa_std_min+pred_betas$b4[i]),# Media = 2
                         plogis(pred_betas$b0[i]+pred_betas$b1[i]*xa_std_min+pred_betas$b6[i]+pred_betas$b4[i]), # Temp = 1, Media = 2
                         plogis(pred_betas$b0[i]+pred_betas$b1[i]*xa_std_min+pred_betas$b7[i]+pred_betas$b4[i]), # Temp = 2, Media = 2
                         plogis(pred_betas$b0[i]+pred_betas$b1[i]*xa_std_min+pred_betas$b8[i]+pred_betas$b4[i]), # Temp = 4, Media = 2
                         plogis(pred_betas$b0[i]+pred_betas$b1[i]*xa_std_min+pred_betas$b9[i]+pred_betas$b4[i]), # Treatment = 1, Media = 2
                         plogis(pred_betas$b0[i]+pred_betas$b1[i]*xa_std_min+pred_betas$b9[i]+pred_betas$b6[i]+pred_betas$b4[i]),# Treatment = 1, Temp = 1, Media = 2
                         plogis(pred_betas$b0[i]+pred_betas$b1[i]*xa_std_min+pred_betas$b9[i]+pred_betas$b7[i]+pred_betas$b4[i]),# Treatment = 1, Temp = 2, Media = 2
                         plogis(pred_betas$b0[i]+pred_betas$b1[i]*xa_std_min+pred_betas$b9[i]+pred_betas$b8[i]+pred_betas$b4[i]),# Treatment = 1, Temp = 4, Media = 2
                         # media = 3
                         plogis(pred_betas$b0[i]+pred_betas$b1[i]*xa_std_min+pred_betas$b5[i]), # Media = 3
                         plogis(pred_betas$b0[i]+pred_betas$b1[i]*xa_std_min+pred_betas$b6[i]+pred_betas$b5[i]), # Temp = 1, Media = 3
                         plogis(pred_betas$b0[i]+pred_betas$b1[i]*xa_std_min+pred_betas$b7[i]+pred_betas$b5[i]),# Temp = 2, Media = 3
                         plogis(pred_betas$b0[i]+pred_betas$b1[i]*xa_std_min+pred_betas$b8[i]+pred_betas$b5[i]),# Temp = 4, Media = 3
                         plogis(pred_betas$b0[i]+pred_betas$b1[i]*xa_std_min+pred_betas$b9[i]+pred_betas$b5[i]),# Treatment = 1, Media = 3
                         plogis(pred_betas$b0[i]+pred_betas$b1[i]*xa_std_min+pred_betas$b9[i]+pred_betas$b6[i]+pred_betas$b5[i]),# Treatment = 1, Temp = 1, Media = 3
                         plogis(pred_betas$b0[i]+pred_betas$b1[i]*xa_std_min+pred_betas$b9[i]+pred_betas$b7[i]+pred_betas$b5[i]),# Treatment = 1, Temp = 2, Media = 3
                         plogis(pred_betas$b0[i]+pred_betas$b1[i]*xa_std_min+pred_betas$b9[i]+pred_betas$b8[i]+pred_betas$b5[i]))# Treatment = 1, Temp = 4, Media = 3
  
  photoperiod1_pred[i] <- mean(plogis(pred_betas$b0[i]+pred_betas$b1[i]*xa_std_min+pred_betas$b10[i]), # Photo = 1
                         plogis(pred_betas$b0[i]+pred_betas$b1[i]*xa_std_min+pred_betas$b6[i]+pred_betas$b10[i]), # Temp = 1, Photo = 1
                         plogis(pred_betas$b0[i]+pred_betas$b1[i]*xa_std_min+pred_betas$b7[i]+pred_betas$b10[i]),# Temp = 2, Photo = 1
                         plogis(pred_betas$b0[i]+pred_betas$b1[i]*xa_std_min+pred_betas$b8[i]+pred_betas$b10[i]),# Temp = 4, Photo = 1
                         plogis(pred_betas$b0[i]+pred_betas$b1[i]*xa_std_min+pred_betas$b9[i]+pred_betas$b10[i]),# Treatment = 1, Photo = 1
                         plogis(pred_betas$b0[i]+pred_betas$b1[i]*xa_std_min+pred_betas$b9[i]+pred_betas$b6[i]+pred_betas$b10[i]),# Treatment = 1, Temp = 1, Photo = 1
                         plogis(pred_betas$b0[i]+pred_betas$b1[i]*xa_std_min+pred_betas$b9[i]+pred_betas$b7[i]+pred_betas$b10[i]),# Treatment = 1, Temp = 2, Photo = 1
                         plogis(pred_betas$b0[i]+pred_betas$b1[i]*xa_std_min+pred_betas$b9[i]+pred_betas$b8[i]+pred_betas$b10[i]),# Treatment = 1, Temp = 4, Photo = 1
                         # media = 2
                         plogis(pred_betas$b0[i]+pred_betas$b1[i]*xa_std_min+pred_betas$b4[i]+pred_betas$b10[i]),# Media = 2, Photo = 1
                         plogis(pred_betas$b0[i]+pred_betas$b1[i]*xa_std_min+pred_betas$b6[i]+pred_betas$b4[i]+pred_betas$b10[i]),# Temp = 1, Media = 2, Photo = 1
                         plogis(pred_betas$b0[i]+pred_betas$b1[i]*xa_std_min+pred_betas$b7[i]+pred_betas$b4[i]+pred_betas$b10[i]),# Temp = 2, Media = 2, Photo = 1
                         plogis(pred_betas$b0[i]+pred_betas$b1[i]*xa_std_min+pred_betas$b8[i]+pred_betas$b4[i]+pred_betas$b10[i]),# Temp = 4, Media = 2, Photo = 1
                         plogis(pred_betas$b0[i]+pred_betas$b1[i]*xa_std_min+pred_betas$b9[i]+pred_betas$b4[i]+pred_betas$b10[i]),# Treatment = 1, Media = 2, Photo = 1
                         plogis(pred_betas$b0[i]+pred_betas$b1[i]*xa_std_min+pred_betas$b9[i]+pred_betas$b6[i]+pred_betas$b4[i]+pred_betas$b10[i]),# Treatment = 1, Temp = 1, Media = 2, Photo = 1
                         plogis(pred_betas$b0[i]+pred_betas$b1[i]*xa_std_min+pred_betas$b9[i]+pred_betas$b7[i]+pred_betas$b4[i]+pred_betas$b10[i]),# Treatment = 1, Temp = 2, Media = 2, Photo = 1
                         plogis(pred_betas$b0[i]+pred_betas$b1[i]*xa_std_min+pred_betas$b9[i]+pred_betas$b8[i]+pred_betas$b4[i]+pred_betas$b10[i]),# Treatment = 1, Temp = 4, Media = 2, Photo = 1
                         # media = 3
                         plogis(pred_betas$b0[i]+pred_betas$b1[i]*xa_std_min+pred_betas$b5[i]+pred_betas$b10[i]), # Media = 3, Photo = 1
                         plogis(pred_betas$b0[i]+pred_betas$b1[i]*xa_std_min+pred_betas$b6[i]+pred_betas$b5[i]+pred_betas$b10[i]),# Temp = 1, Media = 3, Photo = 1
                         plogis(pred_betas$b0[i]+pred_betas$b1[i]*xa_std_min+pred_betas$b7[i]+pred_betas$b5[i]+pred_betas$b10[i]),# Temp = 2, Media = 3, Photo = 1
                         plogis(pred_betas$b0[i]+pred_betas$b1[i]*xa_std_min+pred_betas$b8[i]+pred_betas$b5[i]+pred_betas$b10[i]),# Temp = 4, Media = 3, Photo = 1
                         plogis(pred_betas$b0[i]+pred_betas$b1[i]*xa_std_min+pred_betas$b9[i]+pred_betas$b5[i]+pred_betas$b10[i]),# Treatment = 1, Media = 3, Photo = 1
                         plogis(pred_betas$b0[i]+pred_betas$b1[i]*xa_std_min+pred_betas$b9[i]+pred_betas$b6[i]+pred_betas$b5[i]+pred_betas$b10[i]),# Treatment = 1, Temp = 1, Media = 3, Photo = 1
                         plogis(pred_betas$b0[i]+pred_betas$b1[i]*xa_std_min+pred_betas$b9[i]+pred_betas$b7[i]+pred_betas$b5[i]+pred_betas$b10[i]),# Treatment = 1, Temp = 2, Media = 3, Photo = 1
                         plogis(pred_betas$b0[i]+pred_betas$b1[i]*xa_std_min+pred_betas$b9[i]+pred_betas$b8[i]+pred_betas$b5[i]+pred_betas$b10[i]))# Treatment = 1, Temp = 4, Media = 3, Photo = 1
  
  photoperiod2_pred[i] <- mean(plogis(pred_betas$b0[i]+pred_betas$b1[i]*xa_std_min+pred_betas$b11[i]), # Photo = 2
                         plogis(pred_betas$b0[i]+pred_betas$b1[i]*xa_std_min+pred_betas$b6[i]+pred_betas$b11[i]), # Temp 1, Photo 2
                         plogis(pred_betas$b0[i]+pred_betas$b1[i]*xa_std_min+pred_betas$b7[i]+pred_betas$b11[i]),# Temp 2, Photo 2
                         plogis(pred_betas$b0[i]+pred_betas$b1[i]*xa_std_min+pred_betas$b8[i]+pred_betas$b11[i]),# Temp 4, Photo 2
                         plogis(pred_betas$b0[i]+pred_betas$b1[i]*xa_std_min+pred_betas$b9[i]+pred_betas$b11[i]),# Treatment 1, Photo 2
                         plogis(pred_betas$b0[i]+pred_betas$b1[i]*xa_std_min+pred_betas$b9[i]+pred_betas$b6[i]+pred_betas$b11[i]),# Treatment 1,Temp = 1, Photo 2
                         plogis(pred_betas$b0[i]+pred_betas$b1[i]*xa_std_min+pred_betas$b9[i]+pred_betas$b7[i]+pred_betas$b11[i]),# Treatment 1,Temp = 2, Photo 2
                         plogis(pred_betas$b0[i]+pred_betas$b1[i]*xa_std_min+pred_betas$b9[i]+pred_betas$b8[i]+pred_betas$b11[i]),# Treatment 1,Temp = 4, Photo 2
                         # media = 2
                         plogis(pred_betas$b0[i]+pred_betas$b1[i]*xa_std_min+pred_betas$b4[i]+pred_betas$b11[i]),# Media = 2, Photo = 2
                         plogis(pred_betas$b0[i]+pred_betas$b1[i]*xa_std_min+pred_betas$b6[i]+pred_betas$b4[i]+pred_betas$b11[i]),# Temp = 1, Media = 2, Photo = 2
                         plogis(pred_betas$b0[i]+pred_betas$b1[i]*xa_std_min+pred_betas$b7[i]+pred_betas$b4[i]+pred_betas$b11[i]),# Temp = 2, Media = 2, Photo = 2
                         plogis(pred_betas$b0[i]+pred_betas$b1[i]*xa_std_min+pred_betas$b8[i]+pred_betas$b4[i]+pred_betas$b11[i]),# Temp = 4, Media = 2, Photo = 2
                         plogis(pred_betas$b0[i]+pred_betas$b1[i]*xa_std_min+pred_betas$b9[i]+pred_betas$b4[i]+pred_betas$b11[i]),# Treatment = 1, Media = 2, Photo = 2
                         plogis(pred_betas$b0[i]+pred_betas$b1[i]*xa_std_min+pred_betas$b9[i]+pred_betas$b6[i]+pred_betas$b4[i]+pred_betas$b11[i]),# Treatment = 1,Temp = 1, Media = 2, Photo = 2
                         plogis(pred_betas$b0[i]+pred_betas$b1[i]*xa_std_min+pred_betas$b9[i]+pred_betas$b7[i]+pred_betas$b4[i]+pred_betas$b11[i]),# Treatment = 1,Temp = 2, Media = 2, Photo = 2
                         plogis(pred_betas$b0[i]+pred_betas$b1[i]*xa_std_min+pred_betas$b9[i]+pred_betas$b8[i]+pred_betas$b4[i]+pred_betas$b11[i]),# Treatment = 1,Temp = 4, Media = 2, Photo = 2
                         # media = 3
                         plogis(pred_betas$b0[i]+pred_betas$b1[i]*xa_std_min+pred_betas$b5[i]+pred_betas$b11[i]),# Media = 3, Photo = 2
                         plogis(pred_betas$b0[i]+pred_betas$b1[i]*xa_std_min+pred_betas$b6[i]+pred_betas$b5[i]+pred_betas$b11[i]),# Temp = 1, Media = 3, Photo = 2
                         plogis(pred_betas$b0[i]+pred_betas$b1[i]*xa_std_min+pred_betas$b7[i]+pred_betas$b5[i]+pred_betas$b11[i]),# Temp = 2 Media = 3, Photo = 2
                         plogis(pred_betas$b0[i]+pred_betas$b1[i]*xa_std_min+pred_betas$b8[i]+pred_betas$b5[i]+pred_betas$b11[i]),# Temp = 4, Media = 3, Photo = 2
                         plogis(pred_betas$b0[i]+pred_betas$b1[i]*xa_std_min+pred_betas$b9[i]+pred_betas$b5[i]+pred_betas$b11[i]),# Treatment = 1, Media = 3, Photo = 2
                         plogis(pred_betas$b0[i]+pred_betas$b1[i]*xa_std_min+pred_betas$b9[i]+pred_betas$b6[i]+pred_betas$b5[i]+pred_betas$b11[i]),# Treatment = 1, Temp = 1, Media = 3, Photo = 2
                         plogis(pred_betas$b0[i]+pred_betas$b1[i]*xa_std_min+pred_betas$b9[i]+pred_betas$b7[i]+pred_betas$b5[i]+pred_betas$b11[i]),# Treatment = 1, Temp = 2, Media = 3, Photo = 2
                         plogis(pred_betas$b0[i]+pred_betas$b1[i]*xa_std_min+pred_betas$b9[i]+pred_betas$b8[i]+pred_betas$b5[i]+pred_betas$b11[i]))# Treatment = 1, Temp = 4, Media = 3, Photo = 2
  
  
}

## Create graphics ####

# Create graph for temperature
a <- tibble(value = c(temp1_pred, temp2_pred, temp3_pred, temp4_pred),
       temp = rep(c("30", "20/15", "27/15", "25"), each = length(temp1_pred))) %>%
  ggplot(aes(x = factor(temp), y = value)) +
  geom_violin(draw_quantiles = c(0.025, 0.975), fill = "gray77",  color = "gray27") + 
  stat_summary(fun = median, geom = "point", size = 3, color = "gray27") +
  ylab("P(germination success)") +
  xlab("temperature (°C)") + ylim(0,0.9)

b <- tibble(value = c(media1_pred, media2_pred, media3_pred),
       media = rep(c("sand","growth media","sand + soil"), each = length(media1_pred))) %>%
  ggplot(aes(x = media, y = value)) +
  geom_violin(draw_quantiles = c(0.025, 0.975), fill = "gray77",  color = "gray27") + 
  stat_summary(fun = median, geom = "point", size = 3, color = "gray27") +
  ylab("") +
  xlab("type of media")+ ylim(0,0.9)

c <- tibble(value = c(treatment0_pred, treatment1_pred),
       treatment = rep(c("no", "yes"), each = length(treatment1_pred))) %>%
  ggplot(aes(x = treatment, y = value)) +
  geom_violin(draw_quantiles = c(0.025, 0.975), fill = "gray77",  color = "gray27") + 
  stat_summary(fun = median, geom = "point", size = 3, color = "gray27") + 
  ylab("P(germination success)") +
  xlab("pre-treatment")+ ylim(0,0.9)

d <- tibble(value = c(photoperiod3_pred, photoperiod1_pred, photoperiod2_pred),
            treatment = rep(c("12/12", "15/9", "dark"), each = length(photoperiod3_pred))) %>%
  ggplot(aes(x = treatment, y = value)) +
  geom_violin(draw_quantiles = c(0.025, 0.975), fill = "gray77",  color = "gray27") + 
  stat_summary(fun = median, geom = "point", size = 3, color = "gray27") + 
  ylab("") +
  xlab("photoperiod")+ ylim(0,0.9)

png(here("figs_tables","Fig4.png"), width = 6.5, height = 6.2, res = 300, units = "in")
plot_grid(a,b,c,d,
          nrow = 2,
          align = "h", labels = "auto")
dev.off()

## Calculations for comparisons that are in-text ####
# Comparing temperature regimes (27/15 is best, 25 is worst)
mu_t <- mean(temp3_pred-temp4_pred) 
qt_t <- quantile(temp3_pred - temp4_pred, c(0.025, 0.975))

# Comparing media (sand is the best, growth media is the worst)
mu_m <- mean(media1_pred - media2_pred)
qt_m <- quantile(media1_pred - media2_pred, c(0.025, 0.975))
mean(media1_pred) # 0.3041523 sand
mean(media2_pred) # 0.07290745 growth media

# Comparing pre-treatments
mu_p <-mean(treatment0_pred - treatment1_pred)
qt_p <- quantile(treatment0_pred - treatment1_pred, c(0.025, 0.975))
mean(treatment0_pred) # 0.2591152 no treatment
mean(treatment1_pred) # 0.03845524 treatment

# Comparing photoperiods
mu_ph <-mean(photoperiod1_pred - photoperiod2_pred)
qt_ph <- quantile(photoperiod1_pred - photoperiod2_pred, c(0.025, 0.975))

tibble(condition = c("temp", "media", "pre", "photo"),
       mean = c(mu_t, mu_m, mu_p, mu_ph),
       lower = c(qt_t[1], qt_m[1], qt_p[1], qt_ph[1]),
       upper = c(qt_t[2], qt_m[2], qt_p[2], qt_ph[2]))

# # A tibble: 4 x 4
# condition  mean   lower upper
# <chr>     <dbl>   <dbl> <dbl>
# 1 temp      0.141 -0.0255 0.325
# 2 media     0.231  0.0248 0.452
# 3 pre       0.221  0.125  0.327
# 4 photo     0.110 -0.0346 0.262