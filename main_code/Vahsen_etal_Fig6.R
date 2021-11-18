# Create Figure 6 in Vahsen et al. germination paper
# Also runs the binomial model needed for the tetrazolium test data
# 19 Feb 2021

## Preliminaries ####
# Load libraries
library(tidyverse); library(ggiraphExtra); library(MuMIn); library(here)

# Read in clean tetrazolium data
tetra <- read_csv(here("data","tetrazolium_data.csv"))

## Fit binomial regression model ####

# Fit a full model with Assay_Num (factor) and the interaction between Depth and
# Assay_Num included.

# Set Assay_Num as factor
tetra$Assay_Num <- factor(tetra$Assay_Num)

# Caculate the proportion and store in data frame
tetra$prop <- tetra$Tetra_Viable / tetra$Seeds_Tested

# Fit full model
tetra_model <- glm(prop ~ Assay_Num * Depth, data = tetra, family = "binomial",
                   weights = Seeds_Tested)

# Create alternative models to see which is the best fit for prediction
tetra_model_noint <- update(tetra_model, .~.-Assay_Num:Depth)
tetra_model_noA <- update(tetra_model_noint, .~.-Assay_Num)
tetra_model_noD <- update(tetra_model_noint, .~.-Depth)

MuMIn::AICc(tetra_model)
MuMIn::AICc(tetra_model_noint)
MuMIn::AICc(tetra_model_noA) # Best model -- final model just includes Depth
MuMIn::AICc(tetra_model_noD)

## Create Figure 6 ####
png(here("figs_tables","Fig6.png"), height = 4, width = 2.8, units = "in", res = 300)
tetra %>% 
  ggplot(aes(x = Depth, y = prop)) +
  geom_point() +
  geom_smooth(method = "glm", method.args = list(family = "binomial"),
              aes(weight = Seeds_Tested), fullrange = TRUE, color = "black") + 
  coord_flip() + 
  theme_classic() +
  xlab("seed depth (cm)") +
  ylab("proportion of viable seeds") +
  scale_y_continuous(position = "right") +
  scale_x_reverse()
dev.off()

## Summary statistics for in-text ####  

# Number of total seeds tested
sum(tetra$Seeds_Tested)
# 470

# Total proportion of seeds that had tetrazolium positive
sum(tetra$Tetra_Viable) / sum(tetra$Seeds_Tested)
# 0.1042553