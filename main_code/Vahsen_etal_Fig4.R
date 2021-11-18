# Create Figure 4 in Vahsen et al. germination paper
# Also calculates predicted means and 95% CI for in-text
# 19 Feb 2021

## Preliminaries ####
# Load libraries
library(tidyverse);library(patchwork);library(ggmcmc)
library(EnvStats);library(cowplot); library(here);
library(spatstat)

# Set ggplot theme
theme_set(theme_classic())

# Read in data that was used to fit the models (same as original data file that
# takes out grouping by Assay_Num)
germ_all_final <- read_csv(here("outputs","germ_all_data_forModel.csv"))

# Read in coda objects from all four models
Model3_coda <- read_rds(here("outputs", "Model3_coda_forPlotting.rds"))

## Calculate predicted values for each treatment ####  

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

# Calculate predicted values for each temperature treatment at reference levels
# for other experimental treatments

for(i in 1:nrow(pred_betas)){
  temp3_pred[i] <- pred_betas$b0[i]+pred_betas$b1[i]*xa_std_min
  temp1_pred[i] <- pred_betas$b0[i]+pred_betas$b1[i]*xa_std_min + pred_betas$b6[i]
  temp2_pred[i] <- pred_betas$b0[i]+pred_betas$b1[i]*xa_std_min + pred_betas$b7[i]
  temp4_pred[i] <- pred_betas$b0[i]+pred_betas$b1[i]*xa_std_min + pred_betas$b8[i]
}

# Calculate predicted values for each media treatment at reference levels
# for other experimental treatments
for(i in 1:nrow(pred_betas)){
  media3_pred[i] <- pred_betas$b0[i]+pred_betas$b1[i]*xa_std_min
  media1_pred[i] <- pred_betas$b0[i]+pred_betas$b1[i]*xa_std_min + pred_betas$b4[i]
  media2_pred[i] <- pred_betas$b0[i]+pred_betas$b1[i]*xa_std_min + pred_betas$b5[i]
}

# Calculate predicted values for each pre-treatment at reference levels
# for other experimental treatments
for(i in 1:nrow(pred_betas)){
  treatment0_pred[i] <- pred_betas$b0[i]+pred_betas$b1[i]*xa_std_min
  treatment1_pred[i] <- pred_betas$b0[i]+pred_betas$b1[i]*xa_std_min + pred_betas$b9[i]
}

# Calculate predicted values for each photoperiod treatment at reference levels
# for other experimental treatments
for(i in 1:nrow(pred_betas)){
  photoperiod3_pred[i] <- pred_betas$b0[i]+pred_betas$b1[i]*xa_std_min
  photoperiod1_pred[i] <- pred_betas$b0[i]+pred_betas$b1[i]*xa_std_min + pred_betas$b10[i]
  photoperiod2_pred[i] <- pred_betas$b0[i]+pred_betas$b1[i]*xa_std_min + pred_betas$b11[i]
}

## Create graphics ####

a <- tibble(value = c(plogis(temp1_pred), plogis(temp2_pred),
                      plogis(temp3_pred), plogis(temp4_pred)),
       temp = rep(c("30", "20/15", "27/15", "25"), each = length(temp1_pred))) %>%
  ggplot(aes(x = factor(temp), y = value)) +
  geom_violin(draw_quantiles = c(0.025, 0.975), fill = "gray77",  color = "gray27") + 
  stat_summary(fun = median, geom = "point", size = 3, color = "gray27") +
  ylab("P(germination success)") +
  xlab("temperature (°C)") + ylim(0,0.9)

b <- tibble(value = c(plogis(media1_pred), plogis(media2_pred),
                      plogis(media3_pred)),
       media = rep(c("sand","growth media","sand + soil"), each = length(media1_pred))) %>%
  ggplot(aes(x = media, y = value)) +
  geom_violin(draw_quantiles = c(0.025, 0.975), fill = "gray77",  color = "gray27") + 
  stat_summary(fun = median, geom = "point", size = 3, color = "gray27") +
  ylab("") +
  xlab("type of media")+ ylim(0,0.9)

c <- tibble(value = c(plogis(treatment0_pred), plogis(treatment1_pred)),
       treatment = rep(c("no", "yes"), each = length(treatment1_pred))) %>%
  ggplot(aes(x = treatment, y = value)) +
  geom_violin(draw_quantiles = c(0.025, 0.975), fill = "gray77",  color = "gray27") + 
  stat_summary(fun = median, geom = "point", size = 3, color = "gray27") + 
  ylab("P(germination success)") +
  xlab("pre-treatment")+ ylim(0,0.9)

d <- tibble(value = c(plogis(photoperiod3_pred), plogis(photoperiod1_pred),
                      plogis(photoperiod2_pred)),
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
plogis(median(temp3_pred)) 
plogis(median(temp4_pred)) 
plogis(median(temp3_pred)) - plogis(median(temp4_pred))

# Comparing pre-treatments
plogis(median(treatment0_pred))
plogis(quantile(treatment0_pred, c(0.025, 0.975)))
plogis(median(treatment1_pred))
plogis(quantile(treatment1_pred, c(0.025, 0.975)))

# Comparing photoperiods
median(plogis(photoperiod1_pred)) 
median(plogis(photoperiod2_pred))

# Comparing media
median(plogis(media1_pred)) 
plogis(quantile(media1_pred, c(0.025, 0.975)))
median(plogis(media2_pred)) 
plogis(quantile(media2_pred, c(0.025, 0.975)))
