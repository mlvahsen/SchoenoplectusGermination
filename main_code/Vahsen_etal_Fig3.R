# Create Figure 3 in Vahsen et al. germination paper
# Also calculates in-text values related to average germination at certain depths
# 19 Feb 2021

## Preliminaries ####
# Load libraries
library(tidyverse);library(patchwork);library(ggmcmc)
library(EnvStats);library(cowplot); library(here); library(spatstat)

# Set ggplot theme
theme_set(theme_classic())

# Read in data that was used to fit the models (same as original data file that
# takes out grouping by Assay_Num)
germ_all_final <- read_csv(here("outputs","germ_all_data_forModel.csv"))

# Read in coda objects from all four models
Model1_coda <- read_rds(here("outputs","Model1_coda_forPlotting.rds"))
Model2_coda <- read_rds(here("outputs", "Model2_coda_forPlotting.rds"))
Model3_coda <- read_rds(here("outputs", "Model3_coda_forPlotting.rds"))
Model4_coda <- read_rds(here("outputs", "Model4_coda_forPlotting.rds"))

# Use same settings for the fitting the stan models. Need to know the number of
# iterations and chains to create this plot
n_iter <- 10000 # number of iterations to save
n_chains <- 3 # number of chains
n_burnin <- 2000 # number of warm-up samples
thin_interval <- 3 # sample only every X sample from each chain (reduces autocorrelation)

## Calculate predicted versus observed ####

# Get predicted values and plot against observed
get_pred_obs <- function(coda_object, data, sum_stats, n_iter, n_chains, reps, ylab){
  coda_object %>%  
    ggs(family = "y_new") %>%
    group_by(Parameter) %>%
        summarize(
        predicted = mean(value),
        lower = quantile(value, 0.025),
        upper = quantile(value, 0.975)) %>% 
    mutate(observed = data$Germinated) -> sum_stats
  
  data$pred_mean <- sum_stats$predicted
  
  inset_3a <- sum_stats %>% 
    filter(observed < 6) %>% 
    ggplot(aes(x = observed, y = predicted)) +
    geom_pointrange(aes(ymin = lower, ymax = upper), 
                    position=position_jitter(width=0.3), alpha = 0.3, size = 0.1) +
    geom_abline(aes(intercept = 0, slope = 1)) +
    scale_x_continuous(breaks = 0:10) + ylim(-4,90) +
    geom_vline(aes(xintercept = 0.5), linetype = "dashed", color = "gray87") +
    geom_vline(aes(xintercept = 1.5), linetype = "dashed", color = "gray87") +
    geom_vline(aes(xintercept = 2.5), linetype = "dashed", color = "gray87") +
    geom_vline(aes(xintercept = 3.5), linetype = "dashed", color = "gray87") +
    geom_vline(aes(xintercept = 4.5), linetype = "dashed", color = "gray87") +
    xlab("") + ylab("") + stat_n_text(size = 2, y.pos = 88) + theme(panel.border = element_rect(colour = "gray47", fill=NA, size=1, linetype = "dashed"),
                                                panel.background = element_rect(fill = "transparent", color = NA),
                                                plot.background = element_rect(fill = "transparent", color = NA))
  
  full_3a <- sum_stats %>% 
    ggplot(aes(x = observed, y = predicted)) +
    geom_pointrange(aes(ymin = lower, ymax = upper), alpha = 0.5) +
    geom_abline(aes(intercept = 0, slope = 1)) +
    ylab(ylab) + ylim(-1,250) +
    xlab("observed number of seeds germinated") +
    geom_rect(aes(xmin = -1, xmax = 5.5, ymin = -1, ymax = 90), fill = "NA", color = "gray47",
              linetype = "dashed") +
    geom_segment(aes(x = 0, y = 90, xend = 11, yend = 185), linetype = "dashed", color = "gray47") +
    geom_segment(aes(x = 5, y = 90, xend = 75, yend = 185), linetype = "dashed", color = "gray47")
  
  ggdraw() +
    draw_plot(full_3a) +
    draw_plot(inset_3a, x = 0.12, y = 0.64, width = 0.50, height = 0.36) 
  
}

## Create plots ####

# Figure 3a
Fig3a <- get_pred_obs(Model3_coda, germ_all_final, sum_stats, n_iter, n_chains, 8001,  "predicted number of seeds germinated")
# Supplemental Figure for the rest of the candidate models
a <- get_pred_obs(Model1_coda, germ_all_final, sum_stats, n_iter, n_chains, 8001, "predicted number of seeds germinated")
b <- get_pred_obs(Model2_coda, germ_all_final, sum_stats, n_iter, n_chains, 8001, "")
c <- get_pred_obs(Model4_coda, germ_all_final, sum_stats, n_iter, n_chains, 8001, "")

FigS4 <- plot_grid(a,b,c, nrow = 1, labels = c("M1", "M2", "M4"), label_size = 12, label_x = -0.02)

# Figure 3b
# Fig 3b -- predicted germination probability with CIs
# Predicted means -- temperature treatments
ggs(Model3_coda) %>% 
  filter(Parameter %in% c("b0", "b1", "b4", "b5", "b6", "b7", "b8", "b9", "b10", "b11")) %>% 
  spread(key = Parameter, value = value) -> pred_betas

ggs(Model3_coda, family = "xa_std") %>% 
  group_by(Parameter) %>% 
  summarize(mean = mean(value)) %>% 
  mutate(depth = germ_all_final$Depth_Top) -> pred_xa_std

ggs(Model3_coda, family = "xa") %>% 
  group_by(Parameter) %>%
  filter(substr(Parameter, 1, 3) != "xa_") %>% 
  summarize(mean = mean(value)) %>% 
  filter(Parameter != "bar_xa") -> pred_age

# Add predicted age to the data
germ_all_final$pred_age <- pred_age$mean

xa_range <- range(pred_age$mean)
xa_pred <- seq(xa_range[1], xa_range[2], length.out = nrow(germ_all_final))  
xa_std_range <- range(pred_xa_std$mean)
xa_std_pred <- seq(xa_std_range[1], xa_std_range[2], length.out = nrow(germ_all_final))  

# Need to do weighted mean across treatments
# Calculate how many are in each environmental treatment combination
germ_all_final %>% 
  group_by(Temperature, Media, Treatment, Photoperiod) %>% 
  summarize(n = length(Germinated),
            prop = n / nrow(germ_all_final)) %>% pull(prop) -> weights  # 11 different treatment combinations


out <- array(NA, c(nrow(pred_betas), length(xa_std_pred), 11))
for (i in 1:nrow(pred_betas)){
  out[i,,7] <- (pred_betas$b0[i] + pred_betas$b1[i] * xa_std_pred) # Temp = 3, Media = 1, Treatment = 0, Photo = 3
  out[i,,8] <- (pred_betas$b0[i] + pred_betas$b1[i] * xa_std_pred + pred_betas$b4[i]+ pred_betas$b10[i]) #  Temp = 3, Media = 2, Treatment = 0, Photo = 1
  out[i,,3] <- (pred_betas$b0[i] + pred_betas$b1[i] * xa_std_pred + pred_betas$b4[i] + pred_betas$b6[i] + pred_betas$b10[i]) # Temp = 1, Media = 2, Treatment = 0, Photo = 1
  out[i,,5] <- (pred_betas$b0[i] + pred_betas$b1[i] * xa_std_pred + pred_betas$b4[i] + pred_betas$b7[i] + pred_betas$b10[i]) # Temp = 2, Media = 2, Treatment = 0, Photo = 1
  out[i,,10] <- (pred_betas$b0[i] + pred_betas$b1[i] * xa_std_pred + pred_betas$b4[i] + pred_betas$b8[i] + pred_betas$b10[i]) # Temp = 4, Media = 2, Treatment = 0, Photo = 1
  out[i,,11] <- (pred_betas$b0[i] + pred_betas$b1[i] * xa_std_pred + pred_betas$b5[i] + pred_betas$b8[i] + pred_betas$b10[i]) # Temp = 4, Media = 3, Treatment = 0, Photo = 1
  out[i,,1] <- (pred_betas$b0[i] + pred_betas$b1[i] * xa_std_pred + pred_betas$b6[i] + pred_betas$b10[i]) # Temp = 1, Media = 1, Treatment = 0, Photo = 1
  out[i,,2] <- (pred_betas$b0[i] + pred_betas$b1[i] * xa_std_pred + pred_betas$b6[i] + pred_betas$b9[i] + pred_betas$b10[i]) # Temp = 1, Media = 1, Treatment = 1, Photo = 1
  out[i,,4] <- (pred_betas$b0[i] + pred_betas$b1[i] * xa_std_pred + pred_betas$b4[i] + pred_betas$b6[i] + pred_betas$b11[i]) # Temp = 1, Media = 2, Treatment = 0, Photo = 2
  out[i,,6] <- (pred_betas$b0[i] + pred_betas$b1[i] * xa_std_pred + pred_betas$b4[i] + pred_betas$b7[i] + pred_betas$b11[i]) # Temp = 2, Media = 2, Treatment = 0, Photo = 2
  out[i,,9] <- (pred_betas$b0[i] + pred_betas$b1[i] * xa_std_pred + pred_betas$b4[i]+ pred_betas$b11[i]) #  Temp = 3, Media = 2, Treatment = 0, Photo = 2
}

pred_averages <- apply(out, c(1,2), weighted.mean, w = weights)
pred_quantiles <- apply(pred_averages, 2, quantile, probs = c(0.025,0.5,0.975))

tibble(median = plogis(pred_quantiles[2,]),
       lower = plogis(pred_quantiles[1,]),
       upper = plogis(pred_quantiles[3,]),
       xa_pred = xa_pred,
       depth = seq(range(germ_all_final$Depth_Top)[1], range(germ_all_final$Depth_Top)[2], length.out = nrow(germ_all_final))) %>% 
  ggplot(aes(x = xa_pred, y = median)) +
  geom_line(size = 1.2) +
  geom_line(aes(x = xa_pred, y = lower), linetype = "dashed", size = 1) +
  geom_line(aes(x = xa_pred, y = upper), linetype = "dashed", size = 1) +
  ylab("P(germination success)") +
  geom_point(aes(x = germ_all_final$pred_age, y = germ_all_final$Germinated/germ_all_final$Seeds_Planted), alpha = 0.1, size = 2) +
  scale_x_reverse() + coord_flip() + xlab("seed age (years)") -> Fig3b

Fig3 <- plot_grid(Fig3a, Fig3b, labels = "auto")

# Print Figure 3
png(here("figs_tables","Fig3.png"), res = 300, units = "in", height = 4, width = 8)
Fig3
dev.off()

# Print Supp Figure 2
png(here("figs_tables","FigS4.png"), res = 300, units = "in", height = 4, width = 11.8)
FigS4
dev.off()

## Calculations for in-text ####

# Calculate predicted CIs for different age depths for text
tibble(depth = pred_xa_std$depth,
       age = pred_age$mean,
       age_std = pred_xa_std$mean) %>% 
  filter(depth == 0 | depth == 10 | depth == 20) %>% 
  distinct() %>% 
  pull(age_std) -> age_std_fortext

out_fortext <- array(NA, c(nrow(pred_betas), length(age_std_fortext), 11))
for (i in 1:nrow(pred_betas)){
  out_fortext[i,,7] <- (pred_betas$b0[i] + pred_betas$b1[i] * age_std_fortext) # Temp = 3, Media = 1, Treatment = 0, Photo = 3
  out_fortext[i,,8] <- (pred_betas$b0[i] + pred_betas$b1[i] * age_std_fortext + pred_betas$b4[i]+ pred_betas$b10[i]) #  Temp = 3, Media = 2, Treatment = 0, Photo = 1
  out_fortext[i,,3] <- (pred_betas$b0[i] + pred_betas$b1[i] * age_std_fortext + pred_betas$b4[i] + pred_betas$b6[i] + pred_betas$b10[i]) # Temp = 1, Media = 2, Treatment = 0, Photo = 1
  out_fortext[i,,5] <- (pred_betas$b0[i] + pred_betas$b1[i] * age_std_fortext + pred_betas$b4[i] + pred_betas$b7[i] + pred_betas$b10[i]) # Temp = 2, Media = 2, Treatment = 0, Photo = 1
  out_fortext[i,,10] <- (pred_betas$b0[i] + pred_betas$b1[i] * age_std_fortext + pred_betas$b4[i] + pred_betas$b8[i] + pred_betas$b10[i]) # Temp = 4, Media = 2, Treatment = 0, Photo = 1
  out_fortext[i,,11] <- (pred_betas$b0[i] + pred_betas$b1[i] * age_std_fortext + pred_betas$b5[i] + pred_betas$b8[i] + pred_betas$b10[i]) # Temp = 4, Media = 3, Treatment = 0, Photo = 1
  out_fortext[i,,1] <- (pred_betas$b0[i] + pred_betas$b1[i] * age_std_fortext + pred_betas$b6[i] + pred_betas$b10[i]) # Temp = 1, Media = 1, Treatment = 0, Photo = 1
  out_fortext[i,,2] <- (pred_betas$b0[i] + pred_betas$b1[i] * age_std_fortext + pred_betas$b6[i] + pred_betas$b9[i] + pred_betas$b10[i]) # Temp = 1, Media = 1, Treatment = 1, Photo = 1
  out_fortext[i,,4] <- (pred_betas$b0[i] + pred_betas$b1[i] * age_std_fortext + pred_betas$b4[i] + pred_betas$b6[i] + pred_betas$b11[i]) # Temp = 1, Media = 2, Treatment = 0, Photo = 2
  out_fortext[i,,6] <- (pred_betas$b0[i] + pred_betas$b1[i] * age_std_fortext + pred_betas$b4[i] + pred_betas$b7[i] + pred_betas$b11[i]) # Temp = 2, Media = 2, Treatment = 0, Photo = 2
  out_fortext[i,,9] <- (pred_betas$b0[i] + pred_betas$b1[i] * age_std_fortext + pred_betas$b4[i]+ pred_betas$b11[i]) #  Temp = 3, Media = 2, Treatment = 0, Photo = 2
}

pred_averages_fortext <- apply(out_fortext, c(1,2), weighted.mean, w = weights)
pred_quantiles_fortext <- apply(pred_averages_fortext, 2, quantile, probs = c(0.025,0.975))

# These are the predicted germination probabilities (quantiles and median)
predicted_medians <- plogis(apply(pred_averages_fortext, 2, median))
predicted_medians
# 0.21803969 0.09298733 0.02455019
plogis(pred_quantiles_fortext)
# 2.5%  0.1410249 0.05827305 0.01322760
# 97.5% 0.3128686 0.13911843 0.04262377

# Calculate the CIs for Model 4 parameter to show that there is no benefit of
# adding zero-inflation to the beta-binomial model
ggs(Model4_coda) %>% 
  filter(Parameter == "g1") %>% 
  summarize(mean = mean(value),
            lower = quantile(value, 0.025),
            upper = quantile(value, 0.975))
# mean lower upper
# -0.787 -10.5  7.49

# Calculate predicted CIs for age depths for text (20cm depth)
which(germ_all_final$Depth_Top == 20)
# 221 222 223 224 225 226 227 228 229 230
# These are the observations that have depth = 20

# Get their predicted age from the model
ggs(Model3_coda, family = "xa") %>% 
  filter(Parameter == "xa[221]") %>% 
  summarize(mean = mean(value),
            lower = quantile(value, 0.025),
            upper = quantile(value, 0.975)) 
# 107.  93.9  120.

# Get slope of b1 for in-text
ggs(Model3_coda) %>% 
  filter(Parameter == "b1") %>% 
  summarize(mean = mean(value),
            lower = quantile(value, 0.025),
            upper = quantile(value, 0.975)) 
# mean lower upper
# <dbl> <dbl> <dbl>
# -1.32 -1.62 -1.04

## Calculate R2 for each model ####
give_me_R2 <- function(preds,actual){
  rss <- sum(( preds - actual ) ^ 2)  ## residual sum of squares
  tss <- sum((actual - mean(actual)) ^ 2)  ## total sum of squares
  rsq <- 1 - rss/tss
  return(rsq)
}

Model1_coda %>%  
  ggs(family = "y_new") %>% 
  group_by(Parameter) %>%
  summarize(
    predicted = mean(value),
    lower = quantile(value, 0.025),
    upper = quantile(value, 0.975)
  ) %>%
  mutate(observed = germ_all_final$Germinated) -> Model1_PO

Model2_coda %>%  
  ggs(family = "y_new") %>% 
  group_by(Parameter) %>%
  summarize(
    predicted = mean(value),
    lower = quantile(value, 0.025),
    upper = quantile(value, 0.975)
  ) %>%
  mutate(observed = germ_all_final$Germinated) -> Model2_PO

Model3_coda %>%  
  ggs(family = "y_new") %>% 
  group_by(Parameter) %>%
  summarize(
    predicted = mean(value),
    lower = quantile(value, 0.025),
    upper = quantile(value, 0.975)
  ) %>%
  mutate(observed = germ_all_final$Germinated) -> Model3_PO

Model4_coda %>%  
  ggs(family = "y_new") %>% 
  group_by(Parameter) %>%
  summarize(
    predicted = mean(value),
    lower = quantile(value, 0.025),
    upper = quantile(value, 0.975)
  ) %>%
  mutate(observed = germ_all_final$Germinated) -> Model4_PO

give_me_R2(Model1_PO$predicted, Model1_PO$observed)
# 0.9118568

give_me_R2(Model2_PO$predicted, Model2_PO$observed)
# 0.8899808

give_me_R2(Model3_PO$predicted, Model3_PO$observed)
# 0.8689081

give_me_R2(Model4_PO$predicted, Model4_PO$observed)
# 0.8658285

