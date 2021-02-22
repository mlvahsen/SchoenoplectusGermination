# Create Figure 2 in Vahsen et al. germination paper
# 19 Feb 2021

# Load libraries
library(tidyverse)
library(patchwork)
library(ggmcmc)
library(here)

# Set ggplot theme
theme_set(theme_classic())

# Read in data that was used to fit the models (same as original data file that
# takes out grouping by Assay_Num)
germ_all_final <- read_csv(here("outputs","germ_all_data_forModel.csv"))

# Read in coda objects from all four models
Model1_coda <- read_rds(here("outputs","Model1_coda_forPlotting.rds"))
Model2_coda <- read_rds(here("outputs","Model2_coda_forPlotting.rds"))
Model3_coda <- read_rds(here("outputs","Model3_coda_forPlotting.rds"))
Model4_coda <- read_rds(here("outputs","Model4_coda_forPlotting.rds"))

# Use same settings for the fitting the stan models. Need to know the number of
# iterations and chains to create this plot
n_iter <- 10000 # number of iterations to save
n_chains <- 3 # number of chains
n_burnin <- 2000 # number of warm-up samples
thin_interval <- 3 # sample only every X sample from each chain (reduces autocorrelation)


# Create function that calculates posterior predictive values of the mean,
# standard deviation, and number of zeros in the data set
get_post_preds <- function(coda_object, model_number){
  ggs(coda_object, family = "y_new") %>% 
    group_by(Iteration) %>% 
    summarize(mean = mean(value),
              `standard deviation` = sd(value),
              `number of zeros` = length(which(value == 0))/n_chains) %>% 
    mutate(model = model_number) %>% 
    gather(key = sumstat, value = value, mean:`number of zeros`) -> out
  return(out)}

# Set color palette (colorblind-friendly)
cbPalette <- c("#999999",
               "#E69F00",
               "#56B4E9",
               "#009E73",
               "#F0E442",
               "#0072B2",
               "#D55E00",
               "#CC79A7")

# Apply function to collect posterior predictions for each model
post_pred_all <- rbind(get_post_preds(Model1_coda, "Model 1 - BIN"),
                       get_post_preds(Model2_coda, "Model 2 - ZI-BIN"),
                       get_post_preds(Model3_coda, "Model 3 - BETABIN"),
                       get_post_preds(Model4_coda, "Model 4 - ZI-BETABIN"))

mean <- post_pred_all %>% 
  filter(sumstat == "mean") %>% 
  ggplot(aes(x = value)) +
  geom_density(alpha = 0.2, aes(color = model, fill = model)) +
  xlab("mean number of germinated seeds") +
  scale_color_manual(values = cbPalette[c(2,3,4,8)]) +
  scale_fill_manual(values = cbPalette[c(2,3,4,8)]) +
  geom_vline(aes(xintercept = mean(germ_all_final$Germinated)), linetype = "dashed") +
  theme(legend.title=element_blank())

stdev <- post_pred_all %>% 
  filter(sumstat == "standard deviation") %>% 
  ggplot(aes(x = value)) +
  geom_density(alpha = 0.2, aes(color = model, fill = model)) +
  xlab("std dev number of germinated seeds") +
  scale_color_manual(values = cbPalette[c(2,3,4,8)]) +
  scale_fill_manual(values = cbPalette[c(2,3,4,8)]) +
  geom_vline(aes(xintercept = sd(germ_all_final$Germinated)), linetype = "dashed") +
  ylab("") + theme(legend.title=element_blank())

zeros <- post_pred_all %>% 
  filter(sumstat == "number of zeros") %>% 
  ggplot(aes(x = value)) +
  geom_density(alpha = 0.2, aes(color = model, fill = model)) +
  xlab("number of trials with zero germinated seeds") +
  scale_color_manual(values = cbPalette[c(2,3,4,8)]) +
  scale_fill_manual(values = cbPalette[c(2,3,4,8)]) +
  geom_vline(aes(xintercept = length(which(germ_all_final$Germinated == 0))), linetype = "dashed") +
  ylab("") + theme(legend.title=element_blank())

png(here("figs_tables", "Fig2.png"), res = 300, units = "in", height = 3, width = 11)
mean + labs(title = 'a') + theme(plot.title = element_text(hjust = -0.15, size = 16, face = "bold")) +
  stdev + labs(title = 'b') + theme(plot.title = element_text(hjust = -0.15, size = 16, face = "bold")) +
  zeros + labs(title = 'c') + theme(plot.title = element_text(hjust = -0.15, size = 16, face = "bold")) +
  plot_layout(guides = "collect")
dev.off()             
