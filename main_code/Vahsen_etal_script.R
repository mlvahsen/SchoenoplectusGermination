###
# New analysis of seed germination paper - Gentile et al. 2015
##

# Script created by MLV 7 April 2020
# Updated 22 Feb 2021
# Original models created by RG 2015

# We analyze germination data for resurrected Schoenoplectus americanus seeds
# using Bayesian hierarchical models with four different model likelihoods. 

## Preliminaries & settings ####

# Load libraries
library(tidyverse); library(rstan);
library(ggmcmc); library(cowplot); 
library(rjags); library(loo);
library(EnvStats); library(here)

# Set seed to get the same answers each time
set.seed(1234)

# Source all candidate models from separate file
source(here("main_code", "Vahsen_etal_models.R"))

# Read in germination data
germ_all <- read_csv(here("data", "germination_data.csv"))

# Consolidate data across assays (running models with a random effect of assay
# does not improve model fit)
germ_all %>% 
  group_by(Depth_Top, Core_Location, Temperature, Media, Treatment, Photoperiod,
           Temp1, Temp2, Temp3, Media1, Media2, Photo1, Photo2) %>% 
  summarize(Germinated = sum(Germinated),
            Seeds_Planted = sum(Seeds_Planted)) %>% 
  ungroup()-> germ_all_final

# Set iterations and other convergence specifications
n_iter <- 10000 # number of iterations to save
n_chains <- 3 # number of chains
n_burnin <- 2000 # number of warm-up samples
thin_interval <- 3 # sample only every X sample from each chain (reduces autocorrelation)

# Create function to convert samples to coda format. This makes it easier to
# make plots downstream
stan2coda <- function(germ.fit) {
  mcmc.list(lapply(1:ncol(germ.fit), function(x)
    mcmc(as.array(germ.fit)[, x, ])))
}

# Set plotting theme
theme_set(theme_classic())

# Bring in soil dating calibration analysis to use regression parameters as
# priors
calibration_out <- read_rds(here("outputs","calibration_priors.rds"))

## Summary statistics ####
# Get summary statistics for paper for seed number across experiments,
# locations, and depths

# Total seeds planted
sum(germ_all$Seeds_Planted)
# 10588

# Min and max number of seeds in an experiment
germ_all %>% 
  group_by(factor(Assay_Num)) %>% 
  summarize(total_seeds = sum(Seeds_Planted)) %>% 
  filter(total_seeds == min(total_seeds) | total_seeds == max(total_seeds))
# 1 6                            74
# 2 8                          3868

# Min and max number of seeds from different core locations
germ_all %>% 
  group_by(factor(Core_Location)) %>% 
  summarize(total_seeds = sum(Seeds_Planted)) %>% 
  filter(total_seeds == min(total_seeds) | total_seeds == max(total_seeds))
# 1 1                              3644
# 2 5                                 4

# Number of seeds planted from the top 10cm
germ_all %>% 
  filter(Depth_Top < 10) %>% 
  summarize(sum(Seeds_Planted)) # 3808

# Number of seeds planted from greater than 30 cm
germ_all %>% 
  filter(Depth_Top > 30) %>% 
  summarize(sum(Seeds_Planted)) # 120

# Min and max average germination success by Assay Number
germ_all %>% 
  group_by(factor(Assay_Num)) %>% 
  summarize(germ_success = sum(Germinated)/sum(Seeds_Planted)) %>% 
  filter(germ_success == max(germ_success) | germ_success == min(germ_success))
# 1 3                         0.0132
# 2 13                         0.172 


# Create a histogram to show zero-inflation
png(here("figs_tables","FigS3.png"), height = 3, width = 5.5, units = "in", res = 300)
germ_all_final %>% 
  ggplot(aes(x = Germinated)) +
  geom_histogram(binwidth = 1) +
  xlab("number of seeds germinated")
dev.off()

## Create tibble of true values ####
germ_all_final %>%
  ungroup() %>% 
  summarize(
    `maximum number of germinated seeds` = max(Germinated),
    `mean number of germinated seeds` = mean(Germinated),
    `standard deviation of number of germinated seeds` = sd(Germinated),
    `number of trials with no successful germinants` = length(which(Germinated == 0))
  ) %>%
  gather(key = "parameter", value = "estimate") -> germ_all_sumstats

## Set up data list for stan model (same for all models) ####
germ_all_model <-list(N = nrow(germ_all_final), # total number of observations
                      Nbeta = length(calibration_out$beta_prior_mean), # number of coefficients for the polynomial calibration
                      L = length(unique(germ_all$Core_Location)), # number of core locations
                      y = germ_all_final$Germinated, # number of seeds germinated per depth, trial, temp, media, location
                      s = germ_all_final$Seeds_Planted, # number of seeds planted per depth, trial, temp, media, location
                      x = germ_all_final$Depth_Top, # seed depth
                      t1 = germ_all_final$Temp1, # temperature treatment identifier 1
                      t2 = germ_all_final$Temp2, # temperature treatment identifier 2
                      t3 = germ_all_final$Temp3, # temperature treatment identifier 3
                      m1 = germ_all_final$Media1, # media treatment identifier 1
                      m2 = germ_all_final$Media2, # media treatment identifier 2
                      l = germ_all_final$Core_Location, # core location identifier
                      h = germ_all_final$Treatment, # pre-treatment identifier
                      k1 = germ_all_final$Photo1, # photoperiod identifier 1
                      k2 = germ_all_final$Photo2, #photoperiod identifier 2
                      prior_b_mean = calibration_out$beta_prior_mean,
                      prior_b_covar = calibration_out$beta_prior_covar,
                      prior_sigma_alpha = as.numeric(calibration_out$sigma_prior$alpha),
                      prior_sigma_beta = as.numeric(calibration_out$sigma_prior$beta))

## Fit binomial likelihood with no zero-inflation (Model 1) ####
# Binomial; no zero inflation
Model1_fit <- stan(model_code = model1_stan,
                            data = germ_all_model,
                            iter = n_iter,
                            chains = n_chains,
                            warmup = n_burnin,
                            thin = thin_interval,
                            control = list(max_treedepth = 15, adapt_delta=0.99))

## Fit binomial likelihood with zero-inflation (Model 2) ####

# Binomial with inflation
Model2_fit <- stan(model_code = model2_stan,
                   data = germ_all_model,
                   iter = n_iter,
                   chains = n_chains,
                   warmup = n_burnin,
                   thin = thin_interval,
                   control = list(max_treedepth = 15, adapt_delta=0.99)) 

## Fit beta-binomial likelihood with no zero-inflation (Model 3) ####
Model3_fit <- stan(model_code = model3_stan,
                   data = germ_all_model,
                   iter = n_iter,
                   chains = n_chains,
                   warmup = n_burnin,
                   thin = thin_interval,
                   control = list(max_treedepth = 15, adapt_delta=0.99))

## Fit beta-binomial likelihood with zero-inflation (Model 4) ####
Model4_fit <- stan(model_code = model4_stan,
                   data = germ_all_model,
                   iter = n_iter,
                   chains = n_chains,
                   warmup = n_burnin,
                   thin = thin_interval,
                   control = list(max_treedepth = 15, adapt_delta=0.99))


## Model comparison ####

# Calculate WAIC for each model
waic(extract_log_lik(Model1_fit)) # 1084.3
waic(extract_log_lik(Model2_fit)) # 1030.4
waic(extract_log_lik(Model3_fit)) # 771.9 # Model 3 (BB) is the best model
waic(extract_log_lik(Model4_fit)) # 773.0

# Calculate LOO (leave-one-out cross validation) for each model
loo(Model1_fit) # 1088.1
loo(Model2_fit) # 1034.0
loo(Model3_fit) # 773.2
loo(Model4_fit) # 774.4

## Post-processing of best model ####
samples <- ggs(Model3_coda) %>%
  filter(Parameter %in% c("b0","b1","sigma3",
                     "b2", "b4", "b5", "b6",
                     "b7", "b8", "b9", "b10", "b11",
                     "phi","sigma_pred"))

# Look at traceplots, density plots, and autocorrelation
samples %>%
  ggs_traceplot() + 
  facet_wrap( ~ Parameter, scales = "free")

samples %>%
  ggs_density() + 
  facet_wrap( ~ Parameter, scales = "free")

samples %>%
  ggs_autocorrelation() +
  facet_wrap( ~ Parameter, scales = "free")
# Looks pretty good

## Write rds objects to make plots ####

# Convert from stanfit objects to coda objects
Model1_coda <- stan2coda(Model1_fit)
Model2_coda <- stan2coda(Model2_fit)
Model3_coda <- stan2coda(Model3_fit)
Model4_coda <- stan2coda(Model4_fit)

# Model coda objects
write_rds(Model1_coda, here("outputs","Model1_coda_forPlotting.rds"))
write_rds(Model2_coda, here("outputs","Model2_coda_forPlotting.rds"))
write_rds(Model3_coda, here("outputs","Model3_coda_forPlotting.rds"))
write_rds(Model4_coda, here("outputs","Model4_coda_forPlotting.rds"))

# Raw germination data that takes out assay grouping that is used for plotting
write_csv(germ_all_final, here("outputs","germ_all_data_forModel.csv"))
