# Vahsen et al. "Accounting for variability when resurrecting dormant propagules
# substantiates their use in eco-evolutionary studies"

# Model to calibrate relationship between soil depth and soil age across three
# soil cores collected in 2000 (MX2 & MX4) and 2004 (MX1). 

# Outputs priors for regression coefficients and residual variance to be
# integrated into the hierarchical mixed models. Also, creates Fig S1.

# Load libraries
library(tidyverse); library(lme4); library(mvtnorm); library(cowplot); library(here)

# Read in data
soil <- read_csv(here("supp_data", "Vahsen_etal_Pb210.csv"))
# Only use year estimates for slices that have accretion rates estimated
soil %>% 
  filter(complete.cases(acc_rate)) -> soil

# Calculate adjusted seed year (year before coring date)
soil$year_adj <- soil$core_year - soil$year

# Plot depth by year predictions
soil %>% 
  ggplot(aes(x = depth_cm, y = year_adj, color = sample)) +
  geom_line() + theme_classic()

# Set depth to be the top of the layer to coincide with germination data
soil$depth_cm <- soil$depth_cm -1

##
# Fit with a polynomial with random intercept and slope and compare models
##

fit2_RS <- lmer(year_adj ~ poly(depth_cm,2,raw=TRUE) + (1+depth_cm|sample), data = soil)
fit2_RI <- lmer(year_adj ~ poly(depth_cm,2,raw=TRUE) + (1|sample), data = soil)
fit2_null <- lm(year_adj ~ poly(depth_cm,2,raw=TRUE), data = soil)

AIC(fit2_RS) # 279.2013
AIC(fit2_RI) # 275.7001
AIC(fit2_null) # 268.0081

# Simplest model is the most parsimonious for prediction

# Extract out regression coefficient information 
means <- as.numeric(coef(fit2_null))
covars <- matrix(as.numeric(vcov(fit2_null)), nrow = 3, ncol = 3)

# Extract rse (residual standard deviation)
rse <- summary(fit2_null)$sigma

# Bootstrap RSE to get variance estimate
rse.boot <- function(d){
  di <- d[sample(NROW(d), replace = TRUE), ]
  fit <- lm(year_adj ~ poly(depth_cm,2,raw=TRUE), data = di)
  rse <- summary(fit)$sigma
  rse
}

# 5000 bootstrap samples for the rse
out <- replicate(5000, rse.boot(soil))
hist(out, prob = TRUE, las = 1, xlab = "rse*")

# Calculate mean and standard deviation of rse
mean_sd <- c(mean_rse = rse, sd_rse = sd(out))

# Use moment matching to get to gamma distribution
gamma_moments <- function(mean, tau){
  var <- 1/tau
  alpha <- mean^2/var
  beta <- mean/var
  return(list(alpha = alpha, beta = beta))
}

# Assign priors for hierarchical model
sigma_prior <- gamma_moments(mean_sd[1], mean_sd[2])
beta_prior_mean <- means
beta_prior_covar <- covars

# Generate priors output to read in for source file in hierarchical model code
calibration_out <- list(sigma_prior = sigma_prior,
                        beta_prior_mean = beta_prior_mean,
                        beta_prior_covar = beta_prior_covar)

# Save prior specification for later use in hierarchical model
saveRDS(calibration_out, here("outputs", "calibration_priors.rds"))

# Create plot to predict with confidence intervals
iterations <- 1000
depth = seq(0, 35, length.out = nrow(soil))
pred_value <- matrix(NA, nrow = iterations, ncol = length(depth))
out <-  matrix(NA, nrow = iterations, ncol = length(depth))

for (i in 1:iterations){
  coefs <- mvtnorm::rmvnorm(1, beta_prior_mean, beta_prior_covar)
  for (j in 1:length(depth)){
    pred_value[i,j] <- coefs[1] + coefs[2]*depth[j] + coefs[3]*depth[j]^2 
    out[i,j] <- rnorm(1, pred_value[i,j], rse)
  }
}

# Calculate means and quantiles (both confidence & predictive) and tabulate
pred_tab <- tibble(pred_mean = colMeans(pred_value),
                   conf_low = apply(pred_value, 2, quantile, probs = c(0.025, 0.975))[1,],
                   conf_high = apply(pred_value, 2, quantile, probs = c(0.025, 0.975))[2,],
                   pred_low = apply(out, 2, quantile, probs = c(0.025, 0.975))[1,],
                   pred_high = apply(out, 2, quantile, probs = c(0.025, 0.975))[2,],
                   pred_depth = depth,
                   data_depth = soil$depth_cm,
                   data_year_adj = soil$year_adj)

# Create plot for supplement
envelope <- pred_tab %>%
  ggplot(aes(x = pred_depth, y = pred_mean)) +
  geom_line(size = 1.2) +
  theme_classic() +
  geom_line(aes(x = pred_depth, y = pred_low), linetype = "dashed", col = "orange", size = 1.2) +
  geom_line(aes(x = pred_depth, y = pred_high), linetype = "dashed", col = "orange", size = 1.2) +
  geom_point(aes(x = data_depth, y = data_year_adj), size = 3, col = "black", alpha = 0.5) +
  xlab("depth (cm)") + ylab("soil age (years)")

png(here("figs_tables", "FigS1.png"), height = 4, width = 6, res = 300, units = "in")
envelope
dev.off()


