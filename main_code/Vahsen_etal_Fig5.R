# Create Figure 5 in Vahsen et al. germination paper
# Also calculates in-text values related to average germination for each location
# 19 Feb 2021

## Preliminaries ####
# Load libraries
library(tidyverse);library(patchwork);library(ggmcmc)
library(EnvStats);library(cowplot); library(spatstat);
library(here)

# Set ggplot theme
theme_set(theme_classic())

# Read in data that was used to fit the models (same as original data file that
# takes out grouping by Assay_Num)
germ_all_final <- read_csv(here("outputs","germ_all_data_forModel.csv"))

# Read in coda objects from all four models
Model3_coda <- read_rds(here("outputs", "Model3_coda_forPlotting.rds"))

## Calculate predicted means for each location ####

ggs(Model3_coda, family = "xa_std") %>% 
  group_by(Parameter) %>% 
  summarize(mean = mean(value)) %>% 
  mutate(depth = germ_all_final$Depth_Top) -> pred_xa_std

xa_std_range <- range(pred_xa_std$mean)
xa_std_pred <- seq(xa_std_range[1], xa_std_range[2], length.out = nrow(germ_all_final))  

ggs(Model3_coda, family = "xa") %>% 
  group_by(Parameter) %>%
  filter(substr(Parameter, 1, 3) != "xa_") %>% 
  summarize(mean = mean(value)) %>% 
  filter(Parameter != "bar_xa") -> pred_age

# Store within original data frame
germ_all_final$pred_age <- pred_age$mean

# Create a range of values to smooth across the extend of the data 
xa_range <- range(pred_age$mean)
xa_pred <- seq(xa_range[1], xa_range[2], length.out = nrow(germ_all_final))  

# Extract out regression coefficients
ggs(Model3_coda) %>% 
  group_by(Parameter) %>% 
  filter(Parameter %in% paste("b", 0:11, sep = "")) %>% 
  summarize(mean = mean(value)) %>% 
  spread(key = Parameter, value = mean) -> pred_betas

# Calculate how many are in each environmental treatment combination for weighting
germ_all_final %>% 
  group_by(Temperature, Media, Treatment, Photoperiod) %>% 
  summarize(n = length(Germinated),
            prop = n / nrow(germ_all_final)) %>% pull(prop) -> weights  # 11 different treatment combinations

# Pull betas for each value of seed provenance
ggs(Model3_coda, family = "b3") %>% 
  group_by(Parameter) %>% 
  summarize(mean = mean(value)) %>% 
  pull(mean) -> alpha_prov

pred_prov <- array(NA, c(length(alpha_prov), length(xa_std_pred), 11))
for (i in 1:length(alpha_prov)){
  pred_prov[i,,7] <- (pred_betas$b0 + pred_betas$b1 * xa_std_pred + alpha_prov[i]) # Temp = 3, Media = 1, Treatment = 0, Photo = 3
  pred_prov[i,,8] <- (pred_betas$b0 + pred_betas$b1 * xa_std_pred + pred_betas$b4+ pred_betas$b10+alpha_prov[i]) #  Temp = 3, Media = 2, Treatment = 0, Photo = 1
  pred_prov[i,,3] <- (pred_betas$b0 + pred_betas$b1 * xa_std_pred + pred_betas$b4 + pred_betas$b6+ pred_betas$b10+ alpha_prov[i]) # Temp = 1, Media = 2, Treatment = 0, Photo = 1
  pred_prov[i,,5] <- (pred_betas$b0 + pred_betas$b1 * xa_std_pred + pred_betas$b4 + pred_betas$b7+ pred_betas$b10+ alpha_prov[i]) # Temp = 2, Media = 2, Treatment = 0, Photo = 1
  pred_prov[i,,10] <- (pred_betas$b0 + pred_betas$b1 * xa_std_pred + pred_betas$b4 + pred_betas$b8+ pred_betas$b10+ alpha_prov[i]) # Temp = 4, Media = 2, Treatment = 0, Photo = 1
  pred_prov[i,,11] <- (pred_betas$b0 + pred_betas$b1 * xa_std_pred + pred_betas$b5 + pred_betas$b8+ pred_betas$b10+ alpha_prov[i]) # Temp = 4, Media = 3, Treatment = 0, Photo = 1
  pred_prov[i,,1] <- (pred_betas$b0 + pred_betas$b1 * xa_std_pred + pred_betas$b6+ pred_betas$b10+ alpha_prov[i]) # Temp = 1, Media = 1, Treatment = 0, Photo = 1
  pred_prov[i,,2] <- (pred_betas$b0 + pred_betas$b1 * xa_std_pred + pred_betas$b6+ pred_betas$b10 + pred_betas$b9+ alpha_prov[i]) # Temp = 1, Media = 1, Treatment = 1, Photo = 1
  pred_prov[i,,4] <- (pred_betas$b0 + pred_betas$b1 * xa_std_pred + pred_betas$b4 + pred_betas$b6 + pred_betas$b11+ alpha_prov[i]) # Temp = 1, Media = 2, Treatment = 0, Photo = 2
  pred_prov[i,,6] <- (pred_betas$b0 + pred_betas$b1 * xa_std_pred + pred_betas$b4 + pred_betas$b7 + pred_betas$b11+ alpha_prov[i]) # Temp = 2, Media = 2, Treatment = 0, Photo = 2
  pred_prov[i,,9] <- (pred_betas$b0 + pred_betas$b1 * xa_std_pred + pred_betas$b4+ pred_betas$b11+ alpha_prov[i]) #  Temp = 3, Media = 2, Treatment = 0, Photo = 2
}

# Calculate weighted median across treatment combinations
pred_prov_averages <- plogis(apply(pred_prov, c(1,2), weighted.mean, w = weights))

# Create new label
colnames(pred_prov_averages) <- paste(1:length(xa_std_pred))

# Panel 2 graph -- shows different provenances across seed age
prov <- as_tibble(pred_prov_averages) %>% 
  gather(key = x, value = germ_prob, `1`:`298`) %>% 
  mutate(provenance = factor(rep(1:length(alpha_prov), length(xa_std_pred))),
         age = rep(xa_pred, each = length(alpha_prov)),
         x = as.numeric(x)) %>% 
  ggplot(aes(x = age, y = germ_prob, group = provenance)) +
  geom_line(size = 1.5, alpha = 0.2) + coord_flip() + scale_x_reverse() +
  scale_y_continuous(position = "right", limits = c(0,0.4)) +
  ylab("P(germination success)") + xlab("seed age (years)")

# Panel 1 graph -- shows point estimate and HPD for each location
prov_cat <- ggs(Model3_coda, family = "b3") %>% 
  mutate(Old_Parameter = Parameter) %>% 
  mutate(Parameter = case_when(Old_Parameter == "b3[1]" ~ "01 - Kirkpatrick Marsh",
                               Old_Parameter == "b3[2]" ~ "02 - Corn Island",
                               Old_Parameter == "b3[3]" ~ "03 - Hog Island",
                               Old_Parameter == "b3[4]" ~ "04 - Virginia",
                               Old_Parameter == "b3[5]" ~ "05 - Bay Bridge",
                               Old_Parameter == "b3[6]" ~ "06 - Eastern Shore",
                               Old_Parameter == "b3[7]" ~ "07 - Blackwater Refuge",
                               Old_Parameter == "b3[8]" ~ "08 - Taylor",
                               Old_Parameter == "b3[9]" ~ "09 - Delaware",
                               Old_Parameter == "b3[10]" ~ "NA - Greenhouse",
                               Old_Parameter == "b3[11]" ~ "10 - Sellman Creek")) %>% 
  mutate(Parameter = fct_rev(Parameter)) %>%
  ggs_caterpillar(line = 0, sort = F) +ylab("")

png(here("figs_tables", "Fig5.png"), height = 3.3, width = 6.7, res = 300, units = "in")
cowplot::plot_grid(prov_cat, prov, align = "h",
                   nrow = 1, labels = "auto")
dev.off()

## Calculate effect size for main text ####

# These are at the top layer of the core

ggs(Model3_coda) %>% 
  filter(Parameter %in% c("b0", "b1", "b2", "b3", "b4", "b5", "b6", "b7", "b8", "b9", "b10", "b11")) %>%
  spread(key = Parameter, value = value) -> betas_mc

pred_prov_fortext <- array(NA, c(length(alpha_prov), nrow(betas_mc), 11))
for (i in 1:length(alpha_prov)){
  for (j in 1:nrow(betas_mc)){
    pred_prov_fortext[i,j,7] <- (betas_mc$b0[j] + betas_mc$b1[j] * xa_std_range[1] + alpha_prov[i]) # Temp = 3, Media = 1, Treatment = 0, Photo = 3
    pred_prov_fortext[i,j,8] <- (betas_mc$b0[j] + betas_mc$b1[j] * xa_std_range[1] + betas_mc$b4[j]+ + betas_mc$b10[j]+alpha_prov[i]) #  Temp = 3, Media = 2, Treatment = 0, Photo = 1
    pred_prov_fortext[i,j,3] <- (betas_mc$b0[j] + betas_mc$b1[j] * xa_std_range[1] + betas_mc$b4[j] + betas_mc$b6[j]+ betas_mc$b10[j]+ alpha_prov[i]) # Temp = 1, Media = 2, Treatment = 0, Photo = 1
    pred_prov_fortext[i,j,5] <- (betas_mc$b0[j] + betas_mc$b1[j] * xa_std_range[1] + betas_mc$b4[j] + betas_mc$b7[j]+ betas_mc$b10[j]+ alpha_prov[i]) # Temp = 2, Media = 2, Treatment = 0, Photo = 1
    pred_prov_fortext[i,j,10] <- (betas_mc$b0[j] + betas_mc$b1[j] * xa_std_range[1] + betas_mc$b4[j] + betas_mc$b8[j]+ betas_mc$b10[j]+ alpha_prov[i]) # Temp = 4, Media = 2, Treatment = 0, Photo = 1
    pred_prov_fortext[i,j,11] <- (betas_mc$b0[j] + betas_mc$b1[j] * xa_std_range[1] + betas_mc$b5[j] + betas_mc$b8[j]+ betas_mc$b10[j]+ alpha_prov[i]) # Temp = 4, Media = 3, Treatment = 0, Photo = 1
    pred_prov_fortext[i,j,1] <- (betas_mc$b0[j] + betas_mc$b1[j] * xa_std_range[1] + betas_mc$b6[j]+betas_mc$b10[j]+ alpha_prov[i]) # Temp = 1, Media = 1, Treatment = 0, Photo = 1
    pred_prov_fortext[i,j,2] <- (betas_mc$b0[j] + betas_mc$b1[j] * xa_std_range[1] + betas_mc$b6[j] + betas_mc$b9[j]+ betas_mc$b10[j]+alpha_prov[i]) # Temp = 1, Media = 1, Treatment = 1, Photo = 1
    pred_prov_fortext[i,j,4] <- (betas_mc$b0[j] + betas_mc$b1[j] * xa_std_range[1] + betas_mc$b4[j] + betas_mc$b6[j] + betas_mc$b11[j]+ alpha_prov[i]) # Temp = 1, Media = 2, Treatment = 0, Photo = 2
    pred_prov_fortext[i,j,6] <- (betas_mc$b0[j] + betas_mc$b1[j] * xa_std_range[1] + betas_mc$b4[j] + betas_mc$b7[j] + betas_mc$b11[j]+ alpha_prov[i]) # Temp = 2, Media = 2, Treatment = 0, Photo = 2
    pred_prov_fortext[i,j,9] <- (betas_mc$b0[j] + betas_mc$b1[j] * xa_std_range[1] + betas_mc$b4[j]+ betas_mc$b11[j]+ alpha_prov[i]) #  Temp = 3, Media = 2, Treatment = 0, Photo = 2
  }
 }

pred_prov_averages_fortext <- apply(pred_prov_fortext, c(1,2), weighted.mean, w = weights)
# Calculate difference between Kirkpatrick (1) and Corn (2)
median(plogis(pred_prov_averages_fortext[1,])) - median(plogis(pred_prov_averages_fortext[2,])) 
# 0.1733297
# This is the same as 
plogis(median(pred_prov_averages_fortext[1,])) - plogis(median(pred_prov_averages_fortext[2,])) 
# 0.1733297

# Calculate 95% quantiles around that
quantile(plogis(pred_prov_averages_fortext[1,]) - plogis(pred_prov_averages_fortext[2,]), c(0.025, 0.975)) 
# 2.5%     97.5% 
# 0.1245294 0.2164936 

# Calculate proportion of observations that were 0s
germ_all_final %>% 
  filter(Germinated == 0) %>% 
  nrow() -> zero_trials

zero_trials / nrow(germ_all_final)      
# 0.590604     
      
      