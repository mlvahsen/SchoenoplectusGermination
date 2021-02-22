# Analyze pre-treatment data in assay 3 to see if we can treat all
# pre-treatments as the same

# Load libraries
library(tidyverse); library(brglm2); library(emmeans); library(here)

# Read in data
assay <- read_csv(here("supp_data", "assay3_pretreatment.csv"))

# Filter just assay 3 
assay %>% 
  filter(Assay_Num == 3) %>% 
  group_by(Treatment) %>% 
  summarize(n = length(Assay_Num),
            germ = length(which(Germinated == 1)),
            prop = germ/n)
# Pre-treatments seemed to decrease the efficacy similarly

# Subset data for individual analysis of Assay 3 data
assay3 <- assay %>%
  filter(Assay_Num == 3 & Treatment != "None")  

# There is complete separation so we can fit with brglm()
test <- glm(Germinated ~ Treatment + Depth_Top, data = assay3, family = "binomial")
test_br <- update(test, method = "brglmFit")

test_null <- glm(Germinated ~ Depth_Top, data = assay3, family = "binomial")
test_null_br <- update(test_null, method = "brglmFit")

# Model comparisons via anova() are not recommended. We can compare the residual
# deviance (higher deviance means less is explained by the model). Including
# treatment seems to decrease the residual deviance considerably.
summary(test_br)
summary(test_null_br)

png(here("figs_tables", "FigS2.png"), height = 4, width = 5, res = 300, units = "in")
plot(emmeans(test_br, ~Treatment), type = "response") +
  theme_classic() +
  xlab('germination probability') +
  ylab("pre-treatment for assay 3")
dev.off()

# We could group all of the pre-treatment seeds together
assay3 %>% 
  mutate(treatment_presence = ifelse(Treatment == "None", 0, 1)) %>% 
  group_by(treatment_presence) %>% 
  count()

# There are 1262 seeds that experienced a pre-treatment which is >10% of the
# total number of seeds in the germination manuscript. I think we should group
# by pre-treatment given the large number of seeds that experienced
# pre-treatments.

