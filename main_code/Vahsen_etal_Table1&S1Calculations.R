# Vahsen et al. germination
# Get values for Table 1 and Table S1

## Preliminaries ####
# Load libraries
library(tidyverse);library(ggmcmc); library(here)

# Read in data that was used to fit the models (same as original data file that
# takes out grouping by Assay_Num)
germ_all <- read_csv(here("data","germination_data.csv"))

# Read in Model 3 coda object
Model3_coda <- read_rds(here("outputs","Model3_coda_forPlotting.rds"))

# Calculate a predicted age for each observed value
ggs(Model3_coda, family = "xa") %>% 
  group_by(Parameter) %>%
  filter(substr(Parameter, 1, 3) != "xa_") %>% 
  summarize(mean = mean(value)) %>% 
  filter(Parameter != "bar_xa") -> pred_age

# Add predicted age to final data frame
germ_all_final$pred_age <- pred_age$mean

germ_all %>% 
  select(Depth_Top, Core_Location, Assay_Num, Assay_Date,
         Germinated, Seeds_Planted) %>% 
  mutate(id = paste(Core_Location, Depth_Top, sep = "_"))-> germ_all_sub

germ_all_final %>% 
  select(Depth_Top, Core_Location, pred_age) %>% 
  mutate(id = paste(Core_Location, Depth_Top, sep = "_"))-> germ_all_final_sub

# List of unique depth by location combinations
unique_ids <- unique(germ_all_sub$id)

# Create empty vector to fill in pred_age for the full germ_all_sub
germ_all_sub$pred_age <- NULL

# Merge to original data set so we can get age ranges for each Assay_Num
for (i in 1:length(unique_ids)){
  pred_depth <- unique(subset(germ_all_final_sub, id == unique_ids[i])$pred_age)
  germ_all_sub[which(germ_all_sub$id == unique_ids[i]),"pred_age"] <- pred_depth
}

# Calculate predicted decade
germ_all_sub %>% 
  mutate(Pred_Decade = round(Assay_Date - pred_age, -1)) -> data_for_Table1

## Table 1 ####

# Write a csv file that has predicted decade for each location by depth
# combination

data_for_Table1 %>% 
  group_by(Core_Location, Pred_Decade) %>% 
  summarize(Total_Seeds_Planted = sum(Seeds_Planted),
            Total_Seeds_Germinated = sum(Germinated)) %>% 
  arrange(Core_Location, Pred_Decade) -> summary_data_Table1

# Write csv
write_csv(summary_data_Table1, here("outputs","Table1.csv"))

## Table S1 ####
data_for_Table1 %>% 
  group_by(Assay_Num, Assay_Date) %>% 
  summarize(lower = min(Pred_Decade),
            upper = max(Pred_Decade)) -> summary_data_TableS1

# Write csv
write_csv(summary_data_TableS1, here("outputs","TableS1.csv"))
