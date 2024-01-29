# dust_single_model.R
#
# Code for model selection to define the best model for explaining dust PFAS levels
# using exposure questionnaire responses 

library(dplyr)
library(ggplot2)
library(MASS)
library(car)

### load pre-processed individual-level data with non-detects as DL / sqrt(2)
df = read.csv('output/data_all_sites_questionnaire.csv', stringsAsFactors=TRUE, na.strings= c('NA', ''))
PFAS = c('PFOA', 'PFOS', 'PFNA', 'PFDA', 'PFUnA', 'PFHxS', 'MeFOSAA')

### make unique household id column
df$hid = paste(df$site, df$householdid, sep="")
length(unique(df$hid))  # 114 households

### #convert concs to log
PFAS = c('PFOA', 'PFOS', 'PFNA', 'PFDA', 'PFUnA', 'PFHxS', 'MeFOSAA')
for (substance in PFAS) {
  df[paste('log_', substance, '_dust', sep="")] = log(df[paste(substance, '_dust_ng.g', sep="")])
  df[paste('log_', substance, '_serum', sep="")] = log(df[paste(substance, '_serum_ng.ml', sep="")])
}

### calculate sum of PFAS in serum + dust 
df$sum_PFAS_serum_ng.ml <- df$PFOA_serum_ng.ml + df$PFOS_serum_ng.ml + df$PFHxS_serum_ng.ml + df$PFDA_serum_ng.ml + df$PFNA_serum_ng.ml +
  df$MeFOSAA_serum_ng.ml + df$PFUnA_serum_ng.ml
df$sum_PFAS_dust_ng.g <- df$PFOA_dust_ng.g + df$PFOS_dust_ng.g + df$PFHxS_dust_ng.g + df$PFDA_dust_ng.g + df$PFNA_dust_ng.g +
  df$MeFOSAA_dust_ng.g + df$PFUnA_dust_ng.g
df$log_sum_PFAS_serum = log(df$sum_PFAS_serum_ng.ml)
df$log_sum_PFAS_dust = log(df$sum_PFAS_dust_ng.g)

### calculate new variables from raw questionnaire answers

# years at home
years_at_home_years = df$AQ3_Years#+ df$AQ3_Months/12
years_at_home_years[is.na(years_at_home_years)] = 0
years_at_home_months = df$AQ3_Months/12
years_at_home_months[is.na(years_at_home_months)] = 0
df$years_at_home = years_at_home_years + years_at_home_months

# drinking water source - leaving out because water source highly confounded w/site
#df$AQ10_WaterSource[df$AQ10_WaterSource==""] = NA
#df$AQ10_WaterSource = factor(df$AQ10_WaterSource)

# cleaning frequency
df$cleaning_bin = ifelse(df$AQ14_Cleaning == "Three times per week or more", 1, 0)
df$cleaning_bin[is.na(df$cleaning_bin)] = 0 

# stain resistant product use
df$stain_bin = ifelse(df$AQ15_StainResist == "Never", 0, 1)
df$stain_bin[is.na(df$stain_bin)] = 0 

# occupational exposure
df$work_bin = ifelse(df$AQ28_None == 'False', 1, 0)
df$work_bin[is.na(df$work_bin)] = 0 
df$work_bin[df$AQ28_DontKnow == 'True']

# race & ethnicity
df$white_nonhispanic = ifelse((df$AQ2_White %in% 'True') & 
                                !((df$AQ1_Ethnicity %in% 'True') | (df$CQ2_Ethnicity %in% 'True')), 1,0)

### aggregate questionnaire to household level
dfhh = df[ , grepl( "dust|hid|site|age|years|bin|white_nonhispanic|Floor" , names( df ) ) ]
dfhh= as.data.frame(dfhh %>% 
                           group_by(hid, site) %>% 
                           summarise(across(where(is.numeric), list(mean=mean,max=max, all=min), na.rm=TRUE),
                                     across(where(is.factor), list(first=first))))
nrow(dfhh) # 114
colnames(dfhh)

### format factors at household level
dfhh$stain_bin_max = factor(dfhh$stain_bin_max, labels=c('no','yes'))
dfhh$cleaning_bin_max = factor(dfhh$cleaning_bin_max, labels=c('Low/Medium', 'High'))
dfhh$stain_bin_max = factor(dfhh$stain_bin_max, labels=c('no','yes'))
dfhh$work_bin_max = factor(dfhh$work_bin_max, labels=c('no', 'yes')) # indicates if anyone in hh has worked in those industries
dfhh$white_nonhispanic_all = factor(dfhh$white_nonhispanic_all, labels=c('no', 'yes'))

# aggregate soft vs hard flooring
dfhh$living_room <- factor(dfhh$AQ16_FloorLR)
dfhh$bedroom <- factor(dfhh$AQ18_FloorB)
dfhh$living_room <- recode(dfhh$AQ16_FloorLR, "c('Hardwood', 'Laminate',
                                  'Tile', 'Vinyl')='hard';c('Carpet') = 'soft'")
dfhh$bedroom <- recode(dfhh$AQ18_FloorB_first, "c('Hardwood', 'Laminate',
                                  'Tile', 'Vinyl')='hard';c('Carpet') = 'soft'")


### single model to explain dust levels
### using stepwise regression

#PFOS
PFOS_full_model = lm(log_PFOS_dust_mean ~ site + age_mean + years_at_home_mean + cleaning_bin_max + 
             stain_bin_max + work_bin_max + living_room + 
             bedroom + white_nonhispanic_all, data = dfhh)
PFOS_step_model = step(PFOS_full_model, direction="both", trace=FALSE)
summary(PFOS_step_model)

#PFOA
PFOA_full_model = lm(log_PFOA_dust_mean ~ site + age_mean + years_at_home_mean + cleaning_bin_max + 
                        stain_bin_max + work_bin_max + living_room + 
                       bedroom + white_nonhispanic_all, data = dfhh)
PFOA_step_model = step(PFOA_full_model, direction="both", trace=FALSE)
summary(PFOA_step_model)


#PFHxS
PFHxS_full_model = lm(log_PFHxS_dust_mean ~ site + age_mean + years_at_home_mean + cleaning_bin_max + 
                        stain_bin_max + work_bin_max + living_room + 
                        bedroom + white_nonhispanic_all, data = dfhh)
PFHxS_step_model = step(PFHxS_full_model, direction="both", trace=FALSE)
summary(PFHxS_step_model)

#PFNA
PFNA_full_model = lm(log_PFNA_dust_mean ~ site + age_mean + years_at_home_mean + cleaning_bin_max + 
                       stain_bin_max + work_bin_max + living_room + 
                       bedroom + white_nonhispanic_all, data = dfhh)
PFNA_step_model = step(PFNA_full_model, direction="both", trace=FALSE)
summary(PFNA_step_model)
TukeyHSD(aov(PFNA_step_model), which=c("site")) 


#PFDA
PFDA_full_model = lm(log_PFDA_dust_mean ~ site + age_mean + years_at_home_mean + cleaning_bin_max + 
                       stain_bin_max + work_bin_max + living_room + 
                       bedroom + white_nonhispanic_all, data = dfhh)
PFDA_step_model = step(PFDA_full_model, direction="both", trace=FALSE)
summary(PFDA_step_model)
TukeyHSD(aov(PFDA_step_model), which=c("site")) 


#MeFOSAA
MeFOSAA_full_model = lm(log_MeFOSAA_dust_mean ~ site + age_mean + years_at_home_mean + cleaning_bin_max + 
                          stain_bin_max + work_bin_max + living_room + 
                          bedroom + white_nonhispanic_all, data = dfhh)
MeFOSAA_step_model = step(MeFOSAA_full_model, direction="both", trace=FALSE)
summary(MeFOSAA_step_model)

#PFUnA
PFUnA_full_model = lm(log_PFUnA_dust_mean ~ site + age_mean + years_at_home_mean + cleaning_bin_max + 
                        stain_bin_max + work_bin_max + living_room + 
                        bedroom + white_nonhispanic_all, data = dfhh)
PFUnA_step_model = step(PFUnA_full_model, direction="both", trace=FALSE)
summary(PFUnA_step_model)
TukeyHSD(aov(PFUnA_step_model), which=c("site")) 


# sum PFAS 
sum_full_model = lm(log_sum_PFAS_dust_mean ~ site + age_mean + years_at_home_mean + cleaning_bin_max + 
                      stain_bin_max + work_bin_max + living_room + 
                      bedroom + white_nonhispanic_all, data = dfhh)
sum_step_model = step(sum_full_model, direction="both", trace=FALSE)
summary(sum_step_model)

