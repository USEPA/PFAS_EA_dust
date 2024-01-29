# analyze_questions.R
#
# Code to analyze associations between dust PFAS levels and exposure questionnaire
# responses using ANOVA.

library(psych)
library(dplyr)
library(ggplot2)
library(tidyr)
library(lme4)
library(car)
library(FSA)

# load pre-processed individual-level data with non-detects as DL / sqrt(2)
df = read.csv('output/data_all_sites_questionnaire.csv', stringsAsFactors=TRUE, na.strings= c('NA', ''))
PFAS = c('PFOA', 'PFOS', 'PFNA', 'PFDA', 'PFUnA', 'PFHxS', 'MeFOSAA')

# make unique household id column
df$hid = paste(df$site, df$householdid, sep="")
length(unique(df$hid))  # 114 households

# convert concs to log
PFAS = c('PFOA', 'PFOS', 'PFNA', 'PFDA', 'PFUnA', 'PFHxS', 'MeFOSAA')
for (substance in PFAS) {
  df[paste('log_', substance, '_dust', sep="")] = log(df[paste(substance, '_dust_ng.g', sep="")])
  df[paste('log_', substance, '_serum', sep="")] = log(df[paste(substance, '_serum_ng.ml', sep="")])
}


# calculate sum of PFAS in serum + dust 
df$sum_PFAS_serum_ng.ml <- df$PFOA_serum_ng.ml + df$PFOS_serum_ng.ml + df$PFHxS_serum_ng.ml + df$PFDA_serum_ng.ml + df$PFNA_serum_ng.ml +
  df$MeFOSAA_serum_ng.ml + df$PFUnA_serum_ng.ml
df$sum_PFAS_dust_ng.g <- df$PFOA_dust_ng.g + df$PFOS_dust_ng.g + df$PFHxS_dust_ng.g + df$PFDA_dust_ng.g + df$PFNA_dust_ng.g +
  df$MeFOSAA_dust_ng.g + df$PFUnA_dust_ng.g
df$log_sum_PFAS_serum = log(df$sum_PFAS_serum_ng.ml)
df$log_sum_PFAS_dust = log(df$sum_PFAS_dust_ng.g)

df$log_dust_load <- log(df$dust_load_g_m2)

# calculate dust load variables
for (substance in PFAS) {
  df[paste('log_', substance, '_dust_load', sep="")] = log(df[paste(substance, '_dust_ng.g', sep="")] * df$dust_load_g_m2)
}
df["log_sum_PFAS_dust_load"] = log(df['sum_PFAS_dust_ng.g'] * df$dust_load_g_m2)

###########
### basic demographics for the dataset
###########

summary(df$sex) # Female: 115, Male 87

df_complete = df[complete.cases(df$age),] # remove 1 age non-response
nrow(df_complete[df_complete$age < 18,]) # age < 18 : 17
nrow(df_complete[(df_complete$age >= 18) & (df_complete$age < 50),]) # age >= 18 but < 50: 57
nrow(df_complete[df_complete$age >= 50,]) # age >= 50: 127

# N for white, hispanic
nrow(df[(df$AQ1_Ethnicity %in% 'True') & (df$AQ2_White %in% 'True'),]) # Adults 3
nrow(df[(df$CQ2_Ethnicity%in%'True') & (df$AQ2_White%in%'True'),]) # Children 3

# N for white, non-hispanic
nrow(df[(df$AQ1_Ethnicity %in% 'False') & (df$AQ2_White %in% 'True'),]) # Adults 162
nrow(df[(df$CQ2_Ethnicity%in%'False') & (df$AQ2_White%in%'True'),]) # Children 13

# N for not white
nrow(df[df$AQ2_White %in% 'False',]) # Adults + children 17

#range of participants per house
max(table(df$hid)) # max 6
min(table(df$hid)) # min 1


###########
### analyze mean HH age (and max HH age)
###########

# aggregate within households
dfhh_age = df[ , grepl( "dust|hid|site|age" , names( df ) ) ]
dfhh_age = as.data.frame(dfhh_age %>% 
                             group_by(hid, site) %>% 
                             summarise(across(everything(), list(mean=mean,max=max))))
nrow(dfhh_age) # 114

plot(log_PFOS_dust_mean ~ age_mean, data=dfhh_age)
plot(log_PFOA_dust_mean ~ age_mean, data=dfhh_age)
plot(log_PFHxS_dust_mean ~ age_mean, data=dfhh_age)
plot(log_PFNA_dust_mean ~ age_mean, data=dfhh_age)
plot(log_PFDA_dust_mean ~ age_mean, data=dfhh_age)
plot(log_MeFOSAA_dust_mean ~ age_mean, data=dfhh_age)
plot(log_PFUnA_dust_mean ~ age_mean, data=dfhh_age)
plot(log_sum_PFAS_dust_mean ~ age_mean, data=dfhh_age)

Anova(lm(log_PFOS_dust_mean ~ age_mean + site, data=dfhh_age))  # sig
Anova(lm(log_PFOA_dust_mean ~ age_mean + site, data=dfhh_age))  # sig
Anova(lm(log_PFHxS_dust_mean ~ age_mean + site, data=dfhh_age))  # sig
Anova(lm(log_PFNA_dust_mean ~ age_mean + site, data=dfhh_age))  # sig
Anova(lm(log_PFDA_dust_mean ~ age_mean + site, data=dfhh_age))  # sig
Anova(lm(log_MeFOSAA_dust_mean ~ age_mean + site, data=dfhh_age))  # sig
Anova(lm(log_PFUnA_dust_mean ~ age_mean + site, data=dfhh_age))  # sig
Anova(lm(log_sum_PFAS_dust_mean ~ age_mean + site, data=dfhh_age))  # sig

Anova(lm(log_PFOS_dust_mean ~ age_max + site, data=dfhh_age))  # sig
Anova(lm(log_PFOA_dust_mean ~ age_max + site, data=dfhh_age))  # sig
Anova(lm(log_PFHxS_dust_mean ~ age_max + site, data=dfhh_age))  # sig
Anova(lm(log_PFNA_dust_mean ~ age_max + site, data=dfhh_age))  # sig
Anova(lm(log_PFDA_dust_mean ~ age_max + site, data=dfhh_age))  # sig
Anova(lm(log_MeFOSAA_dust_mean ~ age_max + site, data=dfhh_age))  # sig
Anova(lm(log_PFUnA_dust_mean ~ age_max + site, data=dfhh_age))  # sig
Anova(lm(log_sum_PFAS_dust_mean ~ age_max + site, data=dfhh_age))  # sig


###########
### analyze years at current address
###########

# aggregate within households
years_at_home_years = df$AQ3_Years
years_at_home_years[is.na(years_at_home_years)] = 0
years_at_home_months = df$AQ3_Months/12
years_at_home_months[is.na(years_at_home_months)] = 0
df$years_at_home = years_at_home_years + years_at_home_months
dfhh_years_at_home= df[ , grepl( "dust|hid|site|years|age" , names( df ) ) ]
dfhh_years_at_home = as.data.frame(dfhh_years_at_home %>% 
                        group_by(hid, site) %>% 
                        summarise(across(everything(), mean)))
nrow(dfhh_years_at_home) # 114

# PFOS
mod = lm(log_PFOS_dust ~ years_at_home + site, data=dfhh_years_at_home)
Anova(mod)  # ANOVA sig +;
hist(resid(mod))
plot(log_PFOS_dust ~ years_at_home, data=dfhh_years_at_home)


# PFOA
mod = lm(log_PFOA_dust ~ years_at_home + site, data=dfhh_years_at_home)
Anova(mod)  # ANOVA sig +
hist(resid(mod))
plot(log_PFOA_dust ~ years_at_home, data=dfhh_years_at_home)

# PFHxS
mod = lm(log_PFHxS_dust ~ years_at_home + site, data=dfhh_years_at_home)
Anova(mod)  # ANOVA sig +
hist(resid(mod))
plot(log_PFHxS_dust ~ years_at_home, data=dfhh_years_at_home)

# PFNA
mod = lm(log_PFNA_dust ~ years_at_home + site, data=dfhh_years_at_home)
Anova(mod)  # ANOVA sig +
hist(resid(mod))
plot(log_PFNA_dust ~ years_at_home, data=dfhh_years_at_home)

# PFDA
mod = lm(log_PFDA_dust ~ years_at_home + site, data=dfhh_years_at_home)
Anova(mod)  # ANOVA SIG + 
hist(resid(mod))
plot(log_PFDA_dust ~ years_at_home, data=dfhh_years_at_home)

# MeFOSAA
mod = lm(log_MeFOSAA_dust ~ years_at_home + site, data=dfhh_years_at_home)
Anova(mod)  # ANOVA sig + 
hist(resid(mod))
plot(log_MeFOSAA_dust ~ years_at_home, data=dfhh_years_at_home)

# PFUnA
mod = lm(log_PFUnA_dust ~ years_at_home + site, data=dfhh_years_at_home)
Anova(mod)  # ANOVA NOT sig
hist(resid(mod))
plot(log_PFUnA_dust ~ years_at_home, data=dfhh_years_at_home)

# Sum of PFAS
mod = lm(log_sum_PFAS_dust ~ years_at_home + site, data=dfhh_years_at_home)
Anova(mod)  # ANOVA SIG
hist(resid(mod))

# age vs years at home (individual-level)
plot(years_at_home ~ age, data=df)
cor.test(df$years_at_home, df$age, method=c('pearson')) # SIG r=0.46

# age vs years at home (household-level)
cor.test(dfhh_years_at_home$years_at_home, dfhh_years_at_home$age, method=c('pearson')) # SIG r=0.51

# does years at home relate to dust load in the house
Anova(lm(log_dust_load ~ years_at_home + site, data=dfhh_years_at_home)) # no effect on dust load


### NOTES
# years lived at current home is sig. positively correlated with
# dust levels for 6 of 7 PFAS (All but PFUnA)
# also correlated for sum(PFAS) 
# May correlate with age of home or age of elements inside the home


###########
### analyze cleaning frequency
###########

summary(df$AQ14_Cleaning)
df$cleaning_bin = ifelse(df$AQ14_Cleaning == "Three times per week or more", 1, 0)
df$cleaning_bin[is.na(df$cleaning_bin)] = 0 
dfhh_clean = df[ , grepl( "dust|hid|site|age|cleaning_bin" , names( df ) ) ]
dfhh_clean = as.data.frame(dfhh_clean %>% 
                                 group_by(hid, site) %>% 
                             summarise(across(everything(), list(mean=mean,max=max))))
dfhh_clean$cleaning_bin_max = factor(dfhh_clean$cleaning_bin_max, labels=c('Low/Medium', 'High'))
nrow(dfhh_clean) # 114
summary(dfhh_clean$cleaning_bin_max) # 25 High frequency, 79 Medium/Low frequency

Anova(lm(log_PFOS_dust_mean ~ cleaning_bin_max + site, data=dfhh_clean)) # SIG (p=0.007)
Anova(lm(log_PFOA_dust_mean ~ cleaning_bin_max + site, data=dfhh_clean)) # not sig (p=0.12)
Anova(lm(log_PFHxS_dust_mean ~ cleaning_bin_max + site, data=dfhh_clean)) # SIG (p=0.01)
Anova(lm(log_PFNA_dust_mean ~ cleaning_bin_max + site, data=dfhh_clean)) # not sig (p=0.48)
Anova(lm(log_PFDA_dust_mean ~ cleaning_bin_max + site, data=dfhh_clean)) # not sig (p=0.58)
Anova(lm(log_MeFOSAA_dust_mean ~ cleaning_bin_max + site, data=dfhh_clean)) # SIG (p=0.001)
Anova(lm(log_PFUnA_dust_mean ~ cleaning_bin_max + site, data=dfhh_clean)) # not sig (p=0.97)
Anova(lm(log_sum_PFAS_dust_mean ~ cleaning_bin_max + site, data=dfhh_clean)) # SIG (p=0.04)

# does cleaning frequency affect dust load in the house
Anova(lm(log_dust_load_max ~ cleaning_bin_max + site, data=dfhh_clean)) # no effect on dust load


### NOTES
# 3 of 7 PFAS (PFOS, PFHxS, MeFOSAA) are correlated with household cleaning frequency
# Also sum(PFAS) is correlated with cleaning frequency


# cleaning freq vs age
plot(age_mean ~ factor(cleaning_bin_max), data=dfhh_clean)
summary(aov(age_mean ~ cleaning_bin_max + site, data=dfhh_clean)) # households with older mean age cleaned less

# cleaning frequency effect sizes
(exp(3.07385-0.96652) - exp(3.07385)) / exp(3.07385) * 100 # PFOS
(exp(2.1734-0.9706) - exp(2.1734)) / exp(2.1734) * 100 # PFHxS
(exp(1.77760-1.08013) - exp(1.77760)) / exp(1.77760) * 100 # MeFOSAA
(exp(4.46286-0.66427) - exp(4.46286)) / exp(4.46286) * 100 # SUM PFAS


###########
### analyze stain resistant products
###########

summary(df$AQ15_StainResist)
df$stain_bin = ifelse(df$AQ15_StainResist == "Never", 0, 1)
df$stain_bin[is.na(df$stain_bin)] = 0 
dfhh_stain = df[ , grepl( "dust|hid|site|stain_bin|age" , names( df ) ) ]
dfhh_stain = as.data.frame(dfhh_stain %>% 
                             group_by(hid, site) %>% 
                             summarise(across(everything(), list(mean=mean,max=max))))
dfhh_stain$stain_bin_max = factor(dfhh_stain$stain_bin_max, labels=c('No','Yes'))
nrow(dfhh_stain) # 114
summary(dfhh_stain$stain_bin_max) # yes = 15, no = 99

# ANOVAs
Anova(lm(log_PFOS_dust_mean ~ stain_bin_max + site, data=dfhh_stain))
Anova(lm(log_PFOA_dust_mean ~ stain_bin_max + site, data=dfhh_stain))
Anova(lm(log_PFHxS_dust_mean ~ stain_bin_max + site, data=dfhh_stain))
Anova(lm(log_PFNA_dust_mean ~ stain_bin_max + site, data=dfhh_stain))  # SIG (p=0.002)
Anova(lm(log_PFDA_dust_mean ~ stain_bin_max + site, data=dfhh_stain))  # SIG (p=0.011)
Anova(lm(log_MeFOSAA_dust_mean ~ stain_bin_max + site, data=dfhh_stain))
Anova(lm(log_PFUnA_dust_mean ~ stain_bin_max + site, data=dfhh_stain))  # SIG (p<0.001)
Anova(lm(log_sum_PFAS_dust_mean ~ stain_bin_max + site, data=dfhh_stain)) 


# age vs stain resistant product use
summary(aov(age_mean ~ stain_bin_max, data=dfhh_stain)) # no rel. between age and SR product use

# stain resistant product use effect sizes
(exp(1.21980+0.85496) - exp(1.21980)) / exp(1.21980) * 100 # PFNA
(exp(1.37293+0.67148 ) - exp(1.37293)) / exp(1.37293) * 100 # PFDA
(exp(0.57252+0.99208 ) - exp(0.57252)) / exp(0.57252) * 100 # PFUnA

### NOTES
# 3 of 7 PFAS (PFNA, PFDA, PFUnA) are correlated with stain resistant product use)
# sum of PFAS not correlated


###########
### analyze flooring
###########

# create flooring household df
dfhh_floors = df[ , grepl( "dust|hid|site|Floor|age" , names( df ) ) ]
dfhh_floors = as.data.frame(dfhh_floors %>% 
                             group_by(hid, site) %>% 
                             summarise(across(everything(), list(mean=mean,first=first))))  
nrow(dfhh_floors) # 114
# taking first listed per hh to be the flooring
dfhh_floors$AQ16_FloorLR = factor(dfhh_floors$AQ16_FloorLR_first)
dfhh_floors$AQ17_FloorK = factor(dfhh_floors$AQ17_FloorK_first)
dfhh_floors$AQ18_FloorB = factor(dfhh_floors$AQ18_FloorB_first)
summary(dfhh_floors$AQ16_FloorLR)
summary(dfhh_floors$AQ17_FloorK)
summary(dfhh_floors$AQ18_FloorB)

#  aggregate to hard vs soft flooring (except kitchen since no soft type there)
dfhh_floors$living_room <- recode(dfhh_floors$AQ16_FloorLR, "c('Hardwood', 'Laminate',
                                  'Tile', 'Vinyl')='hard';c('Carpet') = 'soft'")
dfhh_floors$kitchen <- dfhh_floors$AQ17_FloorK_first # kitchen has no soft type
dfhh_floors$bedroom <- recode(dfhh_floors$AQ18_FloorB_first, "c('Hardwood', 'Laminate',
                                  'Tile', 'Vinyl')='hard';c('Carpet') = 'soft'")
dfhh_floors[dfhh_floors$living_room == "Other","living_room"] <- NA
dfhh_floors[dfhh_floors$bedroom == "Other", "bedroom"] <- NA
dfhh_floors$living_room <- droplevels(dfhh_floors$living_room)
dfhh_floors$bedroom <- droplevels(dfhh_floors$bedroom)

summary(dfhh_floors$living_room) # hard = 53, soft = 60, NA = 1
summary(dfhh_floors$bedroom) # hard = 33, soft = 79, NA = 2
        

### Living Room
Anova(lm(log_PFOS_dust_mean ~ living_room + site, data=dfhh_floors))
Anova(lm(log_PFOA_dust_mean ~ living_room + site, data=dfhh_floors))
Anova(lm(log_PFHxS_dust_mean ~ living_room + site, data=dfhh_floors))
Anova(lm(log_PFNA_dust_mean ~ living_room + site, data=dfhh_floors))
Anova(lm(log_PFDA_dust_mean ~ living_room + site, data=dfhh_floors))  # SIG Soft > Hard
TukeyHSD(aov(log_PFDA_dust_mean ~ living_room + site, data=dfhh_floors)) 
Anova(lm(log_MeFOSAA_dust_mean ~ living_room + site, data=dfhh_floors))
Anova(lm(log_PFUnA_dust_mean ~ living_room + site, data=dfhh_floors))  # SIG soft > Hard
TukeyHSD(aov(log_PFUnA_dust_mean ~ living_room + site, data=dfhh_floors))
Anova(lm(log_sum_PFAS_dust_mean ~ living_room + site, data=dfhh_floors))


# LR vs age
summary(aov(age_mean ~ living_room, data=dfhh_floors)) # SIG
TukeyHSD(aov(age_mean ~ living_room, data=dfhh_floors))# older household more likely to have carpet than Vinyl
plot(age_mean ~ AQ16_FloorLR, data=dfhh_floor_LR)

# effect sizes
(exp(1.19458+0.55) - exp(1.19458)) / exp(1.19458) * 100 # PFDA
(exp(0.489912+0.423566 ) - exp(0.489912)) / exp(0.489912) * 100 # PFDA

### Kitchen - not examined because only 2 households had soft type flooring


### Bedroom

Anova(lm(log_PFOS_dust_mean ~ bedroom + site, data=dfhh_floors))
Anova(lm(log_PFOA_dust_mean ~ bedroom + site, data=dfhh_floors))
Anova(lm(log_PFHxS_dust_mean ~ bedroom + site, data=dfhh_floors))
Anova(lm(log_PFNA_dust_mean ~ bedroom + site, data=dfhh_floors))
Anova(lm(log_PFDA_dust_mean ~ bedroom + site, data=dfhh_floors))
Anova(lm(log_MeFOSAA_dust_mean ~ bedroom + site, data=dfhh_floors))
Anova(lm(log_PFUnA_dust_mean ~ bedroom + site, data=dfhh_floors))
Anova(lm(log_sum_PFAS_dust_mean ~ bedroom + site, data=dfhh_floors))


# bedroom vs age
summary(aov(age_mean ~ AQ18_FloorB, data=dfhh_floor_B)) # not sig
plot(age_mean ~ AQ18_FloorB, data=dfhh_floor_B)


# does floor type relate to dust load in the house
Anova(lm(log_dust_load_mean ~ bedroom + site, data=dfhh_floors)) # SIG effect on dust load
Anova(lm(log_dust_load_mean ~ living_room + site, data=dfhh_floors)) # SIG effect on dust load

plot(log_dust_load_mean ~ bedroom, data=dfhh_floors) # SIG effect on dust load
plot(log_dust_load_mean ~ living_room, data=dfhh_floors) # SIG effect on dust load


###########
### analyze occupational exposure
###########

df$work_bin = ifelse(df$AQ28_None == 'False', 1, 0)
df$work_bin[is.na(df$work_bin)] = 0 
df$work_bin[df$AQ28_DontKnow == 'True'] = 0 
dfhh_work = df[ , grepl( "dust|hid|site|work_bin|age" , names( df ) ) ]
dfhh_work = as.data.frame(dfhh_work %>% 
                            group_by(hid, site) %>% 
                            summarise(across(everything(), list(mean=mean,max=max))))
nrow(dfhh_work) # 114
dfhh_work$work_bin_max = factor(dfhh_work$work_bin_max, labels=c('no', 'yes')) # indicates if anyone in hh has worked in those industries
summary(dfhh_work$work_bin_max) # yes = 16, no = 98

plot(log_PFOS_dust_mean ~ work_bin_max, data=dfhh_work)
plot(log_PFOA_dust_mean ~ work_bin_max, data=dfhh_work)
plot(log_PFHxS_dust_mean ~ work_bin_max, data=dfhh_work)
plot(log_PFNA_dust_mean ~ work_bin_max, data=dfhh_work)
plot(log_PFDA_dust_mean ~ work_bin_max, data=dfhh_work)
plot(log_MeFOSAA_dust_mean ~ work_bin_max, data=dfhh_work)
plot(log_PFUnA_dust_mean ~ work_bin_max, data=dfhh_work)

Anova(lm(log_PFOS_dust_mean ~ work_bin_max + site, data=dfhh_work)) # not sig
Anova(lm(log_PFOA_dust_mean ~ work_bin_max + site, data=dfhh_work)) # not sig
Anova(lm(log_PFHxS_dust_mean ~ work_bin_max + site, data=dfhh_work)) # not sig
Anova(lm(log_PFNA_dust_mean ~ work_bin_max + site, data=dfhh_work)) # not sig
Anova(lm(log_PFDA_dust_mean ~ work_bin_max + site, data=dfhh_work)) # not sig  
Anova(lm(log_MeFOSAA_dust_mean ~ work_bin_max + site, data=dfhh_work)) # not sig
Anova(lm(log_PFUnA_dust_mean ~ work_bin_max + site, data=dfhh_work)) # not sig
Anova(lm(log_sum_PFAS_dust_mean ~ work_bin_max + site, data=dfhh_work)) # not sig


### NOTES
# No significant trends between PFAS levels in dust and someone in HH
#   working in industries with possible exposure.

# correlated with age? no
summary(aov(age_mean ~ work_bin_max, data=dfhh_work))


###########
### race&ethnicity 
###########


# define households with at least 1 non-white or white, Hispanic participant

df$white_nonhispanic = ifelse((df$AQ2_White %in% 'True') & 
                            !((df$AQ1_Ethnicity %in% 'True') | (df$CQ2_Ethnicity %in% 'True')), 1,0)
dfhh_race = df[ , grepl( "dust|hid|site|white_nonhispanic|age" , names( df ) ) ]
dfhh_race = as.data.frame(dfhh_race %>% 
                            group_by(hid, site) %>% 
                            summarise(across(everything(), list(mean=mean,all=min))))
nrow(dfhh_race) # n = 114
dfhh_race$white_nonhispanic_all = factor(dfhh_race$white_nonhispanic_all, labels=c('no', 'yes'))
summary(dfhh_race$white_nonhispanic_all) # yes (all participants white, non-hispanic) = 103, no = 13

Anova(lm(log_PFOS_dust_mean ~ white_nonhispanic_all + site, data=dfhh_race)) # not sig
Anova(lm(log_PFOA_dust_mean ~ white_nonhispanic_all + site, data=dfhh_race)) # not sig
Anova(lm(log_PFHxS_dust_mean ~ white_nonhispanic_all + site, data=dfhh_race)) # not sig
Anova(lm(log_PFNA_dust_mean ~ white_nonhispanic_all + site, data=dfhh_race)) # not sig
Anova(lm(log_PFDA_dust_mean ~ white_nonhispanic_all + site, data=dfhh_race)) # not sig
Anova(lm(log_MeFOSAA_dust_mean ~ white_nonhispanic_all + site, data=dfhh_race)) # not sig
Anova(lm(log_PFUnA_dust_mean ~ white_nonhispanic_all + site, data=dfhh_race)) # not sig
Anova(lm(log_sum_PFAS_dust_mean ~ white_nonhispanic_all + site, data=dfhh_race)) # not sig


### NOTES
# No significant trends between PFAS levels in dust and someone in HH
#   identifying as something other than white and non-Hispanic 
#   But there were only 13/114 households with someone identifying as 
#   anything other than white and non-Hispanic.

####################


###########
### Cleaning frequency manuscript plots
###########

# reform dataframe from wide to long
log_dust_cols <- c('log_PFOA_dust_mean', 'log_PFOS_dust_mean', 'log_PFHxS_dust_mean',
                  'log_PFNA_dust_mean', 'log_PFDA_dust_mean', 'log_PFUnA_dust_mean',
                  'log_MeFOSAA_dust_mean', 'log_sum_PFAS_dust_mean')
df_clean_long <- melt(dfhh_clean, id.vars=c('site', 'cleaning_bin_max'), measure.vars=log_dust_cols,
     variable.name='chemical', value.name='dust_ng.g')
chemical_labels <- c('PFOA', 'bold(PFOS)', 'bold(PFHxS)', 'PFNA', 'PFDA', 'PFUnA', 'bold(MeFOSAA)', 'bold("Sum PFAS")')
df_clean_long$chemical_labeled <- factor(df_clean_long$chemical, labels=chemical_labels)
df_clean_long$cleaning_level = factor(df_clean_long$cleaning_bin_max, labels=c('Low/Med', 'High'))
#sig_chemicals <- c('log_PFOS_dust_mean', 'log_PFHxS_dust_mean',
#                   'log_MeFOSAA_dust_mean', 'log_sum_PFAS_dust_mean')


plot_cleaning <- ggplot(df_clean_long[complete.cases(df_clean_long$cleaning_level),], 
            aes(cleaning_level, dust_ng.g, fill=cleaning_level)) + facet_wrap(~chemical_labeled,nrow=1, labeller = label_parsed) +
  geom_boxplot() + scale_fill_manual(values=c("grey90", "grey60")) +#ylim(0,200) 
  xlab("Cleaning frequency") + ylab("Concentration in dust (log ng/g)") +
  theme_classic() + theme(axis.title = element_text(size = 14), 
                          axis.text = element_text(size = 12),
                          plot.title = element_text(size = 18, hjust = 0.5),
                          legend.position = "none",
                          panel.background = element_rect(fill = NA, color = "black"),
                          strip.text.x = element_text(size = 14))


###########
### SR product use manuscript plots
###########

log_dust_cols <- c('log_PFOA_dust_mean', 'log_PFOS_dust_mean', 'log_PFHxS_dust_mean',
                   'log_PFNA_dust_mean', 'log_PFDA_dust_mean', 'log_PFUnA_dust_mean',
                   'log_MeFOSAA_dust_mean', 'log_sum_PFAS_dust_mean')
df_stain_long <- melt(dfhh_stain, id.vars=c('site', 'stain_bin_max'), measure.vars=log_dust_cols,
                      variable.name='chemical', value.name='dust_ng.g')
chemical_labels <- c('PFOA', 'PFOS', 'PFHxS', 'bold(PFNA)', 'bold(PFDA)', 'bold(PFUnA)', 'MeFOSAA', '"Sum PFAS"')
df_stain_long$chemical_labeled <- factor(df_stain_long$chemical, labels=chemical_labels)
#sig_chemicals <- c('log_PFNA_dust_mean', 'log_PFDA_dust_mean',
#                   'log_PFUnA_dust_mean')
label_shade_df <- data.frame(
  var = log_dust_cols,
  var_color = c("white", "white", "white", "grey50", "grey50", "grey50", "white", "white")
)


plot_stain <- ggplot(df_stain_long[complete.cases(df_stain_long$stain_bin_max),]) + 
  facet_wrap(~chemical_labeled,nrow=1, labeller = label_parsed) +
  geom_boxplot(aes(stain_bin_max, dust_ng.g, fill=stain_bin_max)) + scale_fill_manual(values=c("grey90", "grey60")) +#ylim(0,200) 
  xlab("Use of stain-resistant products") + ylab("Concentration in dust (log ng/g)") +
  theme_classic() + theme(axis.title = element_text(size = 14), 
                          axis.text = element_text(size = 12),
                          plot.title = element_text(size = 18, hjust = 0.5),
                          legend.position = "none",
                          panel.background = element_rect(fill = NA, color = "black"),
                          strip.text.x = element_text(size = 14))

###########
### combine cleaning frequency and SR product use plots
########### 

output = plot_grid(plot_cleaning, plot_stain,align='h',nrow=2, labels=c("a)", "b)"),scale=0.93)# + 
  #draw_label("Concentration in dust (log ng/g)", x=  0.03, y=0.5, vjust= 1.5, angle=90)
ggsave("output/figures/Figure_5.jpeg", output, width = 12, height = 7, units = "in", dpi = 600)
