# serum_lme_modeling.R
#
# Code to analyze associations between dust PFAS levels and serum PFAS levels using
# linear mixed effects modeling to account for the grouping of participants into
# households. 

library(psych)
library(dplyr)
library(ggplot2)
library(tidyr)
library(lme4)
library(car)
library(MuMIn)

### Data loading and preprocessing

# load pre-processed individual-level data with non-detects as DL / sqrt(2)
df = read.csv('output/data_all_sites_questionnaire.csv', stringsAsFactors=TRUE, na.strings= c('NA', ''))
PFAS = c('PFOA', 'PFOS', 'PFNA', 'PFDA', 'PFUnA', 'PFHxS', 'MeFOSAA')

# load pre-processed individual-level data with non-detects as NA
df.nd = read.csv('output/data_all_sites_below_DLs.csv', stringsAsFactors=TRUE)

# make unique household id column
df$hid = paste(df$site, df$householdid, sep="")

# calculate log of participant age
df$log_age <- log(df$age)

# serum and dust concentrations are log-normally distributed, so convert to log space
for (substance in PFAS) {
  df[paste('log_', substance, '_serum', sep="")] = log(df[paste(substance, '_serum_ng.ml', sep="")])
  df[paste('log_', substance, '_dust', sep="")] = log(df[paste(substance, '_dust_ng.g', sep="")])
}

# calculate sum of PFAS in serum + dust 
df$sum_PFAS_serum_ng.ml <- df$PFOA_serum_ng.ml + df$PFOS_serum_ng.ml + df$PFHxS_serum_ng.ml + df$PFDA_serum_ng.ml + df$PFNA_serum_ng.ml +
  df$MeFOSAA_serum_ng.ml + df$PFUnA_serum_ng.ml
df$sum_PFAS_dust_ng.g <- df$PFOA_dust_ng.g + df$PFOS_dust_ng.g + df$PFHxS_dust_ng.g + df$PFDA_dust_ng.g + df$PFNA_dust_ng.g +
  df$MeFOSAA_dust_ng.g + df$PFUnA_dust_ng.g
df$log_sum_PFAS_serum = log(df$sum_PFAS_serum_ng.ml)
df$log_sum_PFAS_dust = log(df$sum_PFAS_dust_ng.g)

# stain resistant product use variable
df$stain_bin = ifelse(df$AQ15_StainResist == "Never", 0, 1)
df$stain_bin[is.na(df$stain_bin)] = 0 

# create data set without 1 row with NA age
df_complete = df[complete.cases(df$age),]


### serum % detection

detect_freq_df = as.data.frame(matrix(nrow = length(PFAS), ncol = 2))
colnames(detect_freq_df) <- c("pfas_substance","freq")
for (i in 1:length(PFAS)) {
  substance = PFAS[i]
  data = df.nd[paste(substance, '_serum_ng.ml', sep="")]
  na_count = sum(is.na(data))
  participant_n = nrow(df.nd) #
  detect_freq = 100 - (na_count / participant_n * 100)
  detect_freq_df$pfas_substance[i] = substance
  detect_freq_df$freq[i] = detect_freq
}
detect_freq_df


### plots of raw log serum vs log dust
plot(log_PFOA_serum ~ log_PFOA_dust, data=df, pch=16)
plot(log_PFOS_serum ~ log_PFOS_dust, data=df, pch=16)
plot(log_PFHxS_serum ~ log_PFHxS_dust, data=df, pch=16)
plot(log_PFNA_serum~ log_PFNA_dust, data=df, pch=16)
plot(log_PFDA_serum ~ log_PFDA_dust, data=df, pch=16)
plot(log_MeFOSAA_serum ~ log_MeFOSAA_dust, data=df, pch=16)
plot(log_PFUnA_serum ~ log_PFUnA_dust, data=df, pch=16)


### Are PFAS serum levels correlated with the age of participants?
# note: random effect not needed as age and serum levels are both measured
#       at the individual level

# analyze age ~ serum relationship with Pearson's test
cor.test(df_complete$log_PFOA_serum, df_complete$age, method = 'pearson') # sig (p=0.002)
cor.test(df_complete$log_PFOS_serum, df_complete$age, method = 'pearson') # sig (p<0.001)
cor.test(df_complete$log_PFHxS_serum, df_complete$age, method = 'pearson') # sig (p=0.002)
cor.test(df_complete$log_PFNA_serum, df_complete$age, method = 'pearson') # sig (p<0.001)
cor.test(df_complete$log_PFDA_serum, df_complete$age, method = 'pearson') # sig (p=0.007)
cor.test(df_complete$log_MeFOSAA_serum, df_complete$age, method = 'pearson') # sig (p<0.001)
cor.test(df_complete$log_PFUnA_serum, df_complete$age, method = 'pearson') # sig (p<0.001)


### Are dust PFAS levels in homes correlated with the (mean) age of residents?
# note: random effect not needed as mean household age and dust levels
#       are both measured at the household level

# aggregate data by household, take mean age of residents in house
dfhh_age = df[ , grepl( "dust|hid|site|age" , names( df ) ) ]
dfhh_age = as.data.frame(dfhh_age %>% 
                           group_by(hid, site) %>% 
                           summarise(across(everything(), mean)))
nrow(dfhh_age) # 114

# analyze age ~ dust relationship with Pearson's test
cor.test(dfhh_age$log_PFOA_dust, dfhh_age$age, method = 'pearson') # sig (p<0.001)
cor.test(dfhh_age$log_PFOS_dust, dfhh_age$age, method = 'pearson') # sig (p<0.001)
cor.test(dfhh_age$log_PFHxS_dust, dfhh_age$age, method = 'pearson') # sig (p=0.01)
cor.test(dfhh_age$log_PFNA_dust, dfhh_age$age, method = 'pearson') # sig (p=0.006)
cor.test(dfhh_age$log_PFDA_dust, dfhh_age$age, method = 'pearson') # sig (p<0.001)
cor.test(dfhh_age$log_MeFOSAA_dust, dfhh_age$age, method = 'pearson') # sig (p<0.001)
cor.test(dfhh_age$log_PFUnA_dust, dfhh_age$age, method = 'pearson') # sig (p<0.001)
cor.test(dfhh_age$log_sum_PFAS_dust, dfhh_age$age, method = 'pearson') # sig (p<0.001)



### Analyze correlation between participant serum levels and household dust levels
### using linear mixed effects modeling - accounting for participant age


# PFOA

PFOA_age_mixed_model <- lmer(log_PFOA_serum ~ log_PFOA_dust + log_age + site + (1|hid), data=df_complete)
summary(PFOA_age_mixed_model)
drop1(PFOA_age_mixed_model, test="Chisq") # Non-significant (p=0.22)
hist(resid(PFOA_age_mixed_model)) # residuals roughly normal
plot(PFOA_age_mixed_model)
r.squaredGLMM(PFOA_age_mixed_model)  # marginal r2 = 0.50

# PFOS

PFOS_age_mixed_model <- lmer(log_PFOS_serum ~ log_PFOS_dust + log_age + site + (1|hid), data=df_complete)
summary(PFOS_age_mixed_model)
drop1(PFOS_age_mixed_model, test="Chisq") # non-significant dust effect (p=0.12)
hist(resid(PFOS_age_mixed_model)) # residuals roughly normal
plot(PFOS_age_mixed_model)
r.squaredGLMM(PFOS_age_mixed_model) # marginal r2 = 0.59

# PFHxS

PFHxS_age_mixed_model <- lmer(log_PFHxS_serum ~ log_PFHxS_dust + log_age + site + (1|hid), data=df_complete)
summary(PFHxS_age_mixed_model)
drop1(PFHxS_age_mixed_model, test="Chisq") # non-significant dust effect (p=0.09)
hist(resid(PFHxS_age_mixed_model)) # residuals roughly normal
plot(PFHxS_age_mixed_model)
r.squaredGLMM(PFHxS_age_mixed_model) # marginal r2 = 0.52

# PFNA

PFNA_age_mixed_model <- lmer(log_PFNA_serum ~ log_PFNA_dust + log_age + site +  (1|hid), data=df_complete)
summary(PFNA_age_mixed_model)
drop1(PFNA_age_mixed_model, test="Chisq") # Significant dust effect (p=0.002)
hist(resid(PFNA_age_mixed_model)) # residuals roughly normal
plot(PFNA_age_mixed_model)
r.squaredGLMM(PFNA_age_mixed_model) # marginal r2 = 0.38

# PFDA

PFDA_age_mixed_model <- lmer(log_PFDA_serum ~ log_PFDA_dust + log_age + site + (1|hid), data=df_complete)
summary(PFDA_age_mixed_model)
drop1(PFDA_age_mixed_model, test="Chisq") # non-significant dust effect (p=0.039)
hist(resid(PFDA_age_mixed_model)) # residuals roughly normal
plot(PFDA_age_mixed_model)
r.squaredGLMM(PFDA_age_mixed_model) # marginal r2 = 0.19

# PFUnA

PFUnA_age_mixed_model <- lmer(log_PFUnA_serum ~ log_PFUnA_dust + log_age + site +(1|hid), data=df_complete)
summary(PFUnA_age_mixed_model)
drop1(PFUnA_age_mixed_model, test="Chisq") # non-significant dust effect (p=0.055) - SR product use SIG
hist(resid(PFUnA_age_mixed_model)) # residuals roughly normal
plot(PFUnA_age_mixed_model)
r.squaredGLMM(PFUnA_age_mixed_model) # marginal r2 = 0.21

# MeFOSAA

MeFOSAA_age_mixed_model <- lmer(log_MeFOSAA_serum ~ log_MeFOSAA_dust + log_age + site + (1|hid), data=df_complete)
summary(MeFOSAA_age_mixed_model)
drop1(MeFOSAA_age_mixed_model, test="Chisq") # Significant dust effect (p<0.001)
hist(resid(MeFOSAA_age_mixed_model)) # residuals roughly normal
plot(MeFOSAA_age_mixed_model)
r.squaredGLMM(MeFOSAA_age_mixed_model) # marginal r2 = 0.25

# SUM of PFAS

sum_age_mixed_model <- lmer(log_sum_PFAS_serum ~ log_sum_PFAS_dust  + log_age + site + (1|hid), data=df_complete)
summary(sum_age_mixed_model)
drop1(sum_age_mixed_model, test="Chisq") # non-significant dust effect (p=0.07)
plot(sum_age_mixed_model)
r.squaredGLMM(sum_age_mixed_model) # marginal r2 = 0.57


#####
# PFAS dust load (g of chemical per floor area)
#####

# calculate PFAS dust load variables
for (substance in PFAS) {
  df[paste('log_', substance, '_dust_load', sep="")] = log(df[paste(substance, '_dust_ng.g', sep="")] * df$dust_load_g_m2)
}
df["log_sum_PFAS_dust_load"] = log(df['sum_PFAS_dust_ng.g'] * df$dust_load_g_m2)

# look at distribution of dust loading on floors
df$log_dust_load <- log(df$dust_load_g_m2)
hist(df$log_dust_load)
mean(df$log_dust_load)
sd(df$log_dust_load)
mean(df$dust_load_g_m2)
sd(df$dust_load_g_m2)

# are PFAS concentrations and dust loading correlated?
corr.test(df$log_dust_load, df$log_PFOA_dust) # not sig
corr.test(df$log_dust_load, df$log_PFOS_dust) # not sig
corr.test(df$log_dust_load, df$log_PFHxS_dust) # sig
corr.test(df$log_dust_load, df$log_PFNA_dust) # sig
corr.test(df$log_dust_load, df$log_PFDA_dust) # sig
corr.test(df$log_dust_load, df$log_PFUnA_dust) # sig
corr.test(df$log_dust_load, df$log_MeFOSAA_dust) # sig

plot(df$log_dust_load, df$log_PFHxS_dust) # sig
plot(df$log_dust_load, df$log_PFNA_dust) # sig
plot(df$log_dust_load, df$log_PFDA_dust) # sig
plot(df$log_dust_load, df$log_PFUnA_dust) # sig
plot(df$log_dust_load, df$log_MeFOSAA_dust) # sig

# is dust loading correlated with age?
mm <- lmer(age ~ log_dust_load + (1|hid), data=df) # not sig
drop1(mm, test="Chisq") # not sig

# is PFAS loading correlated with age?
corr.test(df$age, df$log_PFOA_dust_load) # sig
corr.test(df$age, df$log_PFOS_dust_load) # sig
corr.test(df$age, df$log_PFHxS_dust_load) # not sig
corr.test(df$age, df$log_PFNA_dust_load) # not sig
corr.test(df$age, df$log_PFDA_dust_load) # not sig
corr.test(df$age, df$log_PFUnA_dust_load) # not sig
corr.test(df$age, df$log_MeFOSAA_dust_load) # sig


### Analyze correlation between participant serum levels and household dust 
### PFAS load using linear mixed effects modeling - accounting for participant age

# remove the 1 row with missing age
dfa <- subset(df, !is.na(age))

# PFOA

mixed_model <- lmer(log_PFOA_serum ~ log_PFOA_dust_load + site + log_age + (1|hid), data=dfa)
summary(mixed_model)
drop1(mixed_model, test="Chisq") # no significant dust load effect
hist(resid(mixed_model)) # residuals roughly normal
plot(mixed_model)
r.squaredGLMM(mixed_model)  # marginal r2 = 0.50

# PFOS

mixed_model <- lmer(log_PFOS_serum ~ log_PFOS_dust_load + site + log_age + (1|hid), data=dfa)
summary(mixed_model)
drop1(mixed_model, test="Chisq") # no significant dust load effect
hist(resid(mixed_model)) # residuals roughly normal
plot(mixed_model)
r.squaredGLMM(mixed_model) # marginal r2 = 0.58

# PFHxS

mixed_model <- lmer(log_PFHxS_serum ~ log_PFHxS_dust_load + site + log_age + (1|hid), data=dfa)
summary(mixed_model)
drop1(mixed_model, test="Chisq") # no significant dust load effect
hist(resid(mixed_model)) # residuals roughly normal
plot(mixed_model)
r.squaredGLMM(mixed_model) # marginal r2 = 0.51

# PFNA

mixed_model <- lmer(log_PFNA_serum ~ log_PFNA_dust_load + site + log_age + (1|hid), data=dfa)
summary(mixed_model)
drop1(mixed_model, test="Chisq") # Significant dust effect (p=0.01)
hist(resid(mixed_model)) # residuals roughly normal
plot(mixed_model)
r.squaredGLMM(mixed_model) # marginal r2 = 0.37

# PFDA

mixed_model <- lmer(log_PFDA_serum ~ log_PFDA_dust_load + site + log_age + (1|hid), data=dfa)
summary(mixed_model)
drop1(mixed_model, test="Chisq") # no significant dust load effect
hist(resid(mixed_model)) # residuals roughly normal
plot(mixed_model)
r.squaredGLMM(mixed_model) # marginal r2 = 0.16

# PFUnA

mixed_model <- lmer(log_PFUnA_serum ~ log_PFUnA_dust_load + site + log_age + (1|hid), data=dfa)
summary(mixed_model)
drop1(mixed_model, test="Chisq") # no significant dust load effect
hist(resid(mixed_model)) # residuals roughly normal
plot(mixed_model)
r.squaredGLMM(mixed_model) # marginal r2 = 0.20

# MeFOSAA

mixed_model <- lmer(log_MeFOSAA_serum ~ log_MeFOSAA_dust_load + site + log_age + (1|hid), data=dfa)
summary(mixed_model)
drop1(mixed_model, test="Chisq") # Significant dust effect (p<0.001)
hist(resid(mixed_model)) # residuals roughly normal
plot(mixed_model)
r.squaredGLMM(mixed_model) # marginal r2 = 0.25

# SUM PFAS
mixed_model <- lmer(log_sum_PFAS_serum ~ log_sum_PFAS_dust_load  + site + log_age + (1|hid), data=dfa)
summary(mixed_model)
drop1(mixed_model, test="Chisq") # no significant dust load effect
plot(mixed_model)
r.squaredGLMM(mixed_model) # marginal r2 = 0.55



##### Plots for manuscript


# PFNA dust vs serum - figure A1
ggplot(df, aes(y = log_PFNA_serum, x = log_PFNA_dust)) + geom_point(size = 2, color='grey40',alpha=0.5) +
  theme_classic() + ggtitle("") +
  geom_abline(intercept=-1.33, slope=0.19, lwd=1, lty='dashed', color='black') + 
  xlab("PFNA in dust (log ng/g)") + ylab("PFNA in serum (log ng/ml)") +
  theme(axis.title = element_text(size = 14), 
        axis.text = element_text(size = 12),
        plot.title = element_text(size = 18, hjust = 0.5)) +
  annotate("text", x = 4.5, y = 1, label = "p == 0.002", parse = T, hjust = 0, size = 5) + 
  annotate("text", x = 4.5, y = 0.7, label = "R^2 == 0.38", parse = T, hjust = 0, size = 5)
ggsave("output/figures/Figure_A1.jpeg", width = 6, height = 5, units = "in", dpi = 300)

##### Generate Appendix table A1
PFAS_order = c('PFHxS', 'PFOA', 'PFOS', 'PFNA', 'PFDA', 'PFUnA', 'MeFOSAA', 'sum')

table_a1 <- list()
for(i in 1:length(PFAS_order)){
  model_obj <- eval(parse(text=paste(PFAS_order[i],"_age_mixed_model",sep="")))
  estimates <- sprintf("%.2f", round(coef(summary((model_obj)))[,'Estimate'],2))
  ses <- sprintf("%.2f", round(coef(summary((model_obj)))[,'Std. Error'],2))
  table_a1[[i]] <- paste(estimates,ses, sep=" ± ")
}
names(table_a1) <- PFAS_order
table_a1 <- do.call("cbind",table_a1)
rownames(table_a1) <- c("Intercept", "ln(dust ug/g)", "ln(age)","SiteWV", "SiteDE",
                    "SiteWA", "SiteTX","SiteNY", "SiteAK", "SiteCO")
write.csv(table_a1, "output/table_a1.csv")

