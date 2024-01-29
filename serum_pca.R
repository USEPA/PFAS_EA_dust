# serum_pca.R
#
# Code for multivariate analysis of serum PFAS levels using principal component
# analysis (PCA)

library(dplyr)
library(ggplot2)
library(cowplot)
library(ggfortify)
require(ggrepel)
library(tidyr)
library(lme4)
library(car)
library(MuMIn)

### Data loading and preprocessing

# load pre-processed individual-level data with non-detects as DL / sqrt(2)
df = read.csv('output/data_all_sites.csv', stringsAsFactors=TRUE)
PFAS = c('PFOA', 'PFOS', 'PFNA', 'PFDA', 'PFUnA', 'PFHxS', 'MeFOSAA')

# load pre-processed individual-level data with non-detects as NA
df.nd = read.csv('output/data_all_sites_below_DLs.csv', stringsAsFactors=TRUE)

# make unique household id column
df$hid = paste(df$site, df$householdid, sep="")

# create dataset without 1 row with NA age
df_complete = df[complete.cases(df$age),]

# serum and dust concentrations are log-normally distributed, so convert to log space
for (substance in PFAS) {
  df[paste('log_', substance, '_serum', sep="")] = log(df[paste(substance, '_serum_ng.ml', sep="")])
  df[paste('log_', substance, '_dust', sep="")] = log(df[paste(substance, '_dust_ng.g', sep="")])
}
for (substance in PFAS) {
  df_complete[paste('log_', substance, '_serum', sep="")] = log(df_complete[paste(substance, '_serum_ng.ml', sep="")])
  df_complete[paste('log_', substance, '_dust', sep="")] = log(df_complete[paste(substance, '_dust_ng.g', sep="")])
}

# calculate sum of PFAS in serum + dust 
df$sum_PFAS_serum_ng.ml <- df$PFOA_serum_ng.ml + df$PFOS_serum_ng.ml + df$PFHxS_serum_ng.ml + df$PFDA_serum_ng.ml + df$PFNA_serum_ng.ml +
  df$MeFOSAA_serum_ng.ml + df$PFUnA_serum_ng.ml
df$sum_PFAS_dust_ng.g <- df$PFOA_dust_ng.g + df$PFOS_dust_ng.g + df$PFHxS_dust_ng.g + df$PFDA_dust_ng.g + df$PFNA_dust_ng.g +
  df$MeFOSAA_dust_ng.g + df$PFUnA_dust_ng.g
df$log_sum_PFAS_serum = log(df$sum_PFAS_serum_ng.ml)
df$log_sum_PFAS_dust = log(df$sum_PFAS_dust_ng.g)

df_complete$sum_PFAS_serum_ng.ml <- df_complete$PFOA_serum_ng.ml + df_complete$PFOS_serum_ng.ml + df_complete$PFHxS_serum_ng.ml + df_complete$PFDA_serum_ng.ml + df_complete$PFNA_serum_ng.ml +
  df_complete$MeFOSAA_serum_ng.ml + df_complete$PFUnA_serum_ng.ml
df_complete$sum_PFAS_dust_ng.g <- df_complete$PFOA_dust_ng.g + df_complete$PFOS_dust_ng.g + df_complete$PFHxS_dust_ng.g + df_complete$PFDA_dust_ng.g + df_complete$PFNA_dust_ng.g +
  df_complete$MeFOSAA_dust_ng.g + df_complete$PFUnA_dust_ng.g
df_complete$log_sum_PFAS_serum = log(df_complete$sum_PFAS_serum_ng.ml)
df_complete$log_sum_PFAS_dust = log(df_complete$sum_PFAS_dust_ng.g)


#### PCA of serum data

serum.pca <- prcomp(df[,c("log_PFOA_serum", "log_PFOS_serum", "log_PFHxS_serum", "log_PFNA_serum", "log_PFDA_serum",
                          "log_MeFOSAA_serum", "log_PFUnA_serum")], center=T, scale=T)
serum.pca$x[,'PC1']
serum.pca$rotation[,'PC1'] <- serum.pca$rotation[,'PC1'] * -1  # flip direction of PC1 for clarity
rownames(serum.pca$rotation) <- c('PFOA', 'PFOS', 'PFHxS', 'PFNA', 'PFDA', 'MeFOSAA', 'PFUnA')
summary(serum.pca)
df.pca <- cbind(df, serum.pca$x)


### plot serum loadings

#PC1 vs PC2 biplot
autoplot(serum.pca, size=2, colour='grey40',data=df,loadings = TRUE, lwd=5,loadings.colour = 'black',loadings.size=2,loadings.label.colour="black",
         loadings.label = TRUE, loadings.label.size = 5, loadings.label.repel=T, frame=F) + theme_classic(base_size=15)
ggsave("output/figures/serum_biplot.jpeg", width = 6, height = 5, units = "in", dpi = 300)


#PC1 vs PC2 with site labeled
autoplot(serum.pca, size=2, data=df, colour='site',loadings = TRUE, loadings.colour = 'black',loadings.label.colour="black",
         loadings.label = TRUE, loadings.label.size = 5,frame=F) + theme_classic()

#PC1 vs PC3
autoplot(serum.pca, x=2,y=3,size=2, data=df, colour='site',loadings = TRUE, loadings.colour = 'black',loadings.label.colour="black",
         loadings.label = TRUE, loadings.label.size = 5,frame=F) + theme_classic()


### regressions of total dust levels vs PCs

#PC1 - with age
df.pca.complete <- df.pca[complete.cases(df.pca$age),]
mixed_model <- lmer(PC1 ~ log_sum_PFAS_dust + log(age) + site + (1|hid), data=subset(df.pca.complete))
summary(mixed_model)
drop1(mixed_model, test="Chisq") # Significant dust effect (p=0.01), also age effect
r.squaredGLMM(mixed_model) # R2 = 0.48

#PC2 - with age
mixed_model <- lmer(PC2 ~ log_sum_PFAS_dust + log(age) + site + (1|hid), data=subset(df.pca.complete))
summary(mixed_model)
drop1(mixed_model, test="Chisq") # Significant dust effect (p=0.02), and no age effect
r.squaredGLMM(mixed_model) # R2 = 0.33


#### Figures for manuscript

# Figure 4 - (a) biplot of seurm PCA, (b) PC2 vs sum PFAS dust

fig4a <- autoplot(serum.pca, size=2, scale=0, colour='grey40',alpha=0.5,data=df,loadings = TRUE, lwd=5,loadings.colour = 'black',loadings.size=2,loadings.label.colour="black",
                  loadings.label = TRUE, loadings.label.size = 5, loadings.label.repel=T, loadings.label.fontface='bold',frame=F) + 
  theme_classic(base_size=15) +
  annotate("text", x = -5.3, y = 2.8, label="a)", size=5)

fig4b <- ggplot(df.pca, aes(y = PC2, x = log_sum_PFAS_dust)) + geom_point(size = 2,colour='grey40',alpha=0.5) +
  theme_classic(base_size=15) + 
  xlab("Sum PFAS in dust (log ng/g)") + ylab("PC2") +
  geom_abline(intercept=0.89665, slope=-0.23257, lwd=1, lty='dashed', color='black') + 
  theme(axis.title = element_text(size = 14), 
        axis.text = element_text(size = 12),
        plot.title = element_text(size = 18, hjust = 0.5)) +
  annotate("text", x = 8.9, y = 2.5, label = "p == 0.02", parse = T, hjust = 0, size = 5) + 
  annotate("text", x = 8.9, y = 2.1, label = "R^2 == 0.33", parse = T, hjust = 0, size = 5) +
  annotate("text", x = 1.7, y = 2.8, label="b)", size=5)

ggsave("output/figures/figure_4.jpeg", plot_grid(fig4a, fig4b,align='h'),width = 12, height = 5, units = "in", dpi = 300)

