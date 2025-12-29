# Load packages
# To install packages before loading, use: install.package('nameofpackage')
library(GGally)
library(ggplot2)
library(carData)
library(scales)
library(car)
library(tidyverse)
library (aod)
library(knitr)
library(broom)
library(lme4)
library(nlme)
library(dplyr)
library(stats)
library(emmeans)
library(forcats)
library(ggsignif)

# Set directory where you would like this script and outputs of this script to be saved
setwd("~/data/NCfiles/BrainMetabolism/ImagingNeuroscience/")

# Import data (these files must be in the same directory as your present working directory set above)
rawdata <- read.csv(file = 'BrainMetabolismMasterDataSheet_WithAgeAndSex.csv')
gtv_invol <- read.csv(file = 'GTVInvolvedBrainRegions_Patient.csv')
SegID <- read.csv(file='SegID_Volumes.csv')
SegID_class <- read.csv(file='SegID_WMGMid.csv')
SegID_class$SegID <- factor(SegID_class$SegID) 

# Source custom functions (these files must be in the same directory as your present working directory set above)
source("generate_regression_formulas.R")
source("pval_to_asterisk.R")

# Remove ventricles from rawdata based on the SLANT region code
rawdata = rawdata[rawdata$SegID != 4 & rawdata$SegID != 11 & rawdata$SegID != 49 & 
                    rawdata$SegID != 50 & rawdata$SegID != 51 & rawdata$SegID != 52, ]
rownames(rawdata) <- NULL

# Remove brain regions lower than 13C resolution limit from rawdata
res_lim = 1.5*1.5*1.5 # cc
segID_greaterthanres = SegID[SegID$CI_upper >= res_lim, ]
segID_lessthanres = SegID[SegID$CI_upper < res_lim, ]
n_lessthanres = nrow(segID_lessthanres)
rownames(segID_greaterthanres) <- NULL
rownames(segID_lessthanres) <- NULL

rawdata_VThresh = rawdata[!(rawdata$SegID %in% segID_lessthanres$SegID), ] # VThresh = Volume Thresholded
rownames(rawdata_VThresh) <- NULL

# Exclude tumor involved brain regions from patients using the gtv_invol csv file
rawdata_VThresh_copy <- rawdata_VThresh
gtv_invol$Remove=TRUE 
temp1=left_join(rawdata_VThresh_copy, gtv_invol)
temp1$Remove[is.na(temp1$Remove)]<-FALSE 
rawdata_VThresh_notumor=temp1[temp1$Remove==FALSE,] # tumor involved brain regions removed for specific patients
rownames(rawdata_VThresh_notumor) <- NULL

# Set categorical variable coding scheme 
rawdata_VThresh_notumor$ParticipantType <- factor(rawdata_VThresh_notumor$ParticipantType) # Control is reference level, dummy coded 
rawdata_VThresh_notumor$ID <- factor(rawdata_VThresh_notumor$ID) # dummy coded 
rawdata_VThresh_notumor$Sex <- factor(rawdata_VThresh_notumor$Sex) # dummy coded 
rawdata_VThresh_notumor$Sex <- relevel(rawdata_VThresh_notumor$Sex, ref = "M") # set male to reference level
rawdata_VThresh_notumor$SegID <- factor(rawdata_VThresh_notumor$SegID) 
contrasts(rawdata_VThresh_notumor$SegID) <- contr.sum(132-40-6) # deviation code 86 regions (132 minus 40 regions < res lim - 6 ventricles).

# Filter data by grey and white matter regions, tumor regions removed
# Grey matter regions 
segID_GM = SegID_class[SegID_class$WMGMid == 'GM' | SegID_class$WMGMid == 'Both', ]
rownames(segID_GM) <- NULL

# White matter regions
segID_WM = SegID_class[SegID_class$WMGMid == 'WM' | SegID_class$WMGMid == 'Both', ]
rownames(segID_WM) <- NULL

rawdata_VThresh_notumor_GM = rawdata_VThresh_notumor[(rawdata_VThresh_notumor$SegID %in% segID_GM$SegID), ] 
rownames(rawdata_VThresh_notumor_GM) <- NULL
rawdata_VThresh_notumor_WM = rawdata_VThresh_notumor[(rawdata_VThresh_notumor$SegID %in% segID_WM$SegID), ] 
rownames(rawdata_VThresh_notumor_WM) <- NULL

# Set categorical variable coding scheme 
rawdata_VThresh_notumor_GM$ParticipantType <- factor(rawdata_VThresh_notumor_GM$ParticipantType) # Control is reference level, dummy coded 
rawdata_VThresh_notumor_GM$ID <- factor(rawdata_VThresh_notumor_GM$ID) # dummy coded 
rawdata_VThresh_notumor_GM$Sex <- factor(rawdata_VThresh_notumor_GM$Sex) # dummy coded
rawdata_VThresh_notumor_GM$Sex <- relevel(rawdata_VThresh_notumor_GM$Sex, ref = "M") # male is reference level
rawdata_VThresh_notumor_GM$SegID <- factor(rawdata_VThresh_notumor_GM$SegID) 
contrasts(rawdata_VThresh_notumor_GM$SegID) <- contr.sum(82) # deviation code 82 GM regions above res lim

# Dumby code categorical data using factor() function for linear modeling below
rawdata_VThresh_notumor_WM$ParticipantType <- factor(rawdata_VThresh_notumor_WM$ParticipantType) # Control is reference level, dummy coded
rawdata_VThresh_notumor_WM$ID <- factor(rawdata_VThresh_notumor_WM$ID) # dummy coded
rawdata_VThresh_notumor_WM$Sex <- factor(rawdata_VThresh_notumor_WM$Sex) # dummy coded
rawdata_VThresh_notumor_WM$Sex <- relevel(rawdata_VThresh_notumor_WM$Sex, ref = "M") # set male to reference level
rawdata_VThresh_notumor_WM$SegID <- factor(rawdata_VThresh_notumor_WM$SegID) 
contrasts(rawdata_VThresh_notumor_WM$SegID) <- contr.sum(5) # deviation code 12 WM regions above res lim





# ========================
# Raw Lactate Signal 
# ========================

# Histogram of raw lactate signal
ggplot(rawdata_VThresh_notumor, aes(x=Lactate)) + geom_histogram()

# Log transform and plot histogram
rawdata_VThresh_notumor_log <- rawdata_VThresh_notumor
rawdata_VThresh_notumor_log$Lactate <- log10(rawdata_VThresh_notumor$Lactate)
ggplot(rawdata_VThresh_notumor_log, aes(x=Lactate)) + geom_histogram()

# Generate regression formulas and perform model selection via AIC minimization
indep_var= c("ParticipantType","SegID","Sex","Age")
dep_var = "log10(Lactate) ~"
regression_formulas = generate_regression_formulas(indep_var, dep_var)

model_sel_AIC_lac = data.frame(matrix(ncol = 2, nrow = 0))
colnames(model_sel_AIC_lac) <- c("Model", "AIC")
for (f in 1:length(regression_formulas)){
  print(f)
  melr_model_selection = lme(regression_formulas[[f]],
                             random=~1|ID,
                             data = rawdata_VThresh_notumor) # recall contrast coding scheme set above
  summary_melr_model_selection=summary(melr_model_selection)
  model_sel_AIC_lac = rbind(model_sel_AIC_lac, list(Model=paste(deparse(regression_formulas[[f]]), collapse=" "), 
                                                  AIC=summary_melr_model_selection[["AIC"]]))
}
model_sel_AIC_lac <- model_sel_AIC_lac[order(model_sel_AIC_lac$AIC), ]
write.csv(model_sel_AIC_lac,"~/data/NCfiles/BrainMetabolism/ImagingNeuroscience/Regression_Results/lactate/lactate_ModelSelection_AIC.csv", row.names=FALSE)

# Run selected mixed effects model post AIC minimization (view exported csv file above first)
melr_model_lac = lme(log10(Lactate) ~ ParticipantType + SegID + ParticipantType:SegID,
                 random=~1|ID, 
                 data=rawdata_VThresh_notumor) # note the contrast coding set above
summary(melr_model_lac)

# Save regression results 
model_results_lac = data.frame(coef(summary(melr_model_lac)))
model_results_lac <- tibble::rownames_to_column(model_results_lac, "Variables") 
write.table(model_results_lac , file = "~/data/NCfiles/BrainMetabolism/ImagingNeuroscience/Regression_Results/lactate/lactate_MELR_Summary.csv", row.names=FALSE)

# Calculate estimated marginal means
emm_lac = emmeans(melr_model_lac, ~ParticipantType|SegID, type='response') # disregard warning message about contrast dropping - this is normal
summary(emm_lac)
emmdf_lac = as.data.frame(emm_lac)
write.table(emmdf_lac, file = "~/data/NCfiles/BrainMetabolism/ImagingNeuroscience/Regression_Results/lactate/lactate_MELR_EMM.csv", row.names=FALSE)

EMM_lac = data.frame(contrast(emm_lac, method='revpairwise', by='SegID', adjust='none', infer = c(TRUE, TRUE)))
EMM_lac$Ratio_Percent = (EMM_lac$ratio - 1)*100 # convert ratio of means to percent difference
write.table(EMM_lac, file = "~/data/NCfiles/BrainMetabolism/ImagingNeuroscience/Regression_Results/lactate/lactate_MELR_EMM_t-test.csv", row.names=FALSE)

emm_lac_average = emmeans(melr_model_lac, ~ParticipantType,type='response') # disregard warning message about contrast dropping - this is normal
emm_lac_average_df=as.data.frame(emm_lac_average)

emm_lac_segid_class = left_join(EMM_lac, SegID_class) 
emm_lac_segid_class = emm_lac_segid_class[emm_lac_segid_class$Lobes != 'CerebralWhiteMatter',]
emm_lac_segid_class$SegID = fct_reorder(emm_lac_segid_class$SegID, emm_lac_segid_class$Ratio_Percent) # reorder brain regions according to percent change in y-variable

# Run ANOVA and Tukey Post Hoc test to assess differences between brain lobes
# p-values and asterisks for boxplots
anova_lac <- aov(Ratio_Percent ~ Lobes, data = emm_lac_segid_class) # anova between lobes
tukey_lac <- TukeyHSD(anova_lac) # tukey HSD 
tukey_lac_df <- as.data.frame(tukey_lac$Lobes)
tukey_lac_df$comparison <- rownames(tukey_lac_df)
tukey_lac_df$signif <- sapply(tukey_lac_df$`p adj`, pval_to_asterisk)
tukey_lac_df
write.table(tukey_lac_df, file = "~/data/NCfiles/BrainMetabolism/ImagingNeuroscience/Regression_Results/lactate/lactate_lobes_tukey.csv", row.names=FALSE)





# ========================
# Raw Pyruvate Signal 
# ========================

# Histogram of raw pyruvate signal
ggplot(rawdata_VThresh_notumor, aes(x=Pyruvate)) + geom_histogram()

# Log transform and plot histogram
rawdata_VThresh_notumor_log <- rawdata_VThresh_notumor
rawdata_VThresh_notumor_log$Pyruvate <- log10(rawdata_VThresh_notumor$Pyruvate)
ggplot(rawdata_VThresh_notumor_log, aes(x=Pyruvate)) + geom_histogram()

# Generate regression formulas and perform model selection via AIC minimization
indep_var = c("ParticipantType","SegID","Sex","Age")
dep_var = "log10(Pyruvate) ~"
regression_formulas = generate_regression_formulas(indep_var, dep_var)

model_sel_AIC_pyr = data.frame(matrix(ncol = 2, nrow = 0))
colnames(model_sel_AIC_pyr) <- c("Model", "AIC")
for (f in 1:length(regression_formulas)){
  print(f)
  melr_model_selection = lme(regression_formulas[[f]],
                             random=~1|ID,
                             data = rawdata_VThresh_notumor) # recall contrast coding scheme set above
  summary_melr_model_selection=summary(melr_model_selection)
  model_sel_AIC_pyr = rbind(model_sel_AIC_pyr, list(Model=paste(deparse(regression_formulas[[f]]), collapse=" "), 
                                                    AIC=summary_melr_model_selection[["AIC"]]))
}
model_sel_AIC_pyr <- model_sel_AIC_pyr[order(model_sel_AIC_pyr$AIC), ]
write.csv(model_sel_AIC_pyr,"~/data/NCfiles/BrainMetabolism/ImagingNeuroscience/Regression_Results/pyruvate/pyruvate_ModelSelection_AIC.csv", row.names=FALSE)


# Run selected mixed effects model post AIC minimization (view exported csv file above first)
melr_model_pyr = lme(log10(Pyruvate) ~ ParticipantType + SegID + ParticipantType:SegID,
                     random=~1|ID, 
                     data=rawdata_VThresh_notumor) # note the contrast coding set above
summary(melr_model_pyr)

# Save regression results 
model_results_pyr = data.frame(coef(summary(melr_model_pyr)))
model_results_pyr <- tibble::rownames_to_column(model_results_pyr, "Variables") 
write.table(model_results_pyr , file = "~/data/NCfiles/BrainMetabolism/ImagingNeuroscience/Regression_Results/pyruvate/pyruvate_MELR_Summary.csv", row.names=FALSE)

# Calculate estimated marginal means
emm_pyr = emmeans(melr_model_pyr, ~ParticipantType|SegID, type='response') # disregard warning message about contrast dropping - this is normal
summary(emm_pyr)
emmdf_pyr = as.data.frame(emm_pyr)
write.table(emmdf_pyr, file = "~/data/NCfiles/BrainMetabolism/ImagingNeuroscience/Regression_Results/pyruvate/pyruvate_MELR_EMM.csv", row.names=FALSE)

EMM_pyr = data.frame(contrast(emm_pyr, method='revpairwise', by='SegID', adjust='none', infer = c(TRUE, TRUE)))
EMM_pyr$Ratio_Percent = (EMM_pyr$ratio - 1)*100 # convert ratio of means to percent difference
write.table(EMM_pyr, file = "~/data/NCfiles/BrainMetabolism/ImagingNeuroscience/Regression_Results/pyruvate/pyruvate_MELR_EMM_t-test.csv", row.names=FALSE)

emm_pyr_average = emmeans(melr_model_pyr, ~ParticipantType,type='response') # disregard warning message about contrast dropping - this is normal
emm_pyr_average_df=as.data.frame(emm_pyr_average)

emm_pyr_segid_class = left_join(EMM_pyr, SegID_class) 
emm_pyr_segid_class = emm_pyr_segid_class[emm_pyr_segid_class$Lobes != 'CerebralWhiteMatter',]
emm_pyr_segid_class$SegID = fct_reorder(emm_pyr_segid_class$SegID, emm_pyr_segid_class$Ratio_Percent) # reorder brain regions according to percent change in y-variable

# Run ANOVA and Tukey Post Hoc test to assess differences between brain lobes
# p-values and asterisks for boxplots
anova_pyr <- aov(Ratio_Percent ~ Lobes, data = emm_pyr_segid_class) # anova between lobes
tukey_pyr <- TukeyHSD(anova_pyr) # tukey HSD 
tukey_pyr_df <- as.data.frame(tukey_pyr$Lobes)
tukey_pyr_df$comparison <- rownames(tukey_pyr_df)
tukey_pyr_df$signif <- sapply(tukey_pyr_df$`p adj`, pval_to_asterisk)
tukey_pyr_df
write.table(tukey_pyr_df, file = "~/data/NCfiles/BrainMetabolism/ImagingNeuroscience/Regression_Results/pyruvate/pyruvate_lobes_tukey.csv", row.names=FALSE)





# ========================
# Raw Bicarbonate Signal 
# ========================

# Histogram of raw bicarbonate signal
ggplot(rawdata_VThresh_notumor, aes(x=Bicarbonate)) + geom_histogram()

# Log transform and plot histogram
rawdata_VThresh_notumor_log <- rawdata_VThresh_notumor
rawdata_VThresh_notumor_log$Bicarbonate <- log10(rawdata_VThresh_notumor$Bicarbonate)
ggplot(rawdata_VThresh_notumor_log, aes(x=Bicarbonate)) + geom_histogram()

# Generate regression formulas and perform model selection via AIC minimization
indep_var = c("ParticipantType","SegID","Sex","Age")
dep_var = "log10(Bicarbonate) ~"
regression_formulas = generate_regression_formulas(indep_var, dep_var)

model_sel_AIC_bic = data.frame(matrix(ncol = 2, nrow = 0))
colnames(model_sel_AIC_bic) <- c("Model", "AIC")
for (f in 1:length(regression_formulas)){
  print(f)
  melr_model_selection = lme(regression_formulas[[f]],
                             random=~1|ID,
                             data = rawdata_VThresh_notumor) # recall contrast coding scheme set above
  summary_melr_model_selection=summary(melr_model_selection)
  model_sel_AIC_bic = rbind(model_sel_AIC_bic, list(Model=paste(deparse(regression_formulas[[f]]), collapse=" "), 
                                                    AIC=summary_melr_model_selection[["AIC"]]))
}
model_sel_AIC_bic <- model_sel_AIC_bic[order(model_sel_AIC_bic$AIC), ]
write.csv(model_sel_AIC_bic,"~/data/NCfiles/BrainMetabolism/ImagingNeuroscience/Regression_Results/bicarbonate/bicarbonate_ModelSelection_AIC.csv", row.names=FALSE)

# Run selected mixed effects model post AIC minimization (view exported csv file above first)
melr_model_bic = lme(log10(Bicarbonate) ~ ParticipantType + SegID + ParticipantType:SegID,
                     random=~1|ID, 
                     data=rawdata_VThresh_notumor) # note the contrast coding set above
summary(melr_model_bic)

# Save regression results 
model_results_bic = data.frame(coef(summary(melr_model_bic)))
model_results_bic <- tibble::rownames_to_column(model_results_bic, "Variables") 
write.table(model_results_bic , file = "~/data/NCfiles/BrainMetabolism/ImagingNeuroscience/Regression_Results/bicarbonate/bicarbonate_MELR_Summary.csv", row.names=FALSE)

# Calculate estimated marginal means
emm_bic = emmeans(melr_model_bic, ~ParticipantType|SegID, type='response') # disregard warning message about contrast dropping - this is normal
summary(emm_bic)
emmdf_bic = as.data.frame(emm_bic)
write.table(emmdf_bic, file = "~/data/NCfiles/BrainMetabolism/ImagingNeuroscience/Regression_Results/bicarbonate/bicarbonate_MELR_EMM.csv", row.names=FALSE)

EMM_bic = data.frame(contrast(emm_bic, method='revpairwise', by='SegID', adjust='none', infer = c(TRUE, TRUE)))
EMM_bic$Ratio_Percent = (EMM_bic$ratio - 1)*100 # convert ratio of means to percent difference
write.table(EMM_bic, file = "~/data/NCfiles/BrainMetabolism/ImagingNeuroscience/Regression_Results/bicarbonate/bicarbonate_MELR_EMM_t-test.csv", row.names=FALSE)

emm_bic_average = emmeans(melr_model_bic, ~ParticipantType,type='response') # disregard warning message about contrast dropping - this is normal
emm_bic_average_df = as.data.frame(emm_bic_average)

emm_bic_segid_class = left_join(EMM_bic, SegID_class) 
emm_bic_segid_class = emm_bic_segid_class[emm_bic_segid_class$Lobes != 'CerebralWhiteMatter',]
emm_bic_segid_class$SegID = fct_reorder(emm_bic_segid_class$SegID, emm_bic_segid_class$Ratio_Percent) # reorder brain regions according to percent change in y-variable

# Run ANOVA and Tukey Post Hoc test to assess differences between brain lobes
# p-values and asterisks for boxplots
anova_bic <- aov(Ratio_Percent ~ Lobes, data = emm_bic_segid_class) # anova between lobes
tukey_bic <- TukeyHSD(anova_bic) # tukey HSD 
tukey_bic_df <- as.data.frame(tukey_bic$Lobes)
tukey_bic_df$comparison <- rownames(tukey_bic_df)
tukey_bic_df$signif <- sapply(tukey_bic_df$`p adj`, pval_to_asterisk)
tukey_bic_df
write.table(tukey_bic_df, file = "~/data/NCfiles/BrainMetabolism/ImagingNeuroscience/Regression_Results/bicarbonate/bicarbonate_lobes_tukey.csv", row.names=FALSE)





# ====================================
# Lactate-to-Pyruvate Ratio (lac/pyr)
# ====================================

# Histogram of lac/pyr ratio
ggplot(rawdata_VThresh_notumor, aes(x=LacPyr)) + geom_histogram()

# Log transform and plot histogram
rawdata_VThresh_notumor_log <- rawdata_VThresh_notumor
rawdata_VThresh_notumor_log$LacPyr <- log10(rawdata_VThresh_notumor$LacPyr)
ggplot(rawdata_VThresh_notumor_log, aes(x=LacPyr)) + geom_histogram()

# Generate regression formulas and perform model selection via AIC minimization
indep_var = c("ParticipantType","SegID","Sex","Age")
dep_var = "log10(LacPyr) ~"
regression_formulas = generate_regression_formulas(indep_var, dep_var)

model_sel_AIC_lacpyr = data.frame(matrix(ncol = 2, nrow = 0))
colnames(model_sel_AIC_lacpyr) <- c("Model", "AIC")
for (f in 1:length(regression_formulas)){
  print(f)
  melr_model_selection = lme(regression_formulas[[f]],
                             random=~1|ID,
                             data = rawdata_VThresh_notumor) # recall contrast coding scheme set above
  summary_melr_model_selection=summary(melr_model_selection)
  model_sel_AIC_lacpyr = rbind(model_sel_AIC_lacpyr, list(Model=paste(deparse(regression_formulas[[f]]), collapse=" "), 
                                                    AIC=summary_melr_model_selection[["AIC"]]))
}
model_sel_AIC_lacpyr <- model_sel_AIC_lacpyr[order(model_sel_AIC_lacpyr$AIC), ]
write.csv(model_sel_AIC_lacpyr,"~/data/NCfiles/BrainMetabolism/ImagingNeuroscience/Regression_Results/lacpyr/lacpyr_ModelSelection_AIC.csv", row.names=FALSE)

# Run selected mixed effects model post AIC minimization (view exported csv file above first)
melr_model_lacpyr = lme(log10(LacPyr) ~ ParticipantType + SegID + ParticipantType:SegID,
                 random=~1|ID, 
                 data=rawdata_VThresh_notumor) # note the contrast coding set above
summary(melr_model_lacpyr)

# Save regression results
model_results_lacpyr = data.frame(coef(summary(melr_model_lacpyr)))
model_results_lacpyr <- tibble::rownames_to_column(model_results_lacpyr, "Variables") 
write.table(model_results_lacpyr, file = "~/data/NCfiles/BrainMetabolism/ImagingNeuroscience/Regression_Results/lacpyr/lacpyr_MELR_Summary.csv")

# Estimated marginal means
emm_lacpyr = emmeans(melr_model_lacpyr, ~ParticipantType|SegID,type='response') # disregard warning message about contrast dropping - this is normal
emmdf_lacpyr = as.data.frame(emm_lacpyr)
write.table(emmdf_lacpyr, file = "~/data/NCfiles/BrainMetabolism/ImagingNeuroscience/Regression_Results/lacpyr/lacpyr_MELR_EMM.csv")

EMM_lacpyr = data.frame(contrast(emm_lacpyr,method='revpairwise', by='SegID', adjust='none', infer = c(TRUE, TRUE)))
EMM_lacpyr$Ratio_Percent = (EMM_lacpyr$ratio - 1)*100 # convert ratio of means to percent difference
write.table(EMM_lacpyr, file = "~/data/NCfiles/BrainMetabolism/ImagingNeuroscience/Regression_Results/lacpyr/lacpyr_MELR_EMM_t-test.csv")

emm_lacpyr_average = emmeans(melr_model_lacpyr, ~ParticipantType,type='response') # disregard warning message about contrast dropping - this is normal
emm_lacpyr_average_df=as.data.frame(emm_lacpyr_average)

emm_lacpyr_segid_class = left_join(EMM_lacpyr, SegID_class) 
emm_lacpyr_segid_class = emm_lacpyr_segid_class[emm_lacpyr_segid_class$Lobes != 'CerebralWhiteMatter',]
emm_lacpyr_segid_class$SegID = fct_reorder(emm_lacpyr_segid_class$SegID, emm_lacpyr_segid_class$Ratio_Percent) # reorder brain regions according to percent change in y-variable

# Run ANOVA and Tukey Post Hoc test to assess differences between brain lobes
# p-values and asterisks for boxplots
anova_lacpyr <- aov(Ratio_Percent ~ Lobes, data = emm_lacpyr_segid_class) # anova between lobes
tukey_lacpyr <- TukeyHSD(anova_lacpyr)
tukey_lacpyr_df <- as.data.frame(tukey_lacpyr$Lobes)
tukey_lacpyr_df$comparison <- rownames(tukey_lacpyr_df)
tukey_lacpyr_df$signif <- sapply(tukey_lacpyr_df$`p adj`, pval_to_asterisk)
tukey_lacpyr_df
write.table(tukey_lacpyr_df, file = "~/data/NCfiles/BrainMetabolism/ImagingNeuroscience/Regression_Results/lacpyr/lacpyr_lobes_tukey.csv", row.names=FALSE)
tukey_lacpyr_df= tukey_lacpyr_df[tukey_lacpyr_df$'p adj'<=0.05,]
comparisons <- strsplit(tukey_lacpyr_df$comparison, "-") # Create list of comparisons for geom_signif
y_positions <- seq(from = max(emm_lacpyr_segid_class$Ratio_Percent) * 1.15, 
                   by = max(emm_lacpyr_segid_class$Ratio_Percent) * 0.15, length.out = nrow(tukey_lacpyr_df))

# Boxplot of lac/pyr ratio vs brain lobe
lacpyr_VThresh_EMM_lobe <- ggplot(emm_lacpyr_segid_class, aes(x=reorder(Lobes, -Ratio_Percent, FUN=median), y=Ratio_Percent))+
  geom_boxplot(outlier.shape=NA, color='black')+
  geom_signif(comparisons = comparisons, annotations = tukey_lacpyr_df$signif,
              y_position=y_positions, tip_length = 0, vjust=0.5,textsize = 4) +
  theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank(),
        panel.background=element_blank(), axis.line=element_line(color='black'),
        axis.text.x = element_text(angle=45, hjust=1.1, color="black"),
        axis.text.y = element_text(color="black"),
        text=element_text(size=20), legend.position="none")+
  scale_y_continuous(breaks=seq(0,90,10.0))+
  labs(x='', y=expression('EMM' ~ Delta * '% lac/pyr'),fill='')
lacpyr_VThresh_EMM_lobe

# Save figure as .eps and .tiff
ggsave(lacpyr_VThresh_EMM_lobe, file="Figure_3a.eps", width = 9, height = 12, units = "cm", device="eps", 
       path="/home/ncappelletto/data/NCfiles/BrainMetabolism/ImagingNeuroscience/Figure_3/")
ggsave(lacpyr_VThresh_EMM_lobe, file="Figure_3a.tiff", width = 9, height = 12, units = "cm", device="tiff", 
       dpi=300, path="/home/ncappelletto/data/NCfiles/BrainMetabolism/ImagingNeuroscience/Figure_3/")

# Boxplot of average lac/pyr ratio in patients vs controls
lacpyr_VThresh_EMM_allregions <- ggplot(rawdata_VThresh_notumor, aes(x=ParticipantType, y=LacPyr)) +
  geom_boxplot(data=rawdata_VThresh_notumor, aes(x=ParticipantType, y=LacPyr), 
               color='black', position=position_dodge(width=0.85), outlier.shape=NA)+
  geom_signif(comparisons=list(c('Control','Patient')), annotation = '***',
              y_position = 3.8,tip_length = 0,textsize = 5) +
  geom_point(data=emm_lacpyr_average_df, aes(x=ParticipantType, y=response),
             position=position_dodge(width=0.85), size=2, color='red3')+
  geom_errorbar(data=emm_lacpyr_average_df,aes(x=ParticipantType, y=response, ymin=lower.CL, ymax=upper.CL, width=0.2),
                position=position_dodge(width=0.85), color='red3')+
  theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank(),
        panel.background=element_blank(), axis.line=element_line(color='black'),
        axis.text.x = element_text(color="black"), axis.text.y = element_text(color="black"),
        text=element_text(size=20),legend.position = 'top')+
  scale_y_continuous(breaks=seq(0,5,0.4))+
  coord_cartesian(ylim = c(0,4.2))+
  scale_x_discrete(labels=c('Controls','Patients'))+
  labs(x='', y='lac/pyr',fill='')
lacpyr_VThresh_EMM_allregions 

# save figures as .eps and .tiff
ggsave(lacpyr_VThresh_EMM_allregions, file="Figure_2a.eps", width = 9, height = 12, units = "cm", device="eps", 
       path="/home/ncappelletto/data/NCfiles/BrainMetabolism/ImagingNeuroscience/Figure_2/")
ggsave(lacpyr_VThresh_EMM_allregions, file="Figure_2a.tiff", width = 9, height = 12, units = "cm", device="tiff", 
       dpi=300, path="/home/ncappelletto/data/NCfiles/BrainMetabolism/ImagingNeuroscience/Figure_2/")

# Region with max change: R Cun
Cun = rawdata_VThresh_notumor[rawdata_VThresh_notumor$SegID == 114 | rawdata_VThresh_notumor$SegID == 115,]
Cun_emm <- subset(emmdf_lacpyr, SegID==114 | SegID==115)

# Set p-value asterisks for boxplots based on EMM t-tests
asterisks <- c("***", "***")   # SegID 114 → ***, SegID 115 → ***
y_ast     <- c(3.3, 3.3)
y_bracket <- c(3.2, 3.2)

# Boxplot of region with max change: R Cun
lacpyr_VThresh_EMM_Cun <- ggplot(Cun, aes(x=as.factor(SegID), y=LacPyr, fill=ParticipantType)) +
  geom_boxplot(data=Cun, aes(x=as.factor(SegID), y=LacPyr, fill=ParticipantType), 
               color='black', position=position_dodge(width=0.85), outlier.shape=NA)+
  geom_point(data=Cun_emm, aes(x=as.factor(SegID), y=response, fill=ParticipantType),
             position=position_dodge(width=0.85), size=2, color='red3')+
  geom_errorbar(data=Cun_emm, aes(x=as.factor(SegID), y=response, ymin=lower.CL, ymax=upper.CL, width=0.2,
                                 fill=ParticipantType), position=position_dodge(width=0.85), color='red3')+
  annotate("text", x = 1, y = y_ast[1], label = asterisks[1], size = 5) +
  geom_segment(aes(x = 1 - 0.25, xend = 1 + 0.25,
                   y = y_bracket[1], yend = y_bracket[1])) +
  annotate("text", x = 2, y = y_ast[2], label = asterisks[2], size = 5) +
  geom_segment(aes(x = 2 - 0.25, xend = 2 + 0.25,
                   y = y_bracket[2], yend = y_bracket[2])) +
  theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank(),
        panel.background=element_blank(), axis.line=element_line(color='black'),
        axis.text.x = element_text(color="black"), axis.text.y = element_text(color="black"),
        text=element_text(size=20), legend.position = 'top')+
  scale_fill_manual(values=c('Control'="grey",'Patient'="grey45"), labels=c("Controls","Patients"))+
  scale_y_continuous(breaks=seq(0,5,0.3))+
  coord_cartesian(ylim =  c(0, 3.3))+
  scale_x_discrete(labels=c('R Cun','L Cun'))+
  labs(x='', y='lac/pyr',fill='')
lacpyr_VThresh_EMM_Cun

# Save figure as .eps and .tiff
ggsave(lacpyr_VThresh_EMM_Cun, file="Figure_2b.eps", width = 9, height = 12, units = "cm", device="eps", 
       path="/home/ncappelletto/data/NCfiles/BrainMetabolism/ImagingNeuroscience/Figure_2/")
ggsave(lacpyr_VThresh_EMM_Cun, file="Figure_2b.tiff", width = 9, height = 12, units = "cm", device="tiff", 
       dpi=300, path="/home/ncappelletto/data/NCfiles/BrainMetabolism/ImagingNeuroscience/Figure_2/")

# Region with min change: L OCP
OCP = rawdata_VThresh_notumor[rawdata_VThresh_notumor$SegID == 156 | rawdata_VThresh_notumor$SegID == 157,]
OCP_emm <- subset(emmdf_lacpyr, SegID==156 | SegID==157)

# Set p-value asterisks for boxplots based on EMM t-tests
asterisks <- c("***", "n.s.")  # SegID 156 → **, SegID 157 → n.s.
y_ast     <- c(2.05, 2.1)
y_bracket <- c(2.0, 2.0)

# Boxplot of region with max change: L OCP
lacpyr_VThresh_EMM_OCP <- ggplot(OCP, aes(x=as.factor(SegID), y=LacPyr, fill=ParticipantType)) +
  geom_boxplot(data=OCP, aes(x=as.factor(SegID), y=LacPyr, fill=ParticipantType), 
               color='black', position=position_dodge(width=0.85), outlier.shape=NA)+
  geom_point(data=OCP_emm, aes(x=as.factor(SegID), y=response, fill=ParticipantType),
             position=position_dodge(width=0.85), size=2, color='red3')+
  geom_errorbar(data=OCP_emm, aes(x=as.factor(SegID), y=response,ymin=lower.CL, ymax=upper.CL, width=0.2,
                                 fill=ParticipantType), position=position_dodge(width=0.85), color='red3')+
  annotate("text", x = 1, y = y_ast[1], label = asterisks[1], size = 5) +
  geom_segment(aes(x = 1 - 0.25, xend = 1 + 0.25,
                   y = y_bracket[1], yend = y_bracket[1])) +
  annotate("text", x = 2, y = y_ast[2], label = asterisks[2], size = 5) +
  geom_segment(aes(x = 2 - 0.25, xend = 2 + 0.25,
                   y = y_bracket[2], yend = y_bracket[2])) +
  theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank(),
        panel.background=element_blank(), axis.line=element_line(color='black'),
        axis.text.x = element_text(color="black"), axis.text.y = element_text(color="black"),
        text=element_text(size=20), legend.position = 'top')+
  scale_fill_manual(values=c('Control'="grey",'Patient'="grey45"),labels=c("Controls","Patients"))+
  scale_y_continuous(breaks=seq(0,5,0.2))+
  coord_cartesian(ylim =  c(0, 2.2))+
  scale_x_discrete(labels=c('R OCP','L OCP'))+
  labs(x='', y='lac/pyr',fill='')
lacpyr_VThresh_EMM_OCP

# save figures as .eps and .tiff
ggsave(lacpyr_VThresh_EMM_OCP, file="Figure_2c.eps", width = 9, height = 12, units = "cm", device="eps", 
       path="/home/ncappelletto/data/NCfiles/BrainMetabolism/ImagingNeuroscience/Figure_2/")
ggsave(lacpyr_VThresh_EMM_OCP, file="Figure_2c.tiff", width = 9, height = 12, units = "cm", device="tiff", 
       dpi=300, path="/home/ncappelletto/data/NCfiles/BrainMetabolism/ImagingNeuroscience/Figure_2/")





# =========================================
# Bicarbonate-to-pyruvate ratio (bic/pyr)
# =========================================

# Histogram of bic/pyr ratio
ggplot(rawdata_VThresh_notumor, aes(x=BicPyr)) + geom_histogram()

# Log transform and plot histogram
rawdata_VThresh_notumor_log <- rawdata_VThresh_notumor
rawdata_VThresh_notumor_log$BicPyr <- log10(rawdata_VThresh_notumor$BicPyr)
ggplot(rawdata_VThresh_notumor_log, aes(x=BicPyr)) + geom_histogram()

# Generate regression formulas and perform model selection via AIC minimization
indep_var = c("ParticipantType","SegID","Sex","Age")
dep_var = "log10(BicPyr) ~"
regression_formulas = generate_regression_formulas(indep_var, dep_var)

model_sel_AIC_bicpyr = data.frame(matrix(ncol = 2, nrow = 0))
colnames(model_sel_AIC_bicpyr) <- c("Model", "AIC")
for (f in 1:length(regression_formulas)){
  print(f)
  melr_model_selection = lme(regression_formulas[[f]],
                             random=~1|ID,
                             data = rawdata_VThresh_notumor) # recall contrast coding scheme set above
  summary_melr_model_selection=summary(melr_model_selection)
  model_sel_AIC_bicpyr = rbind(model_sel_AIC_bicpyr, list(Model=paste(deparse(regression_formulas[[f]]), collapse=" "), 
                                                          AIC=summary_melr_model_selection[["AIC"]]))
}
model_sel_AIC_bicpyr <- model_sel_AIC_bicpyr[order(model_sel_AIC_bicpyr$AIC), ]
write.csv(model_sel_AIC_bicpyr,"~/data/NCfiles/BrainMetabolism/ImagingNeuroscience/Regression_Results/bicpyr/bicpyr_ModelSelection_AIC.csv", row.names=FALSE)

# Run selected mixed effects model post AIC minimization (view exported csv file above first)
melr_model_bicpyr = lme(log10(BicPyr) ~ ParticipantType + SegID + ParticipantType:SegID,
                 random=~1|ID, 
                 data=rawdata_VThresh_notumor) # note the contrast coding set above
summary(melr_model_bicpyr)

# Save regression results
model_results_bicpyr = data.frame(coef(summary(melr_model_bicpyr)))
model_results_bicpyr <- tibble::rownames_to_column(model_results_bicpyr,"Variables") 
write.table(model_results_bicpyr, file = "~/data/NCfiles/BrainMetabolism/ImagingNeuroscience/Regression_Results/bicpyr/bicpyr_MELR_Summary.csv")

# Estimated marginal means
emm_bicpyr = emmeans(melr_model_bicpyr, ~ParticipantType|SegID,type='response') # disregard warning message about contrast dropping - this is normal
emmdf_bicpyr = as.data.frame(emm_bicpyr)
write.table(emmdf_bicpyr, file = "~/data/NCfiles/BrainMetabolism/ImagingNeuroscience/Regression_Results/bicpyr/bicpyr_MELR_EMM.csv")

EMM_bicpyr = data.frame(contrast(emm_bicpyr, method='revpairwise', by='SegID', adjust='none',infer = c(TRUE, TRUE)))
EMM_bicpyr$Ratio_Percent = (EMM_bicpyr$ratio - 1)*100 # convert ratio of means to percent difference
write.table(EMM_bicpyr , file = "~/data/NCfiles/BrainMetabolism/ImagingNeuroscience/Regression_Results/bicpyr/bicpyr_MELR_EMM_t-test.csv")

emm_bicpyr_average = emmeans(melr_model_bicpyr, ~ParticipantType,type='response') # disregard warning message about contrast dropping - this is normal
emm_bicpyr_average_df=as.data.frame(emm_bicpyr_average)

emm_bicpyr_segid_class = left_join(EMM_bicpyr, SegID_class) 
emm_bicpyr_segid_class = emm_bicpyr_segid_class[emm_bicpyr_segid_class$Lobes != 'CerebralWhiteMatter',]
emm_bicpyr_segid_class$SegID = fct_reorder(emm_bicpyr_segid_class$SegID, emm_bicpyr_segid_class$Ratio_Percent) # reorder brain regions according to percent change in y-variable

# Run ANOVA and Tukey Post Hoc test to assess differences between brain lobes
# p-values and asterisks for boxplots
anova_bicpyr <- aov(Ratio_Percent ~ Lobes, data = emm_bicpyr_segid_class) # anova between lobes
tukey_bicpyr <- TukeyHSD(anova_bicpyr)
tukey_bicpyr_df <- as.data.frame(tukey_bicpyr$Lobes)
tukey_bicpyr_df$comparison <- rownames(tukey_bicpyr_df)
tukey_bicpyr_df$signif <- sapply(tukey_bicpyr_df$`p adj`, pval_to_asterisk)
tukey_bicpyr_df
write.table(tukey_bicpyr_df, file = "~/data/NCfiles/BrainMetabolism/ImagingNeuroscience/Regression_Results/bicpyr/bicpyr_lobes_tukey.csv", row.names=FALSE)
tukey_bicpyr_df= tukey_bicpyr_df[tukey_bicpyr_df$'p adj'<=0.05,]
comparisons <- strsplit(tukey_bicpyr_df$comparison, "-") # Create list of comparisons for geom_signif
y_positions <- seq(from = max(emm_bicpyr_segid_class$Ratio_Percent) * 1.05, 
                   by = max(emm_bicpyr_segid_class$Ratio_Percent) * 0.06, length.out = nrow(tukey_bicpyr_df))

# Boxplot of bic/pyr ratio vs brain lobe
bicpyr_VThresh_EMM_lobe <- ggplot(emm_bicpyr_segid_class, aes(x=reorder(Lobes, -Ratio_Percent, FUN=median), y=Ratio_Percent))+
  geom_boxplot(outlier.shape=NA, color='black')+
  geom_signif(comparisons = comparisons, annotations = tukey_bicpyr_df$signif,
              y_position=y_positions, tip_length = 0, vjust=0.5, textsize = 4) +
  theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank(),
        panel.background=element_blank(), axis.line=element_line(color='black'),
        axis.text.x = element_text(angle=45, hjust=1.1, color="black"),
        axis.text.y = element_text(color="black"),
        text=element_text(size=20), legend.position="none")+
  scale_y_continuous(breaks=seq(-10,180,10))+
  labs(x='', y=expression('EMM' ~ Delta * '% bic/pyr'),fill='')
bicpyr_VThresh_EMM_lobe 

# Save figure as .eps and .tiff
ggsave(bicpyr_VThresh_EMM_lobe, file="Figure_3b.eps", width = 9, height = 12, units = "cm", device="eps", 
       path="/home/ncappelletto/data/NCfiles/BrainMetabolism/ImagingNeuroscience/Figure_3/")
ggsave(bicpyr_VThresh_EMM_lobe, file="Figure_3b.tiff", width = 9, height = 12, units = "cm", device="tiff", 
       dpi=300, path="/home/ncappelletto/data/NCfiles/BrainMetabolism/ImagingNeuroscience/Figure_3/")

# Boxplot of average lac/pyr ratio in patients vs controls
bicpyr_VThresh_EMM_allregions <- ggplot(rawdata_VThresh_notumor, aes(x=ParticipantType, y=BicPyr)) +
  geom_boxplot(data=rawdata_VThresh_notumor, aes(x=ParticipantType, y=BicPyr), 
               color='black', position=position_dodge(width=0.85), outlier.shape=NA)+
  geom_signif(comparisons=list(c('Control','Patient')), annotation = 'n.s.',
              y_position = 0.85,tip_length = 0, textsize = 5) +
  geom_point(data=emm_bicpyr_average_df, aes(x=ParticipantType, y=response),
             position=position_dodge(width=0.85), size=2, color='red3')+
  geom_errorbar(data=emm_bicpyr_average_df,aes(x=ParticipantType, y=response, ymin=lower.CL, ymax=upper.CL, width=0.2),
                position=position_dodge(width=0.85), color='red3')+
  theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank(),
        panel.background=element_blank(), axis.line=element_line(color='black'),
        axis.text.x = element_text(color="black"), axis.text.y = element_text(color="black"),
        text=element_text(size=20), legend.position = 'top')+
  scale_y_continuous(breaks=seq(0,5,0.1))+
  coord_cartesian(ylim = c(0,1))+
  scale_x_discrete(labels=c('Controls','Patients'))+
  labs(x='', y='bic/pyr',fill='')
bicpyr_VThresh_EMM_allregions

# Save figure as .eps and .tiff
ggsave(bicpyr_VThresh_EMM_allregions, file="Figure_2d.eps", width = 9, height = 12, units = "cm", device="eps", 
       path="/home/ncappelletto/data/NCfiles/BrainMetabolism/ImagingNeuroscience/Figure_2/")
ggsave(bicpyr_VThresh_EMM_allregions, file="Figure_2d.tiff", width = 9, height = 12, units = "cm", device="tiff", 
       dpi=300, path="/home/ncappelletto/data/NCfiles/BrainMetabolism/ImagingNeuroscience/Figure_2/")

# Region with max change: R Cun
Cun = rawdata_VThresh_notumor[rawdata_VThresh_notumor$SegID == 114 | rawdata_VThresh_notumor$SegID == 115,]
Cun_emm <- subset(emmdf_bicpyr, SegID==114 | SegID==115)

# Set p-value asterisks for boxplots based on EMM t-tests
asterisks <- c("**", "*")  # SegID 114 → **, SegID 115 → *
y_ast     <- c(0.82, 0.82)
y_bracket <- c(0.8, 0.8)

# Boxplot of region with max change: R Cun
bicpyr_VThresh_EMM_Cun <- ggplot(Cun, aes(x=as.factor(SegID), y=BicPyr, fill=ParticipantType)) +
  geom_boxplot(data=Cun, aes(x=as.factor(SegID), y=BicPyr, fill=ParticipantType), 
               color='black', position=position_dodge(width=0.85), outlier.shape=NA)+
  geom_point(data=Cun_emm, aes(x=as.factor(SegID), y=response, fill=ParticipantType),
             position=position_dodge(width=0.85), size=2, color='red3')+
  geom_errorbar(data=Cun_emm,aes(x=as.factor(SegID), y=response, ymin=lower.CL, ymax=upper.CL, width=0.2,
                                 fill=ParticipantType), position=position_dodge(width=0.85), color='red3')+
  annotate("text", x = 1, y = y_ast[1], label = asterisks[1], size = 5) +
  geom_segment(aes(x = 1 - 0.25, xend = 1 + 0.25,
                   y = y_bracket[1], yend = y_bracket[1])) +
  annotate("text", x = 2, y = y_ast[2], label = asterisks[2], size = 5) +
  geom_segment(aes(x = 2 - 0.25, xend = 2 + 0.25,
                   y = y_bracket[2], yend = y_bracket[2])) +
  theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank(),
        panel.background=element_blank(), axis.line=element_line(color='black'),
        axis.text.x = element_text(color="black"), axis.text.y = element_text(color="black"),
        text=element_text(size=20), legend.position = 'top')+
  scale_fill_manual(values=c('Control'="grey",'Patient'="grey45"),labels=c("Controls","Patients"))+
  scale_y_continuous(breaks=seq(0,5,0.1))+
  coord_cartesian(ylim =  c(0, 0.9))+
  scale_x_discrete(labels=c('R Cun','L Cun'))+
  labs(x='', y='bic/pyr',fill='')
bicpyr_VThresh_EMM_Cun

# Save figure as .eps and .tiff
ggsave(bicpyr_VThresh_EMM_Cun, file="Figure_2e.eps", width = 9, height = 12, units = "cm", device="eps", 
       path="/home/ncappelletto/data/NCfiles/BrainMetabolism/ImagingNeuroscience/Figure_2/")
ggsave(bicpyr_VThresh_EMM_Cun, file="Figure_2e.tiff", width = 9, height = 12, units = "cm", device="tiff", 
       dpi=300, path="/home/ncappelletto/data/NCfiles/BrainMetabolism/ImagingNeuroscience/Figure_2/")

# Region with min change: R VDC
VDC = rawdata_VThresh_notumor[rawdata_VThresh_notumor$SegID == 61 | rawdata_VThresh_notumor$SegID == 62,]
VDC_emm <- subset(emmdf_bicpyr, SegID==61 | SegID==62)

# Set p-value asterisks for boxplots based on EMM t-tests
asterisks <- c("n.s.", "n.s.")  # SegID 61 → n.s., SegID 62 → n.s.
y_ast     <- c(0.53, 0.53)      # Height for text
y_bracket <- c(0.51, 0.51)      # Height for bracket line

# Boxplot of region with min change: R VDC
bicpyr_VThresh_EMM_VDC <- ggplot(VDC, aes(x=as.factor(SegID), y=BicPyr, fill=ParticipantType)) +
  geom_boxplot(data=VDC, aes(x=as.factor(SegID), y=BicPyr, fill=ParticipantType), 
               color='black', position=position_dodge(width=0.85), outlier.shape=NA)+
  geom_point(data=VDC_emm, aes(x=as.factor(SegID), y=response, fill=ParticipantType),
             position=position_dodge(width=0.85), size=2, color='red3')+
  geom_errorbar(data=VDC_emm, aes(x=as.factor(SegID), y=response, ymin=lower.CL, ymax=upper.CL, width=0.2,
                                 fill=ParticipantType), position=position_dodge(width=0.85), color='red3')+
  annotate("text", x = 1, y = y_ast[1], label = asterisks[1], size = 5) +
  geom_segment(aes(x = 1 - 0.25, xend = 1 + 0.25,
                   y = y_bracket[1], yend = y_bracket[1])) +
  annotate("text", x = 2, y = y_ast[2], label = asterisks[2], size = 5) +
  geom_segment(aes(x = 2 - 0.25, xend = 2 + 0.25,
                   y = y_bracket[2], yend = y_bracket[2])) +
  theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank(),
        panel.background=element_blank(), axis.line=element_line(color='black'),
        axis.text.x = element_text(color="black"), axis.text.y = element_text(color="black"),
        text=element_text(size=20), legend.position = 'top')+
  scale_fill_manual(values=c('Control'="grey",'Patient'="grey45"),labels=c("Controls","Patients"))+
  scale_y_continuous(breaks=seq(0,5,0.05))+
  coord_cartesian(ylim =  c(0, 0.55))+
  scale_x_discrete(labels=c('R VDC','L VDC'))+
  labs(x='', y='bic/pyr',fill='')
bicpyr_VThresh_EMM_VDC

# Save figure as .eps and .tiff
ggsave(bicpyr_VThresh_EMM_VDC, file="Figure_2f.eps", width = 9, height = 12, units = "cm", device="eps", 
       path="/home/ncappelletto/data/NCfiles/BrainMetabolism/ImagingNeuroscience/Figure_2/")
ggsave(bicpyr_VThresh_EMM_VDC, file="Figure_2f.tiff", width = 9, height = 12, units = "cm", device="tiff", 
       dpi=300, path="/home/ncappelletto/data/NCfiles/BrainMetabolism/ImagingNeuroscience/Figure_2/")





# =========================================
# Lactate-to-bicarbonate ratio (lac/bic)
# =========================================

# Histogram of signal ratios
ggplot(rawdata_VThresh_notumor, aes(x=LacBic)) + geom_histogram()

# Log transform and plot histogram
rawdata_VThresh_notumor_log <- rawdata_VThresh_notumor
rawdata_VThresh_notumor_log$LacBic <- log10(rawdata_VThresh_notumor$LacBic)
ggplot(rawdata_VThresh_notumor_log, aes(x=LacBic)) + geom_histogram()

# Generate regression formulas and perform model selection via AIC minimization
indep_var = c("ParticipantType","SegID","Sex","Age")
dep_var = "log10(LacBic) ~"
regression_formulas = generate_regression_formulas(indep_var, dep_var)

model_sel_AIC_lacbic = data.frame(matrix(ncol = 2, nrow = 0))
colnames(model_sel_AIC_lacbic) <- c("Model", "AIC")
for (f in 1:length(regression_formulas)){
  print(f)
  melr_model_selection = lme(regression_formulas[[f]],
                             random=~1|ID,
                             data = rawdata_VThresh_notumor) # recall contrast coding scheme set above
  summary_melr_model_selection=summary(melr_model_selection)
  model_sel_AIC_lacbic = rbind(model_sel_AIC_lacbic, list(Model=paste(deparse(regression_formulas[[f]]), collapse=" "), 
                                                          AIC=summary_melr_model_selection[["AIC"]]))
}
model_sel_AIC_lacbic <- model_sel_AIC_lacbic[order(model_sel_AIC_lacbic$AIC), ]
write.csv(model_sel_AIC_lacbic, "~/data/NCfiles/BrainMetabolism/ImagingNeuroscience/Regression_Results/lacbic/lacbic_ModelSelection_AIC.csv", row.names=FALSE)

# Run selected mixed effects model post AIC minimization (view exported csv file above first)
melr_model_lacbic = lme(log10(LacBic) ~ ParticipantType + SegID + ParticipantType:SegID,
                 random=~1|ID, 
                 data=rawdata_VThresh_notumor) # note the contrast coding set above
summary(melr_model_lacbic)

# Save regression results
model_results_lacbic = data.frame(coef(summary(melr_model_lacbic)))
model_results_lacbic <- tibble::rownames_to_column(model_results_lacbic, "Variables") 
write.table(model_results_lacbic, file = "~/data/NCfiles/BrainMetabolism/ImagingNeuroscience/Regression_Results/lacbic/lacbic_MELR_Summary.csv")

# Estimated marginal means
emm_lacbic = emmeans(melr_model_lacbic, ~ParticipantType|SegID,type='response') # disregard warning message about contrast dropping - this is normal
emmdf_lacbic = as.data.frame(emm_lacbic)
write.table(emmdf_lacbic, file = "~/data/NCfiles/BrainMetabolism/ImagingNeuroscience/Regression_Results/lacbic/lacbic_MELR_EMM.csv")

EMM_lacbic = data.frame(contrast(emm_lacbic, method='revpairwise', by='SegID', adjust='none', infer = c(TRUE, TRUE)))
EMM_lacbic$Ratio_Percent = (EMM_lacbic$ratio - 1)*100 # convert ratio of means to percent difference
write.table(EMM_lacbic, file = "~/data/NCfiles/BrainMetabolism/ImagingNeuroscience/Regression_Results/lacbic/lacbic_MELR_EMM_t-test.csv")

emm_lacbic_average = emmeans(melr_model_lacbic, ~ParticipantType,type='response') # disregard warning message about contrast dropping - this is normal
emm_lacbic_average_df = as.data.frame(emm_lacbic_average)

emm_lacbic_segid_class = left_join(EMM_lacbic, SegID_class) 
emm_lacbic_segid_class = emm_lacbic_segid_class[emm_lacbic_segid_class$Lobes != 'CerebralWhiteMatter',]
emm_lacbic_segid_class$SegID = fct_reorder(emm_lacbic_segid_class$SegID, emm_lacbic_segid_class$Ratio_Percent) # reorder brain regions according to percent change in y-variable

# Run ANOVA and Tukey Post Hoc test to assess differences between brain lobes
# p-values and asterisks for boxplots
anova_lacbic <- aov(Ratio_Percent ~ Lobes, data = emm_lacbic_segid_class) # anova between lobes
tukey_lacbic <- TukeyHSD(anova_lacbic) 
tukey_lacbic_df <- as.data.frame(tukey_lacbic$Lobes)
tukey_lacbic_df$comparison <- rownames(tukey_lacbic_df)
tukey_lacbic_df$signif <- sapply(tukey_lacbic_df$`p adj`, pval_to_asterisk)
tukey_lacbic_df
write.table(tukey_lacbic_df, file = "~/data/NCfiles/BrainMetabolism/ImagingNeuroscience/Regression_Results/lacbic/lacbic_lobes_tukey.csv", row.names=FALSE)
tukey_lacbic_df= tukey_lacbic_df[tukey_lacbic_df$'p adj'<=0.05,]
comparisons <- strsplit(tukey_lacbic_df$comparison, "-") # Create list of comparisons for geom_signif
y_positions <- seq(from = max(emm_lacbic_segid_class$Ratio_Percent) * 1.02, 
                   by = max(emm_lacbic_segid_class$Ratio_Percent) * 0.06, length.out = nrow(tukey_lacbic_df))

# Boxplot of lac/bic ratio vs brain lobe
lacbic_VThresh_EMM_lobe <- ggplot(emm_lacbic_segid_class, aes(x=reorder(Lobes, -Ratio_Percent, FUN=median), y=Ratio_Percent))+
  geom_boxplot(outlier.shape=NA, color='black')+
  geom_signif(comparisons = comparisons, annotations = tukey_lacbic_df$signif,
              y_position=y_positions, tip_length = 0, vjust=0.5, textsize = 4) +
  theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank(),
        panel.background=element_blank(), axis.line=element_line(color='black'),
        axis.text.x = element_text(angle=45, hjust=1.1, color="black"),
        axis.text.y = element_text(color="black"),
        text=element_text(size=20), legend.position="none")+
  scale_y_continuous(breaks=seq(0,80,5))+
  coord_cartesian(ylim = c(0,52))+
  labs(x='', y=expression('EMM' ~ Delta * '% lac/bic'),fill='')
lacbic_VThresh_EMM_lobe

# Save figure as .eps and .tiff
ggsave(lacbic_VThresh_EMM_lobe, file="Figure_3c.eps", width = 9, height = 12, units = "cm", device="eps", 
       path="/home/ncappelletto/data/NCfiles/BrainMetabolism/ImagingNeuroscience/Figure_3/")
ggsave(lacbic_VThresh_EMM_lobe, file="Figure_3c.tiff", width = 9, height = 12, units = "cm", device="tiff", 
       dpi=300, path="/home/ncappelletto/data/NCfiles/BrainMetabolism/ImagingNeuroscience/Figure_3/")

# Boxplot of average lac/pyr ratio in patients vs controls
lacbic_VThresh_EMM_allregions <- ggplot(rawdata_VThresh_notumor, aes(x=ParticipantType, y=LacBic)) +
  geom_boxplot(data=rawdata_VThresh_notumor, aes(x=ParticipantType, y=LacBic), 
               color='black', position=position_dodge(width=0.85), outlier.shape=NA)+
  geom_signif(comparisons=list(c('Control','Patient')), annotation = 'n.s.',
              y_position = 7.2, tip_length = 0, textsize = 5) +
  geom_point(data=emm_lacbic_average_df,aes(x=ParticipantType, y=response),
             position=position_dodge(width=0.85),size=2, color='red3')+
  geom_errorbar(data=emm_lacbic_average_df, aes(x=ParticipantType, y=response, ymin=lower.CL, ymax=upper.CL, width=0.2),
                position=position_dodge(width=0.85), color='red3')+
  theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank(),
        panel.background=element_blank(), axis.line=element_line(color='black'),
        axis.text.x = element_text(color="black"), axis.text.y = element_text(color="black"),
        text=element_text(size=20), legend.position = 'top')+
  scale_y_continuous(breaks=seq(0,10,0.6))+
  coord_cartesian(ylim = c(0,8))+
  scale_x_discrete(labels=c('Controls','Patients'))+
  labs(x='', y='lac/bic',fill='')
lacbic_VThresh_EMM_allregions

# Save figure as .eps and .tiff
ggsave(lacbic_VThresh_EMM_allregions, file="Figure_2g.eps", width = 9, height = 12, units = "cm", device="eps", 
       path="/home/ncappelletto/data/NCfiles/BrainMetabolism/ImagingNeuroscience/Figure_2/")
ggsave(lacbic_VThresh_EMM_allregions, file="Figure_2g.tiff", width = 9, height = 12, units = "cm", device="tiff", 
       dpi=300, path="/home/ncappelletto/data/NCfiles/BrainMetabolism/ImagingNeuroscience/Figure_2/")

# Region with max change: R VDC
VDC = rawdata_VThresh_notumor[rawdata_VThresh_notumor$SegID == 61 | rawdata_VThresh_notumor$SegID == 62,]
VDC_emm <- subset(emmdf_lacbic, SegID==61 | SegID==62)

# Set p-value asterisks for boxplots based on EMM t-tests
asterisks <- c("***", "n.s.")  # SegID 61 → ***, SegID 62 → n.s.
y_ast     <- c(11.2, 11.4)
y_bracket <- c(11, 11)

# Boxplot of region with max change: R VDC
lacbic_VThresh_EMM_VDC <- ggplot(VDC, aes(x=as.factor(SegID), y=LacBic, fill=ParticipantType)) +
  geom_boxplot(data=VDC, aes(x=as.factor(SegID), y=LacBic, fill=ParticipantType), 
               color='black', position=position_dodge(width=0.85), outlier.shape=NA)+
  geom_point(data=VDC_emm, aes(x=as.factor(SegID), y=response, fill=ParticipantType),
             position=position_dodge(width=0.85), size=2, color='red3')+
  geom_errorbar(data=VDC_emm, aes(x=as.factor(SegID), y=response, ymin=lower.CL, ymax=upper.CL, width=0.2,
                                 fill=ParticipantType), position=position_dodge(width=0.85), color='red3')+
  annotate("text", x = 1, y = y_ast[1], label = asterisks[1], size = 5) +
  geom_segment(aes(x = 1 - 0.25, xend = 1 + 0.25,
                   y = y_bracket[1], yend = y_bracket[1])) +
  annotate("text", x = 2, y = y_ast[2], label = asterisks[2], size = 5) +
  geom_segment(aes(x = 2 - 0.25, xend = 2 + 0.25,
                   y = y_bracket[2], yend = y_bracket[2])) +
  theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank(),
        panel.background=element_blank(), axis.line=element_line(color='black'),
        axis.text.x = element_text(color="black"), axis.text.y = element_text(color="black"),
        text=element_text(size=20), legend.position = 'top')+
  scale_fill_manual(values=c('Control'="grey",'Patient'="grey45"),labels=c("Controls","Patients"))+
  scale_y_continuous(breaks=seq(0,20,1.5))+
  coord_cartesian(ylim =  c(0, 12))+
  scale_x_discrete(labels=c('R VDC','L VDC'))+
  labs(x='', y='lac/bic',fill='')
lacbic_VThresh_EMM_VDC

# Save figure as .eps and .tiff
ggsave(lacbic_VThresh_EMM_VDC, file="Figure_2h.eps", width = 9, height = 12, units = "cm", device="eps", 
       path="/home/ncappelletto/data/NCfiles/BrainMetabolism/ImagingNeuroscience/Figure_2/")
ggsave(lacbic_VThresh_EMM_VDC, file="Figure_2h.tiff", width = 9, height = 12, units = "cm", device="tiff", 
       dpi=300, path="/home/ncappelletto/data/NCfiles/BrainMetabolism/ImagingNeuroscience/Figure_2/")

# Region with min change:R AnG
AnG = rawdata_VThresh_notumor[rawdata_VThresh_notumor$SegID == 106 | rawdata_VThresh_notumor$SegID == 107,]
AnG_emm <- subset(emmdf_lacbic, SegID==106 | SegID==107)

# Set p-value asterisks for boxplots based on EMM t-tests
asterisks <- c("n.s.", "n.s.")  # SegID 106 → n.s., SegID 107 → n.s.
y_ast     <- c(5.8, 5.8)
y_bracket <- c(5.6, 5.6)

# Boxplot of region with min change: R AnG
lacbic_VThresh_EMM_AnG <- ggplot(AnG, aes(x=as.factor(SegID), y=LacBic, fill=ParticipantType)) +
  geom_boxplot(data=AnG, aes(x=as.factor(SegID), y=LacBic, fill=ParticipantType), 
               color='black', position=position_dodge(width=0.85), outlier.shape=NA)+
  geom_point(data=AnG_emm, aes(x=as.factor(SegID), y=response, fill=ParticipantType),
             position=position_dodge(width=0.85),size=2, color='red3')+
  geom_errorbar(data=AnG_emm, aes(x=as.factor(SegID), y=response, ymin=lower.CL, ymax=upper.CL, width=0.2,
                                 fill=ParticipantType), position=position_dodge(width=0.85), color='red3')+
  annotate("text", x = 1, y = y_ast[1], label = asterisks[1], size = 5) +
  geom_segment(aes(x = 1 - 0.25, xend = 1 + 0.25,
                   y = y_bracket[1], yend = y_bracket[1])) +
  annotate("text", x = 2, y = y_ast[2], label = asterisks[2], size = 5) +
  geom_segment(aes(x = 2 - 0.25, xend = 2 + 0.25,
                   y = y_bracket[2], yend = y_bracket[2])) +
  theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank(),
        panel.background=element_blank(), axis.line=element_line(color='black'),
        axis.text.x = element_text(color="black"), axis.text.y = element_text(color="black"),
        text=element_text(size=20), legend.position = 'top')+
  scale_fill_manual(values=c('Control'="grey",'Patient'="grey45"),labels=c("Controls","Patients"))+
  scale_y_continuous(breaks=seq(0,12,0.4))+
  coord_cartesian(ylim =  c(1.6, 6))+
  scale_x_discrete(labels=c('R AnG','L AnG'))+
  labs(x='', y='lac/bic',fill='')
lacbic_VThresh_EMM_AnG

# Save figure as .eps and .tiff
ggsave(lacbic_VThresh_EMM_AnG, file="Figure_2i.eps", width = 9, height = 12, units = "cm", device="eps", 
       path="/home/ncappelletto/data/NCfiles/BrainMetabolism/ImagingNeuroscience/Figure_2/")
ggsave(lacbic_VThresh_EMM_AnG, file="Figure_2i.tiff", width = 9, height = 12, units = "cm", device="tiff", 
       dpi=300, path="/home/ncappelletto/data/NCfiles/BrainMetabolism/ImagingNeuroscience/Figure_2/")





# ==============================================================
# Grey Matter Density (GMD) / Partial Volume Estimate 1 (PVE1)
# ==============================================================

# Histogram of GMD
ggplot(rawdata_VThresh_notumor_GM, aes(x=PVE1)) + geom_histogram()

# Log transform and plot histogram
rawdata_VThresh_notumor_GM_log <- rawdata_VThresh_notumor_GM
rawdata_VThresh_notumor_GM_log$PVE1 <- log10(rawdata_VThresh_notumor_GM$PVE1)
ggplot(rawdata_VThresh_notumor_GM_log, aes(x=PVE1)) + geom_histogram()

# Generate regression formulas and perform model selection via AIC minimization
indep_var = c("ParticipantType","SegID","Sex","Age")
dep_var = "PVE1 ~"
regression_formulas = generate_regression_formulas(indep_var, dep_var)

model_sel_AIC_GMD = data.frame(matrix(ncol = 2, nrow = 0))
colnames(model_sel_AIC_GMD) <- c("Model", "AIC")
for (f in 1:length(regression_formulas)){
  print(f)
  melr_model_selection = lme(regression_formulas[[f]],
                             random=~1|ID,
                             data = rawdata_VThresh_notumor_GM) # recall contrast coding scheme set above
  summary_melr_model_selection=summary(melr_model_selection)
  model_sel_AIC_GMD = rbind(model_sel_AIC_GMD, list(Model=paste(deparse(regression_formulas[[f]]), collapse=" "), 
                                                    AIC=summary_melr_model_selection[["AIC"]]))
}
model_sel_AIC_GMD <- model_sel_AIC_GMD[order(model_sel_AIC_GMD$AIC), ]
write.csv(model_sel_AIC_GMD,"~/data/NCfiles/BrainMetabolism/ImagingNeuroscience/Regression_Results/GMD/GMD_ModelSelection_AIC.csv", row.names=FALSE)

# Run selected mixed effects model post AIC minimization (view exported csv file above first)
melr_model_GMD = lme(PVE1 ~ ParticipantType + SegID + Age + ParticipantType:SegID,
                 random=~1|ID, 
                 data=rawdata_VThresh_notumor_GM) # note the contrast coding set above
summary(melr_model_GMD)

# Save regression results 
model_results_GMD = data.frame(coef(summary(melr_model_GMD)))
model_results_GMD <- tibble::rownames_to_column(model_results_GMD,"Variables") 
write.table(model_results_GMD, file = "~/data/NCfiles/BrainMetabolism/ImagingNeuroscience/Regression_Results/GMD/GMD_MELR_Summary.csv", row.names=FALSE)

# Calculate estimated marginal means
emm_GMD = emmeans(melr_model_GMD, ~ParticipantType|SegID, type='response') # disregard warning message about contrast dropping - this is normal
emmdf_GMD = as.data.frame(emm_GMD)
write.table(emmdf_GMD, file = "~/data/NCfiles/BrainMetabolism/ImagingNeuroscience/Regression_Results/GMD/GMD_MELR_EMM.csv", row.names=FALSE)

EMM_GMD = data.frame(contrast(emm_GMD, method='revpairwise',by='SegID',adjust='none',infer = c(TRUE, TRUE)))
write.table(EMM_GMD, file = "~/data/NCfiles/BrainMetabolism/ImagingNeuroscience/Regression_Results/GMD/GMD_MELR_EMM_t-test.csv")

emm_GMD_average = emmeans(melr_model_GMD, ~ParticipantType,type='response') # disregard warning message about contrast dropping - this is normal
emm_GMD_average_df = as.data.frame(emm_GMD_average)

emm_GMD_segid_class = left_join(EMM_GMD, SegID_class) 
emm_GMD_segid_class = emm_GMD_segid_class[emm_GMD_segid_class$Lobes != 'CerebralWhiteMatter',]
emm_GMD_segid_class$SegID = fct_reorder(emm_GMD_segid_class$SegID, emm_GMD_segid_class$estimate) # reorder brain regions according to percent change in y-variable





# ========================================================================
# Grey Matter Volume (GMV) / Partial Volume Estimate 1 Volume (PVE1_Vol)
# ========================================================================

# Histogram of GMV
ggplot(rawdata_VThresh_notumor_GM, aes(x=PVE1_Vol)) + geom_histogram()

# Log transform and plot histogram
rawdata_VThresh_notumor_GM_log <- rawdata_VThresh_notumor_GM
rawdata_VThresh_notumor_GM_log$PVE1_Vol <- log10(rawdata_VThresh_notumor_GM$PVE1_Vol)
ggplot(rawdata_VThresh_notumor_GM_log, aes(x=PVE1_Vol)) + geom_histogram()

# Generate regression formulas and perform model selection via AIC minimization
indep_var = c("ParticipantType","SegID","Sex","Age")
dep_var = "log10(PVE1_Vol) ~"
regression_formulas = generate_regression_formulas(indep_var, dep_var)

model_sel_AIC_GMV = data.frame(matrix(ncol = 2, nrow = 0))
colnames(model_sel_AIC_GMV) <- c("Model", "AIC")
for (f in 1:length(regression_formulas)){
  print(f)
  melr_model_selection = lme(regression_formulas[[f]],
                             random=~1|ID,
                             data = rawdata_VThresh_notumor_GM) # recall contrast coding scheme set above
  summary_melr_model_selection = summary(melr_model_selection)
  model_sel_AIC_GMV = rbind(model_sel_AIC_GMV, list(Model=paste(deparse(regression_formulas[[f]]), collapse=" "), 
                                                    AIC=summary_melr_model_selection[["AIC"]]))
}
model_sel_AIC_GMV <- model_sel_AIC_GMV[order(model_sel_AIC_GMV$AIC), ]
write.csv(model_sel_AIC_GMV, "~/data/NCfiles/BrainMetabolism/ImagingNeuroscience/Regression_Results/GMV/GMV_ModelSelection_AIC.csv", row.names=FALSE)

# Run selected mixed effects model post AIC minimization (view exported csv file above first)
melr_model_GMV = lme(log10(PVE1_Vol) ~ ParticipantType + SegID + Sex + Age + ParticipantType:SegID,
                 random=~1|ID, 
                 data=rawdata_VThresh_notumor_GM) # note the contrast coding set above
summary(melr_model_GMV)

# Save regression results
model_results_GMV = data.frame(coef(summary(melr_model_GMV)))
model_results_GMV <- tibble::rownames_to_column(model_results_GMV, "Variables") 
write.table(model_results_GMV, file = "~/data/NCfiles/BrainMetabolism/ImagingNeuroscience/Regression_Results/GMV/GMV_MELR_Summary.csv")

# Estimated marginal means
emm_GMV = emmeans(melr_model_GMV, ~ParticipantType|SegID, type='response') # disregard warning message about contrast dropping - this is normal
emmdf_GMV = as.data.frame(emm_GMV)
write.table(emmdf_GMV, file = "~/data/NCfiles/BrainMetabolism/ImagingNeuroscience/Regression_Results/GMV/GMV_MELR_EMM.csv")

EMM_GMV = data.frame(contrast(emm_GMV, method='revpairwise', by='SegID', adjust='none', infer = c(TRUE, TRUE)))
EMM_GMV$Ratio_Percent = (EMM_GMV$ratio - 1)*100 # convert ratio of means to percent difference
write.table(EMM_GMV, file = "~/data/NCfiles/BrainMetabolism/ImagingNeuroscience/Regression_Results/GMV/GMV_MELR_EMM_t-test.csv")

emm_GMV_average = emmeans(melr_model_GMV, ~ParticipantType,type='response') # disregard warning message about contrast dropping - this is normal
emm_GMV_average_df=as.data.frame(emm_GMV_average)

emm_GMV_segid_class = left_join(EMM_GMV, SegID_class) 
emm_GMV_segid_class = emm_GMV_segid_class[emm_GMV_segid_class$Lobes != 'CerebralWhiteMatter',]
emm_GMV_segid_class$SegID = fct_reorder(emm_GMV_segid_class$SegID, emm_GMV_segid_class$Ratio_Percent) # reorder brain regions according to percent change in y-variable

# Boxplot of average GMV in patients vs controls
GMV_VThresh_EMM_allregions <- ggplot(rawdata_VThresh_notumor_GM, aes(x=ParticipantType, y=PVE1_Vol)) +
  geom_boxplot(data=rawdata_VThresh_notumor_GM, aes(x=ParticipantType, y=PVE1_Vol),
               color='black', position=position_dodge(width=0.85), outlier.shape=NA)+
  geom_point(data=emm_GMV_average_df, aes(x=ParticipantType, y=response),
             position=position_dodge(width=0.85), size=2, color='red3')+
  geom_errorbar(data=emm_GMV_average_df, aes(x=ParticipantType, y=response, ymin=lower.CL, ymax=upper.CL, width=0.2),
                position=position_dodge(width=0.85), color='red3')+
  geom_signif(comparisons=list(c('Control','Patient')), annotation = '*',
              y_position = 15, tip_length = 0, textsize = 5) +
  theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank(),
        panel.background=element_blank(), axis.line=element_line(color='black'),
        axis.text.x = element_text(color="black"), axis.text.y = element_text(color="black"),
        text=element_text(size=20), legend.position = 'top')+
  scale_y_continuous(breaks=seq(0,25,2))+
  coord_cartesian(ylim =  c(0, 18))+
  scale_x_discrete(labels=c('Controls','Patients'))+
  labs(x='', y='GMV',fill='')
GMV_VThresh_EMM_allregions

# Save figure as .eps and .tiff
ggsave(GMV_VThresh_EMM_allregions, file="Figure_5d.eps", width = 9, height = 12, units = "cm", device="eps", 
       path="/home/ncappelletto/data/NCfiles/BrainMetabolism/ImagingNeuroscience/Figure_5/")
ggsave(GMV_VThresh_EMM_allregions, file="Figure_5d.tiff", width = 9, height = 12, units = "cm", device="tiff", 
       dpi=300, path="/home/ncappelletto/data/NCfiles/BrainMetabolism/ImagingNeuroscience/Figure_5/")

# Region with max change: R TrIFG
TrIFG = rawdata_VThresh_notumor_GM[rawdata_VThresh_notumor_GM$SegID == 204 | rawdata_VThresh_notumor_GM$SegID == 205,]
TrIFG_emm <- subset(emmdf_GMV, SegID==204 | SegID==205)

# Set p-value asterisks for boxplots based on EMM t-tests
asterisks <- c("***", "*")     # SegID 204 → ***, SegID 205 → *
y_ast     <- c(4.8, 4.8)
y_bracket <- c(4.7, 4.7) 

# Boxplot of region with max change: R TrIFG
GMV_VThresh_EMM_TrIFG <- ggplot(TrIFG, aes(x=as.factor(SegID), y=PVE1_Vol, fill=ParticipantType)) +
  geom_boxplot(data=TrIFG, aes(x=as.factor(SegID), y=PVE1_Vol, fill=ParticipantType),
               color='black', position=position_dodge(width=0.85), outlier.shape=NA)+
  geom_point(data=TrIFG_emm, aes(x=as.factor(SegID), y=response, fill=ParticipantType),
             position = position_dodge(width=0.85), size=2, color='red3')+
  geom_errorbar(data=TrIFG_emm, aes(x=as.factor(SegID), y=response, ymin=lower.CL, ymax=upper.CL, width=0.2,
                                   fill=ParticipantType), position=position_dodge(width=0.85), color='red3')+
  annotate("text", x = 1, y = y_ast[1], label = asterisks[1], size = 5) +
  geom_segment(aes(x = 1 - 0.25, xend = 1 + 0.25,
                   y = y_bracket[1], yend = y_bracket[1])) +
  annotate("text", x = 2, y = y_ast[2], label = asterisks[2], size = 5) +
  geom_segment(aes(x = 2 - 0.25, xend = 2 + 0.25,
                   y = y_bracket[2], yend = y_bracket[2])) +
  theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank(),
        panel.background=element_blank(), axis.line=element_line(color='black'),
        axis.text.x = element_text(color="black"), axis.text.y = element_text(color="black"),
        text=element_text(size=20), legend.position = 'top')+
  scale_fill_manual(values=c('Control'="grey",'Patient'="grey45"), labels=c("Controls","Patients"))+
  scale_y_continuous(breaks=seq(0,10,0.5))+
  coord_cartesian(ylim =  c(1, 5))+
  scale_x_discrete(labels=c('R TrIFG','L TrIFG'))+
  labs(x='', y='GMV',fill='')
GMV_VThresh_EMM_TrIFG

# Save figure as .eps and .tiff
ggsave(GMV_VThresh_EMM_TrIFG, file="Figure_5e.eps", width = 9, height = 12, units = "cm", device="eps", 
       path="/home/ncappelletto/data/NCfiles/BrainMetabolism/ImagingNeuroscience/Figure_5/")
ggsave(GMV_VThresh_EMM_TrIFG, file="Figure_5e.tiff", width = 9, height = 12, units = "cm", device="tiff", 
       dpi=300, path="/home/ncappelletto/data/NCfiles/BrainMetabolism/ImagingNeuroscience/Figure_5/")

# Region with Min Change: L SPL
SPL = rawdata_VThresh_notumor_GM[rawdata_VThresh_notumor_GM$SegID == 198 | rawdata_VThresh_notumor_GM$SegID == 199,]
SPL_emm <- subset(emmdf_GMV, SegID==198 | SegID==199)

# Set p-value asterisks for boxplots based on EMM t-tests
asterisks <- c("n.s.", "n.s.")     # SegID 198 → n.s., SegID 199 → n.s.
y_ast     <- c(12.7, 12.7)
y_bracket <- c(12.4, 12.4)

# Boxplot of region with min change: L SPL
GMV_VThresh_EMM_SPL <- ggplot(SPL, aes(x=as.factor(SegID), y=PVE1_Vol, fill=ParticipantType)) +
  geom_boxplot(data=SPL, aes(x=as.factor(SegID), y=PVE1_Vol, fill=ParticipantType),
               color='black', position=position_dodge(width=0.85), outlier.shape=NA)+
  geom_point(data=SPL_emm, aes(x=as.factor(SegID), y=response, fill=ParticipantType),
             position=position_dodge(width=0.85), size=2, color='red3')+
  geom_errorbar(data=SPL_emm,aes(x=as.factor(SegID), y=response, ymin=lower.CL, ymax=upper.CL, width=0.2,
                                 fill=ParticipantType), position=position_dodge(width=0.85), color='red3')+
  annotate("text", x = 1, y = y_ast[1], label = asterisks[1], size = 5) +
  geom_segment(aes(x = 1 - 0.25, xend = 1 + 0.25,
                   y = y_bracket[1], yend = y_bracket[1])) +
  annotate("text", x = 2, y = y_ast[2], label = asterisks[2], size = 5) +
  geom_segment(aes(x = 2 - 0.25, xend = 2 + 0.25,
                   y = y_bracket[2], yend = y_bracket[2])) +
  theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank(),
        panel.background=element_blank(), axis.line=element_line(color='black'),
        axis.text.x = element_text(color="black"), axis.text.y = element_text(color="black"),
        text=element_text(size=20), legend.position = 'top')+
  scale_fill_manual(values=c('Control'="grey",'Patient'="grey45"),labels=c("Controls","Patients"))+
  scale_y_continuous(breaks=seq(0,20,0.8))+
  scale_x_discrete(labels=c('R SPL','L SPL'))+
  labs(x='', y='GMV',fill='')
GMV_VThresh_EMM_SPL

# save figures as .eps and .tiff
ggsave(GMV_VThresh_EMM_SPL, file="Figure_5f.eps", width = 9, height = 12, units = "cm", device="eps", 
       path="/home/ncappelletto/data/NCfiles/BrainMetabolism/ImagingNeuroscience/Figure_5/")
ggsave(GMV_VThresh_EMM_SPL, file="Figure_5f.tiff", width = 9, height = 12, units = "cm", device="tiff", 
       dpi=300, path="/home/ncappelletto/data/NCfiles/BrainMetabolism/ImagingNeuroscience/Figure_5/")





# ==============================================================
# White Matter Density (WMD) / Partial Volume Estimate 2 (PVE2)
# ==============================================================

# Histogram of WMD
ggplot(rawdata_VThresh_notumor_WM, aes(x=PVE2)) + geom_histogram()

# Log transform and plot histogram
rawdata_VThresh_notumor_WM_log <- rawdata_VThresh_notumor_WM
rawdata_VThresh_notumor_WM_log$PVE2 <- log10(rawdata_VThresh_notumor_WM$PVE2)
ggplot(rawdata_VThresh_notumor_WM_log, aes(x=PVE2)) + geom_histogram()

# Generate regression formulas and perform model selection via AIC minimization
indep_var = c("ParticipantType","SegID","Sex","Age")
dep_var = "PVE2 ~"
regression_formulas = generate_regression_formulas(indep_var, dep_var)

model_sel_AIC_WMD = data.frame(matrix(ncol = 2, nrow = 0))
colnames(model_sel_AIC_WMD) <- c("Model", "AIC")
for (f in 1:length(regression_formulas)){
  print(f)
  melr_model_selection = lme(regression_formulas[[f]],
                             random=~1|ID,
                             data = rawdata_VThresh_notumor_GM) # recall contrast coding scheme set above
  summary_melr_model_selection = summary(melr_model_selection)
  model_sel_AIC_WMD = rbind(model_sel_AIC_WMD, list(Model=paste(deparse(regression_formulas[[f]]), collapse=" "), 
                                                    AIC=summary_melr_model_selection[["AIC"]]))
}
model_sel_AIC_WMD <- model_sel_AIC_WMD[order(model_sel_AIC_WMD$AIC), ]
write.csv(model_sel_AIC_WMD, "~/data/NCfiles/BrainMetabolism/ImagingNeuroscience/Regression_Results/WMD/WMD_ModelSelection_AIC.csv", row.names=FALSE)

# Run selected mixed effects model post AIC minimization (view exported csv file above first)
melr_model_WMD = lme(PVE2 ~ ParticipantType + SegID + ParticipantType:SegID,
                 random=~1|ID, 
                 data=rawdata_VThresh_notumor_WM) # note the contrast coding set above
summary(melr_model_WMD)

# Save regression results
model_results_WMD = data.frame(coef(summary(melr_model_WMD)))
model_results_WMD <- tibble::rownames_to_column(model_results_WMD, "Variables") 
write.table(model_results_WMD , file = "~/data/NCfiles/BrainMetabolism/ImagingNeuroscience/Regression_Results/WMD/WMD_MELR_Summary.csv")

# Estimated marginal means
emm_WMD = emmeans(melr_model_WMD, ~ParticipantType|SegID, type='response') # disregard warning message about contrast dropping - this is normal
emmdf_WMD = as.data.frame(emm_WMD)
write.table(emmdf_WMD, file = "~/data/NCfiles/BrainMetabolism/ImagingNeuroscience/Regression_Results/WMD/WMD_MELR_EMM.csv")

EMM_WMD = data.frame(contrast(emm_WMD, method='revpairwise', by='SegID', adjust='none', infer = c(TRUE, TRUE)))
write.table(EMM_WMD, file = "~/data/NCfiles/BrainMetabolism/ImagingNeuroscience/Regression_Results/WMD/WMD_MELR_EMM_t-test.csv")

emm_WMD_average = emmeans(melr_model_WMD, ~ParticipantType, type='response') # disregard warning message about contrast dropping - this is normal
emm_WMD_average_df=as.data.frame(emm_WMD_average)

# Boxplot of average WMD in patients vs controls
WMD_VThresh_EMM_boxplot_allregions <- ggplot(rawdata_VThresh_notumor_WM, aes(x=ParticipantType, y=PVE2)) +
  geom_boxplot(data=rawdata_VThresh_notumor_WM, aes(x=ParticipantType, y=PVE2), 
               color='black', position=position_dodge(width=0.85), outlier.shape=NA)+
  geom_point(data=emm_WMD_average_df, aes(x=ParticipantType, y=emmean),
             position=position_dodge(width=0.85), size=2, color='red3')+
  geom_errorbar(data=emm_WMD_average_df, aes(x=ParticipantType, y=emmean, ymin=lower.CL, ymax=upper.CL, width=0.2),
                position=position_dodge(width=0.85), color='red3')+
  geom_signif(comparisons=list(c('Control','Patient')), annotation = 'n.s.',
              y_position = 0.98,tip_length = 0,textsize = 5) +
  theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank(),
        panel.background=element_blank(), axis.line=element_line(color='black'),
        axis.text.x = element_text(color="black"), axis.text.y = element_text(color="black"),
        text=element_text(size=20,family='arial'), legend.position = 'top')+
  scale_y_continuous(breaks=seq(0,1,0.02))+
  coord_cartesian(ylim=c(0.78,1))+
  scale_x_discrete(labels=c('Controls','Patients'))+
  labs(x='', y='WMD',fill='')
WMD_VThresh_EMM_boxplot_allregions

# Save figure as .eps and .tiff
ggsave(WMD_VThresh_EMM_boxplot_allregions, file="Figure_5i.eps", width = 8, height = 11, units = "cm", device="eps", 
       path="/home/ncappelletto/data/NCfiles/BrainMetabolism/ImagingNeuroscience/Figure_5/")
ggsave(WMD_VThresh_EMM_boxplot_allregions, file="Figure_5i.tiff", width = 8, height = 11, units = "cm", device="tiff", 
       dpi = 300, path="/home/ncappelletto/data/NCfiles/BrainMetabolism/ImagingNeuroscience/Figure_5/")





# ========================================================================
# White Matter Volume (WMV) / Partial Volume Estimate 2 Volume (PVE2_Vol)
# ========================================================================

# Histogram of WMV
ggplot(rawdata_VThresh_notumor_WM, aes(x=PVE2_Vol)) + geom_histogram()

# Log transform and plot histogram
rawdata_VThresh_notumor_WM_log <- rawdata_VThresh_notumor_WM
rawdata_VThresh_notumor_WM_log$PVE2_Vol <- log10(rawdata_VThresh_notumor_WM$PVE2_Vol)
ggplot(rawdata_VThresh_notumor_WM_log, aes(x=PVE2_Vol)) + geom_histogram()

# Generate regression formulas and perform model selection via AIC minimization
indep_var = c("ParticipantType","SegID","Sex","Age")
dep_var = "log10(PVE2_Vol) ~"
regression_formulas = generate_regression_formulas(indep_var, dep_var)

model_sel_AIC_WMV = data.frame(matrix(ncol = 2, nrow = 0))
colnames(model_sel_AIC_WMV) <- c("Model", "AIC")
for (f in 1:length(regression_formulas)){
  print(f)
  melr_model_selection = lme(regression_formulas[[f]],
                             random=~1|ID,
                             data = rawdata_VThresh_notumor_GM) # recall contrast coding scheme set above
  summary_melr_model_selection = summary(melr_model_selection)
  model_sel_AIC_WMV = rbind(model_sel_AIC_WMV, list(Model=paste(deparse(regression_formulas[[f]]), collapse=" "), 
                                                    AIC=summary_melr_model_selection[["AIC"]]))
}
model_sel_AIC_WMV <- model_sel_AIC_WMV[order(model_sel_AIC_WMV$AIC), ]
write.csv(model_sel_AIC_WMV, "~/data/NCfiles/BrainMetabolism/ImagingNeuroscience/Regression_Results/WMV/WMV_ModelSelection_AIC.csv", row.names=FALSE)

# Run selected mixed effects model post AIC minimization (view exported csv file above first)
melr_model_WMV = lme(log10(PVE2_Vol) ~ ParticipantType + SegID + ParticipantType:SegID,
                 random=~1|ID, 
                 data=rawdata_VThresh_notumor_WM) # note the contrast coding set above
summary(melr_model_WMV)

# Save regression results
model_results_WMV = data.frame(coef(summary(melr_model_WMV)))
model_results_WMV <- tibble::rownames_to_column(model_results_WMV, "Variables") 
write.table(model_results_WMV , file = "~/data/NCfiles/BrainMetabolism/ImagingNeuroscience/Regression_Results/WMV/WMV_MELR_Summary.csv")

# Estimated marginal means
emm_WMV = emmeans(melr_model_WMV, ~ParticipantType|SegID, type='response') # disregard warning message about contrast dropping - this is normal
emmdf_WMV = as.data.frame(emm_WMV)
write.table(emmdf_WMV, file = "~/data/NCfiles/BrainMetabolism/ImagingNeuroscience/Regression_Results/WMV/WMV_MELR_EMM.csv")

EMM_WMV = data.frame(contrast(emm_WMV, method='revpairwise',by='SegID',adjust='none', infer = c(TRUE, TRUE)))
EMM_WMV$Ratio_Percent = (EMM_WMV$ratio - 1)*100 # convert ratio of means to percent difference
write.table(EMM_WMV, file = "~/data/NCfiles/BrainMetabolism/ImagingNeuroscience/Regression_Results/WMV/WMV_MELR_EMM_t-test.csv")

emm_WMV_average = emmeans(melr_model_WMV, ~ParticipantType,type='response') # disregard warning message about contrast dropping - this is normal
emm_WMV_average_df=as.data.frame(emm_WMV_average)

# Boxplot of average WMV in patients vs controls
WMV_VThresh_EMM_boxplot_allregions <- ggplot(rawdata_VThresh_notumor_WM, aes(x=ParticipantType, y=PVE2_Vol)) +
  geom_boxplot(data=rawdata_VThresh_notumor_WM, aes(x=ParticipantType,y=PVE2_Vol), 
               color='black', position=position_dodge(width=0.85), outlier.shape=NA)+
  geom_point(data=emm_WMV_average_df, aes(x=ParticipantType,y=response),
             position=position_dodge(width=0.85), size=2, color='red3')+
  geom_errorbar(data=emm_WMV_average_df, aes(x=ParticipantType, y=response, ymin=lower.CL, ymax=upper.CL, width=0.2),
                position=position_dodge(width=0.85), color='red3')+
  geom_signif(comparisons=list(c('Control','Patient')), annotation = '*',
              y_position = 250, tip_length = 0, textsize = 5) +
  theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank(),
        panel.background=element_blank(), axis.line=element_line(color='black'),
        axis.text.x = element_text(color="black"), axis.text.y = element_text(color="black"),
        text=element_text(size=20), legend.position = 'top')+
  scale_y_continuous(breaks=seq(0,300,25))+
  coord_cartesian(ylim = c(0,285))+
  scale_x_discrete(labels=c('Controls','Patients'))+
  labs(x='', y='WMV',fill='')
WMV_VThresh_EMM_boxplot_allregions 

# save figures as .eps and .tiff
ggsave(WMV_VThresh_EMM_boxplot_allregions, file="Figure_5g.eps", width = 9, height = 12, units = "cm", device="eps", 
       path="/home/ncappelletto/data/NCfiles/BrainMetabolism/ImagingNeuroscience/Figure_5/")
ggsave(WMV_VThresh_EMM_boxplot_allregions, file="Figure_5g.tiff", width = 9, height = 12, units = "cm", device="tiff", 
       dpi = 300, path="/home/ncappelletto/data/NCfiles/BrainMetabolism/ImagingNeuroscience/Figure_5/")

# White matter regions with max and min change: R Cerebellum WM & L Cerebellum WM
CWM = rawdata_VThresh_notumor_WM[rawdata_VThresh_notumor_WM$SegID == 40 | rawdata_VThresh_notumor_WM$SegID == 41,]
CWM_emm <- subset(emmdf_WMV, SegID==40 | SegID==41)

# Set p-value asterisks for boxplots based on EMM t-tests
asterisks <- c("*", "n.s.")     # SegID 40 → *, SegID 41 → n.s.
y_ast     <- c(16.2, 16.5)
y_bracket <- c(16, 16)

# Boxplot of WM regions with max and min change: R Cerebellum WM & L Cerebellum WM
WMV_VThresh_EMM_CWM <- ggplot(CWM, aes(x=as.factor(SegID), y=PVE2_Vol, fill=ParticipantType)) +
  geom_boxplot(data=CWM, aes(x=as.factor(SegID), y=PVE2_Vol, fill=ParticipantType), 
               color='black', position=position_dodge(width=0.85), outlier.shape=NA)+
  geom_point(data=CWM_emm, aes(x=as.factor(SegID), y=response, fill=ParticipantType),
             position=position_dodge(width=0.85), size=2, color='red3')+
  geom_errorbar(data=CWM_emm, aes(x=as.factor(SegID), y=response, ymin=lower.CL, ymax=upper.CL, width=0.2,
                                 fill=ParticipantType), position=position_dodge(width=0.85), color='red3')+
  annotate("text", x = 1, y = y_ast[1], label = asterisks[1], size = 5) +
  geom_segment(aes(x = 1 - 0.25, xend = 1 + 0.25,
                   y = y_bracket[1], yend = y_bracket[1])) +
  annotate("text", x = 2, y = y_ast[2], label = asterisks[2], size = 5) +
  geom_segment(aes(x = 2 - 0.25, xend = 2 + 0.25,
                   y = y_bracket[2], yend = y_bracket[2])) +
  theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank(),
        panel.background=element_blank(),axis.line=element_line(color='black'),
        axis.text.x = element_text(color="black"),axis.text.y = element_text(color="black"),
        text=element_text(size=20),legend.position = 'top')+
  scale_fill_manual(values=c('Control'="grey",'Patient'="grey45"),labels=c("Controls","Patients"))+
  scale_y_continuous(breaks=seq(0,55,1.5))+
  coord_cartesian(ylim = c(5,18))+
  scale_x_discrete(labels=c('R CWM','L CWM'))+
  labs(x='', y='WMV',fill='')
WMV_VThresh_EMM_CWM

# save figures as .eps and .tiff
ggsave(WMV_VThresh_EMM_CWM, file="Figure_5h.eps", width = 9, height = 12, units = "cm", device="eps", 
       path="/home/ncappelletto/data/NCfiles/BrainMetabolism/ImagingNeuroscience/Figure_5/")
ggsave(WMV_VThresh_EMM_CWM, file="Figure_5h.tiff", width = 9, height = 12, units = "cm", device="tiff", 
       dpi=300, path="/home/ncappelletto/data/NCfiles/BrainMetabolism/ImagingNeuroscience/Figure_5/")





# =======================
# Brain Region Volume
# =======================

# Histogram of volume
ggplot(rawdata_VThresh_notumor, aes(x=VolumeCC)) + geom_histogram()

# Log transform and plot histogram
rawdata_VThresh_notumor_log <- rawdata_VThresh_notumor
rawdata_VThresh_notumor_log$VolumeCC <- log10(rawdata_VThresh_notumor$VolumeCC)
ggplot(rawdata_VThresh_notumor_log, aes(x=VolumeCC)) + geom_histogram()


# Generate regression formulas and perform model selection via AIC minimization
indep_var = c("ParticipantType","SegID","Sex","Age")
dep_var = "log10(VolumeCC) ~"
regression_formulas = generate_regression_formulas(indep_var, dep_var)

model_sel_AIC_Vol = data.frame(matrix(ncol = 2, nrow = 0))
colnames(model_sel_AIC_Vol) <- c("Model", "AIC")
for (f in 1:length(regression_formulas)){
  print(f)
  melr_model_selection = lme(regression_formulas[[f]],
                             random=~1|ID,
                             data = rawdata_VThresh_notumor_GM) # recall contrast coding scheme set above
  summary_melr_model_selection = summary(melr_model_selection)
  model_sel_AIC_Vol = rbind(model_sel_AIC_Vol, list(Model=paste(deparse(regression_formulas[[f]]), collapse=" "), 
                                                    AIC=summary_melr_model_selection[["AIC"]]))
}
model_sel_AIC_Vol <- model_sel_AIC_Vol[order(model_sel_AIC_Vol$AIC), ]
write.csv(model_sel_AIC_Vol, "~/data/NCfiles/BrainMetabolism/ImagingNeuroscience/Regression_Results/volume/volume_ModelSelection_AIC.csv", row.names=FALSE)

# Run selected mixed effects model post AIC minimization (view exported csv file above first)
melr_model_Vol = lme(log10(VolumeCC) ~ ParticipantType + SegID + Sex + ParticipantType:SegID,
                 random=~1|ID, 
                 data=rawdata_VThresh_notumor) # note the contrast coding set above
summary(melr_model_Vol)

# make a dataframe with model coefficients and pvalues 
model_results_Vol = data.frame(coef(summary(melr_model_Vol)))
model_results_Vol <- tibble::rownames_to_column(model_results_Vol,"Variables") 
write.table(model_results_Vol, file = "~/data/NCfiles/BrainMetabolism/ImagingNeuroscience/Regression_Results/volume/volume_MELR_Summary.csv")

# Estimated marginal means
emm_Vol = emmeans(melr_model_Vol, ~ParticipantType|SegID,type='response') # disregard warning message about contrast dropping - this is normal
emmdf_Vol = as.data.frame(emm_Vol)
write.table(emmdf_Vol, file = "~/data/NCfiles/BrainMetabolism/ImagingNeuroscience/Regression_Results/volume/volume_MELR_EMM.csv")

EMM_Vol = data.frame(contrast(emm_Vol,method='revpairwise',by='SegID',adjust='none',infer = c(TRUE, TRUE)))
EMM_Vol$Ratio_Percent = (EMM_Vol$ratio - 1)*100 # convert ratio of means to percent difference
write.table(EMM_Vol, file = "~/data/NCfiles/BrainMetabolism/ImagingNeuroscience/Regression_Results/volume/volume_MELR_EMM_t-test.csv")

emm_Vol_average = emmeans(melr_model_Vol, ~ParticipantType,type='response') # disregard warning message about contrast dropping - this is normal
emm_Vol_average_df=as.data.frame(emm_Vol_average)

emm_Vol_segid_class = left_join(EMM_Vol, SegID_class) 
emm_Vol_segid_class = emm_Vol_segid_class[emm_Vol_segid_class$Lobes != 'CerebralWhiteMatter',]
emm_Vol_segid_class$SegID = fct_reorder(emm_Vol_segid_class$SegID, emm_Vol_segid_class$Ratio_Percent) # reorder brain regions according to percent change in y-variable

# Run ANOVA and Tukey Post Hoc test to assess differences between brain lobes
# p-values and asterisks for boxplots
anova_Vol <- aov(Ratio_Percent ~ Lobes, data = emm_Vol_segid_class) # anova between lobes
tukey_Vol_result <- TukeyHSD(anova_Vol) 
tukey_Vol_df <- as.data.frame(tukey_Vol_result$Lobes)
tukey_Vol_df$comparison <- rownames(tukey_Vol_df)
tukey_Vol_df$signif <- sapply(tukey_Vol_df$`p adj`, pval_to_asterisk)
tukey_Vol_df
write.table(tukey_Vol_df, file = "~/data/NCfiles/BrainMetabolism/ImagingNeuroscience/Regression_Results/volume/volume_lobes_tukey.csv", row.names=FALSE)
tukey_Vol_df= tukey_Vol_df[tukey_Vol_df$'p adj'<=0.05,]
comparisons <- strsplit(tukey_Vol_df$comparison, "-") # Create list of comparisons for geom_signif
y_positions <- seq(from = max(emm_Vol_segid_class$Ratio_Percent) * 1, 
                   by = max(emm_Vol_segid_class$Ratio_Percent) * 0.15, length.out = nrow(tukey_Vol_df))

# Volume average vs brain lobes
volume_VThresh_EMM_lobe <- ggplot(emm_Vol_segid_class, aes(x=reorder(Lobes, -Ratio_Percent, FUN=median), y=Ratio_Percent))+
  geom_boxplot(outlier.shape=NA, color='black')+
  geom_signif(comparisons = comparisons, annotations = tukey_Vol_df$signif,
              y_position=y_positions, tip_length = 0, vjust=0.5, textsize = 4) +
  theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank(),
        panel.background=element_blank(), axis.line=element_line(color='black'),
        axis.text.x = element_text(angle=45, hjust=1.1, color="black"),
        axis.text.y = element_text(color="black"),
        text=element_text(size=20), legend.position="none")+
  scale_y_continuous(breaks=seq(-40,20,5))+
  labs(x='', y=expression('EMM'~ Delta * '% volume'),fill='')
volume_VThresh_EMM_lobe 

# Save figure as .eps and .tiff
ggsave(volume_VThresh_EMM_lobe, file="Figure_3d.eps", width = 9, height = 12, units = "cm", device="eps", 
       path="/home/ncappelletto/data/NCfiles/BrainMetabolism/ImagingNeuroscience/Figure_3/")
ggsave(volume_VThresh_EMM_lobe, file="Figure_3d.tiff", width = 9, height = 12, units = "cm", device="tiff", 
       dpi=300, path="/home/ncappelletto/data/NCfiles/BrainMetabolism/ImagingNeuroscience/Figure_3/")

# Boxplot of average WMV in patients vs controls
volume_VThresh_EMM_allregions <- ggplot(rawdata_VThresh_notumor, aes(x=ParticipantType, y=VolumeCC)) +
  geom_boxplot(data=rawdata_VThresh_notumor, aes(x=ParticipantType, y=VolumeCC), 
               color='black', position=position_dodge(width=0.85), outlier.shape=NA)+
  geom_point(data=emm_Vol_average_df, aes(x=ParticipantType, y=response),
             position=position_dodge(width=0.85), size=2, color='red3')+
  geom_errorbar(data=emm_Vol_average_df, aes(x=ParticipantType, y=response, ymin=lower.CL, ymax=upper.CL, width=0.2),
                position=position_dodge(width=0.85), color='red3')+
  geom_signif(comparisons=list(c('Control','Patient')), annotation = '***',
              y_position = 13,tip_length = 0,textsize = 5) +
  theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank(),
        panel.background=element_blank(), axis.line=element_line(color='black'),
        axis.text.x = element_text(color="black"), axis.text.y = element_text(color="black"),
        text=element_text(size=20), legend.position = 'top')+
  scale_y_continuous(breaks=seq(0,50,3))+
  coord_cartesian(ylim =  c(0, 26))+
  scale_x_discrete(labels=c('Controls','Patients'))+
  labs(x='', y='Volume (cc)',fill='')
volume_VThresh_EMM_allregions

# Save figure as .eps and .tiff
ggsave(volume_VThresh_EMM_allregions, file="Figure_5a.eps", width = 9, height = 12, units = "cm", device="eps", 
       path="/home/ncappelletto/data/NCfiles/BrainMetabolism/ImagingNeuroscience/Figure_5/")
ggsave(volume_VThresh_EMM_allregions, file="Figure_5a.tiff", width = 9, height = 12, units = "cm", device="tiff", 
       dpi=300, path="/home/ncappelletto/data/NCfiles/BrainMetabolism/ImagingNeuroscience/Figure_5/")

# Region with max change: R AnG
AnG = rawdata_VThresh_notumor[rawdata_VThresh_notumor$SegID == 106 | rawdata_VThresh_notumor$SegID == 107,]
AnG_emm <- subset(emmdf_Vol, SegID==106 | SegID==107)

# Set p-value asterisks for boxplots based on EMM t-tests
asterisks <- c("***", "**")   # SegID 106 → ***, SegID 107 → **
y_ast     <- c(18.2, 18.2)
y_bracket <- c(18, 18)

# Boxplot of region with max change: R AnG
volume_VThresh_EMM_AnG <- ggplot(AnG, aes(x=as.factor(SegID), y=VolumeCC, fill=ParticipantType)) +
  geom_boxplot(data=AnG, aes(x=as.factor(SegID), y=VolumeCC, fill=ParticipantType), 
               color='black', position=position_dodge(width=0.85), outlier.shape=NA)+
  geom_point(data=AnG_emm, aes(x=as.factor(SegID), y=response, fill=ParticipantType),
             position=position_dodge(width=0.85), size=2, color='red3')+
  geom_errorbar(data=AnG_emm,aes(x=as.factor(SegID), y=response, ymin=lower.CL, ymax=upper.CL, width=0.2,
                                 fill=ParticipantType), position=position_dodge(width=0.85), color='red3')+
  annotate("text", x = 1, y = y_ast[1], label = asterisks[1], size = 5) +
  geom_segment(aes(x = 1 - 0.25, xend = 1 + 0.25,
                   y = y_bracket[1], yend = y_bracket[1])) +
  annotate("text", x = 2, y = y_ast[2], label = asterisks[2], size = 5) +
  geom_segment(aes(x = 2 - 0.25, xend = 2 + 0.25,
                   y = y_bracket[2], yend = y_bracket[2])) +
  theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank(),
        panel.background=element_blank(), axis.line=element_line(color='black'),
        axis.text.x = element_text(color="black"), axis.text.y = element_text(color="black"),
        text=element_text(size=20),legend.position = 'top')+
  scale_fill_manual(values=c('Control'="grey",'Patient'="grey45"),labels=c("Controls","Patients"))+
  scale_y_continuous(breaks=seq(0,30,1))+
  coord_cartesian(ylim = c(8,19))+
  scale_x_discrete(labels=c('R AnG','L AnG'))+
  labs(x='', y='Volume (cc)',fill='')
volume_VThresh_EMM_AnG

# Save figure as .eps and .tiff
ggsave(volume_VThresh_EMM_AnG, file="Figure_5b.eps", width = 9, height = 12, units = "cm", device="eps", 
       path="/home/ncappelletto/data/NCfiles/BrainMetabolism/ImagingNeuroscience/Figure_5/")
ggsave(volume_VThresh_EMM_AnG, file="Figure_5b.tiff", width = 9, height = 12, units = "cm", device="tiff", 
       dpi=300, path="/home/ncappelletto/data/NCfiles/BrainMetabolism/ImagingNeuroscience/Figure_5/")

# Region with min change: L Calc
Calc = rawdata_VThresh_notumor[rawdata_VThresh_notumor$SegID == 108 | rawdata_VThresh_notumor$SegID == 109,]
Calc_emm <- subset(emmdf_Vol, SegID==108 | SegID==109)

# Set p-value asterisks for boxplots based on EMM t-tests
asterisks <- c("n.s.", "n.s.")   # SegID 108 → ns, SegID 109 → ns
y_ast     <- c(7.2, 7.2)
y_bracket <- c(7, 7)

# Boxplot of region with max change: L Calc
volume_VThresh_EMM_Calc <- ggplot(Calc, aes(x=as.factor(SegID), y=VolumeCC, fill=ParticipantType)) +
  geom_boxplot(data=Calc, aes(x=as.factor(SegID), y=VolumeCC, fill=ParticipantType), 
               color='black', position=position_dodge(width=0.85), outlier.shape=NA)+
  geom_point(data=Calc_emm,aes(x=as.factor(SegID), y=response, fill=ParticipantType),
             position=position_dodge(width=0.85), size=2, color='red3')+
  geom_errorbar(data=Calc_emm,aes(x=as.factor(SegID), y=response, ymin=lower.CL, ymax=upper.CL, width=0.2,
                                  fill=ParticipantType), position=position_dodge(width=0.85), color='red3')+
  annotate("text", x = 1, y = y_ast[1], label = asterisks[1], size = 5) +
  geom_segment(aes(x = 1 - 0.25, xend = 1 + 0.25,
                   y = y_bracket[1], yend = y_bracket[1])) +
  annotate("text", x = 2, y = y_ast[2], label = asterisks[2], size = 5) +
  geom_segment(aes(x = 2 - 0.25, xend = 2 + 0.25,
                   y = y_bracket[2], yend = y_bracket[2])) +
  theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank(),
        panel.background=element_blank(), axis.line=element_line(color='black'),
        axis.text.x = element_text(color="black"), axis.text.y = element_text(color="black"),
        text=element_text(size=20), legend.position = 'top')+
  scale_fill_manual(values=c('Control'="grey",'Patient'="grey45"),labels=c("Controls","Patients"))+
  scale_y_continuous(breaks=seq(0,12,0.5))+
  scale_x_discrete(labels=c('R Calc','L Calc'))+
  labs(x='', y='Volume (cc)',fill='')
volume_VThresh_EMM_Calc

# Save figure as .eps and .tiff
ggsave(volume_VThresh_EMM_Calc, file="Figure_5c.eps", width = 9, height = 12, units = "cm", device="eps", 
       path="/home/ncappelletto/data/NCfiles/BrainMetabolism/ImagingNeuroscience/Figure_5/")
ggsave(volume_VThresh_EMM_Calc, file="Figure_5c.tiff", width = 9, height = 12, units = "cm", device="tiff", 
       dpi=300, path="/home/ncappelletto/data/NCfiles/BrainMetabolism/ImagingNeuroscience/Figure_5/")

