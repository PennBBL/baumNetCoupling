rm(list = ls())

#############################
## Load Relevant Libraries ##
#############################
require(ggplot2)
require(Hmisc)
require(MASS)
require(mgcv)
require(ppcor)
require(stringr)
require(visreg)
require(parallel)
require(multilevel)
require(stats)
require(Formula)
require(lavaan)
require(lme4)
require(nlme)
require(gamm4)
require(car)
require(corrplot)
require(stringr)
require(R.matlab)
require(pracma)
require(gridExtra)
require(dimRed)
require(sna)

###################################################
## Load sample construction data (demographics) ##
##################################################
dti_nbackFC_rest_sample <- read.csv("/data/jux/BBL/projects/pncBaumStructFunc/replication/demographics/n727_dti_nbackFC_rest_sample_demographics.csv")

## Define dimensions of the data
nsub <- dim(dti_nbackFC_rest_sample) [1]
nreg <- 400
nedge <- 79800

## Schaefer400 NodeNames
nodeNames <- read.table("/data/joy/BBL/applications/xcpEngine/atlas/schaefer400/schaefer400NodeNames.txt", header=FALSE)

## Myelin Map
mean_myelin <- read.table("/data/jux/BBL/projects/pncBaumStructFunc/myelin_maps/output/hayashi_mean_Myelin_schaefer400.txt", header=FALSE)
colnames(mean_myelin) <- "mean_myelin"

## Margulies FC Gradient
margulies_gradient <- read.table("/data/jux/BBL/projects/pncBaumStructFunc/cortical_scaling/schaefer400x17_mean_regional_margulies_gradient.txt", header=FALSE)
colnames(margulies_gradient) <- "margulies_gradient"

rescaled_evo_expansion <- read.table("/data/jux/BBL/projects/pncBaumStructFunc/cortical_scaling/schaefer400_rescaled_regional_evoExpansion.txt", header=FALSE)
colnames(rescaled_evo_expansion) <- "rescaled_evo_expansion"

restFC_pca_comp1 <- read.table("/data/jux/BBL/projects/pncBaumStructFunc/replication/network_metrics/node_features/n727_Schaefer400_Zscore_mean_restFC_PCA_comp1_loadings.txt", header=FALSE)
colnames(restFC_pca_comp1) <- "restFC_pca_comp1"

restFC_pca_comp2 <- read.table("/data/jux/BBL/projects/pncBaumStructFunc/replication/network_metrics/node_features/n727_Schaefer400_Zscore_mean_restFC_PCA_comp2_loadings.txt", header=FALSE)
colnames(restFC_pca_comp2) <- "restFC_pca_comp2"

######################################################
## Read in Regional Coupling (calculated in MATLAB) ##
######################################################
rest_reg_coupling <- read.csv("/data/jux/BBL/projects/pncBaumStructFunc/replication/network_metrics/n727_Schaefer400_thresh_norm_probSC_restFC_regional_coupling_Spearman_r.txt", header=FALSE)
rest_reg_coupling <- as.data.frame(rest_reg_coupling)

for(i in 1:nreg) {
  colnames(rest_reg_coupling)[i] <- paste("rest_regCoup_V", i, sep = "")
}


nback_reg_coupling <- read.csv("/data/jux/BBL/projects/pncBaumStructFunc/replication/network_metrics/n727_Schaefer400_thresh_norm_probSC_nbackFC_regional_coupling_Spearman_r.txt", header=FALSE)
nback_reg_coupling <- as.data.frame(nback_reg_coupling)

for(i in 1:nreg) {
  colnames(nback_reg_coupling)[i] <- paste("nback_regCoup_V", i, sep = "")
}
################################################

#############################
## Delta PC (Nback - Rest) ##
#############################
regional_delta_posPC <- read.csv("/data/jux/BBL/projects/pncBaumStructFunc/replication/network_metrics/node_features/n727_Schaefer400_Yeo7_regional_DELTA_nbackPC-restPC.txt", header=FALSE)
for(i in 1:nreg) {
  colnames(regional_delta_posPC)[i] <- paste("delta_posPC_V", i, sep = "")
}

regional_delta_negPC <- read.csv("/data/jux/BBL/projects/pncBaumStructFunc/replication/network_metrics/node_features/n727_Schaefer400_Yeo7_regional_DELTA_Neg_nbackPC-restPC.txt", header=FALSE)
for(i in 1:nreg) {
  colnames(regional_delta_negPC)[i] <- paste("delta_negPC_V", i, sep = "")
}

#####################################
## One-sample t-test for Coupling  ##
#####################################
## Apply one-sample T-test across all regions
regCoupling_pvals <- as.numeric(lapply(nback_reg_coupling[,1:400],function(x) t.test(x, mu=0)$p.value))
length(regCoupling_pvals)

## FDR Correction
FDRcorr_regCoupling_pvals <- as.data.frame(p.adjust(regCoupling_pvals, method="fdr"))
colnames(FDRcorr_regCoupling_pvals) <- "FDRcorr_regCoupling_pvals"
FDRcorr_regCoupling_pvals$reg_idx <- 1:400

coupling_results <- cbind(regCoupling_pvals, FDRcorr_regCoupling_pvals)

sig_FDRcorr_regCoupling_pvals <- as.data.frame(subset(coupling_results, coupling_results$FDRcorr_regCoupling_pvals < 0.05))
dim(sig_FDRcorr_regCoupling_pvals)

###############################################
## Merge regional coupling with demographics ##
###############################################
reg_n727_df <- cbind(nback_reg_coupling, rest_reg_coupling, regional_delta_posPC, regional_delta_negPC, dti_nbackFC_rest_sample)

## Set categorical variables as factors
reg_n727_df$sex <- as.factor(reg_n727_df$sex)

## Create non-rounded Age in years measure
reg_n727_df$age <- reg_n727_df$ageAtScan1 / 12

## Create geometric mean motion measure for combining Rest and Nback
reg_n727_df$mean_nbackRest_motion <- sqrt(reg_n727_df$nbackRelMeanRMSMotion * reg_n727_df$restRelMeanRMSMotion)

###############################################
### Run node-wise GAM estimating Age Effect ###
###############################################
covariates=" ~  s(age, k=4) + dti64MeanRelRMS + nbackRelMeanRMSMotion + sex"    
m <- mclapply(names(reg_n727_df[,1:400]), function(x) {as.formula(paste(x, covariates, sep=""))},mc.cores=2)
regCoupling_Age_pvals <- mclapply(m, function(x) { summary(gam(formula = x, fx=TRUE,data=reg_n727_df, REML=T))$s.table[1,4]}, mc.cores=2)
regCoupling_Age_pvals <- as.data.frame(regCoupling_Age_pvals)
regCoupling_Age_pvals <- t(regCoupling_Age_pvals)
regCoupling_Age_pvals <- as.data.frame(regCoupling_Age_pvals)
colnames(regCoupling_Age_pvals) <- "regCoupling_Age_pvals"
regCoupling_Age_pvals$Node_index <- 1:400

n <-1

for(i in 1:400){
  ## Run Linear Regression Model to get t-statistic and direction of effect
  tmp_lm <- lm(as.formula(paste("nback_regCoup_V", i, "~ age + nbackRelMeanRMSMotion + dti64MeanRelRMS + sex", sep='')), data= reg_n727_df)
  ## Output results
  regCoupling_Age_pvals$lm_tvals[i] <- summary(tmp_lm)$coefficients[2,3]
  n <- n + 1
}

## Calculate Z-scores for spline age p-values
regCoupling_Age_pvals$AgeEffect_Zscore <- 0
regCoupling_Age_pvals$AgeEffect_Zscore <- qnorm(regCoupling_Age_pvals$regCoupling_Age_pvals,lower.tail=FALSE)

## Set Z-score sign to positive/negative based on T-value ###
regCoupling_Age_pvals$AgeEffect_Zscore[which(regCoupling_Age_pvals$lm_tvals < 0)] <- -(regCoupling_Age_pvals$AgeEffect_Zscore[which(regCoupling_Age_pvals$lm_tvals < 0)])

## FDR correction
FDRcorr_regCoupling_Age_pvals <- p.adjust(regCoupling_Age_pvals$regCoupling_Age_pvals, method="fdr")
FDRcorr_regCoupling_Age_pvals  <- cbind(regCoupling_Age_pvals,FDRcorr_regCoupling_Age_pvals, nodeNames)

## Define subset of regions showing significant Age effects
sig_FDRcorr_regCoupling_Age_pvals <- FDRcorr_regCoupling_Age_pvals[which(FDRcorr_regCoupling_Age_pvals$FDRcorr_regCoupling_Age_pvals <.05),]

dim(sig_FDRcorr_regCoupling_Age_pvals)

## Number of significant positive effects
sum(sig_FDRcorr_regCoupling_Age_pvals$AgeEffect_Zscore > 0)
## Number of significant negative effects
sum(sig_FDRcorr_regCoupling_Age_pvals$AgeEffect_Zscore < 0)

sigNeg_results <- subset(sig_FDRcorr_regCoupling_Age_pvals, sig_FDRcorr_regCoupling_Age_pvals$AgeEffect_Zscore < 0)
sigPos_results <- subset(sig_FDRcorr_regCoupling_Age_pvals, sig_FDRcorr_regCoupling_Age_pvals$AgeEffect_Zscore > 0)

## Significance Threshold
min(abs(sig_FDRcorr_regCoupling_Age_pvals$AgeEffect_Zscore))

## Write out results
write.table(FDRcorr_regCoupling_Age_pvals$AgeEffect_Zscore, "/data/jux/BBL/projects/pncBaumStructFunc/replication/group_stats/n727_Schaefer400_threshNormProbSC_nbackFC_Spearman_regCoupling_gam_Age_Zscore.txt", col.names = FALSE, row.names=FALSE)

write.table(FDRcorr_regCoupling_Age_pvals$lm_tvals, "/data/jux/BBL/projects/pncBaumStructFunc/replication/group_stats/n727_Schaefer400_threshNormProbSC_nbackFC_Spearman_regCoupling_gam_Age_lm_tstat.txt", col.names = FALSE, row.names=FALSE)


############################################################
## Rest Coupling: Run node-wise GAM estimating Age Effect ##
############################################################
covariates=" ~  s(age, k=4) + dti64MeanRelRMS + restRelMeanRMSMotion + sex"    
m <- mclapply(names(reg_n727_df[,401:800]), function(x) {as.formula(paste(x, covariates, sep=""))},mc.cores=2)
rest_regCoupling_Age_pvals <- mclapply(m, function(x) { summary(gam(formula = x, fx=TRUE,data=reg_n727_df, REML=T))$s.table[1,4]},mc.cores=2)
rest_regCoupling_Age_pvals <- as.data.frame(rest_regCoupling_Age_pvals)
rest_regCoupling_Age_pvals <- t(rest_regCoupling_Age_pvals)
rest_regCoupling_Age_pvals <- as.data.frame(rest_regCoupling_Age_pvals)
colnames(rest_regCoupling_Age_pvals) <- "rest_regCoupling_Age_pvals"
rest_regCoupling_Age_pvals$Node_index <- 1:400

## Apply linear mixed effects model to get t statistic 
## -- use to determine direction of GAM effects 
n <-1

for(i in 1:400){
  ## Run Linear Regression Model
  tmp_lm <- lm(as.formula(paste("rest_regCoup_V", i, "~ age + restRelMeanRMSMotion + dti64MeanRelRMS + sex", sep='')), data= reg_n727_df)
  ## Output results
  rest_regCoupling_Age_pvals$lm_tvals[i] <- summary(tmp_lm)$coefficients[2,3]
  n <- n + 1
}

## Calculate Z-scores for spline age p-values
rest_regCoupling_Age_pvals$AgeEffect_Zscore <- 0
rest_regCoupling_Age_pvals$AgeEffect_Zscore <- qnorm(rest_regCoupling_Age_pvals$rest_regCoupling_Age_pvals,lower.tail=FALSE)

## Set Z-score sign to positive/negative based on T-value ###
rest_regCoupling_Age_pvals$AgeEffect_Zscore[which(rest_regCoupling_Age_pvals$lm_tvals < 0)] <- -(rest_regCoupling_Age_pvals$AgeEffect_Zscore[which(rest_regCoupling_Age_pvals$lm_tvals < 0)])

## FDR correction
FDRcorr_rest_regCoupling_Age_pvals <- p.adjust(rest_regCoupling_Age_pvals$rest_regCoupling_Age_pvals, method="fdr")
FDRcorr_rest_regCoupling_Age_pvals  <- cbind(rest_regCoupling_Age_pvals,FDRcorr_rest_regCoupling_Age_pvals, nodeNames)

## Define subset of regions showing significant Age effects
sig_FDRcorr_rest_regCoupling_Age_pvals <- FDRcorr_rest_regCoupling_Age_pvals[which(FDRcorr_rest_regCoupling_Age_pvals$FDRcorr_rest_regCoupling_Age_pvals <.05),]

dim(sig_FDRcorr_rest_regCoupling_Age_pvals)

## Number of significant positive effects
sum(sig_FDRcorr_rest_regCoupling_Age_pvals$AgeEffect_Zscore > 0)
## Number of significant negative effects
sum(sig_FDRcorr_rest_regCoupling_Age_pvals$AgeEffect_Zscore < 0)

## Significance Threshold
min(abs(sig_FDRcorr_rest_regCoupling_Age_pvals$AgeEffect_Zscore))

## Write out results
write.table(FDRcorr_rest_regCoupling_Age_pvals$AgeEffect_Zscore, "/data/jux/BBL/projects/pncBaumStructFunc/replication/group_stats/n727_Schaefer400_threshNormProbSC_restFC_Spearman_regCoupling_gam_Age_Zscore.txt", col.names = FALSE, row.names=FALSE)

##########################
## Merge COGNITIVE data ##
##########################
cognitive<- read.csv("/data/joy/BBL/studies/pnc/n1601_dataFreeze/cnb/n1601_cnb_factorscores.csv")

cog_df <- merge(reg_n727_df, cognitive, by=c("bblid","scanid"))

############################################################
## F1_Exec_Comp_Res_Accuracy with Regional Coupling (GAM) ##
############################################################
covariates=" ~ F1_Exec_Comp_Res_Accuracy + s(age, k=4) + dti64MeanRelRMS + nbackRelMeanRMSMotion + sex"    
m <- mclapply(names(cog_df[,3:402]), function(x) {as.formula(paste(x, covariates, sep=""))},mc.cores=2)
regCoupling_nbackFC_F1ExecAcc_pvals <- mclapply(m, function(x) { summary(gam(formula = x, fx=TRUE,data=cog_df, REML=T))$p.pv[2]},mc.cores=2)
regCoupling_nbackFC_F1ExecAcc_pvals <- as.data.frame(regCoupling_nbackFC_F1ExecAcc_pvals)
regCoupling_nbackFC_F1ExecAcc_pvals <- t(regCoupling_nbackFC_F1ExecAcc_pvals )
regCoupling_nbackFC_F1ExecAcc_pvals <- as.data.frame(regCoupling_nbackFC_F1ExecAcc_pvals)
colnames(regCoupling_nbackFC_F1ExecAcc_pvals ) <- "regCoupling_nbackFC_F1ExecAcc_pvals "
regCoupling_nbackFC_F1ExecAcc_pvals$Node_index <- 1:400
regCoupling_nbackFC_F1ExecAcc_pvals$nodeNames <- nodeNames

## Extract test statistic for effect of interest
regCoupling_nbackFC_F1ExecAcc_pvals$F1ExecAcc_gam_tstats <- mclapply(m, function(x) { summary(gam(formula = x, fx=TRUE,data=cog_df, REML=T))$p.t[2]},mc.cores=2)
regCoupling_nbackFC_F1ExecAcc_pvals$F1ExecAcc_gam_tstats <- as.numeric(regCoupling_nbackFC_F1ExecAcc_pvals$F1ExecAcc_gam_tstats)

## FDR correction
FDRcorr_regCoupling_nbackFC_F1ExecAcc_pvals <- p.adjust(regCoupling_nbackFC_F1ExecAcc_pvals$regCoupling_nbackFC_F1ExecAcc_pvals , method="fdr")
FDRcorr_regCoupling_nbackFC_F1ExecAcc_pvals <- cbind(regCoupling_nbackFC_F1ExecAcc_pvals ,FDRcorr_regCoupling_nbackFC_F1ExecAcc_pvals, nodeNames)

## Define subset of regions showing significant Age effects on Yeo PC
sig_FDRcorr_regCoupling_nbackFC_F1ExecAcc_pvals  <- FDRcorr_regCoupling_nbackFC_F1ExecAcc_pvals [which(FDRcorr_regCoupling_nbackFC_F1ExecAcc_pvals$FDRcorr_regCoupling_nbackFC_F1ExecAcc_pvals  <.05),]
dim(sig_FDRcorr_regCoupling_nbackFC_F1ExecAcc_pvals )

# Number of significant positive effects
sum(sig_FDRcorr_regCoupling_nbackFC_F1ExecAcc_pvals$F1ExecAcc_gam_tstats > 0)
# Number of significant negative effects
sum(sig_FDRcorr_regCoupling_nbackFC_F1ExecAcc_pvals$F1ExecAcc_gam_tstats < 0)

min(abs(sig_FDRcorr_regCoupling_nbackFC_F1ExecAcc_pvals$F1ExecAcc_gam_tstats))

sigNeg_results <- subset(sig_FDRcorr_regCoupling_nbackFC_F1ExecAcc_pvals, sig_FDRcorr_regCoupling_nbackFC_F1ExecAcc_pvals$F1ExecAcc_gam_tstats < 0)
sigPos_results <- subset(sig_FDRcorr_regCoupling_nbackFC_F1ExecAcc_pvals, sig_FDRcorr_regCoupling_nbackFC_F1ExecAcc_pvals$F1ExecAcc_gam_tstats > 0)

## Write out results
write.table(FDRcorr_regCoupling_nbackFC_F1ExecAcc_pvals$F1ExecAcc_gam_tstats, "/data/jux/BBL/projects/pncBaumStructFunc/replication/group_stats/n727_Schaefer400_threshNormProbSC_nbackFC_regCoupling_gam_F1Exec_tstat.txt", col.names = FALSE, row.names=FALSE)


###################################################################
## F1_Exec_Comp_Res_Accuracy with RestFC Regional Coupling (GAM) ##
###################################################################
covariates=" ~ F1_Exec_Comp_Res_Accuracy + s(age, k=4) + dti64MeanRelRMS + restRelMeanRMSMotion + sex"    
m <- mclapply(names(cog_df[,403:802]), function(x) {as.formula(paste(x, covariates, sep=""))},mc.cores=2)
regCoupling_restFC_F1ExecAcc_pvals <- mclapply(m, function(x) { summary(gam(formula = x, fx=TRUE,data=cog_df, REML=T))$p.pv[2]},mc.cores=2)
regCoupling_restFC_F1ExecAcc_pvals <- as.data.frame(regCoupling_restFC_F1ExecAcc_pvals)
regCoupling_restFC_F1ExecAcc_pvals <- t(regCoupling_restFC_F1ExecAcc_pvals )
regCoupling_restFC_F1ExecAcc_pvals <- as.data.frame(regCoupling_restFC_F1ExecAcc_pvals)
colnames(regCoupling_restFC_F1ExecAcc_pvals ) <- "regCoupling_restFC_F1ExecAcc_pvals "
regCoupling_restFC_F1ExecAcc_pvals$Node_index <- 1:400
regCoupling_restFC_F1ExecAcc_pvals$nodeNames <- nodeNames

## Extract test statistic for effect of interest
regCoupling_restFC_F1ExecAcc_pvals$F1ExecAcc_gam_tstats <- mclapply(m, function(x) { summary(gam(formula = x, fx=TRUE,data=cog_df, REML=T))$p.t[2]},mc.cores=2)
regCoupling_restFC_F1ExecAcc_pvals$F1ExecAcc_gam_tstats <- as.numeric(regCoupling_restFC_F1ExecAcc_pvals$F1ExecAcc_gam_tstats)

## FDR correction
FDRcorr_regCoupling_restFC_F1ExecAcc_pvals <- p.adjust(regCoupling_restFC_F1ExecAcc_pvals$regCoupling_restFC_F1ExecAcc_pvals , method="fdr")
FDRcorr_regCoupling_restFC_F1ExecAcc_pvals <- cbind(regCoupling_restFC_F1ExecAcc_pvals ,FDRcorr_regCoupling_restFC_F1ExecAcc_pvals, nodeNames)

## Define subset of regions showing significant Age effects on Yeo PC
sig_FDRcorr_regCoupling_restFC_F1ExecAcc_pvals  <- FDRcorr_regCoupling_restFC_F1ExecAcc_pvals [which(FDRcorr_regCoupling_restFC_F1ExecAcc_pvals$FDRcorr_regCoupling_restFC_F1ExecAcc_pvals  <.05),]
dim(sig_FDRcorr_regCoupling_restFC_F1ExecAcc_pvals )

# Number of significant positive effects
sum(sig_FDRcorr_regCoupling_restFC_F1ExecAcc_pvals$F1ExecAcc_gam_tstats > 0)
# Number of significant negative effects
sum(sig_FDRcorr_regCoupling_restFC_F1ExecAcc_pvals$F1ExecAcc_gam_tstats < 0)

min(abs(sig_FDRcorr_regCoupling_restFC_F1ExecAcc_pvals$F1ExecAcc_gam_tstats))

sigNeg_results <- subset(sig_FDRcorr_regCoupling_restFC_F1ExecAcc_pvals, sig_FDRcorr_regCoupling_restFC_F1ExecAcc_pvals$F1ExecAcc_gam_tstats < 0)
sigPos_results <- subset(sig_FDRcorr_regCoupling_restFC_F1ExecAcc_pvals, sig_FDRcorr_regCoupling_restFC_F1ExecAcc_pvals$F1ExecAcc_gam_tstats > 0)

## Write out results
write.table(FDRcorr_regCoupling_restFC_F1ExecAcc_pvals$F1ExecAcc_gam_tstats, "/data/jux/BBL/projects/pncBaumStructFunc/replication/group_stats/n727_Schaefer400_threshNormProbSC_restFC_regCoupling_gam_F1Exec_tstat.txt", col.names = FALSE, row.names=FALSE)

###########################################################
## Nback Dprime Association with Regional Coupling (GAM) ##
###########################################################
# nbackBehTwobackDprime
# nbackBehAllDprime
covariates=" ~ nbackBehAllDprime + s(age, k=4) + dti64MeanRelRMS + nbackRelMeanRMSMotion + sex"    
m <- mclapply(names(cog_df[,3:402]), function(x) {as.formula(paste(x, covariates, sep=""))},mc.cores=2)
regCoupling_nbackFC_dprime_pvals <- mclapply(m, function(x) { summary(gam(formula = x, fx=TRUE,data=cog_df, REML=T))$p.pv[2]},mc.cores=2)
regCoupling_nbackFC_dprime_pvals <- as.data.frame(regCoupling_nbackFC_dprime_pvals)
regCoupling_nbackFC_dprime_pvals <- t(regCoupling_nbackFC_dprime_pvals )
regCoupling_nbackFC_dprime_pvals <- as.data.frame(regCoupling_nbackFC_dprime_pvals)
colnames(regCoupling_nbackFC_dprime_pvals ) <- "regCoupling_nbackFC_dprime_pvals "
regCoupling_nbackFC_dprime_pvals$Node_index <- 1:400

## Extract test statistic for effect of interest
regCoupling_nbackFC_dprime_pvals$dprime_gam_tstats <- mclapply(m, function(x) { summary(gam(formula = x, fx=TRUE,data=cog_df, REML=T))$p.t[2]},mc.cores=2)
regCoupling_nbackFC_dprime_pvals$dprime_gam_tstats <- as.numeric(regCoupling_nbackFC_dprime_pvals$dprime_gam_tstats)

## FDR correction
FDRcorr_regCoupling_nbackFC_dprime_pvals <- p.adjust(regCoupling_nbackFC_dprime_pvals$regCoupling_nbackFC_dprime_pvals , method="fdr")
FDRcorr_regCoupling_nbackFC_dprime_pvals <- cbind(regCoupling_nbackFC_dprime_pvals ,FDRcorr_regCoupling_nbackFC_dprime_pvals )

## Define subset of regions showing significant Age effects on Yeo PC
sig_FDRcorr_regCoupling_nbackFC_dprime_pvals  <- FDRcorr_regCoupling_nbackFC_dprime_pvals [which(FDRcorr_regCoupling_nbackFC_dprime_pvals$FDRcorr_regCoupling_nbackFC_dprime_pvals  <.05),]
dim(sig_FDRcorr_regCoupling_nbackFC_dprime_pvals )

# Number of significant positive effects
sum(sig_FDRcorr_regCoupling_nbackFC_dprime_pvals$dprime_gam_tstats > 0)
# Number of significant negative effects
sum(sig_FDRcorr_regCoupling_nbackFC_dprime_pvals$dprime_gam_tstats < 0)

min(abs(sig_FDRcorr_regCoupling_nbackFC_dprime_pvals$dprime_gam_tstats))

## Write out results
write.table(FDRcorr_regCoupling_nbackFC_dprime_pvals$dprime_gam_tstats,  "/data/jux/BBL/projects/pncBaumStructFunc/replication/group_stats/n727_Schaefer400_threshNormProbSC_nbackFC_regCoupling_nback_AllDprime_tstat.txt", col.names = FALSE, row.names=FALSE)

###########################
## Read in Node Features ##
###########################
mean_rest_coupling <- read.table("/data/jux/BBL/projects/pncBaumStructFunc/replication/network_metrics/n727_Schaefer400_thresh_norm_probSC_restFC_groupAvg_coupling_Spearman_r.txt", header=FALSE)
colnames(mean_rest_coupling) <- "mean_rest_coupling"

mean_nback_coupling <- read.table("/data/jux/BBL/projects/pncBaumStructFunc/replication/network_metrics/n727_Schaefer400_thresh_norm_probSC_nbackFC_groupAvg_coupling_Spearman_r.txt", header=FALSE)
colnames(mean_nback_coupling) <- "mean_nback_coupling"

mean_struct_nodeStrength <- read.table("/data/jux/BBL/projects/pncBaumStructFunc/replication/network_metrics/node_features/n727_Schaefer400_groupAvg_struct_thresh_norm_probSC_nodeStrength.txt", header=FALSE)
colnames(mean_struct_nodeStrength) <- "mean_struct_nodeStrength"

mean_nback_posStrength <- read.table("/data/jux/BBL/projects/pncBaumStructFunc/replication/network_metrics/node_features/n727_Schaefer400_mean_nback_pos_nodeStrength.txt", header=FALSE)
colnames(mean_nback_posStrength) <- "mean_nback_posStrength"

mean_nback_negStrength <- read.table("/data/jux/BBL/projects/pncBaumStructFunc/replication/network_metrics/node_features/n727_Schaefer400_mean_nback_neg_nodeStrength.txt", header=FALSE)
colnames(mean_nback_negStrength) <- "mean_nback_negStrength"

struct_PC <- read.table("/data/jux/BBL/projects/pncBaumStructFunc/replication/network_metrics/node_features/n727_thresh_norm_probSC_Schaefer400_Yeo7_mean_structPC.txt", header=FALSE)
colnames(struct_PC) <- "struct_PC"

nback_posPC <- read.table("/data/jux/BBL/projects/pncBaumStructFunc/replication/network_metrics/node_features/n727_Schaefer400_Yeo7_groupAvg_nbackFC_posPC.txt", header=FALSE)
colnames(nback_posPC) <- "nback_posPC"

nback_negPC <- read.table("/data/jux/BBL/projects/pncBaumStructFunc/replication/network_metrics/node_features/n727_Schaefer400_Yeo7_groupAvg_nbackFC_negPC.txt", header=FALSE)
colnames(nback_negPC) <- "nback_negPC"

rest_posPC <- read.table("/data/jux/BBL/projects/pncBaumStructFunc/replication/network_metrics/node_features/n727_Schaefer400_Yeo7_groupAvg_restFC_posPC.txt", header=FALSE)
colnames(rest_posPC) <- "rest_posPC"

rest_negPC <- read.table("/data/jux/BBL/projects/pncBaumStructFunc/replication/network_metrics/node_features/n727_Schaefer400_Yeo7_groupAvg_restFC_negPC.txt", header=FALSE)
colnames(rest_negPC) <- "rest_negPC"

Delta_posPC <- read.table("/data/jux/BBL/projects/pncBaumStructFunc/replication/network_metrics/node_features/n727_Schaefer400_Yeo7_mean_DELTA_nbackPC-restPC.txt", header=FALSE)
colnames(Delta_posPC) <- "Delta_posPC"

Delta_negPC <- read.table("/data/jux/BBL/projects/pncBaumStructFunc/replication/network_metrics/node_features/n727_Schaefer400_Yeo7_mean_DELTA_Neg_nbackPC-restPC.txt", header=FALSE)
colnames(Delta_negPC) <- "Delta_negPC"

mean_nback_activation <- read.table("/data/jux/BBL/projects/pncBaumStructFunc/network_measures/Schaefer400/nback_restFC_coupling/nback_activation/Schaefer400PNC_mean_n727_2back-0back_sigchange_regional_values.txt", header=FALSE)
colnames(mean_nback_activation) <- "mean_nback_activation"


corr_df <- cbind(FDRcorr_regCoupling_Age_pvals$AgeEffect_Zscore, FDRcorr_rest_regCoupling_Age_pvals$AgeEffect_Zscore, mean_rest_coupling, mean_nback_coupling, struct_PC, nback_posPC, nback_negPC, rest_posPC, rest_negPC, rescaled_evo_expansion, margulies_gradient, Delta_posPC, Delta_negPC, mean_nback_activation, restFC_pca_comp1, restFC_pca_comp2)
colnames(corr_df)[1] <- "nback_ageEffect"
colnames(corr_df)[2] <- "rest_ageEffect"

####################################################
## Create correlation matrix of regional measures ##
####################################################
require(corrplot)

## Create Right Hemisphere corr_df
rh_corr_df <- as.data.frame(corr_df[201:400,])

# M <- cor(corr_df)
# corrplot(M, method = "number", col="black")

# rh_M <- cor(rh_corr_df)
# corrplot(rh_M, method = "number", col="black")

## Yeo Partition
Yeo7_part <- read.table("/data/joy/BBL/applications/xcpEngine/atlas/schaefer400/schaefer400x7CommunityAffiliation.1D", header=FALSE)
Yeo17_part <- read.table("/data/joy/BBL/applications/xcpEngine/atlas/schaefer400/schaefer400x17CommunityAffiliation.1D", header=FALSE)

corr_df <- cbind(corr_df, Yeo7_part)
corr_df$dmn_idx <- 0
corr_df$dmn_idx[which(Yeo7_part == 7) ] <- 1
corr_df$Yeo7_part <- as.factor(Yeo7_part[,1])
# corr_df$Yeo17_part <- as.factor(Yeo17_part[,1])

## Create Right Hemisphere corr_df
rh_corr_df <- as.data.frame(corr_df[201:400,])

####################################################
## Permutation Testing for 3rd-level Correlations ##
####################################################
## Age Effect ~ Nback Pos PC
# nperm <- 10000
# perm_ageEffect_nbackPosPC_rvals <-  as.vector(rep(0, nperm))

# for(i in 1:nperm) {
#   perm_ageEffect <- permute(corr_df$nback_ageEffect)
#   perm_r <- rcorr(perm_ageEffect, corr_df$nback_posPC)
#   perm_ageEffect_nbackPosPC_rvals[i] <- perm_r$r[3]
# }

# orig_corr <-  rcorr(corr_df$nback_ageEffect, corr_df$nback_posPC)
# orig_r <- orig_corr$r[3]

# perm_ageEffect_nbackPosPC_pval <- (sum(abs(perm_ageEffect_nbackPosPC_rvals) > abs(orig_r))) / nperm

# hist(perm_ageEffect_scaling_rvals, col="royalblue") 


##########################################
## Figure 2: Associations with Coupling ##
##########################################
require(gridExtra)
axis_text <- element_text(color = "black", size = 14)

nbackCoupling_over_myelination_plot <- ggplot(data = corr_df, aes(mean_myelin, mean_nback_coupling)) + geom_point(aes(alpha=1.0))+ geom_smooth(method="lm",size=2, col="firebrick4", aes(y = mean_nback_coupling)) + theme(axis.text.x = axis_text, axis.text.y = axis_text ) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank()) + theme(axis.line = element_line(colour = 'black', size = 1.5), axis.ticks.length = unit(.25, "cm")) + scale_alpha(guide = 'none')

nbackCoupling_over_structPC_plot <- ggplot(data = corr_df, aes(struct_PC, mean_nback_coupling)) + geom_point(aes(alpha=1.0))+ geom_smooth(method="lm",size=2, col="firebrick4", aes(y = mean_nback_coupling)) + theme(axis.text.x = axis_text, axis.text.y = axis_text ) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank()) + theme(axis.line = element_line(colour = 'black', size = 1.5), axis.ticks.length = unit(.25, "cm")) + scale_alpha(guide = 'none')

nbackCoupling_over_nbackposPC_plot <- ggplot(data = corr_df, aes(nback_posPC, mean_nback_coupling)) + geom_point(aes(alpha=1.0))+ geom_smooth(method="lm",size=2, col="firebrick4", aes(y = mean_nback_coupling)) + theme(axis.text.x = axis_text, axis.text.y = axis_text ) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank()) + theme(axis.line = element_line(colour = 'black', size = 1.5), axis.ticks.length = unit(.25, "cm")) + scale_alpha(guide = 'none')

nbackCoupling_over_nbacknegPC_plot <- ggplot(data = corr_df, aes(nback_negPC, mean_nback_coupling)) + geom_point(aes(alpha=1.0))+ geom_smooth(method="lm",size=2, col="firebrick4", aes(y = mean_nback_coupling)) + theme(axis.text.x = axis_text, axis.text.y = axis_text ) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank()) + theme(axis.line = element_line(colour = 'black', size = 1.5), axis.ticks.length = unit(.25, "cm")) + scale_alpha(guide = 'none')

nbackCoupling_over_evoExpansion_plot <- ggplot(data = rh_corr_df, aes(rescaled_evo_expansion, mean_nback_coupling)) + geom_point(aes(alpha=1.0))+ geom_smooth(method="lm",size=2, col="firebrick4", aes(y = mean_nback_coupling)) + theme(axis.text.x = axis_text, axis.text.y = axis_text ) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank()) + theme(axis.line = element_line(colour = 'black', size = 1.5), axis.ticks.length = unit(.25, "cm")) + scale_alpha(guide = 'none')

nbackCoupling_over_gradient_plot <- ggplot(data = corr_df, aes(margulies_gradient, mean_nback_coupling)) + geom_point(aes(alpha=1.0)) + geom_smooth(method="lm",size=2, col="firebrick4", aes(y = mean_nback_coupling)) + theme(axis.text.x = axis_text, axis.text.y = axis_text ) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank()) + theme(axis.line = element_line(colour = 'black', size = 1.5), axis.ticks.length = unit(.25, "cm")) + scale_alpha(guide = 'none')

nbackCoupling_over_nbackSigChange <- ggplot(data = corr_df, aes(mean_nback_activation, mean_nback_coupling)) + geom_point(aes(alpha=1.0)) + geom_smooth(method="lm",size=2, col="firebrick4", aes(y = mean_nback_coupling)) + theme(axis.text.x = axis_text, axis.text.y = axis_text ) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank()) + theme(axis.line = element_line(colour = 'black', size = 1.5), axis.ticks.length = unit(.25, "cm")) + scale_alpha(guide = 'none')

nbackCoupling_over_restFC_PCA_gradient_plot <- ggplot(data = corr_df, aes(restFC_pca_comp1, mean_nback_coupling)) + geom_point(aes(alpha=1.0)) + geom_smooth(method="lm",size=2, col="firebrick4", aes(y = mean_nback_coupling)) + theme(axis.text.x = axis_text, axis.text.y = axis_text ) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank()) + theme(axis.line = element_line(colour = 'black', size = 1.5), axis.ticks.length = unit(.25, "cm")) + scale_alpha(guide = 'none')

###################################
## Scatterplots in Figure 2 B-E ##
###################################
grid.arrange(nbackCoupling_over_structPC_plot, nbackCoupling_over_nbackposPC_plot, nbackCoupling_over_gradient_plot, nbackCoupling_over_evoExpansion_plot, nrow = 2)

###################################################
## FIGURE 3B: NBACK COUPLING AGE EFFECTS BOXPLOT ##
###################################################
dev.off()

sig_min <- min(abs(sig_FDRcorr_regCoupling_Age_pvals$AgeEffect_Zscore))
neg_sig_min <- -(sig_min)

dodge <- position_dodge(width = 2)
modular_nbackCoupling_ageEffect_boxplot <- ggplot(data=corr_df, aes(x=factor(Yeo7_part), y=nback_ageEffect, fill=factor(Yeo7_part)))
modular_nbackCoupling_ageEffect_boxplot + geom_hline(yintercept= sig_min, linetype="dashed", color = "firebrick4", alpha=0.4) + geom_hline(yintercept= neg_sig_min, linetype="dashed", color = "red", alpha=0.4) + geom_boxplot(notch = FALSE, outlier.size = -1, color="black", lwd=0.7, alpha = 0.7, outlier.alpha = 0.1) +  geom_point(shape = 21,size=2, position = position_jitterdodge(), color="black", alpha=0.4) + theme(panel.grid.minor = element_blank(), panel.background = element_blank()) + scale_fill_manual(values = c("darkorchid4", "lightskyblue3", "forestgreen", "orchid3", "khaki2", "tan2", "indianred3")) + theme(axis.line = element_line(colour = 'black', size = 1.5), axis.ticks.length = unit(.25, "cm"))

##############################################
## ANOVA FOR AGE EFFECTS IN DMN (Figure 3B) ##
##############################################
corr_df$Yeo7_part <- as.factor(corr_df$Yeo7_part)
corr_df$relevel_Yeo7 <- relevel(corr_df$Yeo7_part, ref = 7)
test_anova <- lm(aov(data=corr_df, nback_ageEffect ~ as.factor(relevel_Yeo7)))
##############################################

##########################################################
## Figure 3C-E: Coupling Age effects over Node Features ##
##########################################################
nbackCouplingAgeEffect_over_nbackPosPC <- ggplot(data = corr_df, aes(nback_posPC, nback_ageEffect)) + geom_point(aes(alpha=abs(nback_ageEffect), col=factor(dmn_idx)))  + scale_colour_manual(values = c("lightskyblue3", "firebrick3")) + geom_smooth(method="lm",size=2, col="gray20", aes(y = nback_ageEffect)) + geom_hline(yintercept= sig_min, linetype="dashed", color = "gray36", alpha=0.6) + geom_hline(yintercept= neg_sig_min, linetype="dashed", color = "gray36", alpha=0.6) + theme(axis.text.x = axis_text, axis.text.y = axis_text ) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank()) + theme(axis.line = element_line(colour = 'black', size = 1.5), axis.ticks.length = unit(.25, "cm")) + scale_alpha(guide = 'none') 

nbackCouplingAgeEffect_over_nback_coupling <- ggplot(data = corr_df, aes(mean_nback_coupling, nback_ageEffect)) + geom_point(aes(alpha=abs(nback_ageEffect), col=factor(dmn_idx)))  + scale_colour_manual(values = c("lightskyblue3", "firebrick3")) + geom_smooth(method="lm",size=2, col="gray20", aes(y = nback_ageEffect)) + geom_hline(yintercept= sig_min, linetype="dashed", color = "gray36", alpha=0.6) + geom_hline(yintercept= neg_sig_min, linetype="dashed", color = "gray36", alpha=0.6) + theme(axis.text.x = axis_text, axis.text.y = axis_text ) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank()) + theme(axis.line = element_line(colour = 'black', size = 1.5), axis.ticks.length = unit(.25, "cm")) + scale_alpha(guide = 'none')

nbackCouplingAgeEffect_over_evoExpansion  <- ggplot(data = rh_corr_df, aes(rescaled_evo_expansion, nback_ageEffect)) + geom_point(aes(alpha=abs(nback_ageEffect), col=factor(dmn_idx)))  + scale_colour_manual(values = c("lightskyblue3", "firebrick3")) + geom_smooth(method="lm",size=2, col="gray20", aes(y = nback_ageEffect)) + geom_hline(yintercept= sig_min, linetype="dashed", color = "gray36", alpha=0.6) + geom_hline(yintercept= neg_sig_min, linetype="dashed", color = "gray36", alpha=0.6) + theme(axis.text.x = axis_text, axis.text.y = axis_text ) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank()) + theme(axis.line = element_line(colour = 'black', size = 1.5), axis.ticks.length = unit(.25, "cm")) + scale_alpha(guide = 'none')

nbackCouplingAgeEffect_over_gradient <- ggplot(data = corr_df, aes(margulies_gradient, nback_ageEffect)) + geom_point(aes(alpha=abs(nback_ageEffect), col=factor(dmn_idx)))  + scale_colour_manual(values = c("lightskyblue3", "firebrick3")) + geom_smooth(method="lm",size=2, col="gray20", aes(y = nback_ageEffect)) + geom_hline(yintercept= sig_min, linetype="dashed", color = "gray36", alpha=0.6) + geom_hline(yintercept= neg_sig_min, linetype="dashed", color = "gray36", alpha=0.6) + theme(axis.text.x = axis_text, axis.text.y = axis_text ) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank()) + theme(axis.line = element_line(colour = 'black', size = 1.5), axis.ticks.length = unit(.25, "cm")) + scale_alpha(guide = 'none')

nbackCouplingAgeEffect_over_restFC_PCA_gradient <- ggplot(data = corr_df, aes(restFC_pca_comp1, nback_ageEffect)) + geom_point(aes(alpha=abs(nback_ageEffect), col=factor(dmn_idx)))  + scale_colour_manual(values = c("lightskyblue3", "firebrick3")) + geom_smooth(method="lm",size=2, col="gray20", aes(y = nback_ageEffect)) + geom_hline(yintercept= sig_min, linetype="dashed", color = "gray36", alpha=0.6) + geom_hline(yintercept= neg_sig_min, linetype="dashed", color = "gray36", alpha=0.6) + theme(axis.text.x = axis_text, axis.text.y = axis_text ) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank()) + theme(axis.line = element_line(colour = 'black', size = 1.5), axis.ticks.length = unit(.25, "cm")) + scale_alpha(guide = 'none') # + theme(legend.position = "none")

grid.arrange(nbackCouplingAgeEffect_over_nbackPosPC, nbackCouplingAgeEffect_over_gradient, nbackCouplingAgeEffect_over_evoExpansion, nrow = 2)
#######################################################

###################################################
## Figure 5A: Coupling and Executive Performance ##
###################################################
covariates=" ~ F1_Exec_Comp_Res_Accuracy + s(age, k=4) + dti64MeanRelRMS + nbackRelMeanRMSMotion + sex"    
m <- mclapply(names(cog_df[,3:402]), function(x) {as.formula(paste(x, covariates, sep=""))},mc.cores=2)
regCoupling_nbackFC_F1ExecAcc_pvals <- mclapply(m, function(x) { summary(gam(formula = x, fx=TRUE,data=cog_df, REML=T))$p.pv[2]},mc.cores=2)
regCoupling_nbackFC_F1ExecAcc_pvals <- as.data.frame(regCoupling_nbackFC_F1ExecAcc_pvals)
regCoupling_nbackFC_F1ExecAcc_pvals <- t(regCoupling_nbackFC_F1ExecAcc_pvals )
regCoupling_nbackFC_F1ExecAcc_pvals <- as.data.frame(regCoupling_nbackFC_F1ExecAcc_pvals)
colnames(regCoupling_nbackFC_F1ExecAcc_pvals ) <- "regCoupling_nbackFC_F1ExecAcc_pvals "
regCoupling_nbackFC_F1ExecAcc_pvals$Node_index <- 1:400
regCoupling_nbackFC_F1ExecAcc_pvals$nodeNames <- nodeNames

## Extract test statistic for effect of interest
regCoupling_nbackFC_F1ExecAcc_pvals$F1ExecAcc_gam_tstats <- mclapply(m, function(x) { summary(gam(formula = x, fx=TRUE,data=cog_df, REML=T))$p.t[2]},mc.cores=2)
regCoupling_nbackFC_F1ExecAcc_pvals$F1ExecAcc_gam_tstats <- as.numeric(regCoupling_nbackFC_F1ExecAcc_pvals$F1ExecAcc_gam_tstats)

## FDR correction
FDRcorr_regCoupling_nbackFC_F1ExecAcc_pvals <- p.adjust(regCoupling_nbackFC_F1ExecAcc_pvals$regCoupling_nbackFC_F1ExecAcc_pvals , method="fdr")
FDRcorr_regCoupling_nbackFC_F1ExecAcc_pvals <- cbind(regCoupling_nbackFC_F1ExecAcc_pvals ,FDRcorr_regCoupling_nbackFC_F1ExecAcc_pvals, nodeNames)

## Define subset of regions showing significant Age effects on Yeo PC
sig_FDRcorr_regCoupling_nbackFC_F1ExecAcc_pvals  <- FDRcorr_regCoupling_nbackFC_F1ExecAcc_pvals [which(FDRcorr_regCoupling_nbackFC_F1ExecAcc_pvals$FDRcorr_regCoupling_nbackFC_F1ExecAcc_pvals  <.05),]
dim(sig_FDRcorr_regCoupling_nbackFC_F1ExecAcc_pvals )

# Number of significant positive effects
sum(sig_FDRcorr_regCoupling_nbackFC_F1ExecAcc_pvals$F1ExecAcc_gam_tstats > 0)
# Number of significant negative effects
sum(sig_FDRcorr_regCoupling_nbackFC_F1ExecAcc_pvals$F1ExecAcc_gam_tstats < 0)

min(abs(sig_FDRcorr_regCoupling_nbackFC_F1ExecAcc_pvals$F1ExecAcc_gam_tstats))

sigNeg_results <- subset(sig_FDRcorr_regCoupling_nbackFC_F1ExecAcc_pvals, sig_FDRcorr_regCoupling_nbackFC_F1ExecAcc_pvals$F1ExecAcc_gam_tstats < 0)
sigPos_results <- subset(sig_FDRcorr_regCoupling_nbackFC_F1ExecAcc_pvals, sig_FDRcorr_regCoupling_nbackFC_F1ExecAcc_pvals$F1ExecAcc_gam_tstats > 0)

## Write out results
write.table(FDRcorr_regCoupling_nbackFC_F1ExecAcc_pvals$F1ExecAcc_gam_tstats, "/data/jux/BBL/projects/pncBaumStructFunc/replication/group_stats/n727_Schaefer400_threshNormProbSC_nbackFC_regCoupling_gam_F1Exec_tstat.txt", col.names = FALSE, row.names=FALSE)

######################################################
## Figure 5B-C: Age/Cognitive Effect Plots for V348 ##
######################################################
V348_F1ExecAcc_gam <- gam(nback_regCoup_V348 ~ F1_Exec_Comp_Res_Accuracy + s(age, k=4) + dti64MeanRelRMS + nbackRelMeanRMSMotion + sex, data= cog_df, REML=TRUE)
V348_age_gam <- gam(nback_regCoup_V348 ~s(age, k=4) + dti64MeanRelRMS + nbackRelMeanRMSMotion + sex, data= cog_df, REML=TRUE)

#########################################
## Right vlPFC Nback Coupling Age Plot ##
#########################################
plotdata <- visreg(V348_age_gam,'age', type = "conditional", scale = "linear", plot = FALSE)
smooths <- data.frame(Variable = plotdata$meta$x, 
                      x=plotdata$fit[[plotdata$meta$x]], 
                      smooth=plotdata$fit$visregFit, 
                      lower=plotdata$fit$visregLwr, 
                      upper=plotdata$fit$visregUpr)
predicts <- data.frame(Variable = "dim1", 
                       x=plotdata$res$age,
                       y=plotdata$res$visregRes)

v348_age_plot <- ggplot() +
  geom_point(data = predicts, aes(x, y), colour = "lightblue3", alpha=0.4, size = 1.6 ) +
  geom_line(data = smooths, aes(x = x, y = smooth), colour = "darkslategrey",size=2) +
  geom_line(data = smooths, aes(x = x, y=lower), linetype="dashed", colour = "darkslategrey", alpha = 0.9, size = 0.9) + 
  geom_line(data = smooths, aes(x = x, y=upper), linetype="dashed",colour = "darkslategrey", alpha = 0.9, size = 0.9) +
  theme(legend.position = "none") +
  labs(x = "Age (years)", y = "Right vlPFC Coupling") +
  theme(axis.title.x = element_text(size = rel(1.6))) +
  theme(axis.title.y = element_text(size = rel(1.6))) + 
  theme(axis.text = element_text(size = rel(1.4))) + theme(axis.line = element_line(colour = 'black', size = 1.5), axis.ticks.length = unit(.25, "cm")) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank())


################################################
## Right vlPFC Nback Coupling F1_ExecAcc Plot ##
################################################
plotdata <- visreg(V348_F1ExecAcc_gam,'F1_Exec_Comp_Res_Accuracy', type = "conditional", scale = "linear", plot = FALSE)
smooths <- data.frame(Variable = plotdata$meta$x, 
                      x=plotdata$fit[[plotdata$meta$x]], 
                      smooth=plotdata$fit$visregFit, 
                      lower=plotdata$fit$visregLwr, 
                      upper=plotdata$fit$visregUpr)
predicts <- data.frame(Variable = "dim1", 
                       x=plotdata$res$F1_Exec_Comp_Res_Accuracy,
                       y=plotdata$res$visregRes)

v348_F1ExecAcc_plot <- ggplot() +
  geom_point(data = predicts, aes(x, y), colour = "lightblue3", alpha=0.4, size = 1.6 ) +
  geom_line(data = smooths, aes(x = x, y = smooth), colour = "darkslategrey",size=2) +
  geom_line(data = smooths, aes(x = x, y=lower), linetype="dashed", colour = "darkslategrey", alpha = 0.9, size = 0.9) + 
  geom_line(data = smooths, aes(x = x, y=upper), linetype="dashed",colour = "darkslategrey", alpha = 0.9, size = 0.9) +
  theme(legend.position = "none") +
  labs(x = "Executive Performance", y = "Right vlPFC Coupling") +
  theme(axis.title.x = element_text(size = rel(1.6))) +
  theme(axis.title.y = element_text(size = rel(1.6))) + 
  theme(axis.text = element_text(size = rel(1.4))) + theme(axis.line = element_line(colour = 'black', size = 1.5), axis.ticks.length = unit(.25, "cm")) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank())

v348_F1ExecAcc_plot


grid.arrange(v348_age_plot, v348_F1ExecAcc_plot, nrow = 1)

#################################
## COGNTIVE MEDIATION ANALYSIS ##
#################################
Age_lm <- lm(age ~ sex + dti64MeanRelRMS + nbackRelMeanRMSMotion, data=cog_df)
Age_resid <- resid(Age_lm)

ExecEff_lm <- lm(F1_Exec_Comp_Res_Accuracy ~ sex + dti64MeanRelRMS + nbackRelMeanRMSMotion, data=cog_df)
ExecEff_resid <- resid(ExecEff_lm)

coupling_lm <- lm(nback_regCoup_V348  ~ sex + dti64MeanRelRMS + nbackRelMeanRMSMotion, data=cog_df)
coupling_resid <- resid(coupling_lm)

X <- as.data.frame(scale(Age_resid))
Y <- as.data.frame(scale(ExecEff_resid)) 
M <- as.data.frame(scale(coupling_resid ))

Data <- data.frame(X=X, Y=Y, M=M)
Data <- data.frame(cbind(X,Y,M))
colnames(Data) <- c("X", "Y", "M")

model <- ' # direct effect
Y ~ c*X
# mediator
M ~ a*X
Y ~ b*M
# indirect effect (a*b)
ab := a*b
# total effect
total := c + (a*b)
'
fit_sem <- sem(model, data = Data, se="bootstrap", bootstrap=10000)
summary(fit_sem, fit.measures=TRUE, standardize=TRUE, rsquare=TRUE)

## Calculate bootstrapped confidence intervals for the indirect (c') effect
boot.fit <- parameterEstimates(fit_sem, boot.ci.type="perc",level=0.95, ci=TRUE,standardized = TRUE)
boot.fit

## Extract boot-strapped mediation effect from results:
mediation_pval <- boot.fit$pvalue[7]
print(mediation_pval)

save.image("/data/jux/BBL/projects/pncBaumStructFunc/replication/group_stats/n727_coupling_groupAnalysis_workspace.Rdata")
