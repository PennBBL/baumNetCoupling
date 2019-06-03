rm(list = ls())

#############################
## Load Relevant Libraries ##
#############################
require(lme4)
require(nlme)
require(gamm4)
require(mgcv)
require(ggplot2)
require(Hmisc)
require(MASS)
require(ppcor)
require(stringr)
require(visreg)
require(parallel)
require(multilevel)
require(stats)
require(Formula)
require(lavaan)
require(car)
require(corrplot)
require(stringr)
require(R.matlab)
require(pracma)
require(wesanderson)


###################################
## LOAD GROUP ANALYSIS WORKSPACE ##
###################################
load("/data/jux/BBL/projects/pncBaumStructFunc/replication/demographics/n294_TD_longitudinal_schaefer400_coupling_demographics.RData")

## Schaefer400 NodeNames
nodeNames <- read.table("/data/joy/BBL/applications/xcpEngine/atlas/schaefer400/schaefer400NodeNames.txt", header=FALSE)

## Yeo Partition
Yeo7_part <- read.table("/data/joy/BBL/applications/xcpEngine/atlas/schaefer400/schaefer400x7CommunityAffiliation.1D", header=FALSE)
Yeo17_part <- read.table("/data/joy/BBL/applications/xcpEngine/atlas/schaefer400/schaefer400x17CommunityAffiliation.1D", header=FALSE)


##################################
## Structural and Functional PC ## 
##################################

mean_structural_PC <- read.table("/data/jux/BBL/projects/pncBaumStructFunc/replication/network_metrics/longitudinal/n294_Schaefer400_threshNormProbSC_nbackFC_Yeo7_regional_structPC.txt", header=FALSE)
colnames(mean_structural_PC) <- "mean_structural_PC"

mean_nback_posPC <- read.table("/data/jux/BBL/projects/pncBaumStructFunc/replication/network_metrics/longitudinal/n294_Schaefer400_Yeo7_groupAvg_nbackFC_posPC.txt", header=FALSE)
colnames(mean_nback_posPC) <- "mean_nback_posPC"

#########################################################################
## Allometric Scaling, Evolutionary Expansion, and PET Metabolism maps ##
#########################################################################
## Myelin Map
mean_myelin <- read.table("/data/jux/BBL/projects/pncBaumStructFunc/myelin_maps/output/hayashi_mean_Myelin_schaefer400.txt", header=FALSE)
colnames(mean_myelin) <- "mean_myelin"

## Allometric Scaling, Evolutionary Expansion, and PET Metabolism maps
mean_CBF <- read.table("/data/jux/BBL/projects/pncBaumStructFunc/cortical_scaling/n1172_AslMean_regional_schaefer400x17.txt", header=FALSE)
colnames(mean_CBF) <- "mean_CBF"

cort_scaling <- read.table("/data/jux/BBL/projects/pncBaumStructFunc/cortical_scaling/schaefer400_regional_cortical_scaling.txt", header=FALSE)
colnames(cort_scaling) <- "cort_scaling"

pet_metabolism <- read.table("/data/jux/BBL/projects/pncBaumStructFunc/cortical_scaling/schaefer400_regional_cortical_scaling.txtschaefer400_regional_PET_metabolism.txt", header=FALSE)
colnames(pet_metabolism) <- "pet_metabolism"

margulies_gradient <- read.table("/data/jux/BBL/projects/pncBaumStructFunc/cortical_scaling/schaefer400x17_mean_regional_margulies_gradient.txt", header=FALSE)
# read.table("/Users/Graham/Documents/projects/pncBaumStructFunc/figures/regional_coupling/n819_nbackCoupling/schaefer400x17_regional_margulies_gradient.txt", header=FALSE)
colnames(margulies_gradient) <- "margulies_gradient"

rescaled_evo_expansion <- read.table("/data/jux/BBL/projects/pncBaumStructFunc/cortical_scaling/schaefer400_rescaled_regional_evoExpansion.txt", header=FALSE)
colnames(rescaled_evo_expansion) <- "rescaled_evo_expansion"

evo_expansion <- read.table("/data/jux/BBL/projects/pncBaumStructFunc/cortical_scaling/schaefer400_regional_evoExpansion.txt", header=FALSE)
colnames(evo_expansion) <- "evo_expansion"

dev_expansion <- read.table("/data/jux/BBL/projects/pncBaumStructFunc/cortical_scaling/schaefer400_regional_developmentalExpansion.txt", header=FALSE)
colnames(dev_expansion) <- "dev_expansion"

restFC_pca_comp1 <- read.table("/data/jux/BBL/projects/pncBaumStructFunc/replication/network_metrics/node_features/n727_Schaefer400_Zscore_mean_restFC_PCA_comp1_loadings.txt", header=FALSE)
colnames(restFC_pca_comp1) <- "restFC_pca_comp1"

restFC_pca_comp2 <- read.table("/data/jux/BBL/projects/pncBaumStructFunc/replication/network_metrics/node_features/n727_Schaefer400_Zscore_mean_restFC_PCA_comp2_loadings.txt", header=FALSE)
colnames(restFC_pca_comp2) <- "restFC_pca_comp2"

mean_nback_activation <- read.table("/data/jux/BBL/projects/pncBaumStructFunc/network_measures/Schaefer400/nback_restFC_coupling/nback_activation/Schaefer400PNC_mean_n727_2back-0back_sigchange_regional_values.txt", header=FALSE)
colnames(mean_nback_activation) <- "mean_nback_activation"

##########################################
## Define Transmodal and Unimodal Index ##
##########################################
transmodal_index <- which(margulies_gradient$margulies_gradient > quantile(margulies_gradient$margulies_gradient, 0.9))
unimodal_index <- which(margulies_gradient$margulies_gradient < quantile(margulies_gradient$margulies_gradient, 0.1))

struct_hub_index <- which(mean_structural_PC$mean_structural_PC > quantile(mean_structural_PC$mean_structural_PC, 0.9))
func_hub_index <- which(mean_nback_posPC$mean_nback_posPC > quantile(mean_nback_posPC$mean_nback_posPC, 0.9))

#################################################
## Read in Network Metrics generated in Matlab ## 
#################################################

## Regional Coupling
regionalCoupling_schaefer400 <- read.csv("/data/jux/BBL/projects/pncBaumStructFunc/replication/network_metrics/longitudinal/n294_TD_long_Schaefer400_threshNormProbSC_nbackFC_regional_coupling_Spearman_r.txt", header=FALSE)
for(i in 1:400) {
    colnames(regionalCoupling_schaefer400)[i] <- paste("regCoup_V", i, sep = "")
}

## Structural Node Strength
struct_schaefer400_nodeStrength <- read.csv("/data/jux/BBL/projects/pncBaumStructFunc/replication/network_metrics/longitudinal/n294_Schaefer400_Yeo7_regional_threshNormProbSC_strength.txt", header=FALSE)
for(i in 1:400) {
    colnames(struct_schaefer400_nodeStrength)[i] <- paste("struct_V", i, sep = "")
}

## Structural PC (Yeo-7)
struct_schaefer400_PC <- read.csv("/data/jux/BBL/projects/pncBaumStructFunc/replication/network_metrics/longitudinal/n294_Schaefer400_threshNormProbSC_nbackFC_Yeo7_regional_structPC.txt", header=FALSE)
for(i in 1:400) {
    colnames(struct_schaefer400_PC)[i] <- paste("structPC_V", i, sep = "")
}

## Nback Pos Node Strength
func_schaefer400_pos_nodeStrength <- read.csv("/data/jux/BBL/projects/pncBaumStructFunc/replication/network_metrics/longitudinal/n294_Schaefer400_Yeo7_regional_nbackFC_pos_strength.txt", header=FALSE)
for(i in 1:400) {
    colnames(func_schaefer400_pos_nodeStrength)[i] <- paste("func_posV", i , sep = "")
}

## Nback Neg Node Strength
func_schaefer400_neg_nodeStrength <- read.csv("/data/jux/BBL/projects/pncBaumStructFunc/replication/network_metrics/longitudinal/n294_Schaefer400_Yeo7_regional_nbackFC_neg_strength.txt", header=FALSE)
for(i in 1:400) {
    colnames(func_schaefer400_neg_nodeStrength)[i] <- paste("func_negV", i , sep = "")
}

## Nback Pos PC
func_schaefer400_posPC <- read.csv("/data/jux/BBL/projects/pncBaumStructFunc/replication/network_metrics/longitudinal/n294_Schaefer400_Yeo7_regional_nbackFC_posPC.txt", header=FALSE)
for(i in 1:400) {
    colnames(func_schaefer400_posPC)[i] <- paste("func_posPC_V", i, sep = "")
}

## Nback Neg PC
func_schaefer400_negPC <- read.csv("/data/jux/BBL/projects/pncBaumStructFunc/replication/network_metrics/longitudinal/n294_Schaefer400_Yeo7_regional_nbackFC_negPC.txt", header=FALSE)
for(i in 1:400) {
    colnames(func_schaefer400_negPC)[i] <- paste("func_negPC_V", i, sep = "")
}

#############################################
## Modularity Quality (Q) of Yeo Partition ##
#############################################
struct_Yeo7_Q <- read.table("/data/jux/BBL/projects/pncBaumStructFunc/replication/network_metrics/longitudinal/n294_Schaefer400_threshNormProbSC_Yeo7_struct_Q.txt", header=FALSE)
struct_Yeo7_Q <- as.data.frame(struct_Yeo7_Q)
colnames(struct_Yeo7_Q) <- "struct_Yeo7_Q"

func_Yeo7_Q <- read.table("/data/jux/BBL/projects/pncBaumStructFunc/replication/network_metrics/longitudinal/n294_Schaefer400_nbackFC_Yeo7_Q", header=FALSE)
func_Yeo7_Q <- as.data.frame(func_Yeo7_Q)
colnames(func_Yeo7_Q) <- "func_Yeo7_Q"

###############################################################
## Merge longitudinal data-frame with brain network measures ##
###############################################################
dim(dti_long)

dti_long <- cbind(regionalCoupling_schaefer400, struct_schaefer400_PC, func_schaefer400_posPC, func_schaefer400_negPC, struct_schaefer400_nodeStrength, func_schaefer400_pos_nodeStrength, func_schaefer400_neg_nodeStrength, dti_long, struct_Yeo7_Q, func_Yeo7_Q)

dti_long$sex <- as.factor(dti_long$sex)

dim(dti_long)
###############################################################


#########################################################
## Apply LME across regions for Age effect on coupling ##
#########################################################
coupling_results <- as.data.frame(matrix(NA, nrow=400, ncol=3))
colnames(coupling_results) <- c("region_index", "age_pval", "age_tstat")
coupling_results$region_index <- 1:400
n <-1

for(i in 1:400){

	## Run Linear Mixed Effects Model for each region
	tmp_lm <- lme(as.formula(paste("regCoup_V", i, "~ age_in_years + nbackRelMeanRMSMotion + dti64MeanRelRMS + sex", sep='')), data= dti_long, random= ~ 1|bblid)

	## Output results
	coupling_results[n,2] <- summary(tmp_lm)$tTable[2,5]
	coupling_results[n,3] <- summary(tmp_lm)$tTable[2,4]

	n <- n + 1
}

hist(coupling_results$age_tstat, col="royalblue3")

## FDR Correction
FDRcorr_age_pval <- as.data.frame(p.adjust(coupling_results$age_pval, method="fdr"))
colnames(FDRcorr_age_pval) <- "FDRcorr_age_pval"

## Merge region-level results with mean baseline coupling/FC/SC
coupling_results <- cbind(coupling_results, FDRcorr_age_pval, mean_structural_PC, mean_nback_posPC, nodeNames)
sig_FDRcorr_age_pval <- as.data.frame(subset(coupling_results, coupling_results$FDRcorr_age_pval < 0.05))
dim(sig_FDRcorr_age_pval)

sum(sig_FDRcorr_age_pval$age_tstat > 0)
sum(sig_FDRcorr_age_pval$age_tstat < 0)

min(abs(sig_FDRcorr_age_pval$age_tstat))

###########################################
## FIGURE 4: LME Age Effects on Coupling ##
###########################################

## Results displayed on second-row brains in Fig 4A
write.table(coupling_results$age_tstat, "/data/jux/BBL/projects/pncBaumStructFunc/replication/group_stats/longitudinal/n294_Schaefer400_threshNormProbSC_nbackFC_Coupling_lme_Age_effects.txt", col.names = FALSE, row.names=FALSE)

## Cross-sect age effects: ##
cs_nback_coupling_age_effects <- read.table("/data/jux/BBL/projects/pncBaumStructFunc/replication/group_stats/n727_Schaefer400_threshNormProbSC_nbackFC_Spearman_regCoupling_gam_Age_Zscore.txt", header=FALSE)
cs_nback_coupling_age_effects
coupling_results$cs_nback_coupling_age_effects <- as.numeric(cs_nback_coupling_age_effects$V1)

##################################################################################
## Correlation between Cross-sectional and Longitudinal Age Effects on Coupling ##
##################################################################################
rcorr(coupling_results$age_tstat, coupling_results$cs_nback_coupling_age_effects)

##################################################################
## Separate Timepoints For Longitudinal Delta Coupling Analysis ##
##################################################################
t1_df <- subset(dti_long, dti_long$timepoint==1)
t2_df <- subset(dti_long, dti_long$timepoint==2)

t2_regCoupling_df <- as.data.frame(t2_df[,1:400])
for(i in 1:400) {
    colnames(t2_regCoupling_df)[i] <- paste("t2", colnames(t2_regCoupling_df)[i] , sep = "_")
}

t2_funcPC_df <- as.data.frame(t2_df [,801:1200])
for(i in 1:400) {
    colnames(t2_funcPC_df)[i] <- paste("t2", colnames(t2_funcPC_df)[i] , sep = "_")
}

###############
## Delta Age ##
###############
t1_df$delta_age <- t2_df$age_in_years - t1_df$age_in_years

###################################
## Mean Motion Across Timepoints ##
###################################
t1_df$mean_struct_motion <- (t2_df$dti64MeanRelRMS + t1_df$dti64MeanRelRMS) / 2
t1_df$mean_func_motion <- (t2_df$nbackRelMeanRMSMotion + t1_df$nbackRelMeanRMSMotion) / 2

## Add T2 Coupling values to t1_df ##
t1_df <- cbind(t1_df, t2_regCoupling_df, t2_funcPC_df)

## Define number of subjects and brain regions
nsub <- dim(t1_df)[1]
nreg <- 400

##################################################
## Allocate Empty Matrices for Coupling Output  ##
##################################################
delta_SC <- as.data.frame(matrix(0, nrow=nsub , ncol=nreg))
delta_posFC <- as.data.frame(matrix(0, nrow=nsub , ncol=nreg))
delta_negFC <- as.data.frame(matrix(0, nrow=nsub , ncol=nreg))
delta_regCoup <- as.data.frame(matrix(0, nrow=nsub , ncol=nreg))
delta_structPC <- as.data.frame(matrix(0, nrow=nsub , ncol=nreg))
delta_func_posPC <- as.data.frame(matrix(0, nrow=nsub , ncol=nreg))
delta_func_negPC <- as.data.frame(matrix(0, nrow=nsub , ncol=nreg))
ARC_coupling <- as.data.frame(matrix(0, nrow=nsub , ncol=nreg))
####################################################################

####################
## Delta Coupling ##
####################
for(i in 1:400) {
	t2_coup <- paste("regCoup_V", i, sep = "")
	t1_coup <- paste("regCoup_V", i, sep = "")
	delta_regCoup[,i] <- (t2_df[[t2_coup]] - t1_df[[t1_coup]]) # / (t1_df[[t1_coup]])
 	colnames(delta_regCoup)[i] <- paste("delta_regCoup_V", i , sep = "")
}

mean_delta_coupling <- as.data.frame(lapply(delta_regCoup[,1:400], mean))
mean_delta_coupling <- t(mean_delta_coupling)
colnames(mean_delta_coupling) <- "mean_delta_coupling"

hist(mean_delta_coupling[,1])

sd_delta_coupling <- as.data.frame(lapply(delta_regCoup[,1:400], sd))
sd_delta_coupling <- t(sd_delta_coupling)
colnames(sd_delta_coupling) <- "sd_delta_coupling"

## Export mean Delta Coupling for brain rendering
write.table(mean_delta_coupling[,1], "/data/jux/BBL/projects/pncBaumStructFunc/replication/network_metrics/longitudinal/n294_mean_raw_delta_coupling_nbackFC_Schaefer400.txt", col.names=FALSE, row.names=FALSE)
write.table(sd_delta_coupling[,1], "/data/jux/BBL/projects/pncBaumStructFunc/replication/network_metrics/longitudinal/n294_SD_raw_delta_coupling_nbackFC_Schaefer400.txt", col.names=FALSE, row.names=FALSE)

##################
## ARC Coupling ##
##################
for(i in 1:400) {
	delta_coup <- paste("delta_regCoup_V", i, sep = "")
	t1_coup <- paste("regCoup_V", i, sep = "")
	ARC_coupling[,i] <-  (delta_regCoup[[delta_coup]] / t1_df[[t1_coup]]) / t1_df$delta_age
 	colnames(ARC_coupling)[i] <- paste("ARC_coupling_V", i , sep = "")
}

#######################
## Mean ARC Coupling ##
#######################
mean_ARC_coupling <- as.data.frame(lapply(ARC_coupling[,1:400], mean))
mean_ARC_coupling <- t(mean_ARC_coupling)
colnames(mean_ARC_coupling) <- "mean_ARC_coupling"

sd_ARC_coupling <- as.data.frame(lapply(ARC_coupling[,1:400], sd))
sd_ARC_coupling <- t(sd_ARC_coupling)
colnames(sd_ARC_coupling) <- "sd_ARC_coupling"

# Export mean ARC Coupling for brain rendering
write.table(mean_ARC_coupling[,1], "/data/jux/BBL/projects/pncBaumStructFunc/replication/network_metrics/longitudinal/n294_mean_ARC_threshNormProbSC_nbackFC_coupling.txt", col.names=FALSE, row.names=FALSE)
write.table(sd_ARC_coupling[,1], "/data/jux/BBL/projects/pncBaumStructFunc/replication/network_metrics/longitudinal/n294_SD_ARC_threshNormProbSC_nbackFC.txt", col.names=FALSE, row.names=FALSE)

############################
## Delta SC Node Strength ##
############################
for(i in 1:400) {
	t2_struct <- paste("struct_V", i, sep = "")
	t1_struct <- paste("struct_V", i, sep = "")

	# delta_formula < paste("t2_df$struct_V", i, "- t1_df$struct_V", i,  sep = "")
	delta_SC[,i] <- (t2_df[[t2_struct]] - t1_df[[t1_struct]]) # / (t1_df[[t1_struct]])

	# delta_SC[,i] <- as.formula(paste("t2_df$struct_V", i, " - ","t1_df$struct_V", i,  sep = ""))
    colnames(delta_SC)[i] <- paste("delta_struct_V", i , sep = "")
}

################################
## Delta FC Pos Node Strength ##
################################
for(i in 1:400) {
	t2_func <- paste("func_posV", i, sep = "")
	t1_func <- paste("func_posV", i, sep = "")
	delta_posFC[,i] <- (t2_df[[t2_func]] - t1_df[[t1_func]]) # / ( t1_df[[t1_func]])
	colnames(delta_posFC)[i] <- paste("delta_func_posV", i, sep = "")
}

################################
## Delta FC Neg Node Strength ##
################################
for(i in 1:400) {
	t2_func <- paste("func_negV", i, sep = "")
	t1_func <- paste("func_negV", i, sep = "")
	delta_negFC[,i] <- (t2_df[[t2_func]] - t1_df[[t1_func]]) # / (t1_df[[t1_func]])
	colnames(delta_negFC)[i] <- paste("delta_func_negV", i, sep = "")
}

#####################
## Delta Struct PC ##
#####################
for(i in 1:400) {
	t2_structPC <- paste("structPC_V", i, sep = "")
	t1_structPC <- paste("structPC_V", i, sep = "")
	delta_structPC[,i] <- (t2_df[[t2_structPC]] - t1_df[[t1_structPC]]) # / (t1_df[[t1_structPC]])
 	colnames(delta_structPC)[i] <- paste("delta_structPC_V", i , sep = "")
}

#######################
## Delta Func Pos PC ##
#######################
for(i in 1:400) {
	t2_funcPC <- paste("func_posPC_V", i, sep = "")
	t1_funcPC <- paste("func_posPC_V", i, sep = "")
	delta_func_posPC[,i] <- (t2_df[[t2_funcPC]] - t1_df[[t1_funcPC]]) # / (t1_df[[t1_funcPC]])
 	colnames(delta_func_posPC)[i] <- paste("delta_func_posPC_V", i , sep = "")
}

#######################
## Delta Func Neg PC ##
#######################
for(i in 1:400) {
	t2_funcPC <- paste("func_negPC_V", i, sep = "")
	t1_funcPC <- paste("func_negPC_V", i, sep = "")
	delta_func_negPC[,i] <- (t2_df[[t2_funcPC]] - t1_df[[t1_funcPC]]) # / t1_df[[t1_funcPC]]
 	colnames(delta_func_negPC)[i] <- paste("delta_func_negPC_V", i , sep = "")
}


###########################################################
## Merge Delta Brain Measures with Timepoint 1 dataframe ##
###########################################################
t1_df <- cbind(delta_regCoup, delta_SC, delta_structPC, delta_posFC, delta_negFC, delta_func_posPC, delta_func_negPC, t1_df, t2_regCoupling_df, ARC_coupling)
dim(t1_df)


test_lm <- lm(delta_regCoup_V348 ~ func_posV348 + regCoup_V138 + mean_func_motion + mean_struct_motion + age_in_years + delta_age + sex, data=t1_df)
test_lm <- lm(delta_regCoup_V138 ~ func_posV138 + regCoup_V138 + mean_func_motion + mean_struct_motion + age_in_years + delta_age + sex, data=t1_df)
test_lm <- lm(delta_regCoup_V138 ~ struct_V138 + func_posV138  + mean_func_motion + mean_struct_motion + age_in_years + delta_age + sex, data=t1_df)

#############################################
## Regress T1-Coupling from Delta Coupling ##
#############################################
deltaCoupling_resid <- as.data.frame(matrix(NA, nrow=nsub, ncol=nreg))

for(i in 1:400){

	colnames(deltaCoupling_resid)[i] <- paste("deltaCoupling_resid_V", i , sep = "")

	## Run Linear Regression Model
	tmp_lm <- lm(as.formula(paste("delta_regCoup_V", i, "~", "regCoup_V", i, sep='')), data=t1_df)
	tmp_resid <- stdres(tmp_lm)
	deltaCoupling_resid[,i] <- tmp_resid
}

t1_df <- cbind(t1_df, deltaCoupling_resid)

mean_resid_delta_coupling <- as.data.frame(lapply(deltaCoupling_resid[,1:400], mean))
mean_resid_delta_coupling <- t(mean_resid_delta_coupling)
colnames(mean_resid_delta_coupling) <- "mean_resid_delta_coupling"

# Export mean Delta Coupling for brain rendering
write.table(mean_resid_delta_coupling[,1], "/data/jux/BBL/projects/pncBaumStructFunc/replication/network_metrics/longitudinal/n294_Schaefer400_groupAvg_scaled_resid_delta_nback_coupling.txt", col.names=FALSE, row.names=FALSE)

########################################
## Regress T1-PC from Delta Struct PC ##
########################################
delta_structPC_resid <- as.data.frame(matrix(NA, nrow=nsub, ncol=nreg))

for(i in 1:400){

	colnames(delta_structPC_resid)[i] <- paste("delta_structPC_resid_V", i , sep = "")

	## Run Linear Regression Model
	tmp_lm <- lm(as.formula(paste("delta_structPC_V", i, "~", "structPC_V", i, sep='')), data=t1_df)
	# tmp_resid <- resid(tmp_lm)
	tmp_resid <- stdres(tmp_lm)
	delta_structPC_resid[,i] <- tmp_resid
}

t1_df <- cbind(t1_df, delta_structPC_resid)


##########################################
## Regress T1-PC from Delta Func Pos PC ##
##########################################
delta_funcPosPC_resid <- as.data.frame(matrix(NA, nrow=nsub, ncol=nreg))

for(i in 1:400){

	colnames(delta_funcPosPC_resid)[i] <- paste("delta_funcPosPC_resid_V", i , sep = "")

	## Run Linear Regression Model
	tmp_lm <- lm(as.formula(paste("delta_func_posPC_V", i, "~", "func_posPC_V", i, sep='')), data=t1_df)
	# tmp_resid <- resid(tmp_lm)
	tmp_resid <- stdres(tmp_lm)
	delta_funcPosPC_resid[,i] <- tmp_resid
}

t1_df <- cbind(t1_df, delta_funcPosPC_resid)

##########################################
## Regress T1-PC from Delta Func Neg PC ##
##########################################
delta_funcNegPC_resid <- as.data.frame(matrix(NA, nrow=nsub, ncol=nreg))

for(i in 1:400){

	colnames(delta_funcNegPC_resid)[i] <- paste("delta_funcNegPC_resid_V", i , sep = "")

	## Run Linear Regression Model
	tmp_lm <- lm(as.formula(paste("delta_func_negPC_V", i, "~", "func_negPC_V", i, sep='')), data=t1_df)
	# tmp_resid <- resid(tmp_lm)
	tmp_resid <- stdres(tmp_lm)
	delta_funcNegPC_resid[,i] <- tmp_resid
}

t1_df <- cbind(t1_df, delta_funcNegPC_resid)


######################################################
## Regress T1-Strength from Delta Func Pos Strength ##
######################################################
delta_funcPosStrength_resid <- as.data.frame(matrix(NA, nrow=nsub, ncol=nreg))

for(i in 1:400){

	colnames(delta_funcPosStrength_resid)[i] <- paste("delta_funcPosStrength_resid_V", i , sep = "")

	## Run Linear Regression Model
	tmp_lm <- lm(as.formula(paste("delta_func_posV", i, "~", "func_posV", i, sep='')), data=t1_df)
	# tmp_resid <- resid(tmp_lm)
	tmp_resid <- stdres(tmp_lm)
	delta_funcPosStrength_resid[,i] <- tmp_resid
}

t1_df <- cbind(t1_df, delta_funcPosStrength_resid)

######################################################
## Regress T1-Strength from Delta Func Neg Strength ##
######################################################
delta_funcNegStrength_resid <- as.data.frame(matrix(NA, nrow=nsub, ncol=nreg))

for(i in 1:400){

	colnames(delta_funcNegStrength_resid)[i] <- paste("delta_funcNegStrength_resid_V", i , sep = "")

	## Run Linear Regression Model
	tmp_lm <- lm(as.formula(paste("delta_func_negV", i, "~", "func_negV", i, sep='')), data=t1_df)
	# tmp_resid <- resid(tmp_lm)
	tmp_resid <- stdres(tmp_lm)
	delta_funcNegStrength_resid[,i] <- tmp_resid
}

t1_df <- cbind(t1_df, delta_funcNegStrength_resid)

######################################################
## Regress T1-Strength from Delta Struct Strength ##
######################################################
delta_structStrength_resid <- as.data.frame(matrix(NA, nrow=nsub, ncol=nreg))

for(i in 1:400){

	colnames(delta_structStrength_resid)[i] <- paste("delta_structStrength_resid_V", i , sep = "")

	## Run Linear Regression Model
	tmp_lm <- lm(as.formula(paste("delta_struct_V", i, "~", "struct_V", i, sep='')), data=t1_df)
	# tmp_resid <- resid(tmp_lm)
	tmp_resid <- stdres(tmp_lm)
	delta_structStrength_resid[,i] <- tmp_resid
}

t1_df <- cbind(t1_df,delta_structStrength_resid)

## Yeo Partition
Yeo7_part <- read.table("/data/joy/BBL/applications/xcpEngine/atlas/schaefer400/schaefer400x7CommunityAffiliation.1D", header=FALSE)
Yeo17_part <- read.table("/data/joy/BBL/applications/xcpEngine/atlas/schaefer400/schaefer400x17CommunityAffiliation.1D", header=FALSE)

####################################
## TEST LONGITUDINAL DELTA MODELS ##
##								  ##
##								  ##
## 		RESULTS IN FIG 4B		  ##
####################################

#############################################
## Delta Func PC Predicting Delta Coupling ##
#############################################
deltaCoupling_results <- as.data.frame(matrix(NA, nrow=400, ncol=3))
colnames(deltaCoupling_results) <- c("region_index", "pvals", "tstats")
deltaCoupling_results$region_index <- 1:400
n <-1
  
for(i in 1:400){
  
  ## Run Linear Regression Model
  tmp_lm <- lm(as.formula(paste("delta_regCoup_V", i, "~", "delta_func_posPC_V", i, "+ mean_func_motion + mean_struct_motion + age_in_years + delta_age + sex", sep='')), data=t1_df)
  
  ## Output results
  deltaCoupling_results[n,2] <- summary(tmp_lm)$coeff[2,4]
  deltaCoupling_results[n,3] <- summary(tmp_lm)$coeff[2,3]

  n <- n + 1
}

## Look at distribution of effects
hist(as.numeric(deltaCoupling_results$tstats), col="royalblue3")

## FDR Correction
FDRcorr_coupling_pval <- as.data.frame(p.adjust(deltaCoupling_results$pvals, method="fdr"))
colnames(FDRcorr_coupling_pval) <- "FDRcorr_coupling_pval"

deltaCoupling_results <- cbind(deltaCoupling_results, FDRcorr_coupling_pval, Yeo7_part, mean_nback_activation)

sig_FDRcorr_coupling_pval <- as.data.frame(subset(deltaCoupling_results, deltaCoupling_results$FDRcorr_coupling_pval < 0.05))
dim(sig_FDRcorr_coupling_pval)
sum(sig_FDRcorr_coupling_pval$tstats > 0)
sum(sig_FDRcorr_coupling_pval$tstats < 0)

min(abs(sig_FDRcorr_coupling_pval$tstats))

pos_results <- subset(sig_FDRcorr_coupling_pval, sig_FDRcorr_coupling_pval$tstats > 0)
neg_results <- subset(sig_FDRcorr_coupling_pval, sig_FDRcorr_coupling_pval$tstats < 0)

rcorr(deltaCoupling_results$tstats, deltaCoupling_results$mean_nback_activation)

#####################################
## Write out T-statistic in Fig 4B ##
#####################################
write.table(deltaCoupling_results$tstats, "/data/jux/BBL/projects/pncBaumStructFunc/replication/group_stats/longitudinal/n294_Schaefer400_delta_func_posPC_predicting_delta_Coupling_tstat.txt", col.names = FALSE, row.names=FALSE)
######################################################################################################################################################

#####################
## Save Rdata File ##
#####################
save.image("/data/jux/BBL/projects/pncBaumStructFunc/replication/group_stats/longitudinal/n294_TD_longitudinal_schaefer400_coupling_analysis.RData")
