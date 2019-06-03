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

############################################################
## Construct subject sample for longitudinal DTI analysis ##
############################################################

#######################
## DTI QA Exclusions ##
#######################
# DTI
dti_2416_qa <- read.csv("/data/joy/BBL/studies/pnc/n2416_dataFreeze/neuroimaging/dti/n2416_DTI64/n2416_dti_qa_20170301.csv")

# B0 acquisition
protVal_2416 <- read.csv("/data/joy/BBL/studies/pnc/n2416_dataFreeze/neuroimaging/n2416_pnc_protocol_validation_params_status_20170103.csv")

dti_long <- dti_2416_qa
dti_long <- merge(dti_long, protVal_2416 ,by=c("bblid","scanid")) 

#####################################################
## Only retain subjects with at least 2 timepoints ##
#####################################################
dti_long <- subset(dti_long, dti_long$TotalNtimepoints > 1) 

##############################################
## DTI Protocal Validation Status Exclusion ##
##############################################
dti_long <- dti_long[which(dti_long$dti64ProtocolValidationStatusExclude==0), ] 

################################################
## B0 acquisition (for distortion correction) ##
################################################
dti_long <- dti_long[which(dti_long$b0ProtocolValidationStatus==1), ] 

############################
## DTI Roalf QA exclusion ##
############################
dti_long <- dti_long[which(dti_long$dti64Exclude==0), ] 

## order data by bblid
dti_long <- dti_long[order(dti_long$bblid),]

#####################
## T1 QA Exclusion ##
#####################
t1_qa <- read.csv("/data/joy/BBL/studies/pnc/n2416_dataFreeze/neuroimaging/t1struct/n2416_t1QaData_20170516.csv")

dti_long <- merge(dti_long, t1_qa, by=c("bblid","scanid"))

dti_long <- dti_long[which(dti_long$t1Exclude==0), ] 

#######################
## Rest QA Exclusion ##
#######################
rest_qa <- read.csv("/data/joy/BBL/studies/pnc/n2416_dataFreeze/neuroimaging/rest/n2416_RestQAData_20170714.csv")

dti_long <- merge(dti_long, rest_qa, by=c("bblid","scanid"))

# dti_long <- dti_long[which(dti_long$restExclude==0), ] 

########################
## Nback QA Exclusion ##
########################
nback_qa <- read.csv("/data/joy/BBL/studies/pnc/n2416_dataFreeze/neuroimaging/nback/n2416_NbackQAData_20170427.csv")

dti_long <- merge(dti_long, nback_qa, by=c("bblid","scanid"))

dti_long <- dti_long[which(dti_long$nbackExclude==0), ] 

########################
## Merge Demographics ##
########################
n2416_demogs <- read.csv("/data/joy/BBL/studies/pnc/n2416_dataFreeze/clinical/n2416_demographics_20170310.csv")

dti_long <- merge(dti_long, n2416_demogs, by=c("bblid","scanid"))

## order data by bblid
dti_long <- dti_long[order(dti_long$bblid),]

########################
## Clinical Exclusion ##
########################
clin_df <- read.csv("/data/joy/BBL/studies/pnc/n2416_dataFreeze/clinical/pnc_diagnosis_categorical_20170526.csv")

dti_long <- merge(dti_long, clin_df, by="bblid")

dti_long <- subset(dti_long, dti_long$pncGrpPsychCl=="TD" | dti_long$pncGrpPsychCl=="Resilient")

dim(dti_long)

##################
## Age in Years ##
##################
dti_long$age_in_years <- dti_long$scanageMonths / 12

################################
## READ IN BRAIN NETWORK DATA ##
################################

## Define dimensions of the data
nsub <- dim(dti_long) [1]
nreg <- 400
nedge <- 79800

## Allocate arrays for network data
bad_struct_bblids <- matrix(rep(0, nsub))
bad_struct_scanids <- matrix(rep(0, nsub))
bad_func_bblids <- matrix(rep(0, nsub))
bad_func_scanids <- matrix(rep(0, nsub))

#############################################
## Check for structural brain network data ##
#############################################
for(i in 1:nsub){
  tmp_bblid <- dti_long$bblid [i]
  tmp_scanid <- dti_long$scanid [i]

  struct_netpath <- paste0('/data/joy/BBL/studies/pnc/processedData/diffusion/probabilistic_20171118/',tmp_bblid,'/*x', tmp_scanid, '/output/connectivity/',tmp_bblid, '*.mat')

 if (identical(Sys.glob(struct_netpath), character(0))) {

    bad_struct_bblids[i,] <- tmp_bblid
    bad_struct_scanids[i,] <- tmp_scanid

  } else  {
  	print(Sys.glob(struct_netpath))
    }
}

###############################################
## Remove subjects with missing network data ##
###############################################
bad_scanid_idx <- which(bad_struct_scanids!=0)
length(bad_scanid_idx)
dti_long <- dti_long[-bad_scanid_idx, ]
dim(dti_long)

#############################################
## Check for functional brain network data ##
#############################################
nsub <- dim(dti_long) [1]
for(i in 1:nsub){
  tmp_bblid <- dti_long$bblid [i]
  tmp_scanid <- dti_long$scanid [i]

  func_netpath <- paste0('/data/joy/BBL/studies/pnc/processedData/nback/nbackConnectNoRegression_201608/', tmp_bblid,'/*x',tmp_scanid, '/net/SchaeferPNC/*_SchaeferPNC_network.txt')

 if (identical(Sys.glob(func_netpath), character(0))) {

    bad_func_bblids[i,] <- tmp_bblid
    bad_func_scanids[i,] <- tmp_scanid

  } #else  {
  	# print(Sys.Glob(func_netpath))
  # }
}

###############################################
## Remove subjects with missing network data ##
###############################################
bad_scanid_idx <- which(bad_func_scanids!=0)
length(bad_scanid_idx)
# dti_long <- dti_long[-bad_scanid_idx, ]
# dim(dti_long)

##########################################################
## Retain subjects with 2 timepoints of data passing QA ##
##########################################################
subj_freq <- as.data.frame(table(dti_long$bblid))
long_df <- subset(subj_freq, subj_freq$Freq==2)
dim(long_df)

colnames(long_df) <- c("bblid", "Freq")

dti_long <- merge(dti_long, long_df, by="bblid")

## Define number of subjects after retaining 2-timepoint data
nsub <- dim(dti_long) [1]

###################################
## Read in Connectivity Matrices ##
###################################

## Allocate arrays for network data
struct_nets <- array(rep(0, nsub*nreg*nreg), dim=c(nsub, nreg, nreg))
struct_edgevec <- array(rep(0, nsub*nedge), dim=c(nsub, nedge))

func_nets <- array(rep(0, nsub*nreg*nreg), dim=c(nsub, nreg, nreg))
func_edgevec <- array(rep(0, nsub*nedge), dim=c(nsub, nedge))

nsub <- dim(dti_long) [1]

for(i in 1:nsub){
	tmp_bblid <- dti_long$bblid [i]
	tmp_scanid <- dti_long$scanid [i]

	struct_netpath <- paste0('/data/joy/BBL/studies/pnc/processedData/diffusion/probabilistic_20171118/',tmp_bblid,'/*x', tmp_scanid, '/output/connectivity/', tmp_bblid, '*.mat')
	func_netpath <- paste0('/data/joy/BBL/studies/pnc/processedData/nback/nbackConnectNoRegression_201608/', tmp_bblid,'/*x',tmp_scanid, '/net/SchaeferPNC/*_SchaeferPNC_network.txt')
	
	###############################
	## Read in structural matrix ##
	###############################
	struct_mat <- readMat(Sys.glob(struct_netpath))
	tmp_struct_net <- struct_mat$streamlineCount.mat
	
	# tmp_volNormSC_net <- struct_mat$volNormSC.mat
 	# tmp_connProb_net <- struct_mat$connProbability.mat
	# diag(tmp_connProb_net) <- 0

	## Normalize streamline count by total network strength
	total_strength <- sum(squareform(tmp_struct_net))
 	sq_mat <- squareform(tmp_struct_net)
	norm_struct_net <- squareform(sq_mat / total_strength)
	
	## Append subject network to array
	struct_nets[i, , ] <-  norm_struct_net

	## Append vectorized network to array
	struct_edgevec[i,  ] <-  squareform(norm_struct_net)
	
	#################################
	## Reead in Functional network ##
	#################################
	tmp_func_net <- as.matrix(read.table(Sys.glob(func_netpath)))

	## Set diagonal to zero
	diag(tmp_func_net) <- 0
	
	## Append array with all subject data
	func_nets[i, , ] <-  tmp_func_net

	## Append vectorized network to array
	func_edgevec[i,  ] <-  squareform(tmp_func_net)
}

##################################################################
## Remove Subjects With Disconnected Nodes in Structural Matrix ##
##################################################################
nsub <- dim(dti_long) [1]
struct_rowSums <- array(rep(0, nsub*nreg), dim=c(nsub, nreg))

for(i in 1:nsub){
	tmp_net <- struct_nets[i, , ]   # Rowsum for norm SC
	# tmp_net <- struct_connProb_nets[i, , ]   # Rowsum for connProb
	tmp_struct_rowSums <- rowSums(tmp_net)
	struct_rowSums[i, ] <- tmp_struct_rowSums
}

zero_rowSum <- matrix(rep(0, nsub))
for(i in 1:nsub){
	tmp_sums <- struct_rowSums[i,]
	zero_rowSum[i,] <- sum(tmp_sums==0)
}

## Define index of subjects with disconnected nodes
bad_struct_node_idx <- which(zero_rowSum>0)
length(bad_struct_node_idx)

## BBLIDs for subjects with missing data:
dti_long$bblid[bad_struct_node_idx]

## Remove subjects with missing network data: 
## This removes 2 timepoints for one subject (88705), and 1 timepoint for others (112633 and 120052)
dti_long <- dti_long[-bad_struct_node_idx, ]
dim(dti_long)

struct_edgevec <- struct_edgevec[-bad_struct_node_idx, ]
struct_nets <- struct_nets[-bad_struct_node_idx, , ]

func_edgevec <- func_edgevec[-bad_struct_node_idx, ]
func_nets <- func_nets[-bad_struct_node_idx, , ]

## Remove subjects with only 1 timepoint remaining
which(dti_long$bblid==112633)
which(dti_long$bblid==120052)

## Define index for subjects with 1 timepoint remaining
single_tp_idx <- c(181, 230)

## Remove data for these subjects
dti_long <- dti_long[-single_tp_idx, ]

struct_edgevec <- struct_edgevec[-single_tp_idx, ]
struct_nets <- struct_nets[-single_tp_idx, ,]

func_edgevec <- func_edgevec[-single_tp_idx, ]
func_nets <- func_nets[-single_tp_idx, , ]

nsub <- dim(dti_long) [1]
nsub

##################################################################
## Remove subjects with disconnected nodes in Functional matrix ##
##################################################################
nsub <- dim(dti_long) [1]
func_rowSums <- array(rep(0, nsub*nreg), dim=c(nsub, nreg))

for(i in 1:nsub){
	tmp_net <- func_nets[i, , ]
	tmp_func_rowSums <- rowSums(tmp_net)
	func_rowSums[i, ] <- tmp_func_rowSums
}

zero_rowSum <- matrix(rep(0, nsub))
for(i in 1:nsub){
	tmp_sums <- func_rowSums[i,]
	zero_rowSum[i,] <- length(which(is.na(tmp_sums)));  # which(is.na(tmp_sums))
}

## Define index of subjects with disconnected nodes
bad_func_node_idx <- which(zero_rowSum>0)
length(bad_func_node_idx)

## (NO SUBJECTS WITH DISCONNECTED NODES)

nsub <- dim(dti_long) [1]

####################################################################
## Final check: Make sure all subjects have 2 time points of data ##
####################################################################
dti_long$Freq

#######################################
## Write out vectorized network data ##
#######################################
write.table(struct_edgevec, "/data/jux/BBL/projects/pncBaumStructFunc/replication/edgevec/n294_norm_probSC_edgevec.txt", col.names = FALSE, row.names=FALSE)
write.table(func_edgevec, "/data/jux/BBL/projects/pncBaumStructFunc/replication/edgevec/n294_nbackFC_edgevec.txt", col.names = FALSE, row.names=FALSE)

# write.table(struct_edgevec, "/data/jux/BBL/projects/pncBaumStructFunc/network_measures/Schaefer400/nback_restFC_coupling/n294_longitudinal/edgevec/n294_norm_probSC_edgevec.txt", col.names = FALSE, row.names=FALSE)
# write.table(func_edgevec, "/data/jux/BBL/projects/pncBaumStructFunc/network_measures/Schaefer400/nback_restFC_coupling/n294_longitudinal/edgevec/n294_nbackFC_edgevec.txt", col.names = FALSE, row.names=FALSE)

#########################
## Export demographics ##
#########################
write.csv(dti_long, "/data/jux/BBL/projects/pncBaumStructFunc/replication/demographics/n294_TD_longitudinal_schaefer400_coupling_demographics.csv")

save.image("/data/jux/BBL/projects/pncBaumStructFunc/replication/demographics/n294_TD_longitudinal_schaefer400_coupling_demographics.RData")
