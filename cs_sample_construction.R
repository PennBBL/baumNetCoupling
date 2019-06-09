rm(list = ls())

####################
## Load libraries ##
####################
require(stringr)
require(R.matlab)
require(pracma)
require(Formula)

####################################################################################
### This script will construct a subject sample for Elis's Nback & Rest Analyses ###
####################################################################################
# T1
t1_qa <- read.csv("/data/joy/BBL/studies/pnc/n1601_dataFreeze/neuroimaging/t1struct/n1601_t1QaData_20170306.csv")

# Nback QA
# nback_qa <- read.csv("/data/joy/BBL/studies/pnc/n1601_dataFreeze/neuroimaging/nback/nbackGlmBlockDesign/n1601_NBACKQAData_20181001.csv")
# read.csv("/data/joy/BBL/studies/pnc/n1601_dataFreeze/neuroimaging/nback/n1601_NbackQAData_20170427.csv")

nbackConn_qa <- read.csv("/data/joy/BBL/studies/pnc/n1601_dataFreeze/neuroimaging/nback/nbackConnectNoRegress/n1601_NbackConnectQAData_20170718.csv")

nback_behavior <- read.csv("/data/joy/BBL/studies/pnc/n1601_dataFreeze/neuroimaging/nback/n1601_nbackBehavior_from_20160207_dataRelease_20161027.csv")

# LTN Status
health <- read.csv("/data/joy/BBL/studies/pnc/n1601_dataFreeze/health/n1601_health_20161214.csv")

# Rest QA
rest_qa <- read.csv("/data/joy/BBL/studies/pnc/n1601_dataFreeze/neuroimaging/rest/n1601_RestQAData_20170714.csv")

# DTI QA
dti_qa <- read.csv("/data/joy/BBL/studies/pnc/n1601_dataFreeze/neuroimaging/dti/n1601_dti_qa_20170301.csv")

# B0 Acquisition
protVal <- read.csv("/data/joy/BBL/studies/pnc/n1601_dataFreeze/neuroimaging/n1601_pnc_protocol_validation_params_status_20161220.csv")

# Demographics
demog <- read.csv("/data/joy/BBL/studies/pnc/n1601_dataFreeze/demographics/n1601_demographics_go1_20161212.csv")

# Subject identifier
tracker <- read.csv("/data/joy/BBL/studies/pnc/n1601_dataFreeze/demographics/n1601_tracker_go1_20161212.csv")
tracker <- tracker[c("bblid","scanid")]

########################
## Merge subject data ##
########################
final_df <- demog
final_df <- merge(final_df, health, by=c("bblid","scanid"))
final_df <- merge(final_df, t1_qa, by=c("bblid","scanid"))
final_df <- merge(final_df, nbackConn_qa, by=c("bblid","scanid"))
final_df <- merge(final_df, nback_behavior, by=c("bblid","scanid"))
final_df <- merge(final_df, rest_qa, by=c("bblid","scanid"))
final_df <- merge(final_df, dti_qa, by=c("bblid","scanid"))
final_df <- merge(final_df, protVal, by=c("bblid","scanid"))

##############################
## Apply subject exclusions ##
##############################
dti_nbackFC_rest_sample <- subset(final_df, ltnExcludev2 == 0 & t1Exclude == 0 & nbackFcExclude == 0 & restExclude == 0 & dti64Exclude == 0 & dti64ProtocolValidationStatusExclude == 0 & b0ProtocolValidationStatus == 1)

################################
## READ IN BRAIN NETWORK DATA ##
################################

## Define dimensions of the data
nsub <- dim(dti_nbackFC_rest_sample) [1]
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
  tmp_bblid <- dti_nbackFC_rest_sample$bblid [i]
  tmp_scanid <- dti_nbackFC_rest_sample$scanid [i]

  struct_netpath <- paste0('/data/joy/BBL/studies/pnc/processedData/diffusion/probabilistic_20171118/',tmp_bblid,'/*x', tmp_scanid, '/output/connectivity/',tmp_bblid,'*.mat')

 if (identical(Sys.glob(struct_netpath), character(0))) {

    bad_struct_bblids[i,] <- tmp_bblid
    bad_struct_scanids[i,] <- tmp_scanid

  } else  {
  	print(struct_netpath)
    }
}

###############################################
## Remove subjects with missing network data ##
###############################################
bad_scanid_idx <- which(bad_struct_scanids!=0)
length(bad_scanid_idx)
dti_nbackFC_rest_sample <- dti_nbackFC_rest_sample[-bad_scanid_idx, ]
dim(dti_nbackFC_rest_sample)

#############################################
## Check for functional brain network data ##
#############################################
nsub <- dim(dti_nbackFC_rest_sample) [1]
for(i in 1:nsub){
  tmp_bblid <- dti_nbackFC_rest_sample$bblid [i]
  tmp_scanid <- dti_nbackFC_rest_sample$scanid [i]

  nback_netpath <- paste0('/data/joy/BBL/studies/pnc/processedData/nback/nbackConnectNoRegression_201608/', tmp_bblid,'/*x',tmp_scanid, '/net/SchaeferPNC/*_SchaeferPNC_network.txt')

 if (identical(Sys.glob(nback_netpath), character(0))) {

    bad_func_bblids[i,] <- tmp_bblid
    bad_func_scanids[i,] <- tmp_scanid

  } else  {
  	print(nback_netpath)
  }
}


###############################################
## Remove subjects with missing network data ##
###############################################
bad_scanid_idx <- which(bad_func_scanids!=0)
# dti_nbackFC_rest_sample <- dti_nbackFC_rest_sample[-bad_scanid_idx, ]
# dim(dti_nbackFC_rest_sample)

## Allocate arrays for network data
struct_nets <- array(rep(0, nsub*nreg*nreg), dim=c(nsub, nreg, nreg))
struct_edgevec <- array(rep(0, nsub*nedge), dim=c(nsub, nedge))

nback_nets <- array(rep(0, nsub*nreg*nreg), dim=c(nsub, nreg, nreg))
nback_edgevec <- array(rep(0, nsub*nedge), dim=c(nsub, nedge))

rest_nets <- array(rep(0, nsub*nreg*nreg), dim=c(nsub, nreg, nreg))
rest_edgevec <- array(rep(0, nsub*nedge), dim=c(nsub, nedge))

# struct_volNormSC_nets <- array(rep(0, nsub*nreg*nreg), dim=c(nsub, nreg, nreg))
# struct_volNormSC_edgevec <- array(rep(0, nsub*nedge), dim=c(nsub, nedge))

# struct_length_nets <- array(rep(0, nsub*nreg*nreg), dim=c(nsub, nreg, nreg))
# struct_length_edgevec <- array(rep(0, nsub*nedge), dim=c(nsub, nedge))

# nback_taskreg_nets <- array(rep(0, nsub*nreg*nreg), dim=c(nsub, nreg, nreg))
# nback_taskreg_edgevec <- array(rep(0, nsub*nedge), dim=c(nsub, nedge))

####################################
## Read in subject brain networks ##
####################################
nsub <- dim(dti_nbackFC_rest_sample) [1]

for(i in 1:nsub) {
	tmp_bblid <- dti_nbackFC_rest_sample$bblid [i]
	tmp_scanid <- dti_nbackFC_rest_sample$scanid [i]

	struct_netpath <- paste0('/data/joy/BBL/studies/pnc/processedData/diffusion/probabilistic_20171118/',tmp_bblid,'/*x', tmp_scanid, '/output/connectivity/', tmp_bblid,'*.mat')
	nback_netpath <- paste0('/data/joy/BBL/studies/pnc/processedData/nback/nbackConnectNoRegression_201608/', tmp_bblid,'/*x',tmp_scanid, '/net/SchaeferPNC/*_SchaeferPNC_network.txt')
	rest_netpath <- paste0('/data/joy/BBL/studies/pnc/processedData/restbold/restbold_201607151621/', tmp_bblid,'/*x',tmp_scanid, '/net/SchaeferPNC/*_SchaeferPNC_network.txt')
	# struct_netpath <- paste0('/data/joy/BBL/studies/pnc/processedData/diffusion/probabilistic_20171118/',tmp_bblid,'/*x', tmp_scanid, '/output/connectivity/new*.mat')
	# nback_taskreg_netpath <- paste0('/data/joy/BBL/studies/pnc/n1601_dataFreeze/neuroimaging/nback/nbackConnectTaskRegress/n1601_nback_Network_schaefer400/schaefer400Network/', tmp_bblid,'_*',tmp_scanid,'_schaefer400_network.txt')

	###############################
	## Read in structural matrix ##
	###############################
	struct_mat <- readMat(Sys.glob(struct_netpath))
	tmp_struct_net <- struct_mat$streamlineCount.mat
	# tmp_length_net <- struct_mat$streamlineLength.mat
	# tmp_volNormSC_net <- struct_mat$volNormSC.mat
 
	## Normalize streamline count by total network strength
	total_strength <- sum(squareform(tmp_struct_net))
 	sq_mat <- squareform(tmp_struct_net)
	norm_struct_net <- squareform(sq_mat / total_strength)
	
	## Append subject network to array
	struct_nets[i, , ] <-  norm_struct_net

	## Append vectorized network to array
	struct_edgevec[i,  ] <-  squareform(norm_struct_net)
 	
	##############################
	## Read in Nback-FC network ##
	##############################
	tmp_nback_net <- as.matrix(read.table(Sys.glob(nback_netpath)))

	## Set diagonal to zero
	diag(tmp_nback_net) <- 0
	
	## Append array with all subject data
	nback_nets[i, , ] <-  tmp_nback_net

	## Append vectorized network to array
	nback_edgevec[i,  ] <-  squareform(tmp_nback_net)

	#############################
	## Read in Rest-FC network ##
	#############################
	tmp_rest_net <- as.matrix(read.table(Sys.glob(rest_netpath)))

	## Set diagonal to zero
	diag(tmp_rest_net) <- 0
	
	## Append array with all subject data
	rest_nets[i, , ] <-  tmp_rest_net

	## Append vectorized network to array
	rest_edgevec[i,  ] <-  squareform(tmp_rest_net)

}

##################################################################
## Remove subjects with disconnected nodes in Structural matrix ##
##################################################################
nsub <- dim(dti_nbackFC_rest_sample) [1]
struct_rowSums <- array(rep(0, nsub*nreg), dim=c(nsub, nreg))

for(i in 1:nsub){
	tmp_net <- struct_nets[i, , ]
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

## Remove subjects with missing network data 
dti_nbackFC_rest_sample <- dti_nbackFC_rest_sample[-bad_struct_node_idx, ]
dim(dti_nbackFC_rest_sample)

struct_edgevec <- struct_edgevec[-bad_struct_node_idx, ]
struct_nets <- struct_nets[-bad_struct_node_idx, , ]

nback_edgevec <- nback_edgevec[-bad_struct_node_idx, ]
nback_nets <- nback_nets[-bad_struct_node_idx, , ]

rest_edgevec <- rest_edgevec[-bad_struct_node_idx, ]
rest_nets <- rest_nets[-bad_struct_node_idx, , ]

# struct_volNormSC_edgevec <- struct_volNormSC_edgevec[-bad_struct_node_idx, ]
# struct_volNormSC_nets <- struct_volNormSC_nets[-bad_struct_node_idx, , ]

# struct_length_edgevec <- struct_length_edgevec[-bad_struct_node_idx, ]
# struct_length_nets <- struct_length_nets[-bad_struct_node_idx, , ]

nsub <- dim(dti_nbackFC_rest_sample) [1]

################################################################
## Remove subjects with disconnected nodes in Nback-FC matrix ##
################################################################
nsub <- dim(dti_nbackFC_rest_sample) [1]
func_rowSums <- array(rep(0, nsub*nreg), dim=c(nsub, nreg))

for(i in 1:nsub){
	tmp_net <- nback_nets[i, , ]
	tmp_func_rowSums <- rowSums(tmp_net, na.rm = FALSE)
	func_rowSums[i, ] <- tmp_func_rowSums
}

zero_rowSum <- matrix(rep(0, nsub))
for(i in 1:nsub){
	tmp_sums <- func_rowSums[i,]
	zero_rowSum[i,] <- length(which(is.na(tmp_sums)));  # which(is.na(tmp_sums))
}

## Define index of subjects with disconnected nodes
bad_nbackFC_node_idx <- which(zero_rowSum>0)

## Remove subjects with missing network data 
# dti_nbackFC_rest_sample <- dti_nbackFC_rest_sample[-bad_nbackFC_node_idx, ]
# dim(dti_nbackFC_rest_sample)

# struct_edgevec <- struct_edgevec[-bad_nbackFC_node_idx, ]
# nback_edgevec <- nback_edgevec[-bad_nbackFC_node_idx, ]
# struct_nets <- struct_nets[-bad_nbackFC_node_idx, , ]
# nback_nets <- nback_nets[-bad_nbackFC_node_idx, , ]

nsub <- dim(dti_nbackFC_rest_sample) [1]

###############################################################
## Remove subjects with disconnected nodes in Rest-FC matrix ##
###############################################################
nsub <- dim(dti_nbackFC_rest_sample) [1]
func_rowSums <- array(rep(0, nsub*nreg), dim=c(nsub, nreg))

for(i in 1:nsub){
	tmp_net <- rest_nets[i, , ]
	tmp_func_rowSums <- rowSums(tmp_net, na.rm = FALSE)
	func_rowSums[i, ] <- tmp_func_rowSums
}

zero_rowSum <- matrix(rep(0, nsub))
for(i in 1:nsub){
	tmp_sums <- func_rowSums[i,]
	zero_rowSum[i,] <- length(which(is.na(tmp_sums)));  # which(is.na(tmp_sums))
}

## Define index of subjects with disconnected nodes
bad_restFC_node_idx <- which(zero_rowSum>0)

## Remove subjects with missing network data 
# dti_nbackFC_rest_sample <- dti_nbackFC_rest_sample[-bad_restFC_node_idx, ]
# dim(dti_nbackFC_rest_sample)

# struct_edgevec <- struct_edgevec[-bad_restFC_node_idx, ]
# nback_edgevec <- nback_edgevec[-bad_restFC_node_idx, ]
# struct_nets <- struct_nets[-bad_restFC_node_idx, , ]
# nback_nets <- nback_nets[-bad_restFC_node_idx, , ]

nsub <- dim(dti_nbackFC_rest_sample) [1]


####################################
## Create NbackFC - RestFC Matrix ##
####################################
contrast_edgevec <- array(rep(0, nsub*nedge), dim=c(nsub, nedge))

for(i in 1:nsub){
	
	tmp_rest_edges <- rest_edgevec[i, ]
	tmp_nback_edges <- nback_edgevec[i, ]

	tmp_contrast_edges <- tmp_nback_edges - tmp_rest_edges
	contrast_edgevec[i, ] <- tmp_contrast_edges

}

#######################################
## Write out vectorized network data ##
#######################################
write.table(struct_edgevec, "/data/jux/BBL/projects/pncBaumStructFunc/replication/edgevec/n727_norm_probSC_edgevec.txt", col.names = FALSE, row.names=FALSE)

write.table(nback_edgevec, "/data/jux/BBL/projects/pncBaumStructFunc/replication/edgevec/n727_nbackFC_edgevec.txt", col.names = FALSE, row.names=FALSE)

write.table(rest_edgevec, "/data/jux/BBL/projects/pncBaumStructFunc/replication/edgevec/n727_restFC_edgevec.txt", col.names = FALSE, row.names=FALSE)

write.table(contrast_edgevec, "/data/jux/BBL/projects/pncBaumStructFunc/replication/edgevec/n727_nbackFC_restFC_contrast_edgevec.txt", col.names = FALSE, row.names=FALSE)

###################################
## Write out subject identifiers ##
###################################
write.table(dti_nbackFC_rest_sample$bblid, "/data/jux/BBL/projects/pncBaumStructFunc/replication/demographics/n727_nbackRest_coupling_bblids.txt", col.names = FALSE, row.names=FALSE)

write.table(dti_nbackFC_rest_sample$scanid, "/data/jux/BBL/projects/pncBaumStructFunc/replication/demographics/n727_nbackRest_coupling_scanids.txt", col.names = FALSE, row.names=FALSE)

#########################
## Export demographics ##
#########################
write.csv(dti_nbackFC_rest_sample, "/data/jux/BBL/projects/pncBaumStructFunc/replication/demographics/n727_dti_nbackFC_rest_sample_demographics.csv")

####################
## Save Workspace ##
####################
save.image("/data/jux/BBL/projects/pncBaumStructFunc/replication/demographics/n727_nback_restFC_NormProbSC_coupling_workspace.Rdata")
