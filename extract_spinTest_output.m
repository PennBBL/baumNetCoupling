clear all
close all
clc

addpath('/data/jux/BBL/projects/pncBaumStructFunc/scripts/rotate_parcellation-master/Matlab')
addpath('/data/jux/BBL/projects/pncBaumStructFunc/spin_test/spin-test-master/scripts')
addpath('/data/jux/BBL/projects/pncBaumStructFunc/scripts/colorbrewer/cbrewer/cbrewer')
addpath('/data/jux/BBL/projects/pncBaumStructFunc/scripts/colorbrewer/MatPlotLib2.0_colormaps')
addpath('/data/jux/BBL/projects/pncBaumStructFunc/scripts/WSBM/WSBM_v1.2')

%%%%%%%%%%%%%%%%%
%% Medial Wall %%
%%%%%%%%%%%%%%%%%
lh_mw_label = dlmread('/data/jux/BBL/projects/pncBaumDti/Schaefer2018_LocalGlobal_Parcellation/CBIG-master/stable_projects/brain_parcellation/Schaefer2018_LocalGlobal/Parcellations/FreeSurfer5.3/fsaverage6/label/lh.Medial_wall.label', '',1, 0);
rh_mw_label = dlmread('/data/jux/BBL/projects/pncBaumDti/Schaefer2018_LocalGlobal_Parcellation/CBIG-master/stable_projects/brain_parcellation/Schaefer2018_LocalGlobal/Parcellations/FreeSurfer5.3/fsaverage6/label/rh.Medial_wall.label', '',1, 0);

lh_mw_vert_idx = lh_mw_label(:,1);
rh_mw_vert_idx = rh_mw_label(:,1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Read in Spin Test Output: Mean Nback Coupling %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Load spin test workspace
load('/data/jux/BBL/projects/pncBaumStructFunc/spin_test/output/n727_schaefer400_groupMean_nbackCoupling_fsaverage6_rotationFS_1000perms_06042019.mat');

nreg=400;
nperm=1000;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Yeo7 Community Affiliation Vector %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ci = dlmread('/data/joy/BBL/applications/xcpEngine/atlas/schaefer400/schaefer400x7CommunityAffiliation.1D');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Assign Regional Values by Hemisphere %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
lh_idx=1:200;
rh_idx=201:400;

lh_ci=ci([lh_idx]);
rh_ci=ci([rh_idx]);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Left Hemisphere Surface Information %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[lh_verts,lh_faces] = read_surf('/data/jux/BBL/projects/pncBaumDti/Schaefer2018_LocalGlobal_Parcellation/CBIG-master/stable_projects/brain_parcellation/Schaefer2018_LocalGlobal/Parcellations/FreeSurfer5.3/fsaverage6/surf/lh.inflated');
[lh_verts2, lh_labels, lh_colortable] = read_annotation('/data/jux/BBL/projects/pncBaumDti/Schaefer2018_LocalGlobal_Parcellation/CBIG-master/stable_projects/brain_parcellation/Schaefer2018_LocalGlobal/Parcellations/FreeSurfer5.3/fsaverage6/label/lh.Schaefer2018_400Parcels_17Networks_order.annot');
names = {lh_colortable.struct_names};
names = names{1};

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Right Hemisphere Surface Information %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[rh_verts, rh_faces] = read_surf('/data/jux/BBL/projects/pncBaumDti/Schaefer2018_LocalGlobal_Parcellation/CBIG-master/stable_projects/brain_parcellation/Schaefer2018_LocalGlobal/Parcellations/FreeSurfer5.3/fsaverage6/surf/rh.inflated');
[rh_verts2, rh_labels, rh_colortable] = read_annotation('/data/jux/BBL/projects/pncBaumDti/Schaefer2018_LocalGlobal_Parcellation/CBIG-master/stable_projects/brain_parcellation/Schaefer2018_LocalGlobal/Parcellations/FreeSurfer5.3/fsaverage6/label/rh.Schaefer2018_400Parcels_17Networks_order.annot');
rh_names = {rh_colortable.struct_names};
rh_names = rh_names{1};

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% For each permuted map, extract regional values to build null distribution %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
perm_regional_coupling = zeros(nperm,nreg);

for p = 1:nperm

	%% Define left and right hemi surface values
	lh_currPerm=bigrotl(p,:)';
	rh_currPerm=bigrotr(p,:)';

	%% Resample values for rotated medial wall vertices
	lh_mw_idx = find(lh_currPerm==0);  %% For Margulies gradient, medial wall values == -100
	% lh_currPerm(lh_mw_idx) = datasample(nonzeros(lh_currPerm), length(lh_mw_idx));

	%% Resample values for rotated medial wall vertices
	rh_mw_idx = find(rh_currPerm==0);
	% rh_currPerm(rh_mw_idx) = datasample(nonzeros(rh_currPerm), length(rh_mw_idx));
	
	%% Allocate regional vectors
	regional_vals = zeros(nreg,1);
	lh_reg_vals = zeros(length(lh_ci),1);
	rh_reg_vals = zeros(length(rh_ci),1);
	
	for i = 1:length(lh_ci)
		% Skip medial wall (first row of colortable)
		n= i + 1;  
		lh_ind = lh_labels == lh_colortable.table(n,5);
		% lh_reg_vals(i) =  mean(lh_currPerm(lh_ind));  %% TAKE MEAN OF VERTICES WITHIN REGION
		lh_reg_vals(i) =  mode(lh_currPerm(lh_ind));  %% TAKE MODE OF VERTICES WITHIN REGION
	end

	for i = 1:length(rh_ci)
		% Skip medial wall (first row of colortable)
		n= i + 1;  
		rh_ind = rh_labels == rh_colortable.table(n,5);
		% rh_reg_vals(i) =  mean(rh_currPerm(rh_ind));   %% TAKE MEAN OF VERTICES WITHIN REGION
		rh_reg_vals(i) =  mode(rh_currPerm(rh_ind));     %% TAKE MODE OF VERTICES WITHIN REGION
	end

	%% Define regional measures
	regional_vals(lh_idx) = lh_reg_vals;
	regional_vals(rh_idx) = rh_reg_vals;

	perm_regional_coupling(p,:)= regional_vals';
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Build Null Distribution for Correlations %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
orig_coupling=dlmread('/data/jux/BBL/projects/pncBaumStructFunc/network_measures/Schaefer400/nback_restFC_coupling/coupling/n727_Schaefer400_thresh_norm_probSC_nbackFC_groupAvg_coupling_Spearman_r.txt');
reg_schaefer400_structPC = dlmread('/data/jux/BBL/projects/pncBaumStructFunc/network_measures/Schaefer400/nback_restFC_coupling/node_features/n727_thresh_norm_probSC_Schaefer400_Yeo7_mean_structPC.txt');
reg_schaefer400_myelin=dlmread('/data/jux/BBL/projects/pncBaumStructFunc/myelin_maps/output/hayashi_mean_Myelin_schaefer400.txt');
reg_schaefer400_nback_posPC=dlmread('/data/jux/BBL/projects/pncBaumStructFunc/network_measures/Schaefer400/nback_restFC_coupling/node_features/n727_Schaefer400_Yeo7_groupAvg_nbackFC_posPC.txt');
reg_schaefer400_margulies_gradient=dlmread('/data/jux/BBL/projects/pncBaumStructFunc/network_measures/Schaefer400/margulies_gradient/schaefer400x17_mean_regional_margulies_gradient.txt');
reg_evo_expansion = dlmread('/data/jux/BBL/projects/pncBaumStructFunc/cortical_scaling/schaefer400_rescaled_regional_evoExpansion.txt');
reg_number_edges=dlmread('/data/jux/BBL/projects/pncBaumStructFunc/network_measures/Schaefer400/nback_restFC_coupling/coupling/n727_Schaefer400_threshNormProbSC_nbackFC_groupAvg_number_effective_edges.txt');

%% Supplement %%
restFC_PCA_1_loading = dlmread('/data/jux/BBL/projects/pncBaumStructFunc/replication/network_metrics/node_features/n727_Schaefer400_Zscore_mean_restFC_PCA_comp1_loadings.txt');

%% Evo expansion
rh_evo_expansion = reg_evo_expansion(201:400);
rh_coupling = orig_coupling(201:400);
lh_coupling = orig_coupling(1:200);
mean_rh_coupling = (rh_coupling + lh_coupling) / 2;


%% Allocate null distributions
null_coupling_restFC_pca = zeros(nperm, 1);
null_coupling_myelin=zeros(nperm,1);
null_coupling_marguliesGradient=zeros(nperm,1);
null_coupling_func_PC=zeros(nperm,1);
null_coupling_struct_PC=zeros(nperm,1);
null_coupling_num_edges=zeros(nperm,1);
num_zeros=zeros(nperm,1);

null_eucDist_ageEffect=zeros(nperm,1);


for i = 1:nperm
	%% Myelin map
	tmp_perm_coupling = perm_regional_coupling(i,:)';
	zero_idx=find(tmp_perm_coupling==0);

	num_zeros(i) = length(zero_idx);

	%% Remove regions that are completely covered by rotated medial wall
	tmp_myelin = reg_schaefer400_myelin; 
	tmp_gradient = reg_schaefer400_margulies_gradient;
	tmp_func_PC = reg_schaefer400_nback_posPC;
	tmp_struct_PC = reg_schaefer400_structPC;
	tmp_num_edges= reg_number_edges;
	tmp_restFC_pca = restFC_PCA_1_loading;

	if length(zero_idx) > 1
		tmp_perm_coupling(zero_idx) = [];
		tmp_myelin(zero_idx) = [];
		tmp_gradient(zero_idx) = [];
		tmp_func_PC(zero_idx) = [];
		tmp_struct_PC(zero_idx) = [];
		tmp_num_edges(zero_idx) = [];
		tmp_restFC_pca(zero_idx) = [];
	
		
		%% RestFC PCA 1 Loading
		tmp_restFC_pca_corr=corr(tmp_perm_coupling, tmp_restFC_pca);
		null_coupling_restFC_pca(i) = tmp_restFC_pca_corr;
		
		%% Myelin
		tmp_myelin_corr = corr(tmp_perm_coupling, tmp_myelin);
		null_coupling_myelin(i)=tmp_myelin_corr;

		%% Margulies Principal Gradient
		tmp_gradient_corr = corr(tmp_perm_coupling, tmp_gradient);
		null_coupling_marguliesGradient(i)= tmp_gradient_corr;

		%% NbackFC Pos PC
		tmp_func_PC_corr = corr(tmp_perm_coupling, tmp_func_PC);
		null_coupling_func_PC(i) = tmp_func_PC_corr;

		%% Struct PC
		tmp_struct_PC_corr = corr(tmp_perm_coupling, tmp_struct_PC);
		null_coupling_struct_PC(i) = tmp_struct_PC_corr;

		%% Number of Effective Edges
		tmp_num_edges_corr = corr(tmp_perm_coupling, tmp_num_edges);
		null_coupling_num_edges(i) = tmp_num_edges_corr;

	else

		tmp_restFC_pca_corr=corr(tmp_perm_coupling, tmp_restFC_pca);
		null_coupling_restFC_pca(i) = tmp_restFC_pca_corr;

		tmp_myelin_corr = corr(tmp_perm_coupling, tmp_myelin);
		tmp_gradient_corr = corr(tmp_perm_coupling, tmp_gradient);
		tmp_func_PC_corr = corr(tmp_perm_coupling, tmp_func_PC);


		null_coupling_myelin(i)=tmp_myelin_corr;
		null_coupling_marguliesGradient(i)= tmp_gradient_corr;
		null_coupling_func_PC(i) = tmp_func_PC_corr;
		null_coupling_struct_PC(i) = tmp_struct_PC_corr;
		null_coupling_num_edges(i) = tmp_num_edges_corr;

	end

end

%%%%%%%%%%%%%%%
%% Figure 2B %%
%%%%%%%%%%%%%%%
[orig_r, orig_p] =corr(orig_coupling, reg_schaefer400_structPC)
structPC_sig_perms = find(null_coupling_struct_PC < orig_r);
structPC_perm_p =  length(structPC_sig_perms) / nperm

%%%%%%%%%%%%%%%
%% Figure 2C %%
%%%%%%%%%%%%%%%
[orig_r, orig_p] =corr(orig_coupling, reg_schaefer400_nback_posPC)
funcPC_sig_perms = find(null_coupling_func_PC < orig_r);
funcPC_perm_p =  length(funcPC_sig_perms) / nperm

%%%%%%%%%%%%%%%
%% Figure 2D %%
%%%%%%%%%%%%%%%
[orig_r, orig_p] =corr(orig_coupling, reg_schaefer400_margulies_gradient)
gradient_sig_perms = find(null_coupling_marguliesGradient < orig_r);
gradient_perm_p =  length(gradient_sig_perms) / nperm


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Supplemental Figure S5C: Coupling and RestFC PCA-1 %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[orig_r, orig_p] =corr(orig_coupling, restFC_PCA_1_loading)
restFC_pca_sig_perms = find(null_coupling_restFC_pca < orig_r);
restFC_pca_perm_p =  length(restFC_pca_sig_perms) / nperm


%% Mean Coupling ~ Myelin Map
[orig_r, orig_p] =corr(orig_coupling, reg_schaefer400_myelin)
myelin_sig_perms = find(null_coupling_myelin > orig_r);
myelin_perm_p =  length(myelin_sig_perms) / nperm


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Evolutionary Expansion Perm Test %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
null_coupling_expansion = zeros(nperm,1);

mean_rh_coupling = (rh_coupling + lh_coupling) / 2;

for i = 1:nperm
	%% Nback Pos PC
	tmp_expansion = reg_evo_expansion;
	tmp_rh_expansion = tmp_expansion(201:400);

	tmp_perm_coupling = perm_regional_coupling(i,:)';
	tmp_rh_perm_coupling = tmp_perm_coupling(201:400);
	tmp_lh_perm_coupling = tmp_perm_coupling(1:200);
	tmp_mean_perm_coupling = (tmp_rh_perm_coupling + tmp_lh_perm_coupling) / 2;

	zero_idx=find(tmp_mean_perm_coupling==0);

	if length(zero_idx) > 1
		tmp_mean_perm_coupling(zero_idx) = [];
		tmp_rh_expansion(zero_idx) = [];
		
		%% Evo Expansion
		tmp_expansion_corr = corr(tmp_mean_perm_coupling, tmp_rh_expansion);
		null_coupling_expansion(i)=tmp_expansion_corr;

	else
		tmp_expansion_corr = corr(tmp_mean_perm_coupling, tmp_rh_expansion);
		null_coupling_expansion(i)=tmp_expansion_corr;
	end
end

%%%%%%%%%%%%%%%
%% Figure 2E %%
%%%%%%%%%%%%%%%
[orig_r, orig_p] =corr(mean_rh_coupling, rh_evo_expansion)
expansion_sig_perms = find(null_coupling_expansion < orig_r);
expansion_perm_p =  length(expansion_sig_perms) / nperm

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Read in Spin Test Output: Margulies Gradient Map %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear bigrotl bigrotr

%% Margulies map spin test %%
load('/data/jux/BBL/projects/pncBaumStructFunc/spin_test/output/n727_schaefer400_margulies_gradient_1000perms_06182019.mat');

nreg=400;
nperm=1000;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% For each permuted map, extract regional values to build null distribution %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
perm_margulies_gradient = zeros(nperm,nreg);

for p = 1:nperm

	%% Define left and right hemi surface values
	lh_currPerm=bigrotl(p,:)';
	rh_currPerm=bigrotr(p,:)';

	%% Resample values for rotated medial wall vertices
	lh_mw_idx = find(lh_currPerm==0);  %% For Margulies gradient, medial wall values == -100
	% lh_currPerm(lh_mw_idx) = datasample(nonzeros(lh_currPerm), length(lh_mw_idx));

	%% Resample values for rotated medial wall vertices
	rh_mw_idx = find(rh_currPerm==0);
	% rh_currPerm(rh_mw_idx) = datasample(nonzeros(rh_currPerm), length(rh_mw_idx));
	
	%% Allocate regional vectors
	regional_vals = zeros(nreg,1);
	lh_reg_vals = zeros(length(lh_ci),1);
	rh_reg_vals = zeros(length(rh_ci),1);
	
	for i = 1:length(lh_ci)
		% Skip medial wall (first row of colortable)
		n= i + 1;  
		lh_ind = lh_labels == lh_colortable.table(n,5);
		% lh_reg_vals(i) =  mean(lh_currPerm(lh_ind));  %% TAKE MEAN OF VERTICES WITHIN REGION
		lh_reg_vals(i) =  mode(lh_currPerm(lh_ind));  %% TAKE MODE OF VERTICES WITHIN REGION
	end

	for i = 1:length(rh_ci)
		% Skip medial wall (first row of colortable)
		n= i + 1;  
		rh_ind = rh_labels == rh_colortable.table(n,5);
		% rh_reg_vals(i) =  mean(rh_currPerm(rh_ind));   %% TAKE MEAN OF VERTICES WITHIN REGION
		rh_reg_vals(i) =  mode(rh_currPerm(rh_ind));     %% TAKE MODE OF VERTICES WITHIN REGION
	end

	%% Define regional measures
	regional_vals(lh_idx) = lh_reg_vals;
	regional_vals(rh_idx) = rh_reg_vals;

	perm_margulies_gradient(p,:)= regional_vals';
end

%% Export output
dlmwrite('/data/jux/BBL/projects/pncBaumStructFunc/spin_test/output/permuted_regional_measures/n727_PERMUTED_schaefer400_groupMode_regional_marguliesGradient_06202019.txt', perm_margulies_gradient)


% mean_perm_coupling = mean(perm_regional_coupling);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Build Null Distribution for Correlations %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

F1Exec_coupling_effects=dlmread('/data/jux/BBL/projects/pncBaumStructFunc/local_results/n727_Schaefer400_F1ExecAcc_nbackCoupling_gam_tstat.txt');

%%%%%%%%%%%%%%%%
%% Supplement %%
%%%%%%%%%%%%%%%%
det_coupling_ageEffects = dlmread('/data/jux/BBL/projects/pncBaumStructFunc/local_results/n686_Schaefer400_normThresh_detSC_Communicability_nback_regCoupling_gam_Age_Zscore.txt');
WM_1back2back_coupling_ageEffects = dlmread('/data/jux/BBL/projects/pncBaumStructFunc/local_results/n727_Schaefer400_threshNormProbSC_1back-2back_nbackFC_regCoupling_gam_Age_Zscore.txt');
restFC_PCA_1_loading = dlmread('/data/jux/BBL/projects/pncBaumStructFunc/replication/network_metrics/node_features/n727_Schaefer400_Zscore_mean_restFC_PCA_comp1_loadings.txt');
EucDist_coupling_ageEffect = dlmread('/data/jux/BBL/projects/pncBaumStructFunc/local_results/n727_Schaefer400_EucDist_partialCorr_nback_regCoupling_gam_Age_Zscore.txt');

reg_schaefer400_margulies_gradient=dlmread('/data/jux/BBL/projects/pncBaumStructFunc/network_measures/Schaefer400/margulies_gradient/schaefer400x17_mean_regional_margulies_gradient.txt');
reg_evo_expansion = dlmread('/data/jux/BBL/projects/pncBaumStructFunc/cortical_scaling/schaefer400_rescaled_regional_evoExpansion.txt');
orig_coupling=dlmread('/data/jux/BBL/projects/pncBaumStructFunc/network_measures/Schaefer400/nback_restFC_coupling/coupling/n727_Schaefer400_thresh_norm_probSC_nbackFC_groupAvg_coupling_Spearman_r.txt');

%% Allocate null distributions
null_detCouplingAge_marguliesGradient=zeros(nperm,1);
null_WMCouplingAge_marguliesGradient=zeros(nperm,1);

null_eucDist_ageEffect_marguliesGradient=zeros(nperm,1);

null_F1Exec_marguliesGradient=zeros(nperm,1);

num_zeros=zeros(nperm,1);

for i = 1:nperm
	%% Myelin map
	tmp_perm_gradient = perm_margulies_gradient(i,:)';
	zero_idx=find(tmp_perm_gradient==0);

	num_zeros(i) = length(zero_idx);

	%% Remove regions that are completely covered by rotated medial wall
	tmp_detAgeEffects = det_coupling_ageEffects; 
	tmp_WMcoupling_ageEffects = WM_1back2back_coupling_ageEffects;
	tmp_eucDist_ageEffects= EucDist_coupling_ageEffect;
	tmp_F1Exec_coupling_effects = F1Exec_coupling_effects;

	if length(zero_idx) > 1
		tmp_perm_gradient(zero_idx) = [];
		tmp_detAgeEffects(zero_idx) = [];
		tmp_WMcoupling_ageEffects(zero_idx) = [];
		tmp_eucDist_ageEffects(zero_idx) = [];
		tmp_F1Exec_coupling_effects(zero_idx) = [];

		
		%% EucDist Coupling Age Effects ~ Gradient
		tmp_eucDist_ageEffects_corr = corr(tmp_perm_gradient, tmp_eucDist_ageEffects);
		null_eucDist_ageEffect_marguliesGradient(i) = tmp_eucDist_ageEffects_corr;

		%% Deterministic Coupling Age Effects ~ Gradient
		tmp_detAge_corr = corr(tmp_perm_gradient, tmp_detAgeEffects);
		null_detCouplingAge_marguliesGradient(i)=tmp_detAge_corr;

		%% Wm Coupling ~ Gradient
		tmp_WmAge_corr = corr(tmp_perm_gradient, tmp_WMcoupling_ageEffects);
		null_WMCouplingAge_marguliesGradient(i)= tmp_WmAge_corr;

		%% F1ExecEffect ~ Gradient
		tmp_F1Exec_marguliesGradient_corr = corr(tmp_perm_gradient, tmp_F1Exec_coupling_effects);
		null_F1Exec_marguliesGradient(i) = tmp_F1Exec_marguliesGradient_corr;

	
	else
		%% Deterministic Coupling Age Effects ~ Gradient
		tmp_detAge_corr = corr(tmp_perm_gradient, tmp_detAgeEffects);
		null_detCouplingAge_marguliesGradient(i)=tmp_detAge_corr;

		%% Wm Coupling ~ Gradient
		tmp_WmAge_corr = corr(tmp_perm_gradient, tmp_WMcoupling_ageEffects);
		null_WMCouplingAge_marguliesGradient(i)= tmp_WmAge_corr;

		%% EucDist Coupling Age Effects ~ Gradient
		tmp_eucDist_ageEffects_corr = corr(tmp_perm_gradient, tmp_eucDist_ageEffects);
		null_eucDist_ageEffect_marguliesGradient(i) = tmp_eucDist_ageEffects_corr;

		%% F1ExecEffect ~ Gradient
		tmp_F1Exec_marguliesGradient_corr = corr(tmp_perm_gradient, tmp_F1Exec_coupling_effects);
		null_F1Exec_marguliesGradient(i) = tmp_F1Exec_marguliesGradient_corr;

	end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% F1ExecCouplingEffect ~ Gradient: Main Result %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[orig_r, orig_p] =corr(F1Exec_coupling_effects, reg_schaefer400_margulies_gradient)
F1Exec_sig_perms = find(null_F1Exec_marguliesGradient > orig_r);
F1Exec_perm_p =  length(F1Exec_sig_perms) / nperm


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Det Age Effects: Figure S2-C %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[orig_r, orig_p] =corr(det_coupling_ageEffects, reg_schaefer400_margulies_gradient)
detAge_sig_perms = find(null_detCouplingAge_marguliesGradient > orig_r);
detAge_perm_p =  length(detAge_sig_perms) / nperm

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% WM 1back-2back Coupling Age Effects: Figure S3-C %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[orig_r, orig_p] =corr(WM_1back2back_coupling_ageEffects, reg_schaefer400_margulies_gradient)
WMcouplingAge_sig_perms = find(null_WMCouplingAge_marguliesGradient > orig_r);
WMcouplingAge_perm_p =  length(WMcouplingAge_sig_perms) / nperm

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% EucDist Coupling Age EffectsL Figure S4-C %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[orig_r, orig_p] =corr(EucDist_coupling_ageEffect, reg_schaefer400_margulies_gradient)
eucDist_ageEffect_sig_perms = find(null_eucDist_ageEffect_marguliesGradient > orig_r);
eucDist_ageEffect_perm_p =  length(eucDist_ageEffect_sig_perms) / nperm



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Spin Test for Age Effects on Coupling %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear bigrotr bigrotl


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Read in Spin Test Output: Coupling Age Effects  %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
load('/data/jux/BBL/projects/pncBaumStructFunc/spin_test/output/n727_schaefer400_threshNormProbSC_nbackCoupling_ageEffects_fsaverage6_rotationFS_1000perms.mat')

nreg=400;
nperm=1000;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% For each permuted map, extract regional values to build null distribution %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
perm_ageEffects = zeros(nperm,nreg);

for p = 1:nperm

	lh_currPerm=bigrotl(p,:)';
	rh_currPerm=bigrotr(p,:)';

	%% Allocate regional vectors
	regional_vals = zeros(nreg,1);
	lh_reg_vals = zeros(length(lh_ci),1);
	rh_reg_vals = zeros(length(rh_ci),1);
	
	for i = 1:length(lh_ci)
		% Skip medial wall (first row of colortable)
		n= i + 1;  
		lh_ind = lh_labels == lh_colortable.table(n,5);
		% lh_reg_vals(i) =  mean(lh_currPerm(lh_ind));
		lh_reg_vals(i) =  mode(lh_currPerm(lh_ind));
	end

	for i = 1:length(rh_ci)
		% Skip medial wall (first row of colortable)
		n= i + 1;  
		rh_ind = rh_labels == rh_colortable.table(n,5);
		% rh_reg_vals(i) =  mean(rh_currPerm(rh_ind));
		rh_reg_vals(i) =  mode(rh_currPerm(rh_ind));
	end

	%% Define regional measures
	regional_vals(lh_idx) = lh_reg_vals;
	regional_vals(rh_idx) = rh_reg_vals;

	perm_ageEffects(p,:)= regional_vals';
end

%% Export output
%dlmwrite('/data/jux/BBL/projects/pncBaumStructFunc/spin_test/output/permuted_regional_measures/n727_PERMUTED_schaefer400_regional_nbackCoupling_ageEffects.txt', perm_ageEffects)
dlmwrite('/data/jux/BBL/projects/pncBaumStructFunc/spin_test/output/permuted_regional_measures/n727_PERMUTED_schaefer400_MODE_regional_nbackCoupling_ageEffects.txt', perm_ageEffects)

% mean_perm_coupling = mean(perm_ageEffects);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Build Null Distribution for Correlations %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
orig_coupling = dlmread('/data/jux/BBL/projects/pncBaumStructFunc/network_measures/Schaefer400/nback_restFC_coupling/coupling/n727_Schaefer400_thresh_norm_probSC_nbackFC_groupAvg_coupling_Spearman_r.txt');
reg_long_age_effects=dlmread('/data/jux/BBL/projects/pncBaumStructFunc/local_results/n294_Schaefer400_threshNormProbSC_nbackFC_Coupling_lme_Age_effects.txt');

%%%%%%%%%%%%%%%%
%% Supplement %%
%%%%%%%%%%%%%%%%
det_coupling_ageEffects = dlmread('/data/jux/BBL/projects/pncBaumStructFunc/local_results/n686_Schaefer400_normThresh_detSC_Communicability_nback_regCoupling_gam_Age_Zscore.txt');
mean_1back_2back_coupling = dlmread('/data/jux/BBL/projects/pncBaumStructFunc/network_measures/Schaefer400/nback_restFC_coupling/coupling/n727_Schaefer400_threshNormProbSC_1back-2backFC_groupAvg_coupling_Spearman_r.txt');
restFC_PCA_1_loading = dlmread('/data/jux/BBL/projects/pncBaumStructFunc/replication/network_metrics/node_features/n727_Schaefer400_Zscore_mean_restFC_PCA_comp1_loadings.txt');
EucDist_coupling_ageEffect = dlmread('');

coupling_ageEffects = dlmread('/data/jux/BBL/projects/pncBaumStructFunc/network_measures/Schaefer400/nback_restFC_coupling/gam_results/n727_Schaefer400_thresh_normProbSC_nbackFC_Spearman_coupling_Age_GAM_Zscores.txt');

reg_schaefer400_margulies_gradient = dlmread('/data/jux/BBL/projects/pncBaumStructFunc/network_measures/Schaefer400/margulies_gradient/schaefer400x17_mean_regional_margulies_gradient.txt');
reg_evo_expansion = dlmread('/data/jux/BBL/projects/pncBaumStructFunc/cortical_scaling/schaefer400_rescaled_regional_evoExpansion.txt');

rh_evo_expansion = reg_evo_expansion(201:400);
rh_coupling_ageEffects = coupling_ageEffects(201:400);


null_couplingAgeEffect_restFC_pca = zeros(nperm, 1);
null_couplingAgeEffect_marguliesGradient = zeros(nperm,1);
null_couplingAgeEffect_funcPC = zeros(nperm,1);
null_coupling_long_ageEffects=zeros(nperm,1);

for i = 1:nperm
	tmp_funcPC = reg_schaefer400_nback_posPC;
	tmp_gradient = reg_schaefer400_margulies_gradient;
	tmp_restFC_pca = restFC_PCA_1_loading;
	tmp_long_ageEffect=reg_long_age_effects;
	
	tmp_perm_ageEffect = perm_ageEffects(i,:)';
	zero_idx=find(tmp_perm_ageEffect==0);

	%% Remove regions that are completely covered by rotated medial wall
	
	if length(zero_idx) > 1
		tmp_perm_ageEffect(zero_idx) = [];
		tmp_funcPC(zero_idx) = [];
		tmp_gradient(zero_idx) = [];
		tmp_restFC_pca(zero_idx) = [];
		tmp_long_ageEffect(zero_idx) = [];

		%% Nback Pos PC
		tmp_funcPC_corr = corr(tmp_perm_ageEffect, tmp_funcPC);
		null_couplingAgeEffect_funcPC(i)=tmp_funcPC_corr;
		%% Margulies Principal Gradient
		tmp_gradient_corr = corr(tmp_perm_ageEffect, tmp_gradient);
		null_couplingAgeEffect_marguliesGradient(i)= tmp_gradient_corr;
	
		%% Longitudinal age effect
		tmp_longAge_corr = corr(tmp_perm_ageEffect, tmp_long_ageEffect);
		null_coupling_long_ageEffects(i)= tmp_longAge_corr;

		%% RestFC PCA 1 Loading
		tmp_restFC_pca_corr = corr(tmp_restFC_pca, tmp_perm_ageEffect);
		null_couplingAgeEffect_restFC_pca(i) = tmp_restFC_pca_corr;

	else
		
		%% Func PC
		tmp_funcPC_corr = corr(tmp_perm_ageEffect, tmp_funcPC);
		null_couplingAgeEffect_funcPC(i) = tmp_funcPC_corr;
		%% Margulies gradient
		tmp_gradient_corr = corr(tmp_perm_ageEffect, tmp_gradient);
		null_couplingAgeEffect_marguliesGradient(i) = tmp_gradient_corr;

		%% Longitudinal age effect
		tmp_longAge_corr = corr(tmp_perm_ageEffect, tmp_long_ageEffect);
		null_coupling_long_ageEffects(i)= tmp_longAge_corr;

		%% RestFC PCA 1 Loading
		tmp_restFC_pca_corr = corr(tmp_restFC_pca, tmp_perm_ageEffect);
		null_couplingAgeEffect_restFC_pca(i) = tmp_restFC_pca_corr;
	end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Figure 3B: Age effects and Func PC %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[orig_r, orig_p] =corr(coupling_ageEffects, reg_schaefer400_nback_posPC)
funcPC_sig_perms = find(null_couplingAgeEffect_funcPC < orig_r);
funcPC_perm_p =  length(funcPC_sig_perms) / nperm

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Figure 3C: Age effects and Margulies gradient %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[orig_r, orig_p] =corr(coupling_ageEffects, reg_schaefer400_margulies_gradient)
gradient_sig_perms = find(null_couplingAgeEffect_marguliesGradient > orig_r);
gradient_perm_p =  length(gradient_sig_perms) / nperm

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Figure 4A: Cross-sectional and longitudinal age effects on coupling %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[orig_r, orig_p] =corr(coupling_ageEffects, reg_long_age_effects)
longAge_sig_perms = find(null_coupling_long_ageEffects > orig_r);
longAge_perm_p =  length(longAge_sig_perms) / nperm


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Supplemental Figure S5: CS Age effects and RestFC PCA 1 %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[orig_r, orig_p] =corr(coupling_ageEffects, restFC_PCA_1_loading)
restFC_pca_sig_perms = find(null_coupling_long_ageEffects > orig_r);
restFC_pca_perm_p =  length(restFC_pca_sig_perms) / nperm

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Evolutionary Expansion Perm Test %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
null_couplingAgeEffect_expansion = zeros(nperm,1);

for i = 1:nperm
	%% Nback Pos PC
	tmp_expansion = reg_evo_expansion;
	tmp_rh_expansion = tmp_expansion(201:400);

	tmp_perm_ageEffect = perm_ageEffects(i,:)';
	tmp_rh_perm_ageEffect = tmp_perm_ageEffect(201:400);
	zero_idx=find(tmp_rh_perm_ageEffect==0);

	if length(zero_idx) > 1
		tmp_rh_perm_ageEffect(zero_idx) = [];
		tmp_rh_expansion(zero_idx) = [];
		
		%% Evo Expansion
		tmp_expansion_corr = corr(tmp_rh_perm_ageEffect, tmp_rh_expansion);
		null_couplingAgeEffect_expansion(i)=tmp_expansion_corr;

	else
		tmp_expansion_corr = corr(tmp_rh_perm_ageEffect, tmp_rh_expansion);
		null_couplingAgeEffect_expansion(i)=tmp_expansion_corr;
	end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Figure 3E: Age effects and Evo Expansion %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[orig_r, orig_p] =corr(rh_coupling_ageEffects, rh_evo_expansion)
expansion_sig_perms = find(null_couplingAgeEffect_expansion > orig_r);
expansion_perm_p =  length(expansion_sig_perms) / nperm

