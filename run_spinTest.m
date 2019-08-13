clear all
close all
clc

addpath('/data/jux/BBL/projects/pncBaumStructFunc/spin_test/spin-test-master/scripts')
addpath('/data/jux/BBL/projects/pncBaumStructFunc/scripts/colorbrewer/cbrewer/cbrewer')
addpath('/data/jux/BBL/projects/pncBaumStructFunc/scripts/colorbrewer/MatPlotLib2.0_colormaps')
addpath('/data/jux/BBL/projects/pncBaumStructFunc/scripts/WSBM/WSBM_v1.2')

nreg=400

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Regional values for surface projection and spin test %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
reg_schaefer400_mean_nbackCoupling=dlmread('/data/jux/BBL/projects/pncBaumStructFunc/network_measures/Schaefer400/nback_restFC_coupling/coupling/n727_Schaefer400_thresh_norm_probSC_nbackFC_groupAvg_coupling_Spearman_r.txt');
reg_schaefer400_myelin=dlmread('/data/jux/BBL/projects/pncBaumStructFunc/myelin_maps/output/hayashi_mean_Myelin_schaefer400.txt');
reg_schaefer400_nback_posPC=dlmread('/data/jux/BBL/projects/pncBaumStructFunc/network_measures/Schaefer400/nback_restFC_coupling/node_features/n727_Schaefer400_Yeo7_wholebrain_mean_nbackFC_posPC.txt');
reg_schaefer400_margulies_gradient=dlmread('/data/jux/BBL/projects/pncBaumStructFunc/network_measures/Schaefer400/margulies_gradient/schaefer400x17_mean_regional_margulies_gradient.txt');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% DEFINE REGIONAL VALUES FOR SURFACE PROJECTION %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
reg_vals = reg_schaefer400_margulies_gradient;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Assign values for each Hemisphere %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
lh_idx=1:200;
rh_idx=201:400;

lh_vals= reg_vals([lh_idx]);
rh_vals= reg_vals([rh_idx]);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Left Hemisphere Surface Annotation %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[lh_verts, lh_faces] = read_surf('/data/jux/BBL/projects/pncBaumDti/Schaefer2018_LocalGlobal_Parcellation/CBIG-master/stable_projects/brain_parcellation/Schaefer2018_LocalGlobal/Parcellations/FreeSurfer5.3/fsaverage6/surf/lh.inflated');
[lh_verts2, lh_labels, lh_colortable] = read_annotation('/data/jux/BBL/projects/pncBaumDti/Schaefer2018_LocalGlobal_Parcellation/CBIG-master/stable_projects/brain_parcellation/Schaefer2018_LocalGlobal/Parcellations/FreeSurfer5.3/fsaverage6/label/lh.Schaefer2018_400Parcels_17Networks_order.annot');
lh_names = {lh_colortable.struct_names};
lh_names = lh_names{1};

lh_mw_idx = find(lh_labels == 65793);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Right Hemisphere Surface Annotation %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[rh_verts, rh_faces] = read_surf('/data/jux/BBL/projects/pncBaumDti/Schaefer2018_LocalGlobal_Parcellation/CBIG-master/stable_projects/brain_parcellation/Schaefer2018_LocalGlobal/Parcellations/FreeSurfer5.3/fsaverage6/surf/rh.inflated');
[rh_verts2, rh_labels, rh_colortable] = read_annotation('/data/jux/BBL/projects/pncBaumDti/Schaefer2018_LocalGlobal_Parcellation/CBIG-master/stable_projects/brain_parcellation/Schaefer2018_LocalGlobal/Parcellations/FreeSurfer5.3/fsaverage6/label/rh.Schaefer2018_400Parcels_17Networks_order.annot');
rh_names = {rh_colortable.struct_names};
rh_names = rh_names{1};

rh_mw_idx = find(rh_labels == 65793);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Define index of brain region labels %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
lh_id = [lh_colortable.table];
rh_id = [rh_colortable.table];

% Assign community affiliation to surface labels
labelsnew_lh = zeros(size(lh_labels));
labelsnew_rh = zeros(size(rh_labels));

% Assign community affiliation to surface labels
labelsnew_ci_lh = zeros(size(lh_labels));
labelsnew_ci_rh = zeros(size(rh_labels));

% Assign RestFC Age Effects to surface labels
labelsnew_restCoupling_ageEffect_lh = zeros(size(lh_labels));
labelsnew_restCoupling_ageEffect_rh = zeros(size(rh_labels));

for i = 1:length(lh_vals)
	% Skip medial wall (first row of colortable)
	n= i + 1;  
	ind = lh_labels == lh_id(n,5);
	labelsnew_lh(ind) = lh_vals(i);
end

for i = 1:length(rh_vals)
	% Skip medial wall (first row of colortable)
	n= i + 1;  
	ind = rh_labels == rh_id(n,5);
	labelsnew_rh(ind) = rh_vals(i);
end

%% Set medial wall vertices to Zero
labelsnew_lh(lh_mw_idx) = 0;
labelsnew_rh(rh_mw_idx) = 0;

%% Export medial wall index
dlmwrite('/data/jux/BBL/projects/pncBaumStructFunc/spin_test/brain_maps/lh.schaefer400_medial_wall_index_fsaverage6.txt', lh_mw_idx)
dlmwrite('/data/jux/BBL/projects/pncBaumStructFunc/spin_test/brain_maps/rh.schaefer400_medial_wall_index_fsaverage6.txt', rh_mw_idx)


%% Export fsaverage6 csv files for brain maps
dlmwrite('/data/jux/BBL/projects/pncBaumStructFunc/spin_test/brain_maps/lh.PNC_schaefer400_regional_margulies_gradient_fsaverage6.csv', labelsnew_lh)
dlmwrite('/data/jux/BBL/projects/pncBaumStructFunc/spin_test/brain_maps/rh.PNC_schaefer400_regional_margulies_gradient_fsaverage6.csv', labelsnew_rh)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% RUN SPIN TEST BETWEEN 2 CORTICAL MAPS %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Define parameters
nperm=1000;
output_path ='/data/jux/BBL/projects/pncBaumStructFunc/spin_test/output/n727_schaefer400_margulies_gradient_1000perms_06182019.mat';
left_input ='/data/jux/BBL/projects/pncBaumStructFunc/spin_test/brain_maps/lh.PNC_schaefer400_regional_margulies_gradient_fsaverage6.csv';
right_input ='/data/jux/BBL/projects/pncBaumStructFunc/spin_test/brain_maps/rh.PNC_schaefer400_regional_margulies_gradient_fsaverage6.csv';

%% Run spin test (custom version in fsaverage6)
GB_SpinPermuFS(left_input, right_input, nperm, output_path)
