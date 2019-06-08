clear all
close all

addpath('/data/jag/gbaum/matlab_scripts/BCT/2017_01_15_BCT')
addpath('/data/jux/BBL/projects/pncBaumStructFunc/scripts/colorbrewer/cbrewer/cbrewer')
addpath('/data/jux/BBL/projects/pncBaumStructFunc/scripts/colorbrewer/MatPlotLib2.0_colormaps')
addpath('/data/jux/BBL/projects/pncBaumStructFunc/scripts/WSBM/WSBM_v1.2')
addpath('/home/gbaum/baumNetCoupling')
% addpath('/data/jux/BBL/projects/pncBaumStructFunc/scripts/coupling')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Load in vectorized network data %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
struct_edgevec = dlmread('/data/jux/BBL/projects/pncBaumStructFunc/replication/edgevec/n727_norm_probSC_edgevec.txt');
nback_edgevec = dlmread('/data/jux/BBL/projects/pncBaumStructFunc/replication/edgevec/n727_nbackFC_edgevec.txt');
rest_edgevec = dlmread('/data/jux/BBL/projects/pncBaumStructFunc/replication/edgevec/n727_restFC_edgevec.txt');

%% Read in Yeo network assignments
Yeo7_part = dlmread('/data/joy/BBL/applications/xcpEngine/atlas/schaefer400/schaefer400x7CommunityAffiliation.1D');
Yeo17_part = dlmread('/data/joy/BBL/applications/xcpEngine/atlas/schaefer400/schaefer400x17CommunityAffiliation.1D');

%% Number of subjects
nsub=size(struct_edgevec,1)
%% Number of edges
nedge=size(struct_edgevec,2)
%% Number of nodes (brain regions)
nreg=size(squareform(struct_edgevec(1,:)),1)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Set NaN in FC matrix to 0 %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
nback_isnan_idx=find(isnan(nback_edgevec));
nback_edgevec([nback_isnan_idx])=0;

rest_isnan_idx=find(isnan(rest_edgevec));
rest_edgevec([rest_isnan_idx])=0;

struct_isnan_idx=find(isnan(struct_edgevec));
struct_edgevec([struct_isnan_idx])=0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Calculate Coefficient of Variation %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
group_mats=zeros(nreg,nreg,nsub);
	for i = 1:nsub
		A=squareform(struct_edgevec(i,:)); 
		group_mats(:, :,i) = A;
	end
[W_thr, Wcv] = GLB_threshold_consistency(group_mats, 0.75);
%% Set diagonal to zero
Wcv =(Wcv - diag(diag(Wcv)));
%% Set NaN to zero
isnan_idx=find(isnan(Wcv));
Wcv(isnan_idx)=0;
probSC_Wcv=squareform(Wcv)';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Create thresholded structral matrices %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
sq_Wthr=squareform(W_thr);
thresh_idx=find(sq_Wthr==0);

thresh_struct_edgevec=struct_edgevec;
thresh_struct_edgevec(:, [thresh_idx]) = 0;

thresh_struct_mat = squareform(thresh_struct_edgevec(1,:));
thresh_dens = density_und(thresh_struct_mat)
struct_mat = squareform(struct_edgevec(1,:));
orig_dens = density_und(struct_mat)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% CALCULATE NBACK COUPLING %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[nback_reg_coupling, mean_nback_reg_coupling, num_struct_edges] = fcn_regional_coupling(thresh_struct_edgevec, nback_edgevec, nsub, nreg);

%% Export coupling measures
dlmwrite('/data/jux/BBL/projects/pncBaumStructFunc/replication/network_metrics/n727_Schaefer400_thresh_norm_probSC_nbackFC_regional_coupling_Spearman_r.txt', nback_reg_coupling)
dlmwrite('/data/jux/BBL/projects/pncBaumStructFunc/replication/network_metrics/n727_Schaefer400_thresh_norm_probSC_nbackFC_groupAvg_coupling_Spearman_r.txt', mean_nback_reg_coupling)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% CALCULATE REST COUPLING %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[rest_reg_coupling, mean_rest_reg_coupling, num_struct_edges] = fcn_regional_coupling(thresh_struct_edgevec, rest_edgevec, nsub, nreg);

%% Export coupling measures
dlmwrite('/data/jux/BBL/projects/pncBaumStructFunc/replication/network_metrics/n727_Schaefer400_thresh_norm_probSC_restFC_regional_coupling_Spearman_r.txt', rest_reg_coupling)
dlmwrite('/data/jux/BBL/projects/pncBaumStructFunc/replication/network_metrics/n727_Schaefer400_thresh_norm_probSC_restFC_groupAvg_coupling_Spearman_r.txt', mean_rest_reg_coupling)

corr(mean_rest_reg_coupling, mean_nback_reg_coupling)

mean_num_struct_edges=mean(num_struct_edges)';
dlmwrite('/data/jux/BBL/projects/pncBaumStructFunc/replication/network_metrics/n727_Schaefer400_coupling_regional_numStructEdges.txt', num_struct_edges)
dlmwrite('/data/jux/BBL/projects/pncBaumStructFunc/replication/network_metrics/n727_Schaefer400_coupling_groupAvg_numStructEdges.txt', mean_num_struct_edges)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% PCA on Mean Rest-FC Matrix %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
mean_restFC = squareform(mean(rest_edgevec));
sq_mean_restFC=squareform(mean_restFC)';
z_sq_mean_restFC = zscore(sq_mean_restFC); %% Z-transform Pearson r values
z_mean_restFC = squareform(z_sq_mean_restFC);
% atanh_sq_mean_restFC = atanh(atanh_sq_mean_restFC); %% Z-transform Pearson r values

[coeff,score,latent,tsquared,explained,mu] = pca(z_mean_restFC);

func_PC1_coeff=coeff(:,1);
func_PC2_coeff=coeff(:,2);

%% Export regional loadings from PCA
dlmwrite('/data/jux/BBL/projects/pncBaumStructFunc/replication/network_metrics/node_features/n727_Schaefer400_Zscore_mean_restFC_PCA_comp1_loadings.txt', func_PC1_coeff);
dlmwrite('/data/jux/BBL/projects/pncBaumStructFunc/replication/network_metrics/node_features/n727_Schaefer400_Zscore_mean_restFC_PCA_comp2_loadings.txt', func_PC2_coeff);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Calculate Node Strength %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
struct_nodeStrength = zeros(nsub, nreg);
nback_neg_nodeStrength = zeros(nsub, nreg);
nback_pos_nodeStrength = zeros(nsub, nreg);
rest_neg_nodeStrength = zeros(nsub, nreg);
rest_pos_nodeStrength = zeros(nsub, nreg);

%% Struct node strength
for i = 1:nsub
	A_struct = squareform(thresh_struct_edgevec(i,:));
	struct_str = strengths_und(A_struct);
	struct_nodeStrength(i, :) =  struct_str;
end	

%% Nback node strength
for i = 1:nsub
	A_func = squareform(nback_edgevec(i,:));
	[Spos Sneg] = strengths_und_sign(A_func);
	func_str = strengths_und_sign(A_func);
	nback_pos_nodeStrength(i, :) =  Spos;
	nback_neg_nodeStrength(i, :) =  Sneg;
end	

%% Rest node strength
for i = 1:nsub
	A_func = squareform(rest_edgevec(i,:));
	[Spos Sneg] = strengths_und_sign(A_func);
	func_str = strengths_und_sign(A_func);
	rest_pos_nodeStrength(i, :) =  Spos;
	rest_neg_nodeStrength(i, :) =  Sneg;
end	

%% Group Mean Coupling
group_mean_struct_nodeStrength = mean(struct_nodeStrength)';
group_mean_nback_pos_nodeStrength = mean(nback_pos_nodeStrength)';
group_mean_nback_neg_nodeStrength = mean(nback_neg_nodeStrength)';
group_mean_rest_pos_nodeStrength = mean(rest_pos_nodeStrength)';
group_mean_rest_neg_nodeStrength = mean(rest_neg_nodeStrength)';


% Full structural node strength
dlmwrite('/data/jux/BBL/projects/pncBaumStructFunc/replication/network_metrics/node_features/n727_Schaefer400_regional_struct_thresh_norm_probSC_nodeStrength.txt', struct_nodeStrength)
dlmwrite('/data/jux/BBL/projects/pncBaumStructFunc/replication/network_metrics/node_features/n727_Schaefer400_groupAvg_struct_thresh_norm_probSC_nodeStrength.txt', group_mean_struct_nodeStrength)

% Full func pos node strength
dlmwrite('/data/jux/BBL/projects/pncBaumStructFunc/replication/network_metrics/node_features/n727_Schaefer400_regional_nbackFC_pos_nodeStrength.txt', nback_pos_nodeStrength)
dlmwrite('/data/jux/BBL/projects/pncBaumStructFunc/replication/network_metrics/node_features/n727_Schaefer400_regional_restFC_pos_nodeStrength.txt', rest_pos_nodeStrength)

% Full func neg node strength
dlmwrite('/data/jux/BBL/projects/pncBaumStructFunc/replication/network_metrics/node_features/n727_Schaefer400_regional_nbackFC_neg_nodeStrength.txt', nback_neg_nodeStrength)
dlmwrite('/data/jux/BBL/projects/pncBaumStructFunc/replication/network_metrics/node_features/n727_Schaefer400_regional_restFC_neg_nodeStrength.txt', rest_neg_nodeStrength)

%% Group average node strength
dlmwrite('/data/jux/BBL/projects/pncBaumStructFunc/replication/network_metrics/node_features/n727_Schaefer400_mean_nback_pos_nodeStrength.txt', group_mean_nback_pos_nodeStrength)
dlmwrite('/data/jux/BBL/projects/pncBaumStructFunc/replication/network_metrics/node_features/n727_Schaefer400_mean_nback_neg_nodeStrength.txt', group_mean_nback_neg_nodeStrength)

dlmwrite('/data/jux/BBL/projects/pncBaumStructFunc/replication/network_metrics/node_features/n727_Schaefer400_mean_rest_pos_nodeStrength.txt', group_mean_rest_pos_nodeStrength)
dlmwrite('/data/jux/BBL/projects/pncBaumStructFunc/replication/network_metrics/node_features/n727_Schaefer400_mean_rest_neg_nodeStrength.txt', group_mean_rest_neg_nodeStrength)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% STRUCT Module-specific Segregation (mean PC) %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ncomms=(length(unique(Yeo7_part)))
mod_struct_mean_PC=zeros(nsub,ncomms);
struct_PC=zeros(nsub,nreg);

for i = 1:nsub
	A=squareform(thresh_struct_edgevec(i,:));
	% Calculate Participation Coefficient for each subject's weighted structural network using Shi's Cognitive System assigments for 234-region Lausanne (Nature Communications)
	P=participation_coef(A, Yeo7_part,0);
	struct_PC(i,:) = P';

	for j = 1:length(unique(Yeo7_part))
		% Make S_index
		S_index=find(Yeo7_part==j);
		mod_nodes=P(S_index);
		mod_struct_mean_PC(i,j)=mean(mod_nodes);
	end
end

mean_struct_PC=mean(struct_PC,1)';
wholebrain_mean_struct_PC=mean(struct_PC,2);

% figure; hist(mean_struct_PC)

dlmwrite('/data/jux/BBL/projects/pncBaumStructFunc/replication/network_metrics/node_features/n727_thresh_norm_probSC_Schaefer400_Yeo7_mean_structPC.txt', mean_struct_PC)
% dlmwrite('/data/jux/BBL/projects/pncBaumStructFunc/replication/network_metrics/node_features/n727_thresh_norm_probSC_Schaefer400_Yeo7_wholebrain_mean_structPC.txt', wholebrain_mean_struct_PC)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% NbackFC Module-specific Segregation (mean PC) %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ncomms=(length(unique(Yeo7_part)))
mod_nback_mean_posPC=zeros(nsub,ncomms);
mod_nback_mean_negPC=zeros(nsub,ncomms);
nback_posPC=zeros(nsub,nreg);
nback_negPC=zeros(nsub,nreg);

for i = 1:nsub
	A=squareform(nback_edgevec(i,:));

	% Calculate Participation Coefficient for each subject's weighted structural network using Shi's Cognitive System assigments for 234-region Lausanne (Nature Communications)
	[Ppos Pneg] = participation_coef_sign(A, Yeo7_part);
	nback_posPC(i,:) = Ppos';
	nback_negPC(i,:) = Pneg';

	for j = 1:length(unique(Yeo7_part))
		% Make S_index
		S_index=find(Yeo7_part==j);
		mod_nodes=Ppos(S_index);
		mod_nback_mean_posPC(i,j)=mean(mod_nodes);
		mod_nodes=Pneg(S_index);
		mod_nback_mean_negPC(i,j)=mean(mod_nodes);
	end
end

mean_nback_posPC = mean(nback_posPC)';
mean_nback_negPC = mean(nback_negPC)';

wholebrain_mean_nback_posPC=mean(nback_posPC,2);
wholebrain_mean_nback_negPC=mean(nback_negPC,2);

mean_modular_nback_posPC = mean(mod_nback_mean_posPC)';
mean_modular_nback_negPC = mean(mod_nback_mean_negPC)';

% figure; hist(mean_nback_posPC)
% figure; hist(mean_nback_negPC)

dlmwrite('/data/jux/BBL/projects/pncBaumStructFunc/replication/network_metrics/node_features/n727_Schaefer400_Yeo7_regional_nbackFC_posPC.txt', nback_posPC)
dlmwrite('/data/jux/BBL/projects/pncBaumStructFunc/replication/network_metrics/node_features/n727_Schaefer400_Yeo7_regional_nbackFC_negPC.txt', nback_negPC)

dlmwrite('/data/jux/BBL/projects/pncBaumStructFunc/replication/network_metrics/node_features/n727_Schaefer400_Yeo7_groupAvg_nbackFC_posPC.txt', mean_nback_posPC)
dlmwrite('/data/jux/BBL/projects/pncBaumStructFunc/replication/network_metrics/node_features/n727_Schaefer400_Yeo7_groupAvg_nbackFC_negPC.txt', mean_nback_negPC)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% RestFC Module-specific Segregation (mean PC) %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ncomms=(length(unique(Yeo7_part)))
mod_rest_mean_posPC=zeros(nsub,ncomms);
mod_rest_mean_negPC=zeros(nsub,ncomms);
rest_posPC=zeros(nsub,nreg);
rest_negPC=zeros(nsub,nreg);

for i = 1:nsub
	A=squareform(rest_edgevec(i,:));

	% Calculate Participation Coefficient for each subject's weighted structural network using Shi's Cognitive System assigments for 234-region Lausanne (Nature Communications)
	[Ppos Pneg] = participation_coef_sign(A, Yeo7_part);
	rest_posPC(i,:) = Ppos';
	rest_negPC(i,:) = Pneg';

	for j = 1:length(unique(Yeo7_part))
		% Make S_index
		S_index=find(Yeo7_part==j);
		mod_nodes=Ppos(S_index);
		mod_rest_mean_posPC(i,j)=mean(mod_nodes);
		mod_nodes=Pneg(S_index);
		mod_rest_mean_negPC(i,j)=mean(mod_nodes);
	end
end

mean_rest_posPC = mean(rest_posPC)';
mean_rest_negPC = mean(rest_negPC)';

wholebrain_mean_rest_posPC=mean(rest_posPC,2);
wholebrain_mean_rest_negPC=mean(rest_negPC,2);

% figure; hist(mean_rest_posPC)
% figure; hist(mean_rest_negPC)

dlmwrite('/data/jux/BBL/projects/pncBaumStructFunc/replication/network_metrics/node_features/n727_Schaefer400_Yeo7_regional_restFC_posPC.txt', rest_posPC)
dlmwrite('/data/jux/BBL/projects/pncBaumStructFunc/replication/network_metrics/node_features/n727_Schaefer400_Yeo7_regional_restFC_negPC.txt', rest_negPC)

dlmwrite('/data/jux/BBL/projects/pncBaumStructFunc/replication/network_metrics/node_features/n727_Schaefer400_Yeo7_groupAvg_restFC_posPC.txt', mean_rest_posPC)
dlmwrite('/data/jux/BBL/projects/pncBaumStructFunc/replication/network_metrics/node_features/n727_Schaefer400_Yeo7_groupAvg_restFC_negPC.txt', mean_rest_negPC)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Calculate Delta PC (Nback - Rest) %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
delta_posPC = zeros(nsub,nreg);
for i = 1:nsub
	rest_PC = rest_posPC(i, :)';
	nback_PC = nback_posPC(i, :)';
	tmp_delta = (nback_PC - rest_PC);
	delta_posPC(i,:) = tmp_delta';
end

mean_delta_posPC = mean(delta_posPC)';

delta_negPC = zeros(nsub,nreg);
for i = 1:nsub
	rest_PC = rest_negPC(i, :)';
	nback_PC = nback_negPC(i, :)';
	tmp_delta = (nback_PC - rest_PC);
	delta_negPC(i,:) = tmp_delta';
end

mean_delta_negPC = mean(delta_negPC)';

%% Export pos delta PC
dlmwrite('/data/jux/BBL/projects/pncBaumStructFunc/replication/network_metrics/node_features/n727_Schaefer400_Yeo7_mean_DELTA_nbackPC-restPC.txt', mean_delta_posPC)
dlmwrite('/data/jux/BBL/projects/pncBaumStructFunc/replication/network_metrics/node_features/n727_Schaefer400_Yeo7_regional_DELTA_nbackPC-restPC.txt', delta_posPC)

%% Export neg delta PC
dlmwrite('/data/jux/BBL/projects/pncBaumStructFunc/replication/network_metrics/node_features/n727_Schaefer400_Yeo7_mean_DELTA_Neg_nbackPC-restPC.txt', mean_delta_negPC)
dlmwrite('/data/jux/BBL/projects/pncBaumStructFunc/replication/network_metrics/node_features/n727_Schaefer400_Yeo7_regional_DELTA_Neg_nbackPC-restPC.txt', delta_negPC)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Compute Structural Modularity Quality of Yeo Partition %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
struct_Yeo_Q=zeros(nsub,1);
S=Yeo7_part;
gamma=1;

for i = 1:nsub
	A=squareform(thresh_struct_edgevec(i,:));
	% A=log(A);
	N = size(A,1);
	twomu = 0;
	for s=1
    	k=sum(A(:,:,s));
    	twom=sum(k);
    	twomu=twomu+twom;
    	indx=[1:N]+(s-1)*N;
    	B(indx,indx)=A(:,:,s)-gamma*k'*k/twom;
	end
	struct_Yeo_Q(i) = sum(B(bsxfun(@eq,S,S.'))) ./ twomu;
end

dlmwrite('/data/jux/BBL/projects/pncBaumStructFunc/replication/network_metrics/n727_Schaefer400_thresh_struct_norm_probSC_Yeo7_Q.txt', struct_Yeo_Q)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Compute RestFC Modularity Quality of Yeo Partition %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
rest_Yeo_Q=zeros(nsub,1);
S=Yeo7_part;
gamma=1;
adjmat = A;
twomu = 0;

for i = 1:nsub
	A=squareform(rest_edgevec(i,:));
	N = size(A,1);
	ntype = 'signed';
	Apos = A;
	Aneg = -A;
	Apos(A<0) = 0;
	Aneg(A>0) = 0;
	kpos=sum(Apos)';
	kneg=sum(Aneg)';
	twompos = sum(kpos);
	twomneg = sum(kneg);
	gpos = gamma;
    gneg = gamma;
	pmat = gpos*kpos*kpos'/twompos-gneg*kneg*kneg'/twomneg;
	B = A - pmat;
	H = dummyvar(S);
	delta = H*H';
	degnorm = twompos + twomneg;
	rest_Yeo_Q(i)= sum(sum(B.*delta)) ./ degnorm;
end

dlmwrite('/data/jux/BBL/projects/pncBaumStructFunc/replication/network_metrics/n727_Schaefer400_restFC_Yeo7_Q.txt', rest_Yeo_Q)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Compute NbackFC Modularity Quality of Yeo Partition %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
nback_Yeo_Q = zeros(nsub,1);
S=Yeo7_part;
gamma=1;
adjmat = A;
twomu = 0;

for i = 1:nsub
	A=squareform(nback_edgevec(i,:));
	N = size(A,1);
	ntype = 'signed';
	Apos = A;
	Aneg = -A;
	Apos(A<0) = 0;
	Aneg(A>0) = 0;
	kpos=sum(Apos)';
	kneg=sum(Aneg)';
	twompos = sum(kpos);
	twomneg = sum(kneg);
	gpos = gamma;
    gneg = gamma;
	pmat = gpos*kpos*kpos'/twompos-gneg*kneg*kneg'/twomneg;
	B = A - pmat;
	H = dummyvar(S);
	delta = H*H';
	degnorm = twompos + twomneg;
	nback_Yeo_Q(i)= sum(sum(B.*delta)) ./ degnorm;
end

dlmwrite('/data/jux/BBL/projects/pncBaumStructFunc/replication/network_metrics/n727_Schaefer400_nbackFC_Yeo7_Q.txt', nback_Yeo_Q)

%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Save Matlab Workspace %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%
save('/data/jux/BBL/projects/pncBaumStructFunc/replication/network_metrics/n727_Schaefer400_normProbSC_nback_restFC_network_metrics.mat')

