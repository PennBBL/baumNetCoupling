function [reg_coupling, mean_reg_coupling, num_struct_edges] = fcn_regional_coupling(struct_edgevec, func_edgevec, nsub, nreg)

% fcn_regional_coupling    Regional structure-function coupling
%
%   [reg_coupling, mean_reg_coupling] = fcn_regional_coupling(struct_edgevec, func_edgevec, nsub, nreg);
%
%   Structure-function coupling captures the degree to which regional functional connectivity profiles 
%	are aligned with regional structural connectivity profiles.
%
%   Inputs:     struct_edgevec,        array of vectorized structural connectivity for all subjects (nsub X nedge)
%
%               func_edgevec,        array of vectorized functional connectivity for all subjects (nsub X nedge)
%
%               nsub,				number of subjects in array 
%
%				nreg,				number of regions in connectivity matrix
%
%   Output:     reg_coupling,     regional coupling metric for each subject (nsub X nreg)
%
%               mean_reg_coupling,     group mean regional coupling (nreg)

nsub
nreg

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  NbackFC - probSC Coupling  %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
reg_corr = zeros(nsub,nreg);
pos_reg_corr = zeros(nsub,nreg);
neg_reg_corr = zeros(nsub,nreg);
num_struct_edges = zeros(nsub,nreg);

for i = 1:nsub
     for j = 1:nreg
          struct_mat=squareform(struct_edgevec(i, :));
          struct_edges = struct_mat(:, j);
          zero_struct_idx = find(struct_edges == 0);
          struct_edges(zero_struct_idx) = []; %% Remove edges where there is no direct SC
          
          func_mat=squareform(func_edgevec(i, :));
          func_edges = func_mat(:, j);

          neg_idx = find(func_edges < 0 );
          pos_idx = find(func_edges > 0 );

          neg_func_edges = func_edges;
          neg_func_edges(pos_idx) = 0;
          pos_func_edges = func_edges;
          pos_func_edges(neg_idx) = 0;

          %% Remove edges where there is no direct SC
          func_edges(zero_struct_idx) = []; 
          pos_func_edges(zero_struct_idx) = [];
          neg_func_edges(zero_struct_idx) = [];

          %% Correlation
          [r,p] = corr(struct_edges, func_edges, 'type', 'Spearman');
          reg_corr(i,j) = r;

          num_struct_edges(i,j)= length(struct_edges);
     end
end
