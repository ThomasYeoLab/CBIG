function CBIG_VK2019_example_wrapper(out_dir)

% Written by Valeria Kebets and CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md

CBIG_CODE_DIR = getenv('CBIG_CODE_DIR');

current_dir = fileparts(mfilename('fullpath'));
pos_v = strfind(current_dir,filesep);
root_dir = fullfile(current_dir(1:pos_v(length(pos_v)) - 1));

files_dir = fullfile(root_dir,'examples','input');
ref_dir = fullfile(root_dir,'examples','correct_output');

if ~exist(out_dir), mkdir(out_dir); end

addpath(fullfile(getenv('CBIG_CODE_DIR'),'external_packages','matlab','non_default_packages','PLS_MIPlab'));
scripts_dir = fullfile(root_dir,'replication','code');
addpath(scripts_dir);

% Options
nPerms = 50; % permutations in PLS analysis
normalization_img = 1; % normalization options for RSFC data
normalization_behav = 1; % normalization options for behavior data
% 0 = no normalization
% 1 = zscore across all subjects
% 2 = zscore within groups (default)
% 3 = std normalization across subjects (no centering)
% 4 = std normalization within groups (no centering)

%% 1. Load data

disp('(1) Loading data');

load(fullfile(files_dir,'example_data_CNP.mat'));

%% 2. PLS analysis

disp('(2) Running PLS');

[U,S,V,Lx,Ly,explCovLC,LC_behav_loadings,LC_RSFC_loadings] = ...
    myPLS_analysis(X0,Y0,normalization_img,normalization_behav);

%% 3. Permutation testing

disp('(3) Permutation testing over LCs');
pvals_LC = myPLS_permut(X0,Y0,U,S,nPerms,diagnosis_grouping,normalization_img,normalization_behav,[]);

% FDR correction over the first 5 LCs
[signif_LC, ~] = FDR(pvals_LC(1:5), 0.05);

% Display significant LCs
disp('Significant LCs (after FDR correction):');
for iter_lc = 1:length(signif_LC)
    this_lc = signif_LC(iter_lc);
    disp(['LC' num2str(this_lc) ' (p = ' num2str(pvals_LC(this_lc),'%0.3f') ') explains ' ...
        num2str(round(100*explCovLC(this_lc))) '% of covariance']);
end

if length(signif_LC)==0
    disp('No significant LCs');
end

%% 4. Plots

which_LC = 1;

% RSFC & behavior loadings
nRois = 419;
CBIG_VK2019_plot_loadings(out_dir,LC_behav_loadings,behavNames,LC_RSFC_loadings,nRois,which_LC); 

% RSFC & behavior composite scores
CBIG_VK2019_plot_subjScores(Lx,Ly,CONST_DIAGNOSIS,diagnosis_grouping,which_LC); 

%% Compare results with reference

these_Lx = Lx;
these_Ly = Ly;
these_RSFC_loadings = LC_RSFC_loadings(1:1000,:); % only check first 1000
these_behav_loadings = LC_behav_loadings;

load(fullfile(ref_dir,'PLSresults_example.mat'));
ref_Lx = Lx;
ref_Ly = Ly;
ref_RSFC_loadings = LC_RSFC_loadings;
ref_behav_loadings = LC_behav_loadings;

diff_RSFC_scores = max(max(ref_Lx - these_Lx));
diff_behav_scores = max(max(ref_Ly - these_Ly));
diff_RSFC_loadings = max(max(ref_RSFC_loadings - these_RSFC_loadings));
diff_behav_loadings = max(max(ref_behav_loadings - these_behav_loadings));

% RSFC scores
if (diff_RSFC_scores > 1e-4)
    fprintf(['[FAILED] RSFC scores are different from reference,' ...
        ' max diff: %f.\n'], diff_RSFC_scores);
else
    fprintf('[PASSED] RSFC scores are the same as reference.\n');
end

  % Behav scores
if (diff_behav_scores > 1e-4)
    fprintf(['[FAILED] Behav scores are different from reference,' ...
        ' max diff: %f.\n'], diff_behav_scores);
else
    fprintf('[PASSED] Behav scores are the same as reference.\n');
end

% RSFC loadings
if (diff_RSFC_loadings > 1e-4)
    fprintf(['[FAILED] RSFC loadings are different from reference,' ...
        ' max diff: %f.\n'], diff_RSFC_loadings);
else
    fprintf('[PASSED] RSFC loadings are the same as reference.\n');
end

% Behav loadings
if (diff_behav_loadings > 1e-4)
    fprintf(['[FAILED] Behav loadings are different from reference,' ...
        ' max diff: %f.\n'], diff_behav_loadings);
else
    fprintf('[PASSED] Behav loadings are the same as reference.\n');
end

clear ref_Lx ref_Ly ref_RSFC_loadings ref_behav_loadings these_Lx these_Ly these_RSFC_loadings these_behav_loadings ...
    diff_RSFC_scores diff_behav_scores diff_RSFC_loadings diff_behav_loadings


%% Save results
save(fullfile(out_dir,'PLSresults_example.mat'),...
    'Lx','Ly','LC_RSFC_loadings','LC_behav_loadings');

rmpath(fullfile(getenv('CBIG_CODE_DIR'),'external_packages','matlab','non_default_packages','PLS_MIPlab'));
rmpath(scripts_dir);
