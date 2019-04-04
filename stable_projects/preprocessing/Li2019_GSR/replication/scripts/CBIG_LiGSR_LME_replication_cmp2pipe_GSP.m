function CBIG_LiGSR_LME_replication_cmp2pipe_GSP( LME_dir )

% CBIG_LiGSR_LME_replication_cmp2pipe_GSP( LME_dir )
% 
% This function integrate the output files of
% CBIG_LiGSR_LME_replication_all_GSP.sh and compare the results of Baseline
% and Baseline+GSR pipelines.
% 
% Inputs:
%   - LME_dir
%     The output directory of CBIG_LiGSR_LME_replication_all_GSP.sh. It is
%     assumed that there are two subfolders call "Baseline" and "GSR" storing
%     the variance component model results of each preprocessing pipeline.
%     A file called [LME_dir '/compare_2pipe/allstats_cmp2pipelines.mat']
%     will be created to store the integrated result and the comparison
%     between the two pipelines. It contains:
%     (1) perc_improv: the percentage improvement of mean variance
%                      explained of pipeline 1 w.r.t. pipeline 2.
%     (2) m_jack: the jackknife mean of the explained variance difference
%                 between the two pipelines (averaged across all traits).
%     (3) v_jack: the jackknife variance of the explained variance
%                 difference between the two pipelines (averaged across all
%                 traits).
%     (4) IQR_pos: the number of traits of which the whole interquartile
%                  range is above zero, i.e. in at least 75% jackknife
%                  samples, pipeline 1 > pipeline 2.
%     (5) IQR_neg: the number of traits of which the whole interquartile
%                  range is below zero, i.e. in at least 75% jackknife
%                  samples, pipeline 1 < pipeline 2.
%     (6) med_pos: the number of traits of which the median across all
%                  jackknife samples > 0.
%     (7) med_neg: the number of traits of which the median across all
%                  jackknife samples < 0.
% 
% Written by Jingwei Li and CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md

addpath(fullfile(getenv('CBIG_CODE_DIR'), 'stable_projects', 'preprocessing', ...
    'Li2019_GSR', 'VarianceComponentModel', 'scripts', 'utilities'))

CBIG_LiGSR_LME_cmp2pipe_allstats( fullfile(getenv('CBIG_CODE_DIR'), ...
    'stable_projects', 'preprocessing', 'Li2019_GSR', 'replication', ...
    'scripts', 'GSP_lists/23behaviors_age_sex.txt'), ...
    fullfile(LME_dir, 'GSR'), fullfile(LME_dir, 'Baseline'), 862, 1000, 431, ...
    fullfile(LME_dir, 'compare_2pipe') )

rmpath(fullfile(getenv('CBIG_CODE_DIR'), 'stable_projects', 'preprocessing', ...
    'Li2019_GSR', 'VarianceComponentModel', 'scripts', 'utilities'))

end

