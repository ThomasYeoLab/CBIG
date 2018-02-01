function CBIG_preproc_plot_QC_RSFC_corr_vs_distance_wrapper( sub_dir, sub_list, QC_type, RSFC_stem, output_name, distance_file )

% CBIG_preproc_plot_QC_RSFC_corr_vs_distance_wrapper( sub_dir, sub_list, QC_type, RSFC_stem, output_name )
% 
% Call CBIG_preproc_plot_QC_RSFC_corr_vs_distance_readdata.m to read in the
% QC (quality control) measure and RSFC (resting-state functional
% connectivity) matrix of each subject. Then read in ROIs to ROIs
% volumetric distance matrix. Pass these matrices to
% CBIG_preproc_plot_QC_RSFC_corr_vs_distance_matrix.m to plot the QC-RSFC
% correlation vs ROIs to ROIs distances figure.
% 
% This code assumes the data went through our CBIG_fMRI_Preproc2016
% preprocessing pipeline. For example, the functional connecitivity matrix
% for each subject is
% <path_to_subject>/<subject_name>/FC_metrics/Pearson_r/<subjact_name>_<RSFC_stem>.mat;
% the mean FD for each subject is
% <path_to_subject>/<subject_name>/bold/<run#>/<subject_name>_bld<run#>_rest_skip4_stc_mc_rel_mean.rms;
% etc.
% 
% Inputs:
%     - sub_dir:
%       The absolute path where the preprocessed subjects are stored.
% 
%     - sub_list:
%       The full name of a text file, in which each line is one subject ID.
% 
%     - QC_type:
%       A string specify what is the QC (quality control) measure. Currenly
%       we only support for 'FD_mean', 'FD_std', and 'censored_frames'.
%       'FD_mean' means the QC measure is the mean framewise displacement.
%       'FD_std' means the QC measure is the standard deviation of
%       framewise displacement. 'censored_frames' means the number of
%       censored frames determined by FD and/or DVARS.
% 
%     - RSFC_stem:
%       The stem of RSFC (resting-state functional connectivity) filenames.
%       For example, if the RSFC filename for one subject is:
%       <path_to_data>/Sub0001_Ses1_rest_skip4_stc_mc_residc_interp_FDRMS0.2_DVARS50_bp_0.009_0.08_fs6_sm6_all2all.mat,
%       then RSFC_stem =
%       '_rest_skip4_stc_mc_residc_interp_FDRMS0.2_DVARS50_bp_0.009_0.08_fs6_sm6_all2all'.
% 
%     - output_name:
%       The full filename of output PNG/EPS image output. The user can
%       either specify the filename with or without the extension '.png'.
% 
%     - distance_file:
%       The full filename of ROIs to ROIs volumetric Euclidean distances,
%       if the user has a pre-computed one. It is assumed that the distance
%       file is computed by using CBIG_preproc_FC_metrics.m (when the file
%       is loaded, there will be a structural variable called "distance",
%       and it has a field "distance.distance" which is a num_ROI x
%       num_ROIs distance matrix.). If the user does not pass in this file,
%       this function will read from 
%       [getenv('CBIG_CODE_DIR') ... 
%       '/stable_projects/preprocessing/CBIG_fMRI_Preproc2016/data/ROIs2ROIs_VolAnatDistance/ROIs2ROIs_distance_Schaefer2018_400Parcels_17Networks_order_19aseg_FSLMNI152_1mm.mat']
% 
% Example 1:
%     CBIG_preproc_plot_QC_RSFC_corr_vs_distance_wrapper( ...
% '/mnt/eql/yeo3/data/GSP2016/CBIG_preproc_global_cen_bp/GSP_single_session/CBIG2016_preproc_global_cen_bp', ...
% '/mnt/eql/yeo3/data/GSP2016/CBIG_preproc_global_cen_bp/GSP_single_session/scripts/surf_exist_list.txt', ...
% 'FD_mean', '_rest_skip4_stc_mc_residc_interp_FDRMS0.2_DVARS50_bp_0.009_0.08_fs6_sm6_all2all', ...
% '/data/users/jingweil/storage/PreprocessingPipeline/QC_plots_2/data/QC-RSFC_distance_plot/plots/GSP/GSP_single_session/CBIG2016_global_cen_bp' )
% 
% 
% Written by Jingwei Li and CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md




%% read QC and RSFC for each subject
[QC, RSFC, y_label] = CBIG_preproc_plot_QC_RSFC_corr_vs_distance_readdata( sub_dir, sub_list, QC_type, RSFC_stem );


%% Read in distance matrix
if(~exist('distance_file', 'var'))
    distance_file = [getenv('CBIG_CODE_DIR') ...
        '/stable_projects/preprocessing/CBIG_fMRI_Preproc2016/data/ROIs2ROIs_VolAnatDistance/ROIs2ROIs_dist_Schaefer2018_400Parcels_17Networks_order_19aseg_FSLMNI152_1mm.mat'];
end
fprintf('ROIs to ROIs distances file is: \n\t%s\n', distance_file);
load(distance_file)
distance_matrix = distance.distance;


%% plot
CBIG_preproc_plot_QC_RSFC_corr_vs_distance_matrix( QC, RSFC, distance_matrix, y_label, output_name )


end
