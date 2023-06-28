function CBIG_preproc_FCmetrics( lh_cortical_ROIs_file, rh_cortical_ROIs_file, subcortical_ROIs_file, lh_cortical_data_list, ...
    rh_cortical_data_list, subcortical_data_list, discard_frames_list, metric_type, output_dir, output_prefix )

% CBIG_preproc_FCmetrics( lh_cortical_ROIs_file, rh_cortical_ROIs_file, subcortical_ROIs_file, lh_cortical_data_list, ...
%     rh_cortical_data_list, subcortical_data_list, discard_frames_list, metric_type, output_dir, output_prefix )
%
% Compute FC (functional connectivity) metrics. This function will consider
% both cortical and subcortical ROIs. For example, if the users want to
% compute static Pearson's correlation (metric_type = 'Pearson_r'), this
% function will call CBIG_ComputeROIs2ROIsCorrelationMatrix.m 6 times to
% compute 6 sub-matrices (left-left, left-right, right-right,
% left-subcortical, right-subcortical, and subcortical-subcortical
% correlation matrices), and combine them together to get the final
% correlation matrix.
%
% Inputs:
%     - lh_cortical_ROIs_file:
%       The full name of cortical ROIs in left hemisphere, e.g. '<path to
%       CBIG_CODE_DIR>/stable_projects/brain_parcellation/Schaefer2018_LocalGlobal/Parcellations/FreeSurfer5.3
%       /fsaverage6/label/lh.Schaefer2018_400Parcels_17Networks_order.annot'.
%
%     - rh_cortical_ROIs_file:
%       The full name of cortical ROIs in right hemisphere, e.g. '<path to
%       CBIG_CODE_DIR>/stable_projects/brain_parcellation/Schaefer2018_LocalGlobal/Parcellations/FreeSurfer5.3
%       /fsaverage6/label/rh.Schaefer2018_400Parcels_17Networks_order.annot'.
%
%     - subcortical_ROIs_file:
%       The full name of subcortical ROIs in volumetric space. For example,
%       one subject Sub0001_Ses1's 19 aseg in functional space: '<path to
%       data>/Sub0001_Ses1.subcortex.parcels.func.nii.gz'.
%
%     - lh_cortical_data_list:
%       The cortical fMRI data list in left hemisphere of one subject. It
%       is a one-line text file contains the left hemisphere surface data
%       of all runs. Different runs are separated by spaces.
%       Each fMRI file in this list should be in the same space as
%       "lh_cortical_ROIs_file".
%
%     - rh_cortical_data_list:
%       The cortical fMRI data list in right hemisphere of one subject. It
%       is a one-line text file contains the right hemisphere surface data
%       of all runs. Different runs are separated by spaces.
%       Each fMRI file in this list should be in the same space as
%       "rh_cortical_ROIs_file".
%
%     - subcortical_data_list:
%       The volumetric fMRI data list for one subject. It is used to
%       extract subcortical timseries. It is a one-line text file contains
%       the volumetric data of all runs. Different runs are separated by
%       spaces.
%       Each fMRI file in this list should be in the same space as
%       "subcortical_ROIs_file".
%
%     - discard_frames_list:
%       The list contains the discarded frames files of all runs for one
%       subject. It is a one-line text file contains the discarded frames
%       files of all runs. Different runs are separated by spaces.
%       Each discarded frames file in this list is a text file with T rows
%       of "0" or "1", where T is the number of timepoints. "O" means this
%       frame should be removed; "1" means this frame is kept.
%       If discarded frame files are not needed, use 'NONE' for this input
%       parameter.
%
%     - metric_type:
%       Currently we only support 'Pearson_r'.
%
%     - output_dir:
%       Full path of output directory.
%
%     - output_prefix:
%       A string used in the beginning of the output relative filenames.
%       If metric_type == 'Pearson_r', this function will generate serveral
%       files named as:
%       (1) [output_dir '/' output_prefix '_lh2lh.mat'] --> lh to lh
%       correlation matrix.
%       (2) [output_dir '/' output_prefix '_lh2rh.mat'] --> lh to rh
%       correlation matrix.
%       (3) [output_dir '/' output_prefix '_rh2rh.mat'] --> rh to rh
%       correlation matrix.
%       (4) [output_dir '/' output_prefix '_lh2subcort.mat'] --> lh to
%       subcortical ROIs correlation matrix.
%       (5) [output_dir '/' output_prefix '_rh2subcort.mat'] --> rh to
%       subcortical ROIs correlation matrix.
%       (6) [output_dir '/' output_prefix '_subcort2subcort.mat'] -->
%       subcortical to subcortical ROIs correlation matrix.
%       (7) [output_dir '/' output_prefix '_all2all.mat'] --> this is the
%       final output correlation matrix.
%
% Written by CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md

if(strcmp(metric_type, 'Pearson_r'))
    if(~exist(output_dir))
        mkdir(output_dir);
    end
    
    %%% compute lh to lh correlation
    lh2lh_out_file = [output_dir '/' output_prefix '_lh2lh.mat'];
    CBIG_ComputeROIs2ROIsCorrelationMatrix(lh2lh_out_file, lh_cortical_data_list, lh_cortical_data_list, discard_frames_list, ...
        lh_cortical_ROIs_file, lh_cortical_ROIs_file, 'NONE', 'NONE', 1, 0);
    
    %%% compute lh to rh correlation
    lh2rh_out_file = [output_dir '/' output_prefix '_lh2rh.mat'];
    CBIG_ComputeROIs2ROIsCorrelationMatrix(lh2rh_out_file, lh_cortical_data_list, rh_cortical_data_list, discard_frames_list, ...
        lh_cortical_ROIs_file, rh_cortical_ROIs_file, 'NONE', 'NONE', 1, 0);
    
    %%% compute rh to rh correlation
    rh2rh_out_file = [output_dir '/' output_prefix '_rh2rh.mat'];
    CBIG_ComputeROIs2ROIsCorrelationMatrix(rh2rh_out_file, rh_cortical_data_list, rh_cortical_data_list, discard_frames_list, ...
        rh_cortical_ROIs_file, rh_cortical_ROIs_file, 'NONE', 'NONE', 1, 0);
    
    %%% compute lh to subcortical correlation
    lh2subcort_out_file = [output_dir '/' output_prefix '_lh2subcort.mat'];
    CBIG_ComputeROIs2ROIsCorrelationMatrix(lh2subcort_out_file, lh_cortical_data_list, subcortical_data_list, discard_frames_list, ...
        lh_cortical_ROIs_file, subcortical_ROIs_file, 'NONE', 'NONE', 1, 0);
    
    %%% compute rh to subcortical correlation
    rh2subcort_out_file = [output_dir '/' output_prefix '_rh2subcort.mat'];
    CBIG_ComputeROIs2ROIsCorrelationMatrix(rh2subcort_out_file, rh_cortical_data_list, subcortical_data_list, discard_frames_list, ...
        rh_cortical_ROIs_file, subcortical_ROIs_file, 'NONE', 'NONE', 1, 0);
    
    %%% compute subcortical to subcortical correlation
    subcort2subcort_out_file = [output_dir '/' output_prefix '_subcort2subcort.mat'];
    CBIG_ComputeROIs2ROIsCorrelationMatrix(subcort2subcort_out_file, subcortical_data_list, subcortical_data_list, discard_frames_list, ...
        subcortical_ROIs_file, subcortical_ROIs_file, 'NONE', 'NONE', 1, 0);
    
    %%% combine
    lh2lh = load(lh2lh_out_file);
    lh2rh = load(lh2rh_out_file);
    rh2rh = load(rh2rh_out_file);
    lh2subcort = load(lh2subcort_out_file);
    rh2subcort = load(rh2subcort_out_file);
    subcort2subcort = load(subcort2subcort_out_file);
    if (~isequal(size(subcort2subcort.corr_mat),[19,19]))
        error('The number of subcortical ROIs is not 19. FC matrix will not be saved.')
    end
    lh2all = [lh2lh.corr_mat lh2rh.corr_mat lh2subcort.corr_mat];
    rh2all = [lh2rh.corr_mat' rh2rh.corr_mat rh2subcort.corr_mat];
    subcort2all = [lh2subcort.corr_mat' rh2subcort.corr_mat' subcort2subcort.corr_mat];
    
    corr_mat = [lh2all; rh2all; subcort2all];
    save([output_dir '/' output_prefix '_all2all.mat'], 'corr_mat');
else
    error('Unknown metric_type.');
end


end

