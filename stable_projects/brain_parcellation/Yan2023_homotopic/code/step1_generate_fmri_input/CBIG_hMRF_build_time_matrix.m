function [lh_time_mat, rh_time_mat] = CBIG_hMRF_build_time_matrix(lh_fmri_fullpath_txt, rh_fmri_fullpath_txt, ...
    mesh_type, lh_output_file, rh_output_file, censor_file_fullpath)
% CBIG_hMRF_build_time_matrix(lh_fmri_fullpath_txt, rh_fmri_fullpath_txt, mesh_type, lh_output_file, rh_output_file)
%
% This function generates the premultiplied matrix with the input concatenated time series data.

% Input:
%   - <lh/rh>_fmri_fullpath_txt: (char string)
%     The path to the txt file that contains the list of fmri files for each subject, of <lh/rh> respectively.
%     Accepted fMRI file format: 'nii.gz'.
%     Each line should contain all fmri files of the current scan, separated by a single space.
%
%   - mesh_type: (char string)
%     By default, 'fsaverage6'.
%
% Optional input:
%
%   - <lh/rh>_output_file: (char string)
%     The path to the output <lh/rh> concatenated time series file; if not provided, output would not be saved.
%
%   - censor_file_fullpath: (char string)
%     The path to the txt file that contains the list of censor files for each subject, of right hemisphere.
%     Each censor file is a binary txt file, of length of total number of frames per run per subject.
%     Each line in this file should contain all censor files of the current scan, separated by a single space.
%
% Output:
%   - <lh/rh>_time_mat (matrix)
%     The concatenated and normalized fMRI time courses of NxT dimension. N is the number
%     of vertices of the given input surface mesh. T is the total number of time points in the resultant 
%     concatenated time courses.   
%
% Example:
%   CBIG_hMRF_build_time_matrix('./lh_fullpath.txt', './rh_fullpath.txt', 'fsaverage6',...
%   './lh_output.mat', './rh_output.mat')
%
% Written by Xiaoxuan Yan and CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md

lh_full_paths = table2cell(readtable(lh_fmri_fullpath_txt, 'Delimiter', ' ', 'ReadVariableNames', false));
rh_full_paths = table2cell(readtable(rh_fmri_fullpath_txt, 'Delimiter', ' ', 'ReadVariableNames', false));
if(~isempty(censor_file_fullpath))
    censor_file_fullpath = table2cell(readtable(censor_file_fullpath, 'Delimiter', ' ',...
        'ReadVariableNames', false));
    perform_censor = 1;
else
    perform_censor = 0;
end

[num_subs, num_scans] = size(lh_full_paths);

% pre-load required surface meshes
lh_avg_mesh = CBIG_ReadNCAvgMesh('lh', mesh_type, 'inflated', 'cortex');
lh_cortex_mask = lh_avg_mesh.MARS_label==2; % 2 indicates cortical vertex
rh_avg_mesh = CBIG_ReadNCAvgMesh('rh', mesh_type, 'inflated', 'cortex');
rh_cortex_mask = rh_avg_mesh.MARS_label==2; % 2 indicates cortical vertex

% initialization
for k = 1:num_subs
    fprintf('It is subject number %g \n',k);
    for i = 1:num_scans
        % check if current scan is available for both hemispheres
        if (~isnan(lh_full_paths{k,i}) & ~isempty(lh_full_paths{k,i})...
             & ~isnan(rh_full_paths{k,i}) & ~isempty(rh_full_paths{k,i})) 
            
            % check if the current paths are pointing to actual files
            if(~exist(lh_full_paths{k,i}, 'file'))
                error('Error. \n Input fMRI path %s is pointing to non-existent file.', lh_full_paths{k,l});
            end
            if(~exist(rh_full_paths{k,i}, 'file'))
                error('Error. \n Input fMRI path %s is pointing to non-existent file.', rh_full_paths{k,l});
            end
            
            % data transformation
            lh_input = MRIread(lh_full_paths{k,i});
            lh_input = reshape(lh_input.vol,[size(lh_input.vol,1)*size(lh_input.vol,2)*size(lh_input.vol,3)...
            size(lh_input.vol,4)]);
            rh_input = MRIread(rh_full_paths{k,i});
            rh_input = reshape(rh_input.vol,[size(rh_input.vol,1)*size(rh_input.vol,2)*size(rh_input.vol,3)...
            size(rh_input.vol,4)]);
            
            if(perform_censor)
                if(~exist(censor_file_fullpath{k,i}, 'file'))
                    error('Error. \n Input censor file path %s is pointing to non-existent file.',...
                     censor_file_fullpath{k,i});
                else
                    censor_mask = logical(load(censor_file_fullpath{k,i}));
                    lh_input = lh_input(:, censor_mask);
                    rh_input = rh_input(:, censor_mask);
                end
            end
            
            lh_normalized_data = row_wise_normalize_data(lh_input, lh_cortex_mask);
            rh_normalized_data = row_wise_normalize_data(rh_input, rh_cortex_mask);
            
            if(~exist('lh_time_mat', 'var'))
                lh_time_mat = lh_normalized_data;
                rh_time_mat = rh_normalized_data;
            else
                lh_time_mat = [lh_time_mat lh_normalized_data];
                rh_time_mat = [rh_time_mat rh_normalized_data];
            end
        end     
    end
end

% filter out the NaN data by column (time points), in case there is any
nan_col_idx = union(find(isnan(sum(lh_time_mat,1))), find(isnan(sum(rh_time_mat,1))));
lh_time_mat(:,nan_col_idx) = [];
rh_time_mat(:,nan_col_idx) = [];

if(exist('lh_output_file', 'var'))
    save(lh_output_file,'lh_time_mat', '-v7.3');
    save(rh_output_file,'rh_time_mat', '-v7.3');
end
end

function time_series_masked_demeaned_standardized = row_wise_normalize_data(time_series, cortical_mask)
% given the time series data, this function first masks out the data corresponding to the medial wall vertices
% then it demeans and standardize the input data.
time_series_masked = single(time_series(cortical_mask,:)); % mask out medial wall vertices    
time_series_masked_demeaned = bsxfun(@minus, time_series_masked, mean(time_series_masked,2)); % row wise demeaning
time_series_masked_demeaned_standardized = bsxfun(@rdivide,time_series_masked_demeaned,... 
    std(time_series_masked_demeaned',1)'); % row wise standardization
end