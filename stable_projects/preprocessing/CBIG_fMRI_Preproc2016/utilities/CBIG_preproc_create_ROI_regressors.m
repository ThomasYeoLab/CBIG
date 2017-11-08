function CBIG_preproc_create_ROI_regressors(fMRI_list, output_list, WB_mask, WM_mask, CSF_mask, diff_flag)

% CBIG_preproc_create_ROI_regressors(fMRI_list, output_list, WB_mask, WM_mask, CSF_mask)
%
% This function creates ROI (region of interest) regressors based on given
% masks: Whole Brain mask, White Matter mask or CSF mask. One application of
% this function is to create regressors for a single subject with multiple runs.
% It will
%
% 1) read in fMRI volumetric data from fMRI data list (fMRI_list) and 
%    reshape each input 4D volume to a NxT matrix, where T is number of time 
%    points of each timecourse, N is the number of voxels.
%
% 2) apply the mask to each column of NxT matrix and average the signal in
%    the mask, then we get first 1xT regressor. Take the derivative of the 
%    regressor and append 0 to the vector, we get second 1xT regressor. 
%    Transpose these two vectors to Tx1 vector and concat them column by 
%    column, then we get Tx2 matrix. 
%
% 3) If input have more than one mask, this function will concatenate different
%    ROI regressors (Tx2 matrix) colume by colume in the order: WB, WM,
%    CSF. Regressors will be written into files in output_list.
%
% Input:
%     - fMRI_list:
%       a text file including fMRI data list of all runs, each line is the 
%       full path of the fMRI file. For example:
%       fMRI_list = '/path_to_data_list/fMRI_data_list.txt', each line is
%       '/path_to_fMRI_data/rest_bld00?.nii.gz'.
%
%     - output_list:
%       a text file including output file list of all runs, each line is the
%       full path of the output file name. For example:
%       output_list = '/path_to_output_list/output_list.txt', 
%       each line is '/path_to_output/rest_bld00?_WB_WM_CSF.txt'.   
%
%     - WB_mask:
%       a full path pointing to the Whole Brain mask. If it's [], Whole
%       Brain regressor will not be computed.
%       For example:
%       WB_mask = '/path_to_mask/brainmask.nii.gz', 
%
%     - WM_mask:
%       a full path pointing to the White Matter mask. If it's [], White
%       Matter regressor will not be computed.
%       For example:
%       WM_mask = '/path_to_mask/func.wm.nii.gz', 
%
%     - CSF_mask:
%       a full path pointing to the CSF mask. If it's [], CSF regressor 
%       will not be computed.
%       For example:
%       CSF_mask = '/path_to_mask/func.ventricles.nii.gz', 
%
% 
% Example:
% CBIG_preproc_create_ROI_regressors('fMRI_list.txt', 'output_list.txt', 'wb_mask.nii.gz', 'wm_mask.nii.gz', 'csf_mask.nii.gz')
%
% Written by CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md


if(ischar(diff_flag))
    diff_flag = str2num(diff_flag);
end

%% get fMRI name
[fMRI_name, num_of_fMRI] = text2cell(fMRI_list);
[output_name, num_of_output] = text2cell(output_list);

if (num_of_output ~= num_of_fMRI)
    error('ERROR: number of output file is not even with number of fMRI file')
end

%% calculate ROI regressors and merge them 
for i = 1:num_of_fMRI
    mri = MRIread(fMRI_name{i});
    vol = mri.vol;
    mri_size = size(vol);
    vol_2d = reshape(vol, [mri_size(1)*mri_size(2)*mri_size(3) mri_size(4)]);
    
    all_regressors = [];
    if (~isempty(WB_mask))
    	WB_regressor = mtx2regressor(vol_2d, WB_mask, diff_flag);
        if(~isempty(find(isnan(WB_regressor))))
            WB_regressor = [];
            fprintf('WARNING: Whole brain mask is empty. The whole brain regressor is empty.\n');
        end
        all_regressors = [all_regressors WB_regressor];
    end
    if (~isempty(WM_mask))
    	WM_regressor = mtx2regressor(vol_2d, WM_mask, diff_flag);
        if(~isempty(find(isnan(WM_regressor))))
            WM_regressor = [];
            fprintf('WARNING: White matter mask is empty. The white matter regressor is empty.\n');
        end
        all_regressors = [all_regressors WM_regressor];
    end
    if (~isempty(CSF_mask))
    	CSF_regressor = mtx2regressor(vol_2d, CSF_mask, diff_flag);
        if(~isempty(find(isnan(CSF_regressor))))
            CSF_regressor = [];
            fprintf('WRANING: Ventricles mask is empty. The ventricles regressor is empty.\n');
        end
        all_regressors = [all_regressors CSF_regressor];
    end
    
    % write into output file
    dlmwrite(output_name{i}, all_regressors, ' ');
end




function [cell_array, length] = text2cell(text)
% read text file line by line and write it into cell_array
length = 0;
fid = fopen(text);
while (~feof(fid))
    length = length + 1;
    cell_array{length} = fgetl(fid);
end
fclose(fid);


function regressor = mtx2regressor(mtx, mask, diff_flag)
% Given mask, calculate ROI regressors
mask = MRIread(mask);
mask_vol = mask.vol;
mask_1d = mask_vol(:);

if(sum(logical(mask_1d)) == 0)
    regressor = nan(size(mtx, 2), 2);
    return;
end

mean_tp = mean(mtx(mask_1d==1, :), 1);
mean_tp_d = [0 diff(mean_tp)];

if(diff_flag == 1)
    regressor = [mean_tp' mean_tp_d'];
else
    regressor = mean_tp';
end

