function CBIG_MM_HCP_extract_TC_419(List, out_dir, roi_400_nii, roi_19_nii)

% CBIG_MM_HCP_extract_TC_419(List, out_dir, roi_400_nii)
% 
% This function extract time series file for HCP S1200
%
% Inputs:
%   - List
%     Path of list for time series for HCP S1200.
%
%   - out_dir
%     Path of the your output time series. You can
%     also change this to any place you want.
%
%   - roi_400_nii
%     Path of the Schaefer 2018 400 Parcels ROI file
%   
%   - roi_19_nii
%     Path of the 19 subcortical Parcels ROI file in CBIG order: 
%     (https://github.com/ThomasYeoLab/CBIG/blob/master/stable_projects/preprocessing/CBIG_fMRI_Preproc2016/CBIG_preproc_FCmetrics_wrapper.csh#L559C1-L577C28)
%
% Written by Tong He and CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md

if ~exist(out_dir, 'dir')
   mkdir(out_dir)
end

%% Read info in input files
% read subject run list
fileID = fopen(List);
subject_run_list = textscan(fileID, '%s');
subject_run_list = subject_run_list{1};
fclose(fileID);
subject_list = cell(size(subject_run_list, 1), 1);
run_list = cell(size(subject_run_list, 1), 1);
for i = 1:size(subject_run_list, 1)
    temp = strsplit(subject_run_list{i},'/');
    subject_list{i} = temp{10};
    temp = temp{13};
    run_list{i} = temp(1:14);
end

%% Load Thomas' labelling
x = ft_read_cifti(roi_400_nii,'mapname','array'); 
for i = 1:400
    ROI_cell{i} = find(x.dlabel==i);
end

x = ft_read_cifti(roi_19_nii,'mapname','array'); 
for i = 401:419
    ROI_cell{i} = find(x.dlabel==i-400);
end

% Compute Time courses for each subject
ign=0;
not_ign=0;
for i = 1:size(subject_run_list, 1)

    nii_file = subject_run_list{i};
    try
        cii_full=ft_read_cifti(nii_file);

        fullTC=cii_full.dtseries;
        % Extract TC in Thomas' parcels
        tc=zeros(length(ROI_cell),size(fullTC,2));
        for j=1:length(ROI_cell)
            tc(j,:)=nanmean(fullTC(ROI_cell{j},:));
        end

        %If no regression, use these two lines
        TC=tc';

        save_name = [out_dir, '/TC_419_ROIs_', subject_list{i}, '_', run_list{i}, '.mat'];
        save(save_name,'TC');
        not_ign=not_ign+1;
        Subj_ok{not_ign,1}=[subject_list{i}, run_list{i}];
    catch
    end
end

end
