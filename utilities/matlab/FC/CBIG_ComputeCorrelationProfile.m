function CBIG_ComputeCorrelationProfile(seed_mesh, target, output_file1, output_file2,...
    threshold, varargin_file1, varargin_file2, outlier_text, split_data)

% CBIG_ComputeCorrelationProfile(seed_mesh, target, output_file1,...
% output_file2, threshold, varargin_file1, varargin_file2, outlier_text)
%
% Compute surface correlation profiles of a single subject. The correlation
% profile will be computed for each run and then averaged across runs.
% 
% Inputs:
%     - seed_mesh
%       resolution of seed regions, e.g. 'fsaverage3', 'fs_LR_900'
% 
%     - target
%       resolution of target regions, e.g. 'fsaverage5', 'fs_LR_32k'
% 
%     - output_file1
%       'fsaverage*' spaces:
%       left hemisphere output correlation profile for input in 'fsaverage*' spaces, e.g.      
%       <path>/lh.Sub0033_Ses1.roifsaverage3.thres0.1.surf2surf_profile_scrub.nii.gz;
%       <path>/lh.Sub0033_Ses1.roifsaverage3.thres0.1.surf2surf_profile_scrub.mat;
%
%       'fs_LR*' spaces:
%       entire cortical output correlation profile for input in 'fs_LR*' spaces, e.g.
%       e.g. <path>/100206.roifs_LR_900.thres0.1.surf2surf_profile_scrub.mat.
% 
%     - output_file2  
%       'fsaverage*' spaces:
%       right hemisphere output correlation profile for input in 'fsaverage*' spaces, e.g.      
%       <path>/rh.Sub0033_Ses1.roifsaverage3.thres0.1.surf2surf_profile_scrub.nii.gz;
%       <path>/rh.Sub0033_Ses1.roifsaverage3.thres0.1.surf2surf_profile_scrub.mat;
%
%       'fs_LR*' spaces:
%       This input is not useful for input in 'fs_LR*' spaces (you can pass in 'NONE' or
%       empty string).
% 
%     - threshold
%       threshold for binarization (string, e.g. '0.1'). "threshold" = 0.1
%       means indices with the highest 10% correlations will be set to 1,
%       others will be set to 0.
% 
%     - varargin_file1
%       'fsaverage*' spaces:
%       the text file containing left hemisphere surface data list,
%       e.g. <path>/lh.Sub0033_Ses1.input;
%       Each line in this file corresponds to a single run of this subject. 
%       <varargin_file1> can also be a path which points to a single run if the input
%       is fMRI data of a single run. 
%       e.g. <path>/lh.Sub0033_Ses1.nii.gz.
%
%       'fs_LR*' spaces:
%       the text file containing entire cortical surface data list,
%       e.g. <path>/100206.dtseries_list.txt.
%       Each line in this file corresponds to a single run of this subject.
%       <varargin_file1> can also be a path which points to a single run if the input
%       is fMRI data of a single run.
%       e.g. <path>/Sub0033_Ses1.dtseries.nii.
% 
%     - varargin_file2
%       'fsaverage*' spaces:
%       the text file containing right hemisphere surface data list,
%       e.g. <path>/rh.Sub0033_Ses1.input;
%       Each line in this file corresponds to a single run of this subject. 
%       <varargin_file1> can also be a path which points to a single run if the input
%       is fMRI data of a single run. 
%       e.g. <path>/rh.Sub0033_Ses1.nii.gz.
%
%       'fs_LR*' spaces:
%       This input is not useful for input in 'fs_LR*' spaces (you can pass in 'NONE' or
%       empty string).
%
%     - outlier_text
%       the text file containing outlier file name, 
%       e.g. <path>/outlier.Sub0033_Ses1.input
%       Each line in this list corresponds to a single run of this subject.
%        <outlier_text> can also be a path which points to the outlier file of a single run 
%       if the input is fMRI data of a single run. 
%       e.g. <path>/outlier_Sub0033_Ses1_run1.txt
%       
%     - split_data
%       string or scalar. It can be '0' or other non-zero integer value.
%       If this flag is a non-zero value K, for subject with only one run, this function
%       will split it into K fake runs with equal length. For subject with multiple runs, 
%       it will not split any run.
%       If this flag is 0, whether the input subject is with one run or
%       multiple runs, this function will not split any run.
%       If this flag is 1, for subject with only one run, this function
%       will split it into 2 fake runs with equal length. For subject with multiple runs, 
%       it will not split any run.
%
%
% Written by CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md

% read target and seed meshs
if(~isempty(strfind(seed_mesh, 'fsaverage')))
    lh_seed_avg_mesh = CBIG_ReadNCAvgMesh('lh', seed_mesh, 'inflated', 'cortex');
    rh_seed_avg_mesh = CBIG_ReadNCAvgMesh('rh', seed_mesh, 'inflated', 'cortex');
    
    lh_avg_mesh = CBIG_ReadNCAvgMesh('lh', target, 'inflated', 'cortex');
    rh_avg_mesh = CBIG_ReadNCAvgMesh('rh', target, 'inflated', 'cortex');
elseif(~isempty(strfind(seed_mesh, 'fs_LR')))
    if(strcmp(seed_mesh, 'fs_LR_900') && strcmp(target, 'fs_LR_32k'))
        s_mesh = ft_read_cifti(fullfile(getenv('CBIG_CODE_DIR'), 'data', 'templates', 'surface', ...
            'fs_LR_32k_downsample_900', 'fslr_downsample_900mesh_parcellation.dlabel.nii'));
    else
        error('For fs_LR space, we only support seed_mesh = ''fs_LR_900'' and target = ''fs_LR_32k''.')
    end
else
    error('Unknown ''seed_mesh''.')
end

% if there is only one single run, the input file is the path to the fmri data instead of a txt file
if(contains(varargin_file1, '.nii') | contains(varargin_file1, '.mgh') | contains(varargin_file1, '.dtseries'))
    varargin1{1} = varargin_file1;
else

    % read in both left and right text files.
    fid = fopen(varargin_file1, 'r');
    i = 0;
    while(1);
    tmp = fscanf(fid, '%s\n', 1);
    if(isempty(tmp))
        break
    else
        i = i + 1;
        varargin1{i} = tmp;
    end
    end
    fclose(fid);
end

% read in both left and right text files.
% it is only applicable for data in fsaverage* spaces because for fs_LR
% space, left and right hemispheres are in the same file
if(~isempty(strfind(target, 'fsaverage')))
    if(contains(varargin_file2, '.nii') | contains(varargin_file2, '.mgh'))
        varargin2{1} = varargin_file2;
    else
        fid = fopen(varargin_file2, 'r');
        i = 0;
        while(1);
            tmp = fscanf(fid, '%s\n', 1);
            if(isempty(tmp))
                break
            else
                i = i + 1;
                varargin2{i} = tmp;
            end
        end
        fclose(fid);
    end
end

% read in outlier text files
if(contains(varargin_file1, '.nii') | contains(varargin_file1, '.mgh') | contains(varargin_file1, '.dtseries'))
    outlierin{1} = outlier_text;
else
    if(exist('outlier_text', 'var') && ~strcmp(outlier_text, 'NONE') && ~isempty(outlier_text))
        fid = fopen(outlier_text, 'r');
        i = 0;
        while(1);
            tmp = fscanf(fid, '%s\n', 1);
            if(isempty(tmp))
                break
            else
                i = i + 1;
                outlierin{i} = tmp;
            end
        end
        fclose(fid);
    end
end
if(~exist('split_data', 'var'))
    split_data = 0;
else
    if(ischar(split_data))
        split_data = str2num(split_data);
    end
    if(split_data == 1)
        fprintf('Split flag is set to be 1, will split data into 2 fake runs. \n');
        split_data = 2;
    end
end

% Compute profile for left hemi
for i = 1:length(varargin1)
    if(exist('outlierin', 'var'))
        % {0,1} vector, where uncensored time points are 0
        outliers = dlmread(outlierin{i}); 
    end
    
    input = varargin1{i};
    
    % read input file
    [input_series, t_series, input_size] = read_fmri(input);
    t_series = t_series';
    if(exist('outliers', 'var'))
        t_series = t_series(outliers==1, :);
    end
    
    lh_seed_series = t_series;
    if(~isempty(strfind(seed_mesh, 'fsaverage')))
        lh_seed_ind = (lh_avg_mesh.MARS_label(1:length(lh_seed_avg_mesh.MARS_label)) == 2);
        rh_seed_ind = (rh_avg_mesh.MARS_label(1:length(rh_seed_avg_mesh.MARS_label)) == 2);
        lh_seed_series = lh_seed_series(:, lh_seed_ind);
        
        hemi_index = strfind(input, basename(input));
        rh_seed_file = input; rh_seed_file(hemi_index:hemi_index+1) = 'rh';
        [~, rh_seed_series, ~] = read_fmri(rh_seed_file);
        rh_seed_series = rh_seed_series';
        rh_seed_series = rh_seed_series(:, rh_seed_ind);
        if(exist('outliers', 'var'))
            rh_seed_series = rh_seed_series(outliers==1, :);
        end
        s_series = [lh_seed_series rh_seed_series];
    else
        seed_ind = find(input_series.brainstructure==1 | input_series.brainstructure==2);
        seed_ind = seed_ind(s_mesh.dlabel~=0);
        s_series = lh_seed_series(:, seed_ind);
    end
    
    % normalize series (note that series are now of dimensions: T x N) and compute correlation
    if(split_data~=0 && length(varargin1)==1)
        for split = 1:split_data
            rang_start = floor(size(s_series,1)/split_data) * (split - 1) + 1;
            rang_end = floor(size(s_series,1)/split_data) * split;
            if(split == split_data)
                rang_end = size(s_series,1);
            end
            curr_s_series = s_series(rang_start:rang_end,:);
            curr_t_series = t_series(rang_start:rang_end,:);
            curr_corr_mat1 = CBIG_corr(curr_s_series, curr_t_series);

            disp(['Fake run ' num2str(split) ', isnan: ' num2str(sum(sum(isnan(curr_corr_mat1)))) ' out of '...
            num2str(numel(curr_corr_mat1))]);
            tmp_corr = curr_corr_mat1;
            tmp_corr(isnan(curr_corr_mat1)) = 0;
            corr_mat1(:,:,split) = tmp_corr;
        end
        if(i == 1)
            output = corr_mat1;
        else
            output = output + corr_mat1;
        end
    else        
        corr_mat1 = CBIG_corr(s_series, t_series);
        disp(['isnan: ' num2str(sum(isnan(corr_mat1(:)))) ' out of ' num2str(numel(corr_mat1))]);
        corr_mat1(isnan(corr_mat1)) = 0;
        if(i == 1)
            output = corr_mat1;
        else
            output = output + corr_mat1;
        end
    end
    clear outliers
end

output = output / length(varargin1);
corr_mat1 = output;


% Compute profile for right hemi
if(~isempty(strfind(target, 'fsaverage'))) 
    % if input data is in fs_LR space, this part is not needed.
    for i = 1:length(varargin2)
        if(exist('outlierin', 'var'))
            % {0,1} vector, where uncensored time points are 0
            outliers = dlmread(outlierin{i});        
        end
        
        input = varargin2{i};
        
        % read input file
        [input_series, t_series, input_size] = read_fmri(input);
        t_series = t_series';
        if(exist('outliers', 'var'))
            t_series = t_series(outliers==1, :);
        end
        
        hemi_index = strfind(input, basename(input));
        lh_seed_file = input; lh_seed_file(hemi_index:hemi_index+1) = 'lh';
        [~, lh_seed_series, ~] = read_fmri(lh_seed_file);
        lh_seed_series = lh_seed_series';
        lh_seed_series = lh_seed_series(:, lh_seed_ind);
        if(exist('outliers', 'var'))
            lh_seed_series = lh_seed_series(outliers==1, :);
        end
        
        rh_seed_series = t_series;
        rh_seed_series = rh_seed_series(:, rh_seed_ind);
        
        s_series = [lh_seed_series rh_seed_series];
        
        % normalize series (note that series are now of dimensions: T x N) and compute correlation
        if(split_data~=0 && length(varargin1)==1)
            for split = 1:split_data
                rang_start = floor(size(s_series,1)/split_data) * (split - 1) + 1;
                rang_end = floor(size(s_series,1)/split_data) * split;
                if(split == split_data)
                    rang_end = size(s_series,1);
                end
                curr_s_series = s_series(rang_start:rang_end,:);
                curr_t_series = t_series(rang_start:rang_end,:);
                curr_corr_mat2 = CBIG_corr(curr_s_series, curr_t_series);
    
                disp(['Fake run ' num2str(split) ', isnan: ' num2str(sum(sum(isnan(curr_corr_mat2)))) ' out of '...
                num2str(numel(curr_corr_mat2))]);
                tmp_corr = curr_corr_mat2;
                tmp_corr(isnan(curr_corr_mat2)) = 0;
                corr_mat2(:,:,split) = tmp_corr;
            end
            if(i == 1)
                output = corr_mat2;
            else
                output = output + corr_mat2;
            end
        else        
            corr_mat2 = CBIG_corr(s_series, t_series);
            disp(['isnan: ' num2str(sum(isnan(corr_mat2(:)))) ' out of ' num2str(numel(corr_mat2))]);
            corr_mat2(isnan(corr_mat2)) = 0;
            if(i == 1)
                output = corr_mat2;
            else
                output = output + corr_mat2;
            end
        end
        clear outliers
    end

    output = output / length(varargin2);
    corr_mat2 = output;
end

% combine both hemisphere and threshold
if(split_data~=0 && length(varargin1)==1)
    for run = 1:split_data
        if(~isempty(strfind(target, 'fsaverage')))
            tmp = [corr_mat1(:,:,run) corr_mat2(:,:,run)];
        else
            tmp = corr_mat1(:,:,run);
        end
        tmp = sort(tmp(:), 'descend');
        t = tmp(round(numel(tmp) * str2num(threshold)));
        disp(['threshold: ' num2str(t)]);
        if(str2num(threshold) < 1)
            tmp_corr = corr_mat1(:, :, run);
            tmp_corr(corr_mat1(:,:,run) < t) = 0;
            tmp_corr(corr_mat1(:,:,run) >= t) = 1;
            corr_mat1(:, :, run) = tmp_corr;
            % output_file1 either contains '.nii.gz' or '.mat'
            ind = [strfind(output_file1, '.nii.gz') strfind(output_file1, '.mat')];         
            output_file1_tmp = [output_file1(1:ind-1) '_' num2str(run) output_file1(ind:end)];
            write_fmri(output_file1_tmp, input_series, corr_mat1(:,:,run)',...
             [input_size(1:end-1) size(corr_mat1, 1)]);
            
            if(~isempty(strfind(target, 'fsaverage')))
                tmp_corr = corr_mat2(:, :, run);
                tmp_corr(corr_mat2(:,:,run) < t) = 0;
                tmp_corr(corr_mat2(:,:,run) >= t) = 1;
                corr_mat2(:, :, run) = tmp_corr;
                ind = [strfind(output_file2, '.nii.gz') strfind(output_file2, '.mat')];                                 
                output_file2_tmp = [output_file2(1:ind-1), '_' num2str(run) output_file2(ind:end)];
                write_fmri(output_file2_tmp, input_series, corr_mat2(:,:,run)',...
                 [input_size(1:end-1) size(corr_mat2, 1)]);
            end
        end
    end
else
    if(~isempty(strfind(target, 'fsaverage')))
        tmp = [corr_mat1 corr_mat2];
    else
        tmp = corr_mat1;
    end
    tmp = sort(tmp(:), 'descend');
    t = tmp(round(numel(tmp) * str2num(threshold)));
    disp(['threshold: ' num2str(t)]);
    if(str2num(threshold) < 1)
        corr_mat1(corr_mat1 <  t) = 0;
        corr_mat1(corr_mat1 >= t) = 1;
        write_fmri(output_file1, input_series, corr_mat1', [input_size(1:end-1) size(corr_mat1, 1)]);
        
        if(~isempty(strfind(target, 'fsaverage')))
            corr_mat2(corr_mat2 <  t) = 0;
            corr_mat2(corr_mat2 >= t) = 1;
            write_fmri(output_file2, input_series, corr_mat2', [input_size(1:end-1) size(corr_mat2, 1)]);
        end
    end
end




function [fmri, vol, vol_size] = read_fmri(fmri_name)

% [fmri, vol] = read_fmri(fmri_name)
% Given the name of functional MRI file (fmri_name), this function read in
% the fmri structure and the content of signals (vol).
% 
% Input:
%     - fmri_name:
%       The full path of input file name.
%
% Output:
%     - fmri:
%       The structure read in by MRIread() or ft_read_cifti(). To save
%       the memory, fmri.vol (for NIFTI) or fmri.dtseries (for CIFTI) is
%       set to be empty after it is transfered to "vol".
%
%     - vol:
%       A num_voxels x num_timepoints matrix which is the content of
%       fmri.vol (for NIFTI) or fmri.dtseries (for CIFTI) after reshape.
%
%     - vol_size:
%       The size of fmri.vol (NIFTI) or fmri.dtseries (CIFTI).

if (isempty(strfind(fmri_name, '.dtseries.nii')))
    % if input file is NIFTI file
    fmri = MRIread(fmri_name);
    vol = fmri.vol;
    vol_size = size(vol);
    if(length(vol_size) < 4)
        vol = reshape(vol, prod(vol_size(1:length(vol_size)-1)), vol_size(length(vol_size)));
    else
        vol = reshape(vol, prod(vol_size(1:3)), vol_size(4));
    end
    fmri.vol = [];
else
    % if input file is CIFTI file
    fmri = ft_read_cifti(fmri_name);
    vol = fmri.dtseries;
    vol = vol(fmri.brainstructure==1|fmri.brainstructure==2, :);
    vol_size = size(vol);
    fmri.dtseries = [];
end





function write_fmri(fmri_name, fmri, vol, vol_size)

% function write_fmri(fmri_name, fmri, vol)
% This function write out a fmri strucure (fmri) with signal content (vol)
% into fmri_name.
% 
% Input:
%     - fmri_name:
%       The output file name (full path).
%
%     - fmri:
%       The structure used for MRIwrite() if the output file is in .nii.gz format.
%       For output file in .mat format, this input will not be used.
%
%     - vol:
%       The content of output signals.
%
%     - vol_size:
%       The size that the volume will be reshaped to.

if(~isempty(strfind(fmri_name, '.nii.gz')))
    % if output file is NIFTI file
    vol = reshape(vol, vol_size);
    fmri.vol = vol;
    MRIwrite(fmri, fmri_name);
else
    % if output file is CIFTI file
    profile_mat = single(reshape(vol, vol_size));
    save(fmri_name, 'profile_mat', '-v7.3');
end


