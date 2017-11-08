function CBIG_preproc_CensorQC(QC_dir, subject, run_num, before_cen, interm_cen, after_cen, whole_mask_file, GM_mask_file, outlier_file)

% CBIG_preproc_CensorQC(QC_dir, subject, run_num, before_cen, interm_cen, after_cen, whole_mask_file, GM_mask_file, outlier_file)
%
% This function conducts QC steps for censoring results for one run
% (run_num) of one subject (subject), and save out QC results to QC_dir.
% The original volume file is given by before_cen. The interpolated volume
% is interm_cen. The final volume after censoring is after_cen. If bandpass
% filtering is performed in censoring, the final volume is the interpolated
% signals within passband. If bandpass filtering is not performed in
% censoring, the final volume is the interpolated signals with replacing
% uncesored frames by the original signals. For detailed information,
% please refer to CBIG_preproc_censor_wrapper.m
%
% It first computes correlation and fractional difference for each voxel
% between volumes before_cen and after_cen within whole brain mask (read
% from whole_mask_file) and within grey matter mask (read from
% GM_mask_file). The correlation and fractional difference are computed by
% CBIG_preproc_CensorCorrAdnFracDiff.m
%
% After that, it gives the 6 subplots of three timeseries of the voxels
% (ranking 0%, 20%, ..., 100%, sorted by correlation or fractional
% difference). The three timeseries are: original signal, intermediate
% signal, and the final signal (do both for whole brain and gray matter).
% This step needs to plot the position of outlier frames, which is read
% from outlier_file.
%
% Last, it saves the correlation and fractional difference (whole brain &
% gray matter) to nifti files. The statistics (mean, median, max, min) are
% saved in the corresponding text files.
%
% Input:
%     - QC_dir: 
%       The full path where all quality control results are
%       stored for current subject.
%       e.g. '<path_to_subject>/Sub0001_Ses1/qc'
%
%     - subject: 
%       A string of subject's id, e.g. 'Sub0001_Ses1'
%
%     - run_num: 
%       A string of run number, e.g. '002'
%
%     - before_cen: 
%       The nifti file name (full path) before censoring
%       e.g.
%       '<path_to_subject>/Sub0001_Ses1/bold/002/Sub0001_Ses1_bld002_rest_skip4_stc_mc.nii.gz'
%
%     - interm_cen:
%       The nifti file name (full path) of the interpolated volume,
%       which is the intermediate result
%       e.g.
%       '<path_to_subject>/Sub0001_Ses1/bold/002/Sub0001_Ses1_bld002_rest_skip4_stc_mc_interp_inter_FDRMS0.2_DVARS50.nii.gz'
%
%     - after_cen: 
%       The nifti file name (full path) of the final volume. 
%       e.g.
%       '<path_to_subject>/Sub0001_Ses1/bold/002/Sub0001_Ses1_bld002_rest_skip4_stc_mc_interp_FDRMS0.2_DVARS50.nii.gz'
%
%     - whole_mask_file: 
%       The file name (full path) of whole brain mask. e.g.
%       '<path_to_subject>/Sub0001_Ses1/bold/mask/Sub0001_Ses1.brainmask.bin.nii.gz'
%
%     - GM_mask_file: 
%       The file name (full path) of gray matter mask. e.g.
%       '<path_to_subject>/Sub0001_Ses1/bold/mask/Sub0001_Ses1.func.gm.nii.gz'
%
%     - outlier_file: 
%       The file name (full path) of outliers vector. e.g.
%       '<path_to_subject>/Sub0001_Ses1/qc/Sub0001_Ses1_bld002_FDRMS0.2_DVARS50_motion_outliers.txt'
%
% Example:
% CBIG_preproc_CensorQC('<path_to_subject>/Sub0001_Ses1/qc', 'Sub0001_Ses1', '002', before_cen, interm_cen, after_cen, whole_mask_file, GM_mask_file, outlier_file)
%
% Date: Jun.14, 2016
%
% Written by Jingwei Li.
% Written by CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Read data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[~, whole_mask_vol, vol_sz] = read_fmri(whole_mask_file);
[~, GM_mask_vol, ~] = read_fmri(GM_mask_file);

[~, vol1, ] = read_fmri(before_cen);
[~, vol_m, ~] = read_fmri(interm_cen);
[fMRI2, vol2, ~] = read_fmri(after_cen);
N = size(vol1, 1);
T = size(vol1, 2);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% compute correlation and fractional difference
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
corr = single(zeros(N, 1));
frac_diff = single(zeros(N, 1));

% divide voxels into branches, to reduce memory usage
if(T <= 150)
    voxbinsize = 4000;
elseif(T > 150 && T <= 500)
    voxbinsize = 1000;
elseif(T > 500 && T <= 1000)
    voxbinsize = 500;
elseif(T > 1000 && T <= 3000)
    voxbinsize = 200;
else
    voxbinsize = 100;
end
voxbin = 1:voxbinsize:N;
voxbin = [voxbin N+1];

for v = 1:length(voxbin)-1
    % compute correlation and fractional difference
    fprintf('Dealing with voxels from %d to %d ...\n', voxbin(v), voxbin(v+1)-1);
    [corr(voxbin(v):(voxbin(v+1)-1)), frac_diff(voxbin(v):(voxbin(v+1)-1))] = ...
        CBIG_preproc_CensorCorrAndFracDiff(vol1(voxbin(v):(voxbin(v+1)-1),:), vol2(voxbin(v):(voxbin(v+1)-1),:));
end

% mask correlation and framctional difference
corr_whole = corr;
corr_whole(whole_mask_vol==0) = 0;
frac_diff_whole = frac_diff;
frac_diff_whole(whole_mask_vol==0) = 0;

corr_GM = corr;
corr_GM(GM_mask_vol==0) = 0;
frac_diff_GM = frac_diff;
frac_diff_GM(GM_mask_vol==0) = 0;



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% plot ranking
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

outliers = dlmread(outlier_file);                % binary vector, 0 means censored (high-motion) frames.
outliers = ~outliers;                            % binary vector, 0 means censored (low-motion) frames.


%%%%%% Within GM mask
N_GM_inmask = length(find(GM_mask_vol~=0));

% pick voxels within GM mask
vol1_GM_inmask = vol1(GM_mask_vol~=0,:);
vol_m_GM_inmask = vol_m(GM_mask_vol~=0,:);
vol2_GM_inmask = vol2(GM_mask_vol~=0,:);

corr_GM_inmask = corr_GM(GM_mask_vol~=0);
frac_diff_GM_inmask = frac_diff_GM(GM_mask_vol~=0);

% rank by fractional difference and plot
CBIG_preproc_CensorQC_sort_plot(vol1_GM_inmask, vol_m_GM_inmask, vol2_GM_inmask, N_GM_inmask, outliers, ...
    frac_diff_GM_inmask, corr_GM_inmask, [QC_dir '/' subject '_bld' run_num '_interp_FracDiff_GM_6plots'], 'FRAC_DIFF');

% rank by correlation and plot
CBIG_preproc_CensorQC_sort_plot(vol1_GM_inmask, vol_m_GM_inmask, vol2_GM_inmask, N_GM_inmask, outliers, ...
    frac_diff_GM_inmask, corr_GM_inmask, [QC_dir '/' subject '_bld' run_num '_interp_corr_GM_6plots'], 'CORR');

clear N_GM_inmask vol1_GM_inmask vol2_GM_inmask corr_GM_inmask frac_diff_GM_inmask frac_diff_percent corr_percent
clear vol_m_GM_inmask



%%%%%% Within whole brain mask
N_whole_inmask = length(find(whole_mask_vol~=0));

% pick voxels within whole brain mask
vol1_whole_inmask = vol1(whole_mask_vol~=0,:);
vol_m_whole_inmask = vol_m(whole_mask_vol~=0,:);
vol2_whole_inmask = vol2(whole_mask_vol~=0,:);

corr_whole_inmask = corr_whole(whole_mask_vol~=0);
frac_diff_whole_inmask = frac_diff_whole(whole_mask_vol~=0);

% rank by fractional difference and plot
CBIG_preproc_CensorQC_sort_plot(vol1_whole_inmask, vol_m_whole_inmask, vol2_whole_inmask, N_whole_inmask, outliers, ...
    frac_diff_whole_inmask, corr_whole_inmask, [QC_dir '/' subject '_bld' run_num '_interp_FracDiff_whole_6plots'], 'FRAC_DIFF');

% rank by correlation and plot
CBIG_preproc_CensorQC_sort_plot(vol1_whole_inmask, vol_m_whole_inmask, vol2_whole_inmask, N_whole_inmask, outliers, ...
    frac_diff_whole_inmask, corr_whole_inmask, [QC_dir '/' subject '_bld' run_num '_interp_corr_whole_6plots'], 'CORR');

clear N_whole_inmask vol1_whole_inmask vol2_whole_inmask corr_whole_inmask frac_diff_whole_inmask frac_diff_percent corr_percent
clear vol_m_whole_inmask



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% reshape correlation and fractional difference and save them to nifiti
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% correlation, whole brain mask
median_corr = median(corr_whole(corr_whole~=0));
mean_corr = mean(corr_whole(corr_whole~=0));
max_corr = max(corr_whole(corr_whole~=0));
min_corr = min(corr_whole(corr_whole~=0));

% write into nifti file
write_fmri([QC_dir '/' subject '_bld' run_num '_interp_corr_whole.nii.gz'], fMRI2, corr_whole, vol_sz);

% write into text file
fid = fopen([QC_dir '/' subject '_bld' run_num '_interp_corr_whole.txt'], 'w');
fprintf(fid, '%s  %f\n', 'max', max_corr);
fprintf(fid, '%s  %f\n', 'min', min_corr);
fprintf(fid, '%s  %f\n', 'mean', mean_corr);
fprintf(fid, '%s  %f', 'median', median_corr);
fclose(fid);
clear corr_whole


%%% fractional difference, whole brain mask
median_frac = median(frac_diff_whole(frac_diff_whole~=0));
mean_frac = mean(frac_diff_whole(frac_diff_whole~=0));
max_frac = max(frac_diff_whole(frac_diff_whole~=0));
min_frac = min(frac_diff_whole(frac_diff_whole~=0));

% write into nifti file
write_fmri([QC_dir '/' subject '_bld' run_num '_interp_FracDiff_whole.nii.gz'], fMRI2, frac_diff_whole, vol_sz);

% write into text file
fid = fopen([QC_dir '/' subject '_bld' run_num '_interp_FracDiff_whole.txt'], 'w');
fprintf(fid, '%s  %f\n', 'max', max_frac);
fprintf(fid, '%s  %f\n', 'min', min_frac);
fprintf(fid, '%s  %f\n', 'mean', mean_frac);
fprintf(fid, '%s  %f', 'median', median_frac);
fclose(fid);
clear frac_diff_whole


%%% correlation, GM mask
median_corr = median(corr_GM(corr_GM~=0));
mean_corr = mean(corr_GM(corr_GM~=0));
max_corr = max(corr_GM(corr_GM~=0));
min_corr = min(corr_GM(corr_GM~=0));

% write into nifti file
write_fmri([QC_dir '/' subject '_bld' run_num '_interp_corr_GM.nii.gz'], fMRI2, corr_GM, vol_sz);

% write into text file
fid = fopen([QC_dir '/' subject '_bld' run_num '_interp_corr_GM.txt'], 'w');
fprintf(fid, '%s  %f\n', 'max', max_corr);
fprintf(fid, '%s  %f\n', 'min', min_corr);
fprintf(fid, '%s  %f\n', 'mean', mean_corr);
fprintf(fid, '%s  %f', 'median', median_corr);
fclose(fid);
clear corr_GM


%%% fractional difference, GM mask
median_frac = median(frac_diff_GM(frac_diff_GM~=0));
mean_frac = mean(frac_diff_GM(frac_diff_GM~=0));
max_frac = max(frac_diff_GM(frac_diff_GM~=0));
min_frac = min(frac_diff_GM(frac_diff_GM~=0));

% write into nifti file
write_fmri([QC_dir '/' subject '_bld' run_num '_interp_FracDiff_GM.nii.gz'], fMRI2, frac_diff_GM, vol_sz);

% write into text file
fid = fopen([QC_dir '/' subject '_bld' run_num '_interp_FracDiff_GM.txt'], 'w');
fprintf(fid, '%s  %f\n', 'max', max_frac);
fprintf(fid, '%s  %f\n', 'min', min_frac);
fprintf(fid, '%s  %f\n', 'mean', mean_frac);
fprintf(fid, '%s  %f', 'median', median_frac);
fclose(fid);
clear frac_diff_GM

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
    vol = single(fmri.vol);
    vol_size = size(vol);
    vol = reshape(vol, prod(vol_size(1:3)), prod(vol_size)/prod(vol_size(1:3)));
    fmri.vol = [];
else
    % if input file is CIFTI file
    fmri = ft_read_cifti(fmri_name);
    vol = single(fmri.dtseries);
    vol_size = size(vol);
    fmri.dtseries = [];
end


end



function write_fmri(fmri_name, fmri, vol, vol_size)

% function write_fmri(fmri_name, fmri, vol)
% This function write out a fmri strucure (fmri) with signal content (vol)
% into fmri_name.
% 
% Input:
%     - fmri_name:
%       The output fMRI file name (full path).
%
%     - fmri:
%       The structure for MRIwrite() or ft_write_cifti() to write out.
%
%     - vol:
%       The content of fMRI signals that need to be assgined to fmri.vol
%       (for NIFTI) or fmri.dtseries (for CIFTI).
%
%     - vol_size:
%       The size of volume when it was initially read in.

if(isempty(strfind(fmri_name, '.dtseries.nii')))
    % if output file is NIFTI file
    vol = reshape(vol, vol_size);
    fmri.vol = single(vol);
    MRIwrite(fmri, fmri_name);
else
    % if output file is CIFTI file
    vol = reshape(vol, vol_size);
    fmri.dtseries = single(vol);
    fmri_name = regexprep(fmri_name, '.dtseries.nii', '');
    ft_write_cifti(fmri_name, fmri, 'parameter', 'dtseries');
end


end
