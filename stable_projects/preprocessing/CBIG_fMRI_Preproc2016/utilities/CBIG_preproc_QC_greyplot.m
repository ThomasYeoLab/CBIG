function CBIG_preproc_QC_greyplot( fmri_file, FD_file, DV_file, output, varargin )

% CBIG_preproc_QC_greyplot( fmri_file, gm_mask, FD_file, DV_file, output, varargin )
%
% Inputs:
%     - fmri_file:
%       The filename of fmri image (full path).
%       E.g. '<path-to-image>/Sub0001_Ses1_bld002_rest_skip4_stc_mc_residc_interp_FDRMS0.2_DVARS50_bp_0.009_0.08.nii.gz'
% 
%     - FD_file:
%       The filename of framewise displacement (full path, text file). E.g.
%       '<path-to-file>/Sub0001_Ses1_bld002_rest_skip4_stc_motion_outliers_FDRMS'
% 
%     - DV_file:
%       The filename of DVARS vector (full path, text file). E.g.
%       '<path-to-file>/Sub0001_Ses1_bld002_rest_skip4_stc_motion_outliers_DVARS'
% 
%     - output:
%       The filename of output figure (full path). E.g.
%       '<output_dir>/Sub0001_Ses1_bld0002_rest_skip4_stc_mc_greyplot.png'
% 
%     - varargin:
%       Some optional inputs can be passed in as:
%       'propertyname1', 'propertyvalue1', 'propertyname2',
%       'propertyvalue2', ...
%       See the table below:
%       ====================================================================
%        propertyname      | propertyvalue
%       --------------------------------------------------------------------
%        'GM_mask'         | If 'fmri_file' is a NIFTI volume, then the 
%                          | user needs to pass in a grey matter mask. If
%                          | 'fmri_file' is a CIFTI greyordinate file, then
%                          | GM_mask is not needed. As other mask inputs, 
%                          | it is the full path to the grey mattter mask.
%                          | E.g. '<path-to-mask>/Sub0001_Ses1.gm.func.nii.gz'
%       --------------------------------------------------------------------
%        'WB_mask'         | The filename of whole brain mask (full name).
%                          | If fmri_file is a CIFTI file, the user should
%                          | use 'WBS_text' option. 
%                          | E.g. '<path-to-mask>/Sub0001_Ses1.brainmask.bin.func.nii.gz'
%       --------------------------------------------------------------------
%        'WBS_text'        | The filename of global signal (full path, text 
%                          | file). This file contains one column vector 
%                          | with length of number of frames.
%                          | E.g. '<path-to-file>/Sub0001_Ses1_bld002_rest_skip4_stc_mc_residc_interp_FDRMS0.2_DVARS50_bp_0.009_0.08_WBS.txt'
%                          | This option is useful especially when
%                          | fmri_file is a CIFTI image, because the global
%                          | signal is usually computed outside this 
%                          | function by Workbench wb_command.
%                          | When fmri_file is a NIFTI image, users can
%                          | choose to pass in either WB_mask or WBS_text.
%                          | But if both are passed in, then global signal
%                          | will be computed from NIFTI voxels within the
%                          | mask and this text file will be ignored. 
%       --------------------------------------------------------------------
%        'grey_vox_factor' | Scalar or string. It is the factor to squeeze
%                          | all grey matter voxels in the plot. The height
%                          | of grey matter timeseries subfigure is the
%                          | number of voxels in grey matter divided by
%                          | grey_vox_factor. Default is 200. For Human
%                          | Connectome Project data, we suggest to use 400.
%       --------------------------------------------------------------------
%        'tp_factor'       | Scalar or string. It is the width of grey
%                          | matter timeseries subfigure is tp_factor
%                          | divided by the number of timepoints. Default
%                          | is 0.3. For the Human Connectome Project data,
%                          | we suggest to use 2.
%       --------------------------------------------------------------------
%        'FD_thres'        | Scalar or string. It is the threshold for FD
%                          | in censoring step. Default is 0.2.
%       --------------------------------------------------------------------
%        'DV_thres'        | Scalar or string. It is the threshold for
%                          | DVARS in censoring step. Default is 50.
%       ====================================================================
% 
% Example 1:
%     CBIG_preproc_QC_greyplot( ...
% 'xxx/Sub0001_Ses1/bold/002/Sub0001_Ses1_bld002_rest_skip4_stc_mc_residc_interp_FDRMS0.2_DVARS50_bp_0.009_0.08.nii.gz', ...
% 'xxx/Sub0001_Ses1/bold/mc/Sub0001_Ses1_bld002_rest_skip4_stc_motion_outliers_FDRMS', ...
% 'xxx/Sub0001_Ses1/bold/mc/Sub0001_Ses1_bld002_rest_skip4_stc_motion_outliers_DVARS', ...
% 'xxx/Sub0001_Ses1/qc/Sub0001_Ses1_bld002_rest_skip4_stc_mc_residc_interp_FDRMS0.2_DVARS50_bp_0.009_0.08_greyplot.png', ...
% 'GM_mask', 'xxx/Sub0001_Ses1/bold/mask/Sub0001_Ses1.func.gm.nii.gz'
% 'WB_mask', 'xxx/Sub0001_Ses1/bold/mask/Sub0001_Ses1.brainmask.bin.nii.gz')
% 
%     This example is suitable when the <fmri_file> is NIFTI format. In
%     this case, the grey matter mask is necessary.
% 
% Example 2:
%     CBIG_preproc_QC_greyplot( ...
% 'xxx/HCP_S1200_postprocessing/100206/MNINonLinear/Results/rfMRI_REST1_LR/postprocessing/MSM_reg_wbsgrayordinatecortex/rfMRI_REST1_LR_Atlas_MSMAll_hp2000_clean_regress.dtseries.nii', ...
% 'xxx/HCP/S1200/individuals/100206/MNINonLinear/Results/rfMRI_REST1_LR/Movement_RelativeRMS.txt', ...
% 'xxx/HCP_S1200_postprocessing/100206/MNINonLinear/Results/rfMRI_REST1_LR/postprocessing/QC/DVARS.txt', ...
% 'xxx/HCP/MSM_reg_wbsgrayordinatecortex/100206.rfMRI_REST1_LR.greyplot.png', ...
% 'WBS_text', 'xxx/HCP_S1200_postprocessing/100206/MNINonLinear/Results/rfMRI_REST1_LR/postprocessing/MSM_reg_wbsgrayordinatecortex/scripts/rfMRI_REST1_LR_WBS.txt', ...
% 'grey_vox_factor', 400, 'tp_factor', 0.5);
% 
%     This example is suitable when the <fmri_file> is CIFTI format. In this
%     case, Global signal is pre-computed and passed in by a text file.
% 
% Written by Jingwei Li and CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md

%% Parse input arguments
grey_vox_factor = 200;
tp_factor = 0.3;
FD_thres = 0.2;
DV_thres = 50;

for i = 1:2:length(varargin)
    if(strcmpi(varargin{i}, 'GM_mask'))
        gm_mask = varargin{i+1};
    elseif(strcmpi(varargin{i}, 'WBS_text'))
        wbs_text = varargin{i+1};
    elseif(strcmpi(varargin{i}, 'WB_mask'))
        wb_mask = varargin{i+1};
    elseif(strcmp(varargin{i}, 'grey_vox_factor'))
        grey_vox_factor = varargin{i+1};
        if(ischar(grey_vox_factor))
            grey_vox_factor = str2num(grey_vox_factor);
        end
    elseif(strcmp(varargin{i}, 'tp_factor'))
        tp_factor = varargin{i+1};
        if(ischar(tp_factor))
            tp_factor = str2num(tp_factor);
        end
    elseif(strcmpi(varargin{i}, 'FD_thres'))
        FD_thres = varargin{i+1};
        if(ischar(FD_thres))
            FD_thres = str2num(FD_thres);
        end
    elseif(strcmpi(varargin{i}, 'DV_thres'))
        DV_thres = varargin{i+1};
        if(ischar(DV_thres))
            DV_thres = str2num(DV_thres);
        end
    else
        fprintf(['ERROR: Unknown property name ' varargin{i} '.\n']);
        return
    end
end

if(isempty(strfind(fmri_file, '.dtseries.nii')))
    % When BOLD file is NIFTI, GM mask is needed.
    if(~exist('gm_mask', 'var'))
        fprintf('ERROR: Input BOLD file is a NIFTI file. Grey matter mask is needed.\n');
        return
    end
    
    % Whole brain signal: either a text file or a mask is needed
    if(~exist('wbs_text', 'var') && ~exist('wb_mask', 'var'))
        fprintf('ERROR: Neither whole brain signal nor whole brain mask is passed in.\n');
        return
    end
else
    % When BOLD file is CIFTI, WB signal need to be precomputed and passed in.
    if(~exist('wbs_text', 'var'))
        fprintf('ERROR: Input BOLD file is a CIFTI file. Whole brain signal should be precomputed, saved in a text file, and passed into this function.\n');
        return
    end
end



%% load FD and DV
FD = dlmread(FD_file);
DV = dlmread(DV_file);


%% read in BOLD data, extract gray matter voxels
if(isempty(strfind(fmri_file, '.dtseries.nii')))
    fmri = MRIread(fmri_file);
    fmri_gm = MRIread(gm_mask);
    gm_vol = fmri_gm.vol(:);
    vol = single(reshape(fmri.vol, size(fmri.vol, 1)*size(fmri.vol, 2)*size(fmri.vol, 3), size(fmri.vol, 4)));
    vol_gm = vol(logical(gm_vol)==1, :);
    fmri.vol = [];
    clear fmri_gm gm_vol
else
    fmri = ft_read_cifti(fmri_file);
    vol_gm = single(fmri.dtseries);
    fmri.dtseries = [];
    nan_ind = find(isnan(sum(vol_gm, 2)));
    vol_gm(nan_ind, :) = [];
end


%% Extract global signal
if(exist('wb_mask', 'var'))
    fmri_wb = MRIread(wb_mask);
    wb_vol = fmri_wb.vol(:);
    vol_wb = vol(logical(wb_vol)==1, :);
    wbs_mean = mean(vol_wb', 2);
else
    wbs_mean = dlmread(wbs_text);
end
avg_wbs_mean = mean(wbs_mean);
wbs_mean = wbs_mean - avg_wbs_mean;


%% Normalize grey matter signals
STD_gm = mean(std(vol_gm, [], 2));
vol_gm_zscore = bsxfun(@minus, vol_gm, mean(vol_gm, 2));
vol_gm_zscore = bsxfun(@rdivide, vol_gm_zscore, std(vol_gm, [], 2));
norm_gm = vol_gm_zscore * STD_gm;


%% Start plotting
tp = length(DV);
frames = 1:tp;
xtick_vec = int16(linspace(1, tp, 6));
grey_h = int16(size(norm_gm, 1)/grey_vox_factor);

fig = figure;
set(gcf, 'Visible', 'off')
set(gcf, 'Units', 'points');
set(gcf, 'Position', [0 0 tp/tp_factor+200 50+grey_h+25+55+25+55+25+55+25]);

% first panel: FD
subplot(4,1,1);
hold on
plot(frames, FD, 'r', 'LineWidth', 2);
plot(frames, FD_thres*ones(size(frames)), 'k', 'LineWidth', 2);
hold off
ax = subplot(4,1,1);
set(ax, 'Units', 'points', 'YTick', [0 0.3 0.6], 'YTickLabel', {'0', '0.3', '0.6'}, 'YLim', [0 0.6]);
set(ax, 'FontSize', 18);
set(ax, 'Position', [130 50+grey_h+25+55+25+55+25 tp/tp_factor 55]);
set(ax, 'Xtick', xtick_vec, 'XTickLabel', ''); 
xlim([1 tp]);
set(get(gca, 'YLabel'), 'String', '\color{red}FD');
set(get(gca, 'YLabel'), 'FontSize', 18);
set(gca, 'LineWidth', 2);
set(gca, 'TickDir', 'out', 'box', 'off');

% second panel: DVARS
subplot(4,1,2);
hold on
plot(frames, DV, 'b', 'LineWidth', 2);
plot(frames, DV_thres*ones(size(frames)), 'k', 'LineWidth', 2);
hold off
ax = subplot(4,1,2);
set(ax, 'Units', 'points', 'YTick', [0 40 80], 'YTickLabel', {'0', '40', '80'}, 'YLim', [0 80]);
set(ax, 'FontSize', 18);
set(ax, 'Position', [130 50+grey_h+25+55+25 tp/tp_factor 55]);
set(ax, 'Xtick', xtick_vec, 'XTickLabel', ''); 
xlim([1 tp]);
set(get(gca, 'YLabel'), 'String', '\color{blue}DV');
set(get(gca, 'YLabel'), 'FontSize', 18);
set(gca, 'LineWidth', 2);
set(gca, 'TickDir', 'out', 'box', 'off');

% third panel: GS
subplot(4,1,3);
box on
hold on
plot(frames, wbs_mean, 'k', 'LineWidth', 2);
hold off
ax = subplot(4,1,3);
set(ax, 'Units', 'points')
set(ax, 'FontSize', 18);
set(ax, 'Position', [130 50+grey_h+25 tp/tp_factor 55]);
set(ax, 'Xtick', xtick_vec, 'XTickLabel', ''); 
xlim([1 tp]);
set(get(gca, 'YLabel'), 'String', '\color{black}GS');
set(get(gca, 'YLabel'), 'FontSize', 18);
set(gca, 'LineWidth', 2);
set(gca, 'TickDir', 'out', 'box', 'off');

% grey matter
subplot(4,1,4);
imagesc(norm_gm, [-20 20]); col = colorbar;
colormap gray
set(col, 'FontSize', 18);
ax = subplot(4,1,4);
set(ax, 'Units', 'points')
set(ax, 'FontSize', 18);
set(ax, 'Position', [130 50 tp/tp_factor grey_h]);     % +52 is the offset of colobar
set(ax, 'Xtick', xtick_vec); 
set(ax, 'YTick', [])
xlim([1 tp]);
set(get(gca, 'YLabel'), 'String', 'Grey matter');
col_pos = get(col, 'Position');
col_pos(3) = col_pos(3) / 2;
set(col, 'Position', col_pos);
set(get(gca, 'YLabel'), 'FontSize', 18);
set(gca, 'TickDir', 'out');


% saving
set(gcf, 'PaperPositionMode', 'auto')

[output_dir, ~, ~] = fileparts(output);
if(~exist(output_dir, 'dir'))
    mkdir(output_dir);
end
print(fig, output, '-r300', '-dpng')
close(fig)

end

