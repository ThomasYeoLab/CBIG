function CBIG_preproc_CensorQC_sort_plot(vol1_inmask, vol_m_inmask, vol2_inmask, outliers,...
    FRAC_DIFF_inmask, CORR_inmask, outname_noext, measure_name)

% CBIG_preproc_CensorQC_sort_plot(vol1_inmask, vol_m_inmask, vol2_inmask,
% N_inmask, outliers, FRAC_DIFF_inmask, CORR_inmask, outname_noext, measure_name)
% 
% This function sorts voxels by the values of measure_name ('CORR' -
% correlation or 'FRAC_DIFF' - fractional difference). FRAC_DIFF_inmask is
% a matrix contains the unsorted values of fractional difference within a
% centain mask, like whole brain mask or grey matter mask. CORR_inmask is a
% matrix contains the unsorted values of correlation within the same mask.
% If measure_name == 'FRAC_DIFF', this function will sort FRAC_DIFF_inmask
% in ascending direction; if measure_name == 'CORR', this function will
% sort CORR_inmask in ascending direction. Then it draws the subplots of
% three timeseries of 6 voxels (ranking 0%, 20%, ..., 100%, based on the
% sorted measure). The three timeseries are: original signal, intermediate
% signal, and the final signal. For the description of these three types of
% signal, please refer to CBIG_preproc_CensorQC.m and
% CBIG_preproc_censor_wrapper.m. This step needs to plot the position of
% outlier frames, which is indicated by a vector outliers.
%
% Input:
%     - vol1_inmask: 
%       A num_voxels x num_timepoints matrix. It is the original signal
%       within a mask (whole brain mask or gray matter mask).
%
%     - vol_m_inmask: 
%       A num_voxels x num_timepoints matrix. It is the interpolated
%       intermediate volume within a mask (whole brain mask or gray matter
%       mask).
%
%     - vol2_inmask: 
%       A num_voxels x num_timepoints matrix. It is the final volume within
%       a mask (whole brain mask or gray matter mask) after censoring.
%
%     - outliers: 
%       A num_timepoints x 1 binary vector, where 1 indicates "bad" time
%       points (outliers) and 0 indicates "good" time points.
%
%     - FRAC_DIFF_inmask: 
%       a vector within length of num_voxels. It contains the unsorted
%       measurement of fractional difference within the mask.
%
%     - CORR_inmask: 
%       a vector with length of num_voxels. It contains the unsorted value
%       of correlation within the mask.
%
%     - outname_noext: 
%       The output name (full path) for jpg files without extension. e.g.
%       'subject_dir/Sub001_Ses1/qc/Sub0001_Ses1_bld002_interp_FracDiff_whole_6plots'
%
%     - measure_name: 
%       'CORR' or 'FRAC_DIFF', determines if measure1 is correlation or
%       fractional difference. The voxels are sorted according to
%       measure_name.
%
% Example:
% CBIG_preproc_CensorQC_sort_plot(vol1_inmask, vol_m_inmask, vol2_inmask, N_inmask, outliers, ...
%     FRAC_DIFF_inmask, CORR_inmask,
%     'subject_dir/Sub001_Ses1/qc/Sub0001_Ses1_bld002_interp_FracDiff_whole_6plots', 'FRAC_DIFF') 
% CBIG_preproc_CensorQC_sort_plot(vol1_inmask, vol_m_inmask, vol2_inmask, N_inmask, outliers, ...
%     FRAC_DIFF_inmask, CORR_inmask,
%     'subject_dir/Sub001_Ses1/qc/Sub0001_Ses1_bld002_interp_FracDiff_whole_6plots', 'CORR') 
%
% Date: Jun.14, 2016
%
% Written by Jingwei Li.
% Written by CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md

%% Check input variables
if size(outliers,2) > 1
    error('Input argument ''outliers'' should be a column vector');
end

% Remove NaN
nan_index = isnan(FRAC_DIFF_inmask) | isnan(CORR_inmask);
N_inmask = sum(nan_index == 0);
FRAC_DIFF_inmask(nan_index) = [];
CORR_inmask(nan_index) = [];
vol1_inmask(nan_index, :) = [];
vol_m_inmask(nan_index, :) = [];
vol2_inmask(nan_index, :) = [];

if(strcmp(measure_name, 'FRAC_DIFF'))
    [measure_sorted, measure_ind] = sort(FRAC_DIFF_inmask, 'ascend');
elseif(strcmp(measure_name, 'CORR'))
    [measure_sorted, measure_ind] = sort(CORR_inmask, 'ascend');
else
    error('ERROR: Unrecognized measure name. Please use FRAC_DIFF or CORR as measure name.\n')
end

% find the position of outliers
outlier_1 = find(outliers==1);

% position for 6 subplots
pos_set = [0.05 0.70 0.42 0.25; ...
    0.55 0.70 0.42 0.25; ...
    0.05 0.37 0.42 0.25; ...
    0.55 0.37 0.42 0.25; ...
    0.05 0.04 0.42 0.25; ...
    0.55 0.04 0.42 0.25];

% start to plot
measure1_percent = zeros(11,1);
figure;
set(gcf, 'Visible', 'off')
for i = 1:6
    % find corresponding voxel for 0%, 20%, 40%, 60%, 80%, 100% ranking
    percentage = (i-1)/5;
    n = ceil(percentage * N_inmask);
    if(n==0)
        n = 1;
    end
    
    measure1_percent(i) = measure_sorted(n);
    
    % plot the three BOLD signals: original, interpolated, final
    g = subplot(3,2,i);
    hold on
    f1 = plot(1:length(outliers), vol1_inmask(measure_ind(n),:), 'b', 'LineWidth', 2);
    f2 = plot(1:length(outliers), vol_m_inmask(measure_ind(n),:), 'g', 'LineWidth', 2);
    f3 = plot(1:length(outliers), vol2_inmask(measure_ind(n),:), 'k--', 'LineWidth', 2);
    set(gca, 'XTick', 0:floor(length(outliers)/10):length(outliers));
    set(gca, 'FontSize', 15)
    set(gcf, 'Position', [0,0,1600,1000]);
    set(g, 'Position', pos_set(i,:));
    set(g, 'XLim', [0, length(outliers)]);
    
    pause(1)
    yl = ylim;
    set(gca, 'YLim', [yl(1) yl(2)+0.5*(yl(2)-yl(1))]);
    MinorTicks = 0:floor(length(outliers)/50):(length(outliers));
    plot(repmat(MinorTicks, 2, 1), repmat([yl(1); (yl(2)-yl(1))*1.5/50 + yl(1)], 1, length(MinorTicks)), 'k');
    
    % plot red vertical lines to indicate removed frames
    plot(repmat(outlier_1', 2, 1), repmat([yl(1); yl(2)+0.5*(yl(2)-yl(1))], 1, length(outlier_1)), 'r', 'LineWidth', 1);
    hold off
    
    l = legend([f1, f2, f3], 'Before Censoring', 'Intermediate', 'Censoring Final', 'Location', 'northwest');
    set(l, 'FontSize', 15);
    
    % title
    if(strcmp(measure_name, 'FRAC_DIFF'))
        ti = title(['Percentile:' num2str(percentage*100) '%, Vox:' num2str(measure_ind(n)) ...
            ', Corr:' num2str(CORR_inmask(measure_ind(n))) ', FracDiff:' ...
            num2str(FRAC_DIFF_inmask(measure_ind(n)), '%.3e')]);
    elseif(strcmp(measure_name, 'CORR'))
        ti = title(['Percentile:' num2str(percentage*100) '%, Vox:' num2str(measure_ind(n)) ...
            ', Corr:' num2str(CORR_inmask(measure_ind(n))) ', FracDiff:' ...
            num2str(FRAC_DIFF_inmask(measure_ind(n)), '%.3e')]);
    else
        warning('No matched measure name!')
    end
    set(ti, 'FontSize', 15)
    
end
% hgexport(gcf, outname_noext)
set(gcf, 'PaperPositionMode', 'auto')
saveas(gcf, [outname_noext '.png'])
close(gcf)


