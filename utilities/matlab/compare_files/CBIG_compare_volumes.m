function stats = CBIG_compare_volumes(vol_file1, vol_file2, stat_dir, mask)

% CBIG_compare_volumes(vol_file1, vol_file2, stat_dir, mask)
%
% Compare two NIfTI volumes voxel-wise. Auto-detects dimensionality:
%   3D (T=1): returns absolute difference stats only.
%   4D (T>1): additionally returns per-voxel temporal correlation stats.
% Optionally saves histograms, stat text files, and a correlation volume.
%
% Input:
%   vol_file1 : path to NIfTI volume 1
%   vol_file2 : path to NIfTI volume 2
%   stat_dir  : (optional) output directory for figures and stat files.
%               Pass [] to skip saving.
%   mask      : (optional) brain mask. Can be:
%                - omitted or empty: use all voxels
%                - char: path to a binary NIfTI mask file (non-zero = in mask)
%                - logical vector: pre-flattened voxel indices
%
% Output:
%   stats : struct with fields
%             max_diff, min_diff, median_diff, mean_diff
%             (4D only) max_corr, min_corr, median_corr, mean_corr,
%                       pct_99, pct_95, pct_90, n_nan
%
% Written by XUE Aihuiping and CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md

vol1 = MRIread(vol_file1);
vol2 = MRIread(vol_file2);
[X, Y, Z, T] = size(vol1.vol);

assert(isequal(size(vol1.vol), size(vol2.vol)), ...
    'Volume size mismatch: vol_file1 is [%s] but vol_file2 is [%s]', ...
    num2str(size(vol1.vol)), num2str(size(vol2.vol)));

vol1_m = reshape(vol1.vol, X*Y*Z, T);
vol2_m = reshape(vol2.vol, X*Y*Z, T);

mask_idx = [];
if nargin >= 4 && ~isempty(mask)
    if ischar(mask)
        mask_vol = MRIread(mask);
        mask_idx = logical(mask_vol.vol(:));
    else
        mask_idx = mask;
    end
    assert(numel(mask_idx) == X*Y*Z, ...
        'Mask length (%d) does not match volume spatial size (%d)', ...
        numel(mask_idx), X*Y*Z);
    vol1_m = vol1_m(mask_idx, :);
    vol2_m = vol2_m(mask_idx, :);
end

save_stats = nargin >= 3 && ~isempty(stat_dir);
if save_stats
    [~, fname, fext] = fileparts(vol_file1);
    if strcmp(fext, '.gz'); [~, fname] = fileparts(fname); end
    if ~exist(stat_dir, 'dir'); mkdir(stat_dir); end
end

valid_diff = abs(vol1_m(:) - vol2_m(:));
valid_diff = valid_diff(~isnan(valid_diff));
stats.max_diff    = max(valid_diff);
stats.min_diff    = min(valid_diff);
stats.median_diff = median(valid_diff);
stats.mean_diff   = mean(valid_diff);

disp(['  max diff: '    num2str(stats.max_diff) ...
      '  min diff: '    num2str(stats.min_diff) ...
      '  median diff: ' num2str(stats.median_diff) ...
      '  mean diff: '   num2str(stats.mean_diff)]);

if T == 1
    if save_stats
        fig = figure('Visible', 'off');
        hist(valid_diff);
        xlabel('|diff|'); ylabel('Count'); title('Absolute difference');
        saveas(fig, fullfile(stat_dir, [fname '_diff_hist.png']));
        close(fig);

        fid = fopen(fullfile(stat_dir, [fname '_diff_stat.txt']), 'w');
        fprintf(fid, 'max_diff: %f\nmin_diff: %f\nmedian_diff: %f\nmean_diff: %f\n', ...
            stats.max_diff, stats.min_diff, stats.median_diff, stats.mean_diff);
        fclose(fid);
    end
    return;
end

% Per-voxel temporal correlation: [Nbrain, 1]
% Equivalent to diag(CBIG_corr(vol1_m', vol2_m')) without Nbrain x Nbrain matrix.
vol1_dc  = vol1_m - mean(vol1_m, 2);
vol2_dc  = vol2_m - mean(vol2_m, 2);
vox_corr = sum(vol1_dc .* vol2_dc, 2) ./ ...
           sqrt(sum(vol1_dc.^2, 2) .* sum(vol2_dc.^2, 2));

valid_corr = vox_corr(~isnan(vox_corr));
stats.max_corr    = max(valid_corr);
stats.min_corr    = min(valid_corr);
stats.median_corr = median(valid_corr);
stats.mean_corr   = mean(valid_corr);
stats.pct_99      = mean(valid_corr > 0.99) * 100;
stats.pct_95      = mean(valid_corr > 0.95) * 100;
stats.pct_90      = mean(valid_corr > 0.90) * 100;
stats.n_nan       = sum(isnan(vox_corr));

disp(['  max corr: '    num2str(stats.max_corr) ...
      '  min corr: '    num2str(stats.min_corr) ...
      '  median corr: ' num2str(stats.median_corr) ...
      '  mean corr: '   num2str(stats.mean_corr)]);
disp(['  corr > 0.99: ' num2str(stats.pct_99, '%.1f') '%' ...
      '  corr > 0.95: ' num2str(stats.pct_95, '%.1f') '%' ...
      '  corr > 0.90: ' num2str(stats.pct_90, '%.1f') '%']);

if save_stats
    fig = figure('Visible', 'off');
    hist(valid_diff);
    xlabel('|diff|'); ylabel('Count'); title('Absolute difference');
    saveas(fig, fullfile(stat_dir, [fname '_diff_hist.png']));
    close(fig);

    fig = figure('Visible', 'off');
    hist(valid_corr, linspace(-1, 1, 50));
    xlim([-1 1]);
    xlabel('Correlation'); ylabel('Count'); title('Per-voxel temporal correlation');
    saveas(fig, fullfile(stat_dir, [fname '_corr_hist.png']));
    close(fig);

    % Correlation volume: 0 outside mask, corr values inside.
    corr_flat = zeros(X*Y*Z, 1);
    if ~isempty(mask_idx)
        corr_flat(mask_idx) = vox_corr;
    else
        corr_flat = vox_corr;
    end
    corr_mri = vol1;
    corr_mri.vol = reshape(corr_flat, X, Y, Z);
    MRIwrite(corr_mri, fullfile(stat_dir, [fname '_corr.nii.gz']));

    fid = fopen(fullfile(stat_dir, [fname '_diff_stat.txt']), 'w');
    fprintf(fid, 'max_diff: %f\nmin_diff: %f\nmedian_diff: %f\nmean_diff: %f\n', ...
        stats.max_diff, stats.min_diff, stats.median_diff, stats.mean_diff);
    fclose(fid);

    fid = fopen(fullfile(stat_dir, [fname '_corr_stat.txt']), 'w');
    fprintf(fid, 'max_corr: %f\nmin_corr: %f\nmedian_corr: %f\nmean_corr: %f\n', ...
        stats.max_corr, stats.min_corr, stats.median_corr, stats.mean_corr);
    fprintf(fid, 'pct_99: %.1f%%\npct_95: %.1f%%\npct_90: %.1f%%\nn_nan: %d\n', ...
        stats.pct_99, stats.pct_95, stats.pct_90, stats.n_nan);
    fclose(fid);
end
