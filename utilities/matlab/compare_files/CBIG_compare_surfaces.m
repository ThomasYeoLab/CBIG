function stats = CBIG_compare_surfaces(lh_file1, rh_file1, lh_file2, rh_file2, stat_dir, mesh)

% CBIG_compare_surfaces(lh_file1, rh_file1, lh_file2, rh_file2, stat_dir, mesh)
%
% Compare two surface time series (lh + rh NIfTI files) vertex-wise.
% Auto-detects dimensionality:
%   3D (T=1): returns absolute difference stats only.
%   4D (T>1): additionally returns per-vertex temporal correlation stats.
% Optionally saves surface correlation maps via CBIG_DrawSurfaceMaps,
% histograms, and stat text files.
%
% Input:
%   lh_file1, rh_file1 : left/right hemisphere surface NIfTI for set 1
%   lh_file2, rh_file2 : left/right hemisphere surface NIfTI for set 2
%   stat_dir           : (optional) output directory for figures and stat files.
%                        Pass [] to skip saving.
%   mesh               : (optional) surface mesh name for CBIG_DrawSurfaceMaps
%                        (default: 'fsaverage5')
%
% Output:
%   stats : struct with fields (combined lh+rh)
%             max_diff, min_diff, median_diff, mean_diff
%             (4D only) max_corr, min_corr, median_corr, mean_corr,
%                       pct_99, pct_95, pct_90, n_nan
%
% Written by XUE Aihuiping and CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md

if nargin < 6 || isempty(mesh)
    mesh = 'fsaverage5';
end

lh1 = MRIread(lh_file1); rh1 = MRIread(rh_file1);
lh2 = MRIread(lh_file2); rh2 = MRIread(rh_file2);

% Flatten to [N_vertices, T]
lh1_m = reshape(lh1.vol, [], size(lh1.vol, 4));
rh1_m = reshape(rh1.vol, [], size(rh1.vol, 4));
lh2_m = reshape(lh2.vol, [], size(lh2.vol, 4));
rh2_m = reshape(rh2.vol, [], size(rh2.vol, 4));

assert(isequal(size(lh1_m), size(lh2_m)), ...
    'LH surface size mismatch: file1 is [%s] but file2 is [%s]', ...
    num2str(size(lh1_m)), num2str(size(lh2_m)));
assert(isequal(size(rh1_m), size(rh2_m)), ...
    'RH surface size mismatch: file1 is [%s] but file2 is [%s]', ...
    num2str(size(rh1_m)), num2str(size(rh2_m)));
T = size(lh1_m, 2);

save_stats = nargin >= 5 && ~isempty(stat_dir);
if save_stats
    [~, fname, fext] = fileparts(lh_file1);
    if strcmp(fext, '.gz'); [~, fname] = fileparts(fname); end
    if strncmp(fname, 'lh.', 3); fname = fname(4:end); end
    if ~exist(stat_dir, 'dir'); mkdir(stat_dir); end
end

all1_m = [lh1_m; rh1_m];
all2_m = [lh2_m; rh2_m];
valid_diff = abs(all1_m(:) - all2_m(:));
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

% Per-vertex temporal correlation
lh1_dc = lh1_m - mean(lh1_m, 2);
rh1_dc = rh1_m - mean(rh1_m, 2);
lh2_dc = lh2_m - mean(lh2_m, 2);
rh2_dc = rh2_m - mean(rh2_m, 2);
corr_lh = sum(lh1_dc .* lh2_dc, 2) ./ sqrt(sum(lh1_dc.^2, 2) .* sum(lh2_dc.^2, 2));
corr_rh = sum(rh1_dc .* rh2_dc, 2) ./ sqrt(sum(rh1_dc.^2, 2) .* sum(rh2_dc.^2, 2));

all_corr   = [corr_lh; corr_rh];
valid_corr = all_corr(~isnan(all_corr));
stats.max_corr    = max(valid_corr);
stats.min_corr    = min(valid_corr);
stats.median_corr = median(valid_corr);
stats.mean_corr   = mean(valid_corr);
stats.pct_99      = mean(valid_corr > 0.99) * 100;
stats.pct_95      = mean(valid_corr > 0.95) * 100;
stats.pct_90      = mean(valid_corr > 0.90) * 100;
stats.n_nan       = sum(isnan(all_corr));

disp(['  max corr: '    num2str(stats.max_corr) ...
      '  min corr: '    num2str(stats.min_corr) ...
      '  median corr: ' num2str(stats.median_corr) ...
      '  mean corr: '   num2str(stats.mean_corr)]);
disp(['  corr > 0.99: ' num2str(stats.pct_99, '%.1f') '%' ...
      '  corr > 0.95: ' num2str(stats.pct_95, '%.1f') '%' ...
      '  corr > 0.90: ' num2str(stats.pct_90, '%.1f') '%']);

if save_stats
    corr_lh(isnan(corr_lh)) = 0;
    corr_rh(isnan(corr_rh)) = 0;
    CBIG_DrawSurfaceMaps(corr_lh, corr_rh, mesh, 'inflated', 0, 1);
    set(gcf, 'Visible', 'off');
    saveas(gcf, fullfile(stat_dir, [fname '_corr_surf.png']));
    close(gcf);

    fig = figure('Visible', 'off');
    hist(valid_corr, linspace(-1, 1, 50));
    xlim([-1 1]);
    xlabel('Correlation'); ylabel('Count'); title('Per-vertex temporal correlation');
    saveas(fig, fullfile(stat_dir, [fname '_corr_hist.png']));
    close(fig);

    fig = figure('Visible', 'off');
    hist(valid_diff);
    xlabel('|diff|'); ylabel('Count'); title('Absolute difference');
    saveas(fig, fullfile(stat_dir, [fname '_diff_hist.png']));
    close(fig);

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
