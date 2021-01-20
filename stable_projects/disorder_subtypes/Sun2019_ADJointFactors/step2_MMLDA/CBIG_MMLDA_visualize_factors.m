function CBIG_MMLDA_visualize_factors(in_dir, out_dir, k, min_thresh, max_thresh, ...
    mask, behavior_name, behavior_domain, mni_space)
% CBIG_MMLDA_visualize_factors(in_dir, out_dir, k, min_thresh, max_thresh, mask)
%
% This function will visualize atrophy maps and cognitive deficits of estimated factors. 
% The atrophy maps will be visualized in MNI and fsaverage surface space.
% The cognitive deficits will be visualized by sorted bar plot. All steps are as following
% 1. Pick the initialization with maximum log-likelihood
% 2. Copy final results from MMLDA estimation folder to visualization output folder
% 3. Plot the climbing of log-likelihood as the best run iterates. The log-likelihood 
%    should "go flat" if MMLDA has converged
% 4. Plot the correlation of each run with best run and the runs are sorted based on
%    likelihood. If many initializations yield similar (high correlations) results, then it is enough.
% 5. Write Pr(Factor | Subject), i.e., normalized gamma in MMLDA, to txt, where each row is a subject.
% 6. Write Pr(Brain Atrophy | Factor), i.e., beta1 in MMLDA estimation, to nifti file.
% 7. Project brain from MNI space to fsaverage space and visualize it with Freesurfer.
% 8. Overlay input volume on underlay volumn (e.g., MNI template) with Hot colormap.
%    Then, take screeshots of different slices of the volume and save it as images.
%    Please make sure that your VNC resolution is 2000x1000.
% 9. Visualize behavior deficits of each factor with bar plot. 
%
% Input:
%   - in_dir            : the output directory you specified in CBIG_MMLDA_est.sh
%   - out_dir           : output directory
%   - k                 : number of factors/topics
%   - min_thresh        : minimum of the colorscale for visualizing atrophy maps 
%   - max_thresh        : maximum of the colorscale for visualizing atrophy maps
%   - mask              : path of grey matter mask generated from SPM VBM 
%   - behavior_name     : cell array of behavior names shown in the x tick label
%   - behavior_domain   : cell array of behavior domains, e.g., MEM, EF or NONE.
%   - mni_space         : (optional) Default 'MNI2mm'. MNI space of mask, e.g. MNI2mm, MNI1.5mm
%
% Examples:
% behavior_name = {'MMSE: recall', 'TMT: Part A', 'ADAS: Naming'};
% behavior_domain = {'MEM', 'EF', 'NONE'};
% CBIG_MMLDA_visualize_factors('~/storage/MMLDA_est_output', ...
% '~/storage/MMLDA_visualize', 3, 7.5e-6, 1.5e-5, '~/storage/SPM_VBM/gm_mask.nii.gz', ...
% behavior_name, behavior_domain)
%
% Written by Nanbo Sun and CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md

if nargin < 9
    mni_space = 'MNI2mm';
end

% Pick the initialization with maximum log-likelihood
r = CBIG_MMLDA_find_max_like_run(in_dir, out_dir, k);

% copy final results from MMLDA output folder to 
% visualization output folder
CBIG_MMLDA_copy_best_run(in_dir, out_dir, k, r);

% Plot the climbing of log-likelihood as this run runs
% The log-likelihood should "go flat" if MMLDA has converged 
CBIG_MMLDA_plot_like(out_dir, k, r);

% Check if the current number of initializations are enough
% If many initializations yield similar (high correlations) results, enough
CBIG_MMLDA_plot_corr_with_best(in_dir, out_dir, k, r);

% Write Pr(Factor | Subject), i.e., normalized gamma in MMLDA, to txt, where each row is a subject
CBIG_MMLDA_gamma2table(out_dir, k, r);

% Write Pr(Brain Atrophy | Factor), i.e., beta1 in MMLDA, to nifti file
CBIG_MMDLA_beta2brain(out_dir, k, r, mask);

% Overlay input volume on underlay volumn (e.g., MNI template) with colormap.
% Then, take screeshots of different slices of the volume and save it as images.
CBIG_MMLDA_visualize_atrophy_MNI(out_dir, k, r, min_thresh, max_thresh, mni_space);

% Visualize behavior deficits for each factor
CBIG_MMLDA_visualize_behavior(out_dir, k, r, behavior_name, behavior_domain)

% Project brain from MNI space to fsaverage space and
% visualize it with Freesurfer
CBIG_MMLDA_visualize_atrophy_fsaverage(out_dir, k, r, min_thresh, max_thresh);

close all
end


function r = CBIG_MMLDA_find_max_like_run(in_dir, out_dir, k)
% Find the run with maximum likelihood 
dir_list = dir([in_dir sprintf('/k%s/r*', num2str(k))]);
max_log_like = -inf;
% draw log likelihood - random initialization figure
figure;
hold on;
for idx = 1:numel(dir_list)
    run_name = dir_list(idx).name;
    try
        log_likes = load([in_dir sprintf('/k%s/', num2str(k)) run_name '/likelihood.dat']);
        log_like = log_likes(end, 1);
        run_id = str2num(run_name(2:end));
        scatter(run_id, log_like, 'b.');
        if log_like > max_log_like && log_like ~= 0
            max_log_like = log_like;
            max_run_id = run_id;
        end
    catch
        warning(['k=' num2str(k) ' ' run_name ': no likelihood.dat']);
    end
end
scatter(max_run_id, max_log_like, 'ro');
hold off;
ylabel('Log-likelihood');
xlabel('Random Initialization');
r = max_run_id;
mkdir([out_dir '/k' num2str(k) '/r' num2str(r)])
saveas(gcf, [out_dir '/k' num2str(k) '/r' num2str(r) '/Loglikelihood_Runs.png']);
end

function CBIG_MMLDA_copy_best_run(src_dir, target_dir, k, r)
% Copy best run from source dir to target dir.
mkdir([target_dir '/k' num2str(k) '/r' num2str(r) '/']);
cmd = ['rsync -az ' src_dir '/k' num2str(k) '/r' num2str(r) '/likelihood.dat ' ...
target_dir '/k' num2str(k) '/r' num2str(r) '/'];
system(cmd);
cmd = ['rsync -az ' src_dir '/k' num2str(k) '/r' num2str(r) '/*.log ' target_dir '/k' num2str(k) '/r' num2str(r) '/'];
system(cmd);
cmd = ['rsync -az ' src_dir '/k' num2str(k) '/r' ...
num2str(r) '/final.* ' target_dir '/k' num2str(k) '/r' num2str(r) '/'];
system(cmd);
end

function CBIG_MMLDA_plot_like(out_dir, k, r)
% Plot log likelihood vs. interations.
figure
log_likelihood = load([out_dir '/k' num2str(k) '/r' num2str(r) '/likelihood.dat']);
plot(log_likelihood(:, 1), 'o-');
xlabel('Iteration');
ylabel('Log-likelihood');
grid on;
saveas(gcf, [out_dir '/k' num2str(k) '/r' num2str(r) '/LogLikelihood_Iteration.png']);
end

function CBIG_MMLDA_plot_corr_with_best(in_dir, out_dir, k, r_bestLike)
% Plot correlation of each run with best run, and the runs are
% ordered from high likelihood to low likelihood
dir_list = dir([in_dir sprintf('/k%s/r*', num2str(k))]);
nRuns = numel(dir_list);
% Rank by likelihood
log_like = zeros(nRuns, 1);
for idx = 1:nRuns
    log_likes = load([in_dir sprintf('/k%s/r%s', num2str(k), num2str(idx)) '/likelihood.dat']);
    log_like(idx) = log_likes(end, 1);
end
[~, run_order] = sort(log_like, 'descend');
% Compute average correlation
ave_corr = zeros(nRuns, 1);
for idx = 1:nRuns
    ave_corr(idx) = CBIG_MMLDA_corr_hun_match(in_dir, k, r_bestLike, run_order(idx));
end
% Plot
figure;
plot(1:nRuns, ave_corr, '-o');
box off;
xlabel('Random Initializations');
ylabel('Correlation with the Best Run');
xlim([1, nRuns]);
ylim([0, 1]);
set(gca, 'YTick', 0:0.05:1);
grid on;
saveas(gcf, [out_dir '/k' num2str(k) '/r' num2str(r_bestLike) '/CorrWithBest.png']);
end

function ave_corr = CBIG_MMLDA_corr_hun_match(in_dir, k, r_best, r)
% First do Hungarian matching between factors of different runs.
% Compute average correlation across factors and modalities.
beta_best1 = exp(load([in_dir sprintf('/k%s/r%s', num2str(k), num2str(r_best)) '/final.beta1']));
beta1 = exp(load([in_dir sprintf('/k%s/r%s', num2str(k), num2str(r)) '/final.beta1']));
beta_best2 = exp(load([in_dir sprintf('/k%s/r%s', num2str(k), num2str(r_best)) '/final.beta2']));
beta2 = exp(load([in_dir sprintf('/k%s/r%s', num2str(k), num2str(r)) '/final.beta2']));
%%% Reorder subtypes to obtain the maximal correlation coefficients
% Construct the COST matrix
%         pos1 pos2 ...
% topic1
% topic2
% ...
costMat = zeros(k, k);
for rowIdx = 1:k
    for colIdx = 1:k
        % Assign beta (jobs, column) to beta_best (workers, row)
        corrMat = corrcoef(beta_best1(rowIdx, :)', beta1(colIdx, :)');
        costMat(rowIdx, colIdx) = 1-corrMat(1, 2);
    end
end
% Run the Hungarian matching algorithm
% order: each row (worker)'s matched column (job)
[order, ~] = munkres(costMat);
% Recompute the avergae correlation with sorted topics
corr1 = zeros(k, 1);
for idx = 1:k  
    corrMat = corrcoef(beta_best1(idx, :)', beta1(order(idx), :)');
    corr1(idx) = corrMat(1, 2); 
end
corr2 = zeros(k, 1);
for idx = 1:k  
    corrMat = corrcoef(beta_best2(idx, :)', beta2(order(idx), :)');
    corr2(idx) = corrMat(1, 2); 
end
ave_corr = mean([corr1; corr2]);
end

function CBIG_MMLDA_gamma2table(out_dir, k, r)
% Normalize gamma in MMLDA
gamma = load([out_dir '/k' num2str(k) '/r' num2str(r) '/final.gamma']);
gamma_norm = bsxfun(@times, gamma, 1./(sum(gamma, 2)));
dlmwrite([out_dir '/k' num2str(k) '/r' num2str(r) '/FactorComp.txt'], gamma_norm, 'delimiter', ' ');
end

function CBIG_MMDLA_beta2brain(out_dir, k, r, mask)
% Convert beta file to brain nifti file.
gm_mask = MRIread(mask);
gm_mask = gm_mask.vol;
mask_size = size(gm_mask);
gm_mask(gm_mask~=0) = 1;%binarize the mask
gm_mask_1d = reshape(gm_mask, [1 mask_size(1)*mask_size(2)*mask_size(3)]);
beta = load([out_dir '/k' num2str(k) '/r' num2str(r) '/final.beta1']);
% For each topic
for topic_idx = 1:size(beta, 1)
    beta_row = exp(beta(topic_idx, :));
    % Convert back to 3D
    beta_row_1d = gm_mask_1d; % insert 0's
    beta_row_1d(beta_row_1d==1) = beta_row; % 1's to real values
    beta_row_3d = reshape(beta_row_1d, [mask_size(1) mask_size(2) mask_size(3)]);
    % Write values to MRI
    dummy_mri = MRIread(mask);
    dummy_mri.vol = beta_row_3d;
    mri_filename = [out_dir '/k' num2str(k) '/r' num2str(r) '/topic' num2str(topic_idx) '.nii.gz'];
    MRIwrite(dummy_mri, mri_filename);
end
end

function CBIG_MMLDA_visualize_atrophy_fsaverage(out_dir, k, r, min_thresh, max_thresh)
% Project data from MNI to fsaverage space, overlay it with Thomas 17 network boundary 
% and visualize it with Freeview
CBIG_CODE_DIR = getenv('CBIG_CODE_DIR');
surface_template = 'fsaverage';
% absolute path to output directory
abs_path_to_output_dir = [out_dir '/k' num2str(k) '/r' num2str(r) '/fsaverage'];
mkdir(abs_path_to_output_dir)

for idx = 1:k
    topicMRI = [out_dir '/k' num2str(k) '/r' num2str(r) '/topic' num2str(idx) '.nii.gz'];
    % project data from MNI to fsaverage surface
    mri = MRIread(topicMRI);
    [lh_data, rh_data] = CBIG_ProjectMNI2fsaverage_Ants(mri, surface_template);
    % use fsaverage aparc.annot as reference annotation
    abs_path_to_lh_ref_annot = fullfile(getenv('FREESURFER_HOME'), 'subjects', ...
        surface_template, 'label', 'lh.aparc.annot');
    abs_path_to_rh_ref_annot = fullfile(getenv('FREESURFER_HOME'), 'subjects', ...
        surface_template, 'label', 'rh.aparc.annot');
    % the medial wall in fsaverage aparc.annot is labeled as 0. vertices with
    % this label in the visualization is painted as black.
    ref_medialwall_label = 0;
    % label used in all output files
    label = ['factor' num2str(idx)];
    % use matlab 'parula' for the visualization
    colorscheme = 'parula';
    % number of discrete values the original data is converted to
    discretization_res = 28;
    % visualize the data
    CBIG_DrawSurfaceDataAsAnnotation(lh_data, rh_data, ...
        abs_path_to_lh_ref_annot, abs_path_to_rh_ref_annot, ref_medialwall_label, ...
        surface_template, abs_path_to_output_dir, label, ...
        colorscheme, discretization_res, min_thresh, max_thresh);
    % Annotation file saved in previous example
    underlay_lh_annot_file = fullfile(abs_path_to_output_dir, 'tmp', ['lh.' label '.annot']);
    underlay_rh_annot_file = fullfile(abs_path_to_output_dir, 'tmp', ['rh.' label '.annot']);
    % Overlaying annotation files
    overlay_lh_annot_file = fullfile(CBIG_CODE_DIR, 'utilities', 'matlab', ...
        'figure_utilities', 'draw_surface_data_as_annotation', ...
        'fsaverage_parcel_outlines', 'lh.Yeo2011_17Networks_N1000_Boundary.annot');
    overlay_rh_annot_file = fullfile(CBIG_CODE_DIR, 'utilities', 'matlab', ...
        'figure_utilities', 'draw_surface_data_as_annotation', ...
        'fsaverage_parcel_outlines', 'rh.Yeo2011_17Networks_N1000_Boundary.annot');
    % Output annotation files
    combined_lh_annot_file = fullfile(abs_path_to_output_dir, 'tmp', ['lh.' label '.combined.annot']);
    combined_rh_annot_file = fullfile(abs_path_to_output_dir, 'tmp', ['rh.' label '.combined.annot']);
    % Combine the annotation files
    % This is the label assigned to the medial wall in underlay_lh_annot_file and underlay_rh_annot_file.
    % The medial wall's label was assigned in CBIG_AnnotateSingleHemiMedialWall in the previous example.
    % We want to the medial wall not to be overwritten
    ref_medialwall_label = 50 + 50 * 2^8 + 50 *2^16;
    CBIG_CombineSurfaceAnnotations(underlay_lh_annot_file, overlay_lh_annot_file, ...
        combined_lh_annot_file, ref_medialwall_label);
    CBIG_CombineSurfaceAnnotations(underlay_rh_annot_file, overlay_rh_annot_file, ...
        combined_rh_annot_file, ref_medialwall_label);
    % Visualize the combined annotation
    CBIG_VisualizeSurfaceAnnotationInFreeview(combined_lh_annot_file, combined_rh_annot_file, ...
        surface_template, [label '_with_17network_parcel_outline'], abs_path_to_output_dir);
end
end

function CBIG_MMLDA_visualize_atrophy_MNI(out_dir, k, r, min_thresh, max_thresh, mni_space)
% Visualize atrophy factors with different coronal slices.
CBIG_CODE_DIR = getenv('CBIG_CODE_DIR');
addpath([CBIG_CODE_DIR '/stable_projects/disorder_subtypes/Sun2019_ADJointFactors/utilities'])

if strcmp(mni_space, 'MNI2mm')
    SPM_MNI = [CBIG_CODE_DIR ...
    '/stable_projects/disorder_subtypes/Sun2019_ADJointFactors/step2_MMLDA/' ...
    'Template_T1_IXI555_MNI152_brain_MNI2mm.nii'];
elseif strcmp(mni_space, 'MNI1.5mm')
    SPM_MNI = [CBIG_CODE_DIR ...
    '/stable_projects/disorder_subtypes/Sun2019_ADJointFactors/step2_MMLDA/' ...
    'Template_T1_IXI555_MNI152_brain.nii.gz'];
end

output_dir = [out_dir '/k' num2str(k) '/r' num2str(r) '/MNI'];
mkdir(output_dir)
for idx = 1:k
    in_vol = [out_dir '/k' num2str(k) '/r' num2str(r) '/topic' num2str(idx) '.nii.gz'];
    underlay_vol = SPM_MNI;
    color_map_name = 'MingHot';
    plane = 'coronal';
    out_name = ['topic' num2str(idx)];
    CBIG_MMLDA_view_vol_slice_with_underlay(in_vol, underlay_vol, ...
        mni_space, color_map_name, plane, min_thresh, max_thresh, output_dir, out_name)
end
rmpath([CBIG_CODE_DIR '/stable_projects/disorder_subtypes/Sun2019_ADJointFactors/utilities'])
end

function CBIG_MMLDA_visualize_behavior(out_dir, k, r, behavior_name, behavior_domain)
% Visualize behavior deficits for each factor.
beta2 = exp(load([out_dir '/k' num2str(k) '/r' num2str(r) '/final.beta2']));
behavior_factors = beta2';
% plot all behaivor scors for experiment
for i = 1:size(behavior_factors, 2)
    [behavior_factor_sort, ind] = sort(behavior_factors(:, i), 'descend');
    behavior_domain_sort = behavior_domain(ind);
    figure   
    set(gcf,'Position',[0, 0, 1800, 1000])
    hold on
    for j = 1:size(behavior_factors, 1)
        h = bar(j, behavior_factor_sort(j), 'barwidth', 0.8);
        if strcmp(behavior_domain_sort{j}, 'MEM')
            set(h, 'FaceColor', 'r')
        elseif strcmp(behavior_domain_sort{j}, 'EF')
            set(h, 'FaceColor', 'b')
        else
            set(h, 'FaceColor', 'k')
        end
    end
    set(gca, 'xtick', 1:1:1*size(behavior_factors, 1), 'xticklabel', behavior_name(ind), 'fontsize', 20)
    rotateXLabels(gca(), 45)
    ylim([min(abs(behavior_factors(:))) max(abs(behavior_factors(:)))])
    ylabel('Probability')

    box off
    hgexport(gcf, [out_dir '/k' num2str(k) '/r' num2str(r) '/beta2_F' num2str(i)])
    eps2xxx([out_dir '/k' num2str(k) '/r' num2str(r) '/beta2_F' num2str(i) '.eps'], {'png'})
    hold off
end
% plot top 15 scores for paper
N = 15;
for i = 1:size(behavior_factors, 2)
    [behavior_factor_sort, ind] = sort(behavior_factors(:, i), 'descend');
    behavior_domain_sort = behavior_domain(ind);
    behavior_factor_sort(behavior_factor_sort < 0.03) = 0.030000000001; 
    %% Draw figures
    figure   
    set(gcf,'Position',[0, 0, 1500, 600]);
    hold on
    for j = 1:N
        h(j) = bar(j, behavior_factor_sort(j), 'barwidth', 0.8);
        if strcmp(behavior_domain_sort{j}, 'MEM')
            set(h(j), 'FaceColor', 'r')
        elseif strcmp(behavior_domain_sort{j}, 'EF')
            set(h(j), 'FaceColor', 'b')
        else
            set(h(j), 'FaceColor', 'k')
        end
    end
    
    set(gca, 'xtick', 1:1:1*N, 'xticklabel', behavior_name(ind), 'fontsize', 26)
    set(gca, 'TickDir', 'out')
    if ~isempty(strfind(out_dir, 'ADNI1'))
        ylim([0.025 0.08501])
        set(gca, 'ytick', [0.025 0.045 0.065 0.085]) 
    else
        switch k
            case 3
                ylim([0.03 0.12])
                set(gca, 'ytick', [0.03 0.06 0.09 0.12])
            case 2 
                ylim([0.03 0.09])
                set(gca, 'ytick', [0.03 0.05 0.07 0.09])
        end
    end
    rotateXLabels(gca(), 45)
    set(gca, 'xtick', [])
    ylabel('Pr(Score | Factor)')
    
    ind_mem = find(strcmp(behavior_domain_sort, 'MEM')==1);
    ind_ef = find(strcmp(behavior_domain_sort, 'EF')==1);
    ind_un = find(strcmp(behavior_domain_sort, 'NONE')==1);
    legend([h(ind_mem(1)) h(ind_ef(1)) h(ind_un(1))], ... 
        'MEM-related scores (Crane2012)', ...
        'EF-related scores (Gibbons2012)', ...
        'Others')
    legend('boxoff')

    box off
    hgexport(gcf, [out_dir '/k' num2str(k) '/r' num2str(r) '/beta2_F' num2str(i) '_Top' num2str(N)])
    eps2xxx([out_dir '/k' num2str(k) '/r' num2str(r) '/beta2_F' num2str(i) '_Top' num2str(N) '.eps'], {'png'})
    hold off
end
close all
end
