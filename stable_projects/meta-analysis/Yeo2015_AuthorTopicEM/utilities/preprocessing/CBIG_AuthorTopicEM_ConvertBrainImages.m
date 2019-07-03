function CBIG_AuthorTopicEM_ConvertBrainImages(data_dir, ...
  gibbs_data_file, em_data_file, brain_mask2mm)
% CBIG_AuthorTopicEM_ConvertBrainImages(data_dir, brain_mask, ...
%   gibbs_data_file, em_data_file)
%
% Convert binary-smoothed brain images of activation foci to data input
% for the Gibbs sampler and Expectation-Maximization (EM) algorithm.
% The binary-smoothed brain images are produced by
% CBIG_AuthorTopicEM_PreprocessExpDataFromText.m
%
% Input:
%  - data_dir: path to the base directory containing the preprocessed data.
%      "<processedDir>/BinarySmoothedVolumes" contains the brain images of
%      activation foci smoothed with a binary smoothing kernel.
%      All brain images have 1mm resolution.
%  - gibbs_data_file: path to the input data file for the Gibbs sampler.
%  - em_data_file: path to the input data file for the EM algorithm.
%  - brain_mask2mm: a brain mask at 2mm resolution read from MRIread.
%
% Example:
%   data_dir = '/Work/data';
%   brain_mask2mm = MRIread(fullfile(getenv('CBIG_CODE_DIR'), 'stable_projects', ...
%     'meta-analysis', 'Yeo2015_AuthorTopicEM', 'utilities', ...
%     'mask', 'MNI_mask_conformed.2mm.0.1.nii.gz'));
%   gibbs_data_file = fullfile(data_dir, 'SelfGeneratedThought_GibbsData.mat')
%   em_data_file = fullfile(data_dir, 'SelfGeneratedThought_EMData.mat')
%   CBIG_AuthorTopicEM_ConvertBrainImagesToEMGibbsData(data_dir, ...
%     gibbs_data_file, em_data_file, brain_mask2mm)
%
%   Convert the smoothed brain images of activation foci under
%   '/Work/data/BinarySmoothedVolumes' to input data for Gibbs sampler,
%   saved at '/Work/data/SelfGeneratedThought_GibbsData.mat' and EM algorithm,
%   saved at '/Work/data/SelfGeneratedThought_EMData.mat'.
%
% Written by BT Thomas Yeo and CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md

  disp('Convert brain images to EM/Gibbs data format');

  % Read brain mask
  brain_index = brain_mask2mm.vol(:) == 1;

  experiments_file = fullfile(data_dir, 'ExperimentsData.mat');

  load(experiments_file)
  num_studies = length(experiments);

  % Extract unique tasks
  unique_tasks = {};
  unique_task_count = 0;
  for i = 1:length(exp_tasks)
    match = strcmp(exp_tasks{i}, unique_tasks);
    if (~any(match))
      unique_task_count = unique_task_count + 1;
      unique_tasks{unique_task_count} = cell2mat(exp_tasks{i});
    end
  end

  % Allocate memory
  act = zeros(num_studies, sum(brain_index));
  corpus = cell(num_studies, 1);
  task_by_exp = zeros(unique_task_count, num_studies);
  exp_mask = brain_mask2mm;
  exp_mask.vol = zeros(size(exp_mask.vol));

  % read task activation
  bin_volumes_dir = fullfile(data_dir, 'BinarySmoothedVolumes');
  for i = 1:num_studies
    if (rem(i, 100) == 1)
      disp(['  Exp #', num2str(i)]);
    end
    x = MRIread(fullfile(bin_volumes_dir, ['BinVolume' num2str(i, '%06d') '.nii.gz']));
    x = x.vol(1:2:end, 1:2:end, 1:2:end);
    act(i, :) = x(brain_index)';

    exp_task_indices = strcmp(unique_tasks, exp_tasks{i});
    task_by_exp(exp_task_indices, i) = 1;

    exp_word_count = act(i, act(i, :) ~= 0);
    corpus{i} = exp_word_count';

    exp_mask.vol(x ~= 0) = 1;
  end

  act = sparse(act);

  disp('Transforming to Gibbs Input');
  [WS, DS] = CBIG_AuthorTopicEM_ConvertActToGibbsFormat(act);
  paradigm_by_exp = sparse(double(task_by_exp));
  save(gibbs_data_file, 'WS', 'DS', 'paradigm_by_exp', '-v6');
  
  disp('Transforming to Collapsed Input for EM algorithm');
  [collapsed_act, collapsed_paradigm_by_exp] = CBIG_AuthorTopicEM_CollapseAct(act, paradigm_by_exp);
  corpus = CBIG_AuthorTopicEM_ConvertActToEMFormat(collapsed_act);
  paradigm_by_exp = full(logical(collapsed_paradigm_by_exp));
  corpus1 = corpus(:, 1);
  corpus2 = corpus(:, 2);
  corpus3 = corpus(:, 3);
  save(em_data_file, 'corpus1', 'corpus2', 'corpus3', 'paradigm_by_exp', '-v6');

  exp_mask_dir = fullfile(data_dir, 'mask');
  system(['mkdir -p ' exp_mask_dir]);
  exp_mask_file = fullfile(exp_mask_dir, 'exp_mask.nii.gz');
  MRIwrite(exp_mask, exp_mask_file);
