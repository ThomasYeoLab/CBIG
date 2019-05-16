function CBIG_AuthorTopic_ConvertBrainImagesToCVBData(dataDir, cvbDataFile)
% CBIG_AuthorTopic_ConvertBrainImagesToCVBData(dataDir, cvbDataFile)
%
% Convert binary-smoothed brain images of activation foci to data input
% for the Collapsed Variational Bayes (CVB) algorithm.
% The binary-smoothed brain images are produced by
% CBIG_AuthorTopic_PreprocessExpDataFromText.m
%
% Input:
%  - dataDir: path to the base directory containing the preprocessed data.
%      "<processedDir>/BinarySmoothedVolumes" contains the brain images of
%      activation foci smoothed with a binary smoothing kernel.
%      All brain images have 1mm resolution.
%  - cvbDataFile: path to the input data file for the CVB algorithm.
%
% Example:
%   dataDir = '/Work/data';
%   cvbDataFile = fullfile(dataDir, 'SelfGeneratedThought_CVBData.mat')
%   CBIG_AuthorTopic_ConvertBrainImagesToCVBData(dataDir, cvbDataFile)
%
%   Convert the smoothed brain images of activation foci under
%   '/Work/data/BinarySmoothedVolumes' to input data for CVB algorithm,
%   saved at '/Work/data/SelfGeneratedThought_CVBData.mat'.
%
% Written by Gia H. Ngo and CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md

  disp('Converting brain images to CVB input format');

  % Read brain mask
  baseBrainMaskPath = fullfile(getenv('CBIG_CODE_DIR'), 'stable_projects', 'meta-analysis', 'Ngo2019_AuthorTopic', 'utilities', 'mask', 'MNI_mask_conformed.2mm.0.1.nii.gz');
  brainMask = MRIread(baseBrainMaskPath);
  brainIndex = brainMask.vol(:) == 1;

  experimentsFile = fullfile(dataDir, 'ExperimentsData.mat');

  load(experimentsFile)
  numStudies = length(Experiments);

  % Extract unique tasks
  uniqueTasks = {};
  uniqueTaskCount = 0;
  for i = 1:length(expTasks)
    match = strcmp(expTasks{i}, uniqueTasks);
    if (~any(match))
      uniqueTaskCount = uniqueTaskCount + 1;
      uniqueTasks{uniqueTaskCount} = cell2mat(expTasks{i});
    end
  end

  % Allocate memory
  act = zeros(numStudies, sum(brainIndex));
  corpus = cell(numStudies, 1);
  taskByExp = zeros(uniqueTaskCount, numStudies);
  expMask = brainMask;
  expMask.vol = zeros(size(expMask.vol));

  % read task activation
  binVolumesDir = fullfile(dataDir, 'BinarySmoothedVolumes');
  for i = 1:numStudies
    x = MRIread(fullfile(binVolumesDir, ['BinVolume' num2str(i, '%06d') '.nii.gz']));
    x = x.vol(1:2:end, 1:2:end, 1:2:end);
    act(i, :) = x(brainIndex)';

    expTaskIndices = strcmp(uniqueTasks, expTasks{i});
    taskByExp(expTaskIndices, i) = 1;

    expWordCount = act(i, act(i, :) ~= 0);
    corpus{i} = expWordCount';

    expMask.vol(x ~= 0) = 1;
  end

  outlyingVoxels = ((brainMask.vol == 0) & (expMask.vol ~= 0));
  assert(sum(outlyingVoxels(:)) < 10, ['Activated voxels lie outside of ', ...
      baseBrainMaskPath, ...
      '.Please check that the 1-mm mask used in CBIG_AuthorTopic_PreprocessExpDataFromText.m ', ...
      'corresponds to the 2-mm mask in CBIG_AuthorTopic_ConvertBrainImagesToCVBData.m and CBIG_AuthorTopic_SetupParameters']);
  expMask.vol(brainMask.vol == 0) = 0;

  act = logical(act);
  save(cvbDataFile, 'corpus', 'act', 'taskByExp', 'uniqueTasks', 'uniqueTaskCount', '-v7.3');

  expMaskDir = fullfile(dataDir, 'mask');
  system(['mkdir -p ' expMaskDir]);
  expMaskFile = fullfile(expMaskDir, 'expMask.nii.gz');
  MRIwrite(expMask, expMaskFile);
