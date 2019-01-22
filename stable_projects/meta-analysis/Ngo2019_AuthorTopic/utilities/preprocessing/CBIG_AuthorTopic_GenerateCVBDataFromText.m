function CBIG_AuthorTopic_GenerateCVBDataFromText(textDataPath, ...
    dataDir, cvbDataFileName, binarySmoothKernelDim)
% CBIG_AuthorTopic_GenerateCVBDataFromText(textDataPath, ...
%   dataDir, cvbDataFileName, binarySmoothKernelDim)
%
% Generate input data for Collapsed Variational Bayes (CVB) algorithm
% using coordinate-based meta-analytic data saved in a text file.
%
% Input:
%  - textDataPath: path to the text file containing the data of activation
%    foci in MNI152 space.
%  - dataDir: path to the base directory containing the preprocessed data.
%      "<dataDir>/ExperimentsData.mat" contains the activation foci read
%      into a .mat file.
%      "<dataDir>/ActivationVolumes" contains the activated brain images.
%      "<dataDir>/BinarySmoothedVolumes" contains the activated brain
%      images smoothed with a binary smoothing kernel.
%      All brain images have 1mm resolution.
%  - cvbDataFileName: name of the mat file containing the input data for
%      CVB algorithm.
%      Full path to this file would be "<dataDir>/<cvbDataFileName>"
%  - binarySmoothKernelRad: radius of the binary smoothing kernel in mm.
%      All voxels within binarySmoothKernelRad-mm radius are assigned the
%      value of 1, and 0 otherwise.
%      Default: 10
%
% Example:
%   textDataPath = fullfile(getenv('CBIG_CODE_DIR'), 'stable_projects', ...
%     'meta-analysis', 'Ngo2019_AuthorTopic', 'SelfGeneratedThought', ...
%     'MNI152_activation_coordinates', ...
%     'SelfGeneratedThought_AllCoordinates.txt');
%   dataDir = '/Work/data';
%   CBIG_AuthorTopic_GenerateCVBDataFromText(textDataPath, dataDir, ...
%     "SelfGenearatedThought_CVBData.mat", 10)
%   Process activation foci data saved at `textDataPath`. The preprocessed
%   files are saved under '/Work/data' and the input data file for the CVB
%   algorithm is saved at '/Work/data/SelfGenearatedThought_CVBData.mat'.
%   The smoothed brain images are produced with a 10mm-radius binary
%   smoothing kernel.
%
% Written by Gia H. Ngo and CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md

  if (nargin < 4)
    binarySmoothKernelDim = 10;
  end

  if (~exist(dataDir, 'dir'))
    mkdir(dataDir);
  end

  cvbDataFilePath = fullfile(dataDir, cvbDataFileName);

  CBIG_AuthorTopic_PreprocessExpDataFromText(textDataPath, dataDir, binarySmoothKernelDim);
  CBIG_AuthorTopic_ConvertBrainImagesToCVBData(dataDir, cvbDataFilePath);
