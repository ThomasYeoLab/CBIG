function CBIG_AuthorTopicEM_GenerateDataFromText(text_data_path, ...
    data_dir, gibbs_data_file_name, em_data_file_name, brain_mask1mm, ...
    brain_mask2mm, binary_smooth_kernel_rad)
% CBIG_AuthorTopicEM_GenerateDataFromText(text_data_path, ...
%   data_dir, gibbs_data_file_name, em_data_file_name, brain_mask1mm, ...
%   brain_mask2mm, binary_smooth_kernel_rad)
%
% Generate input data for Gibbs sampler and Expectation-Maximization (EM)
% algorithm using coordinate-based meta-analytic data saved in a text file.
%
% Input:
%  - text_data_path: path to the text file containing the data of activation
%    foci in MNI152 space.
%  - data_dir: path to the output base directory containing the preprocessed data.
%      The following file and directories will be generated:
%      "<data_dir>/ExperimentsData.mat" contains the activation foci read
%      into a .mat file.
%      "<data_dir>/ActivationVolumes" contains the activated brain images.
%      "<data_dir>/BinarySmoothedVolumes" contains the activated brain
%      images smoothed with a binary smoothing kernel.
%      All brain images have 1mm resolution.
%  - gibbs_data_file_name: name of the mat file containing the input data for
%      Gibbs sampler.
%      Full path to this file would be "<data_dir>/<gibbs_data_file_name>"
%  - em_data_file_name: name of the mat fie containing the input data for the
%    EM algorithm.
%      Full path to this file would be "<data_dir>/<em_data_file_name>"
%  _ brain_mask1mm, brain_mask2mm: struct from MRIread containing the mask for
%    brain voxels of the MNI152 template in 1mm and 2mm resolution. brain_mask1mm
%    and brain_mask2mm need to correspond to each other.
%  - binary_smooth_kernel_rad: radius of the binary smoothing kernel in mm.
%      All voxels within binary_smooth_kernel_rad-mm radius are assigned the
%      value of 1, and 0 otherwise.
%      Default: 10
%
% Example:
%   text_data_path = fullfile(getenv('CBIG_CODE_DIR'), 'stable_projects', ...
%     'meta-analysis', 'Yeo2015_AuthorTopicEM', 'examples', ...
%     'SelfGeneratedThought_AllCoordinates.txt');
%   data_dir = '/Work/data';
%   brain_mask1mm = MRIread(fullfile(getenv('CBIG_CODE_DIR'), 'data', ...
%     'templates', 'volume','FSL_MNI152_FS4.5.0', 'mri', 'brain.mgz'));
%   brain_mask2mm = MRIread(fullfile(getenv('CBIG_CODE_DIR'), 'stable_projects', ...
%     'meta-analysis', 'Yeo2015_AuthorTopicEM', 'utilities', ...
%     'mask', 'MNI_mask_conformed.2mm.0.1.nii.gz'));
%   CBIG_AuthorTopicEM_GenerateGibbsEMDataFromText(text_data_path, data_dir, ...
%     "SelfGenearatedThought_GibbsData.mat", "SelfGenearatedThought_EMData.mat", ...
%     brain_mask, 10)
%   Process activation foci data saved at `text_data_path`. The preprocessed
%   files are saved under '/Work/data' and the input data file for the Gibbs sampler
%   is saved at '/Work/data/SelfGenearatedThought_GibbsData.mat', and EM
%   algorithm is saved at '/Work/data/SelfGenearatedThought_EMData.mat'.
%   The smoothed brain images are produced with a 10mm-radius binary
%   smoothing kernel.
%
% Written by B.T.Thomas Yeo and CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md

  if (nargin < 7)
    binary_smooth_kernel_rad = 10;
  end

  if (~exist(data_dir, 'dir'))
    mkdir(data_dir);
  end

  gibbs_data_file_path = fullfile(data_dir, gibbs_data_file_name);
  em_data_file_path = fullfile(data_dir, em_data_file_name);

  CBIG_AuthorTopicEM_PreprocessExpDataFromText(text_data_path, data_dir, binary_smooth_kernel_rad, brain_mask1mm);
  CBIG_AuthorTopicEM_ConvertBrainImages(data_dir, gibbs_data_file_path, em_data_file_path, brain_mask2mm);
