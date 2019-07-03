function CBIG_AuthorTopicEM_PreprocessExpDataFromText(text_data_path, ...
    data_dir, binary_smooth_kernel_rad, brain_mask)
% CBIG_AuthorTopicEM_PreprocessExpDataFromText(text_data_path, ...
%   data_dir, binary_smooth_kernel_rad, brain_mask)
%
% Preprocess data of activation foci's locations in text form: write the
% activation foci into 3D brain images and smooth the brain images with
% a binary kernel
%
% Input:
%  - text_data_path = path to the text file containing the data of activation
%                  foci in MNI152 space.
%  - data_dir       = path to the base directory containing the preprocessed
%                     data.
%      "<data_dir>/ActivationVolumes" contains the activated brain images.
%      "<data_dir>/BinarySmoothedVolumes" contains the activated brain
%      images smoothed with a binary smoothing kernel.
%      All brain images have 1mm resolution.
%  - binary_smooth_kernel_rad = radius of the binary smoothing kernel in mm.
%      All voxels within binary_smooth_kernel_rad-mm radius are assigned the
%      value of 1, and 0 otherwise.
%  - brain_mask1mm = a brain mask at 1mm resolution read by MRIread.
%
% Example:
%   text_data_path = fullfile(getenv('CBIG_CODE_DIR'), 'stable_projects', ...
%     'meta-analysis', 'Ngo2019_AuthorTopic', 'selfGenerated_thought', ...
%     'MNI152_activation_coordinates', ...
%     'SelfGeneratedThought_AllCoordinates.txt');
%   brain_mask1mm = MRIread(fullfile(getenv('CBIG_CODE_DIR'), 'data', ...
%     'templates', 'volume','FSL_MNI152_FS4.5.0', 'mri', 'brain.mgz'));
%   data_dir = 'Work/data';
%   CBIG_AuthorTopicEM_PreprocessExpDataFromText(text_data_path, data_dir, ...
%     brain_mask1mm, 10)
%   Process activation foci data saved at `text_data_path` and save the
%   activated brain images under '/Work/data/ActivationVolumes' and the
%   smoothed brain images under '/Work/data/BinarySmoothedVolumes'.
%   The smoothed brain images are produced with a 10mm-radius binary
%   smoothing kernel.
%
% Written by Gia H. Ngo and CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md

  % load brain mask
  brain_mask.vol(brain_mask.vol ~= 0) = 1;

  % Read into matrix format: 4 x N where N is the number of coordinates
  % First row corresponds to experiment number
  fid = fopen(text_data_path, 'r');
  disp('Reading MNI coordinates from text file');
  exp_count = 0;
  coord_count = 0;
  coord_begin = 0;
  exp_tasks = cell(1, 1);
  while(1)
    line = fgetl(fid);
    if(line == -1)
      break;
    end

    if(strfind(line, '//Task: ')) % contain task information
      tasks = strtrim(strsplit(line(9:end), ','));

      exp_count = exp_count + 1;
      exp_tasks{exp_count} = tasks;

      if(mod(exp_count, 10) == 1)
        disp(['  - Finished reading experiment ' num2str(exp_count)]);
      end
      coord_begin = 1;
      continue
    elseif (isempty(line)) % end of experiment
      coord_begin = 0;
      continue
      elseif (coord_begin == 1)
        coord_count = coord_count + 1;
        coordinates = str2num(line)';
        if (isempty(coordinates))
          error('Extra spaces in text file');
        end

        MNI_coordinates(:, coord_count) = [exp_count; coordinates];
      else
        % ignore line
      end
  end
  disp(' ');

  % Convert to voxel and matrix coordinates
  disp('Transforming MNI coordinates to Matrix Coordinates');
  vox_coord = MNI_coordinates;
  vox_coord(2:end, :) = CBIG_ConvertRas2Vox(MNI_coordinates(2:end, :), brain_mask.vox2ras);
  mat_coord(1, :) = vox_coord(1, :);
  mat_coord(2, :) = round(vox_coord(3, :) + 1);
  mat_coord(3, :) = round(vox_coord(2, :) + 1);
  mat_coord(4, :) = round(vox_coord(4, :) + 1);
  disp(' ');

  experiments = {};
  for e = 1:exp_count
    exp_indices = mat_coord(1, :) == e;

    exp_coords = mat_coord(2:4, exp_indices);

    experiments{e}.coords = exp_coords;
    experiments{e}.task = exp_tasks{e};
  end

  experiments_file = fullfile(data_dir, 'ExperimentsData.mat');
  save(experiments_file, 'experiments', 'exp_tasks');

  % Write activation volumes
  disp('Writing Activation Volumes');
  activation_volumes_dir = fullfile(data_dir, 'ActivationVolumes');
  system(['mkdir -p ', activation_volumes_dir]);
  for i = 1:exp_count
    if(mod(i, 10) == 1)
      disp(['  - Finished processing experiment ' num2str(i)]);
    end

    output_file = fullfile(activation_volumes_dir, ['ActVolume' num2str(i, '%06d') '.nii.gz']);

    if(~exist(output_file, 'file'))
      act_vol = brain_mask;
      act_vol.vol(:) = 0;

      exp_coords = mat_coord(2:4, mat_coord(1, :) == i);
      ind = sub2ind(size(act_vol.vol), exp_coords(1, :), exp_coords(2, :), exp_coords(3, :));
      act_vol.vol(ind) = 1;

      MRIwrite(act_vol, output_file);
    end
  end
  disp(' ');

  disp('Smoothed by binary smoothing kernel');
  bin_volumes_dir = fullfile(data_dir, 'BinarySmoothedVolumes');
  system(['mkdir -p ', bin_volumes_dir]);
  for i = 1:exp_count
    if(mod(i, 10) == 1)
      disp(['  - Finished processing experiment ' num2str(i)]);
    end

    inputFile = fullfile(activation_volumes_dir, ['ActVolume' num2str(i, '%06d') '.nii.gz']);
    output_file = fullfile(bin_volumes_dir, ['BinVolume' num2str(i, '%06d') '.nii.gz']);

    if(~exist(output_file, 'file'))
      bin_vol = MRIread(inputFile);

      bin_vol.vol = bwdist(bin_vol.vol);
      bin_vol.vol(bin_vol.vol <= binary_smooth_kernel_rad) = 1;
      bin_vol.vol(bin_vol.vol > binary_smooth_kernel_rad) = 0;

      bin_vol.vol(brain_mask.vol == 0) = 0; %mask it

      MRIwrite(bin_vol, output_file);
    end
  end
  disp(' ');
