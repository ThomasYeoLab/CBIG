function CBIG_AuthorTopic_PreprocessExpDataFromText(textDataPath, ...
    dataDir, binarySmoothKernelRad)
% CBIG_AuthorTopic_PreprocessExpDataFromText(textDataPath, ...
%   dataDir, binarySmoothKernelRad)
%
% Preprocess data of activation foci's locations in text form: write the
% activation foci into 3D brain images and smooth the brain images with
% a binary kernel
%
% Input:
%  - textDataPath: path to the text file containing the data of activation
%                  foci in MNI152 space.
%  - dataDir     : path to the base directory containing the preprocessed
%                  data.
%      "<dataDir>/ActivationVolumes" contains the activated brain images.
%      "<dataDir>/BinarySmoothedVolumes" contains the activated brain
%      images smoothed with a binary smoothing kernel.
%      All brain images have 1mm resolution.
%  - binarySmoothKernelRad: radius of the binary smoothing kernel in mm.
%      All voxels within binarySmoothKernelRad-mm radius are assigned the
%      value of 1, and 0 otherwise.
%
% Example:
%   textDataPath = fullfile(getenv('CBIG_CODE_DIR'), 'stable_projects', ...
%     'meta-analysis', 'Ngo2019_AuthorTopic', 'selfGenerated_thought', ...
%     'MNI152_activation_coordinates', ...
%     'SelfGeneratedThought_AllCoordinates.txt');
%   dataDir = 'Work/data';
%   CBIG_AuthorTopic_PreprocessExpDataFromText(textDataPath, dataDir, 10)
%   Process activation foci data saved at `textDataPath` and save the
%   activated brain images under '/Work/data/ActivationVolumes' and the
%   smoothed brain images under '/Work/data/BinarySmoothedVolumes'.
%   The smoothed brain images are produced with a 10mm-radius binary
%   smoothing kernel.
%
% Written by Gia H. Ngo and CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md

  % load brain mask
  brain = MRIread(fullfile(getenv('CBIG_CODE_DIR'), 'data', 'templates', 'volume','FSL_MNI152_FS4.5.0', 'mri', 'brain.mgz'));
  brain.vol(brain.vol ~= 0) = 1;

  % Read into matrix format: 4 x N where N is the number of coordinates
  % First row corresponds to experiment number
  fid = fopen(textDataPath, 'r');
  disp('Reading MNI coordinates from text file');
  expCount = 0;
  coordCount = 0;
  coordBegin = 0;
  expTasks = cell(1, 1);
  while(1)
    line = fgetl(fid);
    if(line == -1)
      break;
    end

    if(strfind(line, '//Task: ')) % contain task information
      tasks = strtrim(strsplit(line(9:end), ','));

      expCount = expCount + 1;
      expTasks{expCount} = tasks;

      if(mod(expCount, 10) == 1)
        disp(['  - Finished reading experiment ' num2str(expCount)]);
      end
      coordBegin = 1;
      continue
    elseif (isempty(line)) % end of experiment
      coordBegin = 0;
      continue
      elseif (coordBegin == 1)
        coordCount = coordCount + 1;
        coordinates = str2num(line)';
        if (isempty(coordinates))
          error('Extra spaces in text file');
        end

        MNI_coordinates(:, coordCount) = [expCount; coordinates];
      else
        % ignore line
      end
  end
  disp(' ');

  % Convert to voxel and matrix coordinates
  disp('Transforming MNI coordinates to Matrix Coordinates');
  vox_coord = MNI_coordinates;
  vox_coord(2:end, :) = CBIG_ConvertRas2Vox(MNI_coordinates(2:end, :), brain.vox2ras);
  matCoord(1, :) = vox_coord(1, :);
  matCoord(2, :) = round(vox_coord(3, :) + 1);
  matCoord(3, :) = round(vox_coord(2, :) + 1);
  matCoord(4, :) = round(vox_coord(4, :) + 1);
  disp(' ');

  Experiments = {};
  for e = 1:expCount
    expIndices = matCoord(1, :) == e;

    expCoords = matCoord(2:4, expIndices);

    Experiments{e}.coords = expCoords;
    Experiments{e}.task = expTasks{e};
  end

  experimentsFile = fullfile(dataDir, 'ExperimentsData.mat');
  save(experimentsFile, 'Experiments', 'expTasks');

  % Write activation volumes
  disp('Writing Activation Volumes');
  activationVolumesDir = fullfile(dataDir, 'ActivationVolumes');
  system(['mkdir -p ', activationVolumesDir]);
  for i = 1:expCount
    if(mod(i, 10) == 1)
      disp(['  - Finished processing experiment #' num2str(i)]);
    end

    outputFile = fullfile(activationVolumesDir, ['ActVolume' num2str(i, '%06d') '.nii.gz']);

    if(~exist(outputFile, 'file'))
      actVol = brain;
      actVol.vol(:) = 0;

      expCoords = matCoord(2:4, matCoord(1, :) == i);
      ind = sub2ind(size(actVol.vol), expCoords(1, :), expCoords(2, :), expCoords(3, :));
      actVol.vol(ind) = 1;

      MRIwrite(actVol, outputFile);
    end
  end
  disp(' ');

  disp('Binary Smoothing the Activation Volumes');
  binVolumesDir = fullfile(dataDir, 'BinarySmoothedVolumes');
  system(['mkdir -p ', binVolumesDir]);
  for i = 1:expCount
    if(mod(i, 10) == 1)
      disp(['  - Finished processing experiment #' num2str(i)]);
    end

    inputFile = fullfile(activationVolumesDir, ['ActVolume' num2str(i, '%06d') '.nii.gz']);
    outputFile = fullfile(binVolumesDir, ['BinVolume' num2str(i, '%06d') '.nii.gz']);

    if(~exist(outputFile, 'file'))
      binVol = MRIread(inputFile);

      binVol.vol = bwdist(binVol.vol);
      binVol.vol(binVol.vol <= binarySmoothKernelRad) = 1;
      binVol.vol(binVol.vol > binarySmoothKernelRad) = 0;

      binVol.vol(brain.vol == 0) = 0; %mask it

      MRIwrite(binVol, outputFile);
    end
  end
  disp(' ');
