function [label_structName_avgWinProb_tem, label_structName_avgWinProb_sub, label_structName_avgWinProb_cor] = CBIG_assignFactorsToStructures(tem_sub_cor)

% [label_structName_avgWinProb_tem, label_structName_avgWinProb_sub, label_structName_avgWinProb_cor] = CBIG_assignFactorsToStructures(tem_sub_cor)
%% aparc+aseg.mgz of MNI152
%
% Written by CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md


labelMRIFile = '~/thomas/code/templates/volume/FSL_MNI152_FS4.5.0/mri/aparc+aseg.mgz';
labelMRI = MRIread(labelMRIFile);
labelVol = labelMRI.vol;

% Count how many voxels for each structure (label)
labels = unique(labelVol(:));

%% Find FreeSurfer structure names for these labels

fileID = fopen('/apps/arch/Linux_x86_64/freesurfer/5.3.0/FreeSurferColorLUT.txt');
tbl = textscan(fileID, '%d %s %d %d %d %d', 'CommentStyle', '#');
fclose(fileID);
labels_LUT = tbl{1};
names_LUT = tbl{2};

%% Probability atrophy map volume in MNI152

r = dir('../../lda/postprocess/188blAD/k3/r*');

% Load all probabilistic atrophy maps
probVols = cell(3, 1); % temporal, subcortical, cortical
for idx = 1:3
    whichFactor = tem_sub_cor(idx);
    probMRIFile = ['../../lda/postprocess/188blAD/k3/' r.name '/topic' num2str(whichFactor) '.nii.gz'];
    % Resample this atrophy map from 109x91x91 to 256x256x256
    system(sprintf(...
        'mri_vol2vol --regheader --mov %s --targ %s --o ./atrophyMaps/topic%d_resampled.nii.gz > ./atrophyMaps/topic%d.log',...
        probMRIFile, labelMRIFile, whichFactor, whichFactor));
    % Load resampled MRI
    probMRI = MRIread(['./atrophyMaps/topic' num2str(whichFactor) '_resampled.nii.gz']);
    probVols{idx} = probMRI.vol;
end

%% Assign a winning factor to each structure/label

noStructs = numel(labels);

label_structName_factor_avgWinProb = cell(noStructs, 4);

% For each structure
for idx = 1:noStructs
    label = labels(idx);
    % Label
    label_structName_factor_avgWinProb{idx, 1} = label;
    % Structure name
    label_structName_factor_avgWinProb{idx, 2} = names_LUT{labels_LUT==label};
    % Winning factor
    indVol = labelVol==label;
    noVox = sum(indVol(:));
    sumTemProb = sum(probVols{1}(indVol)); % temporal
    sumSubProb = sum(probVols{2}(indVol)); % subcortical
    sumCorProb = sum(probVols{3}(indVol)); % cortical
    maxSum = max([sumTemProb sumSubProb sumCorProb]);
    switch maxSum
        case sumTemProb
            label_structName_factor_avgWinProb{idx, 3} = 'temporal';
            label_structName_factor_avgWinProb{idx, 4} = sumTemProb/noVox;
        case sumSubProb
            label_structName_factor_avgWinProb{idx, 3} = 'subcortical';
            label_structName_factor_avgWinProb{idx, 4} = sumSubProb/noVox;
        otherwise
            label_structName_factor_avgWinProb{idx, 3} = 'cortical';
            label_structName_factor_avgWinProb{idx, 4} = sumCorProb/noVox;
    end
end

%% Keyword filtering

% Temporal
label_structName_avgWinProb_tem = label_structName_factor_avgWinProb(...
    strcmp(label_structName_factor_avgWinProb(:, 3), 'temporal'), [1 2 4]);
indKeep = ones(size(label_structName_avgWinProb_tem, 1), 1);
keywordList = {'unknown', 'Unknown', 'CSF', 'Lat-Vent', 'WM', 'choroid-plexus', 'CC_', 'corpuscallosum'};
for idx = 1:numel(keywordList) % for each keyword
    keyword = keywordList{idx};
    indKeep_thisKeyword = cellfun(@isempty, strfind(label_structName_avgWinProb_tem(:, 2), keyword));
    indKeep = indKeep&indKeep_thisKeyword; % as long as a cell contains one keyword, it's out
end
fprintf('\n************** Temporal **************\n\n');
fprintf('Removed by keyword filtering:\n');
disp(label_structName_avgWinProb_tem(~indKeep, :));
fprintf('Included:\n');
label_structName_avgWinProb_tem = label_structName_avgWinProb_tem(indKeep, :);
disp(label_structName_avgWinProb_tem);

% Subcortical
label_structName_avgWinProb_sub = label_structName_factor_avgWinProb(...
    strcmp(label_structName_factor_avgWinProb(:, 3), 'subcortical'), [1 2 4]);
indKeep = ones(size(label_structName_avgWinProb_sub, 1), 1);
keywordList = {'unknown', 'Unknown', 'white', 'White', 'vessel', 'Vessel', 'optic', 'Optic', 'CSF', 'Ventricle', 'WM', 'wm', 'CC_'};
for idx = 1:numel(keywordList) % for each keyword
    keyword = keywordList{idx};
    indKeep_thisKeyword = cellfun(@isempty, strfind(label_structName_avgWinProb_sub(:, 2), keyword));
    indKeep = indKeep&indKeep_thisKeyword; % as long as a cell contains one keyword, it's out
end
fprintf('\n************** Subcortical **************\n\n');
fprintf('Removed by keyword filtering:\n');
disp(label_structName_avgWinProb_sub(~indKeep, :));
fprintf('Included:\n');
label_structName_avgWinProb_sub = label_structName_avgWinProb_sub(indKeep, :);
disp(label_structName_avgWinProb_sub);

% Cortical
label_structName_avgWinProb_cor = label_structName_factor_avgWinProb(...
    strcmp(label_structName_factor_avgWinProb(:, 3), 'cortical'), [1 2 4]);
indKeep = ones(size(label_structName_avgWinProb_cor, 1), 1);
keywordList = {'unknown', 'Unknown', 'white', 'White'};
for idx = 1:numel(keywordList) % for each keyword
    keyword = keywordList{idx};
    indKeep_thisKeyword = cellfun(@isempty, strfind(label_structName_avgWinProb_cor(:, 2), keyword));
    indKeep = indKeep&indKeep_thisKeyword; % as long as a cell contains one keyword, it's out
end
fprintf('\n************** Cortical **************\n\n');
fprintf('Removed by keyword filtering:\n');
disp(label_structName_avgWinProb_cor(~indKeep, :));
fprintf('Included:\n');
label_structName_avgWinProb_cor = label_structName_avgWinProb_cor(indKeep, :);
disp(label_structName_avgWinProb_cor);
