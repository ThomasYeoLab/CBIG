
function CBIG_MM_create_FC419_MNI2mm(outputDir, fc400_nii)

% CBIG_MM_create_FC419_MNI2mm(outputDir, fc400_nii)
% 
% This function create functional connectivity nii for MNI2mm space
%
% Inputs:
%   - outputDir
%     Path of the your output data. You can also change this
%     to any place you want.
%
%   - fc400_nii
%     Path of the Schaefer2018 400 Parcels mask for MNI 152 2mm.
%
% Written by Tong He and CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md

    if ~exist(outputDir, 'dir')
        mkdir(outputDir);
    end

    label = [26 58 18 54 16 11 50 8 47 28 60 17 53 13 52 12 51 10 49];
    create_subcorMask(outputDir, label');

    MNI_2mm=fullfile(getenv('FSL_DIR'), 'data', 'standard', 'MNI152_T1_2mm_brain.nii.gz');

    command = strcat("mri_vol2vol --mov ", outputDir, "/subcor_mask.nii.gz --targ ",...
        MNI_2mm, " --o ", outputDir,...
        "/subcor_mask_MNI2mm.nii.gz --regheader --no-save-reg --nearest");
    system(command);

    x = MRIread(fc400_nii);
    y = MRIread(fullfile(outputDir, 'subcor_mask_MNI2mm.nii.gz'));

    temp = y.vol;
    labels = sort(label);
    for i = 1:length(labels)
        temp(temp == labels(i)) = 400 + i;
    end
    x.vol = x.vol + temp;
    MRIwrite(x,fullfile(outputDir, 'FC419_MNI2mm.nii.gz'));
end

function create_subcorMask(outputDir, labels)

% create_subcorMask(outputDir, labels)
% 
% This function create subcortical mask for 419 FC
%
% Inputs:
%   - outputDir
%     Path of the your output data. You can also change this
%     to any place you want.
%
%   - labels
%     Array of subcortical label used for 419 FC

    CBIG_CODE_DIR = getenv('CBIG_CODE_DIR');
    inputDir = [CBIG_CODE_DIR '/data/templates/volume'];

    if nargin < 2 || (isempty(labels))
        labels = [26 58 18 54 16 11 50 8 47 28 60 17 53 13 52 12 51 10 49]';
    end

    %% read aseg file in the subject's own space
    aparc_aseg = MRIread([inputDir '/FSL_MNI152_FS4.5.0/mri/aparc+aseg.mgz']);
    aparc_aseg_vol = aparc_aseg.vol;
    aparc_aseg_1d = reshape(aparc_aseg_vol, [size(aparc_aseg_vol,1)*size(aparc_aseg_vol,2)*size(aparc_aseg_vol,3) 1]);
    mask = ismember(aparc_aseg_1d,labels); % only subcortical regions

    %% write mask into .nii.gz file
    aparc_aseg_1d(~mask,:) = 0; % set other regions to 0

    dummy = aparc_aseg;
    dummy.vol = reshape(aparc_aseg_1d, size(aparc_aseg_vol));
    file_name = [outputDir '/subcor_mask.nii.gz'];
    MRIwrite(dummy,file_name);
end