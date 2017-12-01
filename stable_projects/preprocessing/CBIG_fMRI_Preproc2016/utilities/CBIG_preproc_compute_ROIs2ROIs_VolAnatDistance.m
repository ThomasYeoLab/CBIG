function CBIG_preproc_compute_ROIs2ROIs_VolAnatDistance( cortical_ROIs_file, subcortical_ROIs_file, output_name )

% CBIG_preproc_compute_ROIs2ROIs_VolAnatDistance( cortical_ROIs_file, subcortical_ROIs_file, output_name )
% 
% Given the ROI files "cortical_ROIs_file" and "subcortical_ROIs_file",
% this function compute the ROIs to ROIs volumetric Euclidean distances in
% anatomical space and save out the distance and centroids matrices into
% "output_name".
% 
% IMPORTANT: this function cannot resolve the overlap between cortical ROIs
% and subcortical ROIs. Pleasse make sure your input ROIs are not
% overlapped.
% 
% Make sure your "cortical_ROIs_file" and "subcortical_ROIs_file" are in
% the same space and with the same vox2ras matrix. If both of them are
% passed in, the code will use the vox2ras matrix of
% "subcortical_ROIs_file".
% 
% Inputs:
%     - cortical_ROIs_file:
%       The full name of the cortical parcellation. For example, the
%       400-ROIs parcellation of Schaefer et al. 2018 in MNI 1mm space:
%       'xxxx/Schaefer2018_LocalGlobal/Parcellations/MNI/Schaefer2018_400Parcels_17Networks_order_FSLMNI152_1mm.nii.gz'.
%       If your ROIs do not include any cortical part, you can use 'NONE'.
% 
%     - subcortical_ROIs_file:
%       The full name of the subcortical parcellation. For example, the
%       aseg file in MNI 1mm space:
%       'xxxx/data/templates/volume/FSL_MNI152_FS4.5.0/mri/aparc+aseg_182x218x182.nii.gz'.
%       If your ROIs do not include any subcortical part, you can use
%       'NONE'.
% 
%     - output_name:
%       The full name of the output distance matrix. The output '.mat' file
%       is a structure called "distance" which contains two fields:
%       (1) "distance.distance" is an M x M distance matrix, where M is the
%       total number of ROIs (cortical + subcortical).
%       (2) "distance.centroids" is an M x 3 matrix contains the centroids
%       of all ROIs.
% 
% Written by Jingwei Li and CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md

%% Read in the ROIs
% We only care about these 19 subcortical ROIs (aseg labels)
% 8   Left-Cerebellum-Cortex
% 10  Left-Thalamus-Proper
% 11  Left-Caudate
% 12  Left-Putamen
% 13  Left-Pallidum
% 16  Brain-Stem
% 17  Left-Hippocampus
% 18  Left-Amygdala
% 26  Left-Accumbens-area
% 28  Left-VentralDC
% 47  Right-Cerebellum-Cortex
% 49  Right-Thalamus-Proper
% 50  Right-Caudate
% 51  Right-Putamen
% 52  Right-Pallidum
% 53  Right-Hippocampus
% 54  Right-Amygdala
% 58  Right-Accumbens-area 
% 60  Right-VentralDC
subcort_labels = [8 10 11 12 13 16 17 18 26 28 47 49 50 51 52 53 54 58 60];

if(~strcmp(cortical_ROIs_file, 'NONE'))
    cort_ROIs = MRIread(cortical_ROIs_file);
    cort_ROIs_vol = reshape(cort_ROIs.vol, numel(cort_ROIs.vol), 1);
    vox2ras = cort_ROIs.vox2ras;
    vol_size = size(cort_ROIs.vol);
end
if(~strcmp(subcortical_ROIs_file, 'NONE'))
    subcort_ROIs = MRIread(subcortical_ROIs_file);
    subcort_ROIs_vol = reshape(subcort_ROIs.vol, numel(subcort_ROIs.vol), 1);
    vox2ras = subcort_ROIs.vox2ras;
    vol_size = size(subcort_ROIs.vol);
    
    subcort_ind = [];
    subcort_ROIs_vol_masked = zeros(size(subcort_ROIs_vol));
    for i = 1:length(subcort_labels)
        tmp_ind = find(subcort_ROIs_vol == subcort_labels(i));
        subcort_ind = [subcort_ind; tmp_ind];
        subcort_ROIs_vol_masked(tmp_ind) = i;   % relabel the ROIs from 1 to 19
    end
    clear subcort_ROIs_vol
end
if(strcmp(cortical_ROIs_file, 'NONE') && strcmp(subcortical_ROIs_file, 'NONE'))
    error('Cortical ROIs and subcortical ROIs cannot be both ''NONE''.');
end
if(exist('cort_ROIs', 'var') && exist('subcort_ROIs', 'var'))
    vox2ras_diff = sum(abs(cort_ROIs.vox2ras(:) - subcort_ROIs.vox2ras(:)));
    if(vox2ras_diff ~= 0)
        error('Cortical ROIs vox2ras is not equal to subcortical ROIs vox2ras.');
    end
end


%% Combine 
if(~strcmp(cortical_ROIs_file, 'NONE') && ~strcmp(subcortical_ROIs_file, 'NONE'))
    subcort_ROIs_vol_masked(subcort_ROIs_vol_masked~=0) = max(cort_ROIs_vol) + subcort_ROIs_vol_masked(subcort_ROIs_vol_masked~=0);
    ROIs_vol = cort_ROIs_vol + subcort_ROIs_vol_masked;
    clear cort_ROIs_vol subcort_ROIs_vol_masked
elseif(strcmp(cortical_ROIs_file, 'NONE'))
    ROIs_vol = subcort_ROIs_vol_masked;
    clear subcort_ROIs_vol_masked
else
    ROIs_vol = cort_ROIs_vol;
    clear cort_ROIs_vol
end

%% Compute centroids
ROIs_unique = unique(ROIs_vol);
ROIs_unique = setdiff(ROIs_unique, 0);
centroids = zeros(length(ROIs_unique), 3);
for i = 1:length(ROIs_unique)
    parcel_ind = find(ROIs_vol==ROIs_unique(i));
    [parcel_x, parcel_y, parcel_z] = ind2sub(vol_size, parcel_ind);
    parcel_vox = [parcel_y'; parcel_x'; parcel_z'];
    parcel_ras = CBIG_ConvertVox2Ras(parcel_vox - 1, vox2ras);
    centroids(i, :) = mean(parcel_ras, 2)';
end


%% Compute distances between centroids
% distance_tmp = zeros(length(ROIs_unique) * (length(ROIs_unique)-1) / 2, 1);
distance_tmp = zeros(length(ROIs_unique));
count = 0;
for i = 1:length(ROIs_unique)-1
    for j = i+1:length(ROIs_unique)
        count = count + 1;
        distance_tmp(i, j) = sqrt(sum( (centroids(i, :) - centroids(j, :)).^2 ));
    end
end
distance_tmp = (distance_tmp + distance_tmp') ;

distance.centroids = centroids;
distance.distance = distance_tmp;


%% save
[output_dir, ~, ~] = fileparts(output_name);
if(~exist(output_dir, 'dir'))
    mkdir(output_dir);
end
save(output_name, 'distance');


end

