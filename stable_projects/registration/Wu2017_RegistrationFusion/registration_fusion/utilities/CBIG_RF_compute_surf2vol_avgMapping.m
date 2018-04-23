function CBIG_RF_compute_surf2vol_avgMapping(sub_list, input_dir, input_prefix, output_dir, output_prefix)
% CBIG_RF_compute_surf2vol_avgMapping(sub_list, input_dir, input_prefix, output_dir, output_prefix)
%
% This function computes surface-to-volume avgerage mapping of x/y/z index 
% across a given list of subjects.
%
% Input:
%     - sub_list    :
%                      absolute/relative path to subject list file,
%                        which contains a subject name each line
%     - input_dir   :
%                      absolute/relative path to directory where input can be found
%                      inputs are index files projected to a volumetric atlas space, 
%                        which can be read by MRIread()
%     - input_prefix:
%                      a fixed prefix for each subject's mapping,
%                        which may specify the approach and template used
%                      For example, an index file projected from left hemisphere of fsaverage 
%                        through Subject 0001 to the volumetric atlas space is named: 
%                        lh.xIndex_fsaverage_to_Sub0001_Ses1_FS_to_input_prefix.nii.gz
%     - output_dir:
%                      absolute/relative path to directory where output should be stored
%     - output_prefix:
%                      desired prefix for the outputs
%
% Output:
%     - There is no function output.
%     - a mapping file is created in output_dir:
%           [output_prefix]_avgMapping.prop.mat
%     - a count map file is created in output_dir:
%           [output_prefix]_count.mat
%
% Example:
% CBIG_RF_compute_vol2surf_avgMapping('~/data/GSP_subjectnames.csv', '../results/index_MNI152/', 
%               'RF_M3Z_MNI', '../results/mappings/', 'RF_M3Z_MNI2fsaverage_avg1490Sub')
% This command reads in subject mappings from  ../results/index_fsaverage directory with 
% prefix of 'RF_M3Z_MNI'. An average mapping is generated across the subjects provided in 
% GSP_subjectnames.csv and put into ../results/mappings directory.
%
% Written by Wu Jianxiao and CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md

if nargin < 5
    disp('usage: CBIG_RF_compute_surf2vol_avgMapping(sub_list, input_dir, input_prefix, output_dir, output_prefix)');
    return
end

%Read subject names
fid = fopen(sub_list, 'r');
i = 1;
while ~feof(fid)
    sub_names{i} = fgetl(fid);
    i = i + 1;
end
fclose(fid);

%Loop through each hemisphere
for hemis = {'lh', 'rh'}
    hemi = hemis{1};
    
    %Loop through each subject
    for i = 1:length(sub_names)
        disp([num2str(i) ':' sub_names{i} ' in ' hemi]);
        
        %Read x, y, and z index mapping results
        x = MRIread([input_dir '/' hemi '.xIndex_fsaverage_to_' sub_names{i} '_to_' input_prefix '.nii.gz']);
        y = MRIread([input_dir '/' hemi '.yIndex_fsaverage_to_' sub_names{i} '_to_' input_prefix '.nii.gz']);
        z = MRIread([input_dir '/' hemi '.zIndex_fsaverage_to_' sub_names{i} '_to_' input_prefix '.nii.gz']);
        
        %Track sum of xyz-coordinates at each voxel
        %For each voxel, number of subjects that has a vertex mapped to it is also tracked.
        %A subject is considered to have a vertex mapped to a voxel, if the xyz-coordiante
        % at that voxel is not (0, 0, 0)
        if i == 1
            coord = [x.vol(:)'; y.vol(:)'; z.vol(:)'];
            count = double(sum(abs(coord))~=0);
        else
            coord_sub = [x.vol(:)'; y.vol(:)'; z.vol(:)'];
            coord = coord + coord_sub;
            count = count + double(sum(abs(coord_sub))~=0);
        end
    end
    
    %Compute average mapping and save results
    coord = bsxfun(@rdivide, coord, count+eps);
    save([output_dir '/' hemi '.' output_prefix '_temp.mat'], 'coord', 'count');
end

load([output_dir '/lh.' output_prefix '_temp.mat']);
lh_count = count;
lh_coord = coord;
load([output_dir '/rh.' output_prefix '_temp.mat']);
rh_count = count;
rh_coord = coord;

%Resolve difference near midline
lh_count(lh_count~=0 & rh_count~=0 & rh_count>=lh_count) = 0;
rh_count(lh_count~=0 & rh_count~=0 & lh_count>rh_count) = 0;
lh_coord(:, lh_count==0) = 0;
rh_coord(:, rh_count==0) = 0;

%Normalise the coordinates to fsaverage surface (radius = 100)
lh_norm = sqrt(sum(lh_coord.^2, 1));
lh_coord = lh_coord ./ (repmat(lh_norm, 3, 1) + eps) * 100;
rh_norm = sqrt(sum(rh_coord.^2, 1));
rh_coord = rh_coord ./ (repmat(rh_norm, 3, 1) + eps) * 100;

%Propagate to whole brain at threshold of 50% of total number of subjects
thresh = round(0.5 * length(sub_names));
count_map = double((lh_count + rh_count) >= thresh);
[~, id] = bwdist(reshape(count_map, size(x.vol)));
lh_coord = lh_coord(:, id);
rh_coord = rh_coord(:, id);
save([output_dir '/' output_prefix '_avgMapping.prop.mat'], 'lh_coord', 'rh_coord');
save([output_dir '/' output_prefix '_count.mat'], 'lh_count', 'rh_count');

%Delete temporary files
delete([output_dir '/lh.' output_prefix '_temp.mat'], [output_dir '/rh.' output_prefix '_temp.mat']);
end

