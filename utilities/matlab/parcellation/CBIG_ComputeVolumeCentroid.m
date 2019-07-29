function [center_ras, center_vox, label_list, name_list] = ...
    CBIG_ComputeVolumeCentroid(input_file, output_file, lut_file, min_label, max_label)

% [center_ras, center_vox, label_list, name_list] = ...
%       CBIG_ComputeVolumeCentroid(input_file, output_file, lut_file, min_label, max_label)
%
% This function computes ROI centroids in volume.
% 
% Please note that this function can be slow when an ROI includes too many 
% voxels (e.g. more than 10000).
%
% Input:
%      -input_file: 
%       Path of the volume data. It can be a NIFTI file (eg, .nii) or a MGH
%       file (eg, .mgz).
%
%      -output_file:
%       Path of the output text file. It will contain the RAS coordinates 
%       of each ROI. 
%       If you do not wish to create a text file, please pass [].
%       If the lut_file is used, then the format of the output file will
%       be: ROI Label,ROI Name,R,A,S. Otherwise it will be: ROI Label,R,A,S
% 
%      -lut_file:
%       Path of the lookup table file. This file shows the ROI name of each
%       ROI label . It should be a text file with freesurfer lookup table  
%       format. e.g. ROI_label_number ROI_label_name R G B A
%       If you do not need the name of the ROIs, please pass [].
%
%      -min_label, max_label:
%       First and last ROI label numbers. If min_label and max_label are
%       not given, then compute all the labels in the volume file.
%       If the input volume file is the output of mri_aparc2aseg, the 
%       cortex label will start from 1001 for the left hemisphere and from  
%       2001 for the right hemisphere. 
%
% Output:
%       -center_ras:
%        RAS coordinates of ROI centroids.
% 
%       -center_vox:
%        VOX coordinates of ROI centroids.
% 
%       -label_list:
%        ROI label list.
%
%       -name_list:
%        ROI name list. If lut_file is not used, this will be [].
%
% Example:
% center_ras = CBIG_ComputeVolumeCentroid(input_file, [], [])
%       
% Written by XUE Aihuiping and CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md

if nargin < 4
    min_label = [];
    max_label = [];
elseif nargin < 5
    error('ERROR: Please provide both min_label and max_label.');
end
if ischar(min_label)
    min_label = str2num(min_label);
end
if ischar(max_label)
    max_label = str2num(max_label);
end
if isempty(output_file)
    save_flag = 0;
else
    save_flag = 1;
end
if isempty(lut_file)
    name_flag = 0;
    name_list = [];
else
    name_flag = 1;
end

% Read input volume file
disp(['Reading input file: ' input_file]);
MRI_struct = MRIread(input_file);
vol = MRI_struct.vol;
vol = permute(vol,[2,1,3]);
sizevol = size(vol);
vox2ras = MRI_struct.vox2ras;

% Find ROI labels
label_list = unique(vol);
label_list = setdiff(label_list, 0);
N_ROI = length(label_list);
disp([num2str(N_ROI), ' ROIs found.']);

% Pick ROI labels
if isempty(min_label)
    disp(['Label range is empty. Compute all ' num2str(N_ROI) ' ROI(s).']);
else
    label_list = intersect(min_label:max_label, label_list);
    N_ROI = length(label_list);
    disp(['Label range is from ' num2str(min_label) ' to ' num2str(max_label) '.']);
    disp(['Compute ' num2str(N_ROI) ' ROI(s).']);
end

% Find voxels for each ROI
ROI_cell = cell(1, N_ROI);
center_vox = zeros(N_ROI, 3);
L = zeros(N_ROI,1);
N_com = cell(N_ROI,1);
for i = 1:N_ROI
    [temp(:,1), temp(:,2), temp(:,3)] = ind2sub(sizevol, find(label_list(i) == vol));
    ROI_cell{i} = temp;
    clear temp
    L(i) = size(ROI_cell{i},1);
end

% Compute the center of each ROI
for i = 1:N_ROI
    disp(['Computing center for ROI ' num2str(label_list(i)) ', ' num2str(L(i)) ' voxels included. ']);
    G = zeros(L(i));
    for j = 1:L(i)
        for k = (j+1):L(i)
            diff = abs(ROI_cell{i}(j,:)-ROI_cell{i}(k,:));
            if max(diff) <= 1
                if sum(diff) == 1
                    G(j,k) = 1;
                    G(k,j) = 1;
                elseif sum(diff) == 2
                    G(j,k) = sqrt(2);
                    G(k,j) = sqrt(2);   
                end
            end
        end
    end
    Gs = sparse(G);
    [Gs,N_com{i}] = MergeComponents(Gs,ROI_cell{i});

    d = zeros(L(i));
    for j = 1:L(i)
        d(:,j) = shortest_paths(Gs,j);
    end
    temp_sum=sum(d);
    [~, temp_min_i] = min(temp_sum);
    center_vox(i,:) = ROI_cell{i}(temp_min_i,:);
end

% Convert the VOX coordinates to RAS coordinates
center_vox = center_vox - 1;
center_ras = CBIG_ConvertVox2Ras(center_vox', vox2ras);
center_ras = center_ras';

% Read lookup table
if name_flag
    fid = fopen(lut_file,'r');
    LUT = textscan(fid,'%s%s%s%s%s%s','MultipleDelimsAsOne', 1);
    fclose(fid);
    lut_label = zeros(size(LUT{1}, 1), 1);
    for i = 1:size(LUT{1}, 1)
        temp = str2num(char(LUT{1}(i, :)));
        if isempty(temp)
            lut_label(i) = nan;
        else
            lut_label(i) = temp;
        end
    end
    name_list = cell(size(label_list));
    for i = 1:length(label_list)
        index = find(lut_label == label_list(i));
        if isempty(index)
            name_list{i} = 'Not_found_in_LUT';
        else
            name_list{i} = LUT{2}(index, :);
        end
    end
end

% Save RAS coordinates
if save_flag
    if exist(output_file, 'file')
        delete(output_file);
    end
    fid = fopen(output_file, 'w');
    if name_flag
        line = 'ROI Label,ROI Name,R,A,S\n';
    else
        line = 'ROI Label,R,A,S\n';
    end
    fprintf(fid, line);
    for i = 1:1:N_ROI
        line = [num2str(label_list(i))];
        if name_flag
            line = [line ',' char(name_list{i, :})];
        end
        for j = 1:3
            line = [line ',' num2str(round(center_ras(i, j) * 100) / 100)];
        end
        line = [line '\n'];
        fprintf(fid, line);
    end
    fclose(fid);
    disp(['Results saved to ' output_file]);
end
end

function [Gout, con] = MergeComponents(Gin, cor)

if (~issparse(Gin))
    Gin = sparse(Gin);
end
L = size(Gin,1);
[ci, csize] = components(Gin);
con = csize;
while size(csize,1) > 1
    [csize, sorti] = sort(csize);
    d = zeros(csize(1), L-csize(1));
    c1 = find(ci == sorti(1));
    c2 = find(ci ~= sorti(1));
    for j = 1:csize(1)
        for k = 1:(L-csize(1))
            d(j,k) = pdist([cor(c1(j),:); cor(c2(k), :)], 'euclidean');
        end
    end
    temp_min = min(min(d));
    [p1, p2] = find(d == temp_min);
    for j = 1:size(p1, 1)
        Gin(c1(p1(j)), c2(p2(j))) = temp_min;
        Gin(c2(p2(j)), c1(p1(j))) = temp_min;
    end
    [ci, csize] = components(Gin);
end
Gout = Gin;

end
