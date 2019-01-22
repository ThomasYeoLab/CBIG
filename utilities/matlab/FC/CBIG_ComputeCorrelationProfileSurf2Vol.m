function CBIG_ComputeCorrelationProfileSurf2Vol(seed_mesh, output_file, threshold, mask_file, input_file_txt, seed_file_txt)

% CBIG_ComputeCorrelationProfileSurf2Vol(seed_mesh, output_file, threshold, mask_file, input_file_txt, seed_file_txt)
%
% seed_mesh:       resolution of seed regions, e.g. 'fsaverage3'
% threshold:       correlation matrix absolute threshold for binarization
% mask_file:       volumetric mask
% input_file_text: text file of input data (volume)
% seed_file_text:  text file of seed data (surface)
%
% Compute surface to volume correlation profiles.
%
% Written by CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md


% read in seed files.
fid = fopen(seed_file_txt, 'r');
i = 0;
while(1);
   tmp = fscanf(fid, '%s\n', 1);
   if(isempty(tmp))
       break
   else
       i = i + 1;
       seed_files{i} = tmp;
   end
end
fclose(fid);

% read in input files.
fid = fopen(input_file_txt, 'r');
i = 0;
while(1);
   tmp = fscanf(fid, '%s\n', 1);
   if(isempty(tmp))
       break
   else
       i = i + 1;
       input_files{i} = tmp;
   end
end
fclose(fid);

mask = MRIread(mask_file);
mask_index = find(mask.vol == 1);

lh_avg_mesh = CBIG_ReadNCAvgMesh('lh', seed_mesh, 'inflated', 'cortex');
rh_avg_mesh = CBIG_ReadNCAvgMesh('rh', seed_mesh, 'inflated', 'cortex');

for i = 1:length(input_files)
    input_file = input_files{i};
    input_series = MRIread(input_file);
    t_series = transpose(reshape(input_series.vol, size(input_series.vol, 1) * size(input_series.vol, 2) * size(input_series.vol, 3), size(input_series.vol, 4)));
    t_series = t_series(:, mask_index);
    
    seed_file = seed_files{i};
    hemi_index = strfind(seed_file, basename(seed_file));
    lh_seed_file = seed_file; lh_seed_file(hemi_index:hemi_index+1) = 'lh';
    lh_seed_series = MRIread(lh_seed_file);
    lh_seed_series = transpose(reshape(lh_seed_series.vol, size(lh_seed_series.vol, 1) * size(lh_seed_series.vol, 2) * size(lh_seed_series.vol, 3), size(lh_seed_series.vol, 4)));
    lh_seed_series = lh_seed_series(:, lh_avg_mesh.MARS_label == 2);
    
    rh_seed_file = seed_file; rh_seed_file(hemi_index:hemi_index+1) = 'rh';
    rh_seed_series = MRIread(rh_seed_file);
    rh_seed_series = transpose(reshape(rh_seed_series.vol, size(rh_seed_series.vol, 1) * size(rh_seed_series.vol, 2) * size(rh_seed_series.vol, 3), size(rh_seed_series.vol, 4)));
    rh_seed_series = rh_seed_series(:, rh_avg_mesh.MARS_label == 2);
    
    s_series = [lh_seed_series rh_seed_series];

    % normalize series (note that series are now of dimensions: T x N)
    s_series = s_series - repmat(mean(s_series, 1), size(s_series, 1), 1);
    s_series = s_series./repmat(sqrt(sum(s_series.^2, 1)), size(s_series, 1), 1);

    t_series = t_series - repmat(mean(t_series, 1), size(t_series, 1), 1);
    t_series = t_series./repmat(sqrt(sum(t_series.^2, 1)), size(t_series, 1), 1);

    corr_mat = zeros(size(s_series, 2), size(t_series, 2));
    for k = 1:size(s_series, 2)
        corr_mat(k, :) = sum(repmat(s_series(:, k), 1, size(t_series, 2)) .* t_series, 1);
    end
    
    if(i == 1)
        output = corr_mat;
    else
        output = output + corr_mat;
    end
end
output = output / length(input_files);
corr_mat = output;
disp(['isnan: ' num2str(sum(isnan(corr_mat(:)))) ' out of ' num2str(numel(corr_mat))]);
corr_mat(isnan(corr_mat)) = 0;

tmp = sort(corr_mat(:), 'descend');
t = tmp(round(numel(corr_mat) * str2num(threshold)));
disp(['threshold: ' num2str(t)]);
if(str2num(threshold) < 1)
  corr_mat(corr_mat <  t) = 0;
  corr_mat(corr_mat >= t) = 1;
end

output_series = input_series;
size_input = size(input_series.vol);
size_input(4) = size(corr_mat, 1);
output_series.vol = zeros(size_input);
frame = zeros(size_input(1:3));
for i = 1:size(corr_mat, 1)
    frame(mask_index) = corr_mat(i, :);
    output_series.vol(:, :, :, i) = frame;
end
MRIwrite(output_series, output_file); 

% exit
