function CBIG_ComputeCorrelationProfileSurf2Vol_noregress(seed_mesh, target, output_file, threshold, mask_file, input_file_txt, seed_file_txt)

% CBIG_ComputeCorrelationProfileSurf2Vol_noregress(seed_mesh, target, output_file, threshold, mask_file, input_file_txt, seed_file_txt)
%
% seed_mesh:      resolution of seed regions
% target:         resolution of target regions.
% output_file:    output file name
% threshold:      relative threshold of correlation matrix for binarization
% mask_file:      volumetric mask
% input_file_txt: a text file containing input volume file list
% seed_file_txt:  a text file containing seed surface file list
%
% Compute surface to volume correlation profiles without regression.
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

lh_seed_avg_mesh = CBIG_ReadNCAvgMesh('lh', seed_mesh, 'inflated', 'cortex');
rh_seed_avg_mesh = CBIG_ReadNCAvgMesh('rh', seed_mesh, 'inflated', 'cortex');

lh_avg_mesh = CBIG_ReadNCAvgMesh('lh', target, 'inflated', 'cortex');
rh_avg_mesh = CBIG_ReadNCAvgMesh('rh', target, 'inflated', 'cortex');

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
    lh_seed_series = lh_seed_series(:, lh_avg_mesh.MARS_label(1:length(lh_seed_avg_mesh.MARS_label)) == 2);
    
    rh_seed_file = seed_file; rh_seed_file(hemi_index:hemi_index+1) = 'rh';
    rh_seed_series = MRIread(rh_seed_file);
    rh_seed_series = transpose(reshape(rh_seed_series.vol, size(rh_seed_series.vol, 1) * size(rh_seed_series.vol, 2) * size(rh_seed_series.vol, 3), size(rh_seed_series.vol, 4)));
    rh_seed_series = rh_seed_series(:, rh_avg_mesh.MARS_label(1:length(rh_seed_avg_mesh.MARS_label)) == 2);
    
    s_series = [lh_seed_series rh_seed_series];

    if(nargin > 8)
        s_series = s_series(:, index);
    end
    
    % normalize series (note that series are now of dimensions: T x N)
    disp('Computing Surf2Vol Profiles');
    s_series = bsxfun(@minus, s_series, mean(s_series, 1));
    s_series = bsxfun(@times, s_series, 1./sqrt(sum(s_series.^2, 1)));

    t_series = bsxfun(@minus, t_series, mean(t_series, 1));
    t_series = bsxfun(@times, t_series, 1./sqrt(sum(t_series.^2, 1)));
    
    corr_mat = s_series' * t_series;
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

surf2vol_correlation_profile = transpose(corr_mat); % result is N voxels x # ROI 
save(output_file, 'surf2vol_correlation_profile', '-v7.3');

