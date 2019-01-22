function CBIG_ComputeCorrelationProfileVol2Vol2mm(output_file, threshold, mask_file, input_file_txt, roi_file, regress_bool)

% CBIG_ComputeCorrelationProfileVol2Vol2mm(output_file, threshold, mask_file, input_file_txt, roi_file, regress_bool)
%
% output_file:    output file name
% threshold:      relative threshold of correlation matrix for binarization
% mask_file:      volumetric mask
% input_file_txt: a text file containing input file list
% roi_file:       seed regions in volume
% regress_bool (not used)
% 
% Compute volume to volume correlation profiles.
%
% Written by CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md


% read in roi mask.
roi_mask = MRIread(roi_file);
roi_index = find(roi_mask.vol == 1);
mask = MRIread(mask_file);
mask_index = find(mask.vol == 1);

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

for i = 1:length(input_files)
    input_file = input_files{i};
    input_series = MRIread(input_file);
    series = transpose(reshape(input_series.vol, size(input_series.vol, 1) * size(input_series.vol, 2) * size(input_series.vol, 3), size(input_series.vol, 4)));
    t_series = series(:, mask_index);
    s_series = series(:, roi_index);
    
    % normalize series (note that series are now of dimensions: T x N)
    disp('Computing Vol2Vol Profiles');
    s_series = s_series - repmat(mean(s_series, 1), size(s_series, 1), 1);
    s_series = s_series./repmat(sqrt(sum(s_series.^2, 1)), size(s_series, 1), 1);

    t_series = t_series - repmat(mean(t_series, 1), size(t_series, 1), 1);
    t_series = t_series./repmat(sqrt(sum(t_series.^2, 1)), size(t_series, 1), 1);

    corr_mat = zeros(size(s_series, 2), size(t_series, 2)); % ROI x N voxels
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

vol2vol_correlation_profile = transpose(corr_mat); % result is N voxels x # ROI 
save(output_file, 'vol2vol_correlation_profile','-v7.3');

% exit
