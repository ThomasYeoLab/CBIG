function CBIG_ComputeCorrelationProfileVol2Surf(roi_file, output_file, threshold, surf_input_files_txt, vol_input_files_txt, regress_bool)

% CBIG_ComputeCorrelationProfileVol2Surf(roi_file, output_file, threshold, surf_input_files_txt, vol_input_files_txt, regress_bool) 
%
% roi_file:            seeds in volume
% output_file:         output file name
% threshold:           relative threshold of correlation matrix for binarization
% surf_input_file_txt: a text file containing input surface file list
% vol_input_file_txt:  a text file containing input volume file list
% regress_bool (not used)
%
% Compute volume to surface correlation profiles.
%
% Written by CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md


% read in roi mask.
roi_mask = MRIread(roi_file);
roi_index = find(roi_mask.vol == 1);

% read in surf input files.
fid = fopen(surf_input_files_txt, 'r');
i = 0;
while(1);
   tmp = fscanf(fid, '%s\n', 1);
   if(isempty(tmp))
       break
   else
       i = i + 1;
       surf_input_files{i} = tmp;
   end
end
fclose(fid);

% read in surf input files.
fid = fopen(vol_input_files_txt, 'r');
i = 0;
while(1);
   tmp = fscanf(fid, '%s\n', 1);
   if(isempty(tmp))
       break
   else
       i = i + 1;
       vol_input_files{i} = tmp;
   end
end
fclose(fid);

for i = 1:length(surf_input_files)
    surf_input_file = surf_input_files{i};
    input_series = MRIread(surf_input_file);
    t_series = transpose(reshape(input_series.vol, size(input_series.vol, 1) * size(input_series.vol, 2) * size(input_series.vol, 3), size(input_series.vol, 4)));
    
    vol_input_file = vol_input_files{i};
    roi_series = MRIread(vol_input_file);
    s_series = transpose(reshape(roi_series.vol, size(roi_series.vol, 1) * size(roi_series.vol, 2) * size(roi_series.vol, 3), size(roi_series.vol, 4)));
    s_series = s_series(:, roi_index);
    
    % normalize series (note that series are now of dimensions: T x N)
    disp('Computing Vol2Surf Profiles');
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
output = output / length(surf_input_files);
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
output_series.vol = reshape(transpose(corr_mat), size_input);
MRIwrite(output_series, output_file); 

% exit


