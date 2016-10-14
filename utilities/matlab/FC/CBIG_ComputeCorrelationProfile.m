function CBIG_ComputeCorrelationProfile(seed_mesh, target, output_file1, output_file2, threshold, varargin_text1, varargin_text2, outlier_text)

% CBIG_ComputeCorrelationProfile(seed_mesh, target, output_file1, output_file2, threshold, varargin_text1, varargin_text2, outlier_text)
%
% seed_mesh:      resolution of seed regions, e.g. 'fsaverage5'
% target:         resolution of target resions, e.g. 'fsaverage3'
% output_file1:   left hemisphere output correlation profile
% output_file2:   right hemisphere output correlation profile
% threshold:      relative correlation threshold for binarization
% varargin_text1: the text file containing left hemisphere surface data list
% varargin_text2: the text file containing right hemisphere surface data list
% outlier_text:   the text file containing outlier file naame
%
% Compute surface correlation profiles.
%
% Written by CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md


lh_seed_avg_mesh = CBIG_ReadNCAvgMesh('lh', seed_mesh, 'inflated', 'cortex');
rh_seed_avg_mesh = CBIG_ReadNCAvgMesh('rh', seed_mesh, 'inflated', 'cortex');

lh_avg_mesh = CBIG_ReadNCAvgMesh('lh', target, 'inflated', 'cortex');
rh_avg_mesh = CBIG_ReadNCAvgMesh('rh', target, 'inflated', 'cortex');

% read in both left and right text files.
fid = fopen(varargin_text1, 'r');
i = 0;
while(1);
   tmp = fscanf(fid, '%s\n', 1);
   if(isempty(tmp))
       break
   else
       i = i + 1;
       varargin1{i} = tmp;
   end
end
fclose(fid);

% read in both left and right text files.
fid = fopen(varargin_text2, 'r');
i = 0;
while(1);
   tmp = fscanf(fid, '%s\n', 1);
   if(isempty(tmp))
       break
   else
       i = i + 1;
       varargin2{i} = tmp;
   end
end
fclose(fid);

% read in outlier text files
if(nargin == 8)
    fid = fopen(outlier_text, 'r');
    i = 0;
    while(1);
        tmp = fscanf(fid, '%s\n', 1);
        if(isempty(tmp))
            break
        else
            i = i + 1;
            outlierin{i} = tmp;
        end
    end
    fclose(fid);
end

% Compute profile for left hemi
for i = 1:length(varargin1)
    if(exist('outlierin', 'var'))
        outliers = dlmread(outlierin{i});                                   % {0,1} vector, where uncensored time points are 0
    end
    
    input = varargin1{i};
    input_series = MRIread(input);
    t_series = transpose(reshape(input_series.vol, size(input_series.vol, 1) * size(input_series.vol, 2) * size(input_series.vol, 3), size(input_series.vol, 4)));
    if(exist('outliers', 'var'))
        t_series = t_series(outliers==1, :);
    end
    
    hemi_index = strfind(input, basename(input));
    lh_seed_file = input; lh_seed_file(hemi_index:hemi_index+1) = 'lh';
    lh_seed_series = MRIread(lh_seed_file);
    lh_seed_series = transpose(reshape(lh_seed_series.vol, size(lh_seed_series.vol, 1) * size(lh_seed_series.vol, 2) * size(lh_seed_series.vol, 3), size(lh_seed_series.vol, 4)));
    lh_seed_series = lh_seed_series(:, lh_avg_mesh.MARS_label(1:length(lh_seed_avg_mesh.MARS_label)) == 2);
    if(exist('outliers', 'var'))
        lh_seed_series = lh_seed_series(outliers==1, :);
    end
    
    rh_seed_file = input; rh_seed_file(hemi_index:hemi_index+1) = 'rh';
    rh_seed_series = MRIread(rh_seed_file);
    rh_seed_series = transpose(reshape(rh_seed_series.vol, size(rh_seed_series.vol, 1) * size(rh_seed_series.vol, 2) * size(rh_seed_series.vol, 3), size(rh_seed_series.vol, 4)));
    rh_seed_series = rh_seed_series(:, rh_avg_mesh.MARS_label(1:length(rh_seed_avg_mesh.MARS_label)) == 2);
    if(exist('outliers', 'var'))
        rh_seed_series = rh_seed_series(outliers==1, :);
    end
    
    s_series = [lh_seed_series rh_seed_series];

    % normalize series (note that series are now of dimensions: T x N)
    s_series = bsxfun(@minus, s_series, mean(s_series, 1));
    s_series = bsxfun(@times, s_series, 1./sqrt(sum(s_series.^2, 1)));

    t_series = bsxfun(@minus, t_series, mean(t_series, 1));
    t_series = bsxfun(@times, t_series, 1./sqrt(sum(t_series.^2, 1)));

    corr_mat1 = s_series' * t_series;
    if(i == 1)
        output = corr_mat1;
    else
        output = output + corr_mat1;
    end
    clear outliers
end
output = output / length(varargin1);
corr_mat1 = output;
disp(['isnan: ' num2str(sum(isnan(corr_mat1(:)))) ' out of ' num2str(numel(corr_mat1))]);
corr_mat1(isnan(corr_mat1)) = 0;

% Compute profile for right hemi
for i = 1:length(varargin2)
    if(exist('outlierin', 'var'))
        outliers = dlmread(outlierin{i});                                   % {0,1} vector, where uncensored time points are 0
    end
    
    input = varargin2{i};
    input_series = MRIread(input);
    t_series = transpose(reshape(input_series.vol, size(input_series.vol, 1) * size(input_series.vol, 2) * size(input_series.vol, 3), size(input_series.vol, 4)));
    if(exist('outliers', 'var'))
        t_series = t_series(outliers==1, :);
    end
    
    hemi_index = strfind(input, basename(input));
    lh_seed_file = input; lh_seed_file(hemi_index:hemi_index+1) = 'lh';
    lh_seed_series = MRIread(lh_seed_file);
    lh_seed_series = transpose(reshape(lh_seed_series.vol, size(lh_seed_series.vol, 1) * size(lh_seed_series.vol, 2) * size(lh_seed_series.vol, 3), size(lh_seed_series.vol, 4)));
    lh_seed_series = lh_seed_series(:, lh_avg_mesh.MARS_label(1:length(lh_seed_avg_mesh.MARS_label)) == 2);
    if(exist('outliers', 'var'))
        lh_seed_series = lh_seed_series(outliers==1, :);
    end
    
    rh_seed_file = input; rh_seed_file(hemi_index:hemi_index+1) = 'rh';
    rh_seed_series = MRIread(rh_seed_file);
    rh_seed_series = transpose(reshape(rh_seed_series.vol, size(rh_seed_series.vol, 1) * size(rh_seed_series.vol, 2) * size(rh_seed_series.vol, 3), size(rh_seed_series.vol, 4)));
    rh_seed_series = rh_seed_series(:, rh_avg_mesh.MARS_label(1:length(rh_seed_avg_mesh.MARS_label)) == 2);
    if(exist('outliers', 'var'))
        rh_seed_series = rh_seed_series(outliers==1, :);
    end
    
    s_series = [lh_seed_series rh_seed_series];

    % normalize series (note that series are now of dimensions: T x N)
    s_series = bsxfun(@minus, s_series, mean(s_series, 1));
    s_series = bsxfun(@times, s_series, 1./sqrt(sum(s_series.^2, 1)));

    t_series = bsxfun(@minus, t_series, mean(t_series, 1));
    t_series = bsxfun(@times, t_series, 1./sqrt(sum(t_series.^2, 1)));

    corr_mat2 = s_series' * t_series;
    if(i == 1)
        output = corr_mat2;
    else
        output = output + corr_mat2;
    end
    clear outliers
end
output = output / length(varargin2);
corr_mat2 = output;
disp(['isnan: ' num2str(sum(isnan(corr_mat2(:)))) ' out of ' num2str(numel(corr_mat2))]);
corr_mat2(isnan(corr_mat2)) = 0;

% combine both hemisphere and threshold
tmp = [corr_mat1 corr_mat2];
tmp = sort(tmp(:), 'descend');
t = tmp(round(numel(tmp) * str2num(threshold)));
disp(['threshold: ' num2str(t)]);
if(str2num(threshold) < 1)
  corr_mat1(corr_mat1 <  t) = 0;
  corr_mat1(corr_mat1 >= t) = 1;
  corr_mat2(corr_mat2 <  t) = 0;
  corr_mat2(corr_mat2 >= t) = 1;
end

% tmp = sort(corr_mat1(:), 'descend');
% t1 = tmp(round(numel(tmp) * str2num(threshold)));
% disp(['threshold: ' num2str(t1)]);
% 
% tmp = sort(corr_mat2(:), 'descend');
% t2 = tmp(round(numel(tmp) * str2num(threshold)));
% disp(['threshold: ' num2str(t2)]);
% if(str2num(threshold) < 1)
%     corr_mat1(corr_mat1 <  t1) = 0;
%     corr_mat1(corr_mat1 >= t1) = 1;
%     corr_mat2(corr_mat2 <  t2) = 0;
%     corr_mat2(corr_mat2 >= t2) = 1;
% end

output_series = input_series;
size_input = size(input_series.vol);
size_input(4) = size(corr_mat1, 1);
output_series.vol = reshape(transpose(corr_mat1), size_input);
MRIwrite(output_series, output_file1); 

output_series = input_series;
size_input = size(input_series.vol);
size_input(4) = size(corr_mat2, 1);
output_series.vol = reshape(transpose(corr_mat2), size_input);
MRIwrite(output_series, output_file2);
