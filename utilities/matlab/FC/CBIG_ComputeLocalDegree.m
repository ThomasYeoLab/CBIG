function CBIG_ComputeLocalDegree(input_file, output_file, neighbors_file, threshold)

% CBIG_ComputeLocalDegree(input_file, output_file, neighbors_file, threshold)
% This function is used to compute the local degree with threshold for each
% vertex.
% Input:
%      -input_file:
%       Can be either a .mat file or a .nii/nii.gz file contains fMRI data.
%      -output_file:
%       output .mat file
%      -neighbors_file:
%       a .mat file contains neighborhood_cell, which is a Nx1 cell, each
%       row corresponds to the indices of neighborhood vertices.
%      -threshold:
%       threshold used to compute degree 
%
% Written by CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md


load(neighbors_file, 'neighborhood_cell');
if(strfind(input_file, '.mat'))
   load(input_file, 'fmri'); 
else
   x = MRIread(input_file);
   fmri = reshape(x.vol, length(neighborhood_cell), size(x.vol, 4));
   fmri = fmri';
end

if(ischar(threshold))
   threshold = str2num(threshold); 
end

% normalize series (note that series are now of dimensions: T x N)
fmri = fmri - repmat(mean(fmri, 1), size(fmri, 1), 1);
fmri = fmri./repmat(sqrt(sum(fmri.^2, 1)), size(fmri, 1), 1);

num_vertices = size(fmri, 2);
frac_degree = zeros(1, num_vertices);
for i = 1:num_vertices
    t_series = fmri(:, neighborhood_cell{i});
    corr_vec = sum(repmat(fmri(:, i), 1, size(t_series, 2)) .* t_series, 1);
    frac_degree(i) = sum(corr_vec > threshold)/length(corr_vec);
end

if(strfind(output_file, '.mat'))
   save(output_file, 'frac_degree');
else
   if(length(frac_degree) == 10242)
       write_curv(output_file, frac_degree, int32(20480));
   else
       error(['Does not handle ' num2str(length(frac_degree)) ' vertices']);
   end  
end
exit


