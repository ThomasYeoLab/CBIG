function CBIG_DownsampleDataFSaverage(hemi, input_file, input_mesh_file, output_mesh_file, output_file, var_name, exit_flag)

% CBIG_DownsampleDataFSaverage(hemi, input_file, input_mesh_file, output_mesh_file, output_file, var_name, exit_flag)
% CBIG_DownsampleDataFSaverage('lh', './lh.090722_ZC47CH.014.resid.mgz', 'fsaverage', 'fsaverage5', 'lh.blah.mat', 'fmri');
%
% Written by CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md



if(nargin < 7)
  exit_flag = 0;
else
  if(ischar(exit_flag))
    exit_flag = str2num(exit_flag);
  end
end

input = MRIread(input_file);
input_mesh = CBIG_ReadNCAvgMesh(hemi, input_mesh_file, 'sphere', 'cortex');
num_vertices = size(input_mesh.vertices, 2);
input_data = transpose(reshape(input.vol, num_vertices, numel(input.vol)/num_vertices));

output_mesh = CBIG_ReadNCAvgMesh(hemi, output_mesh_file, 'sphere', 'cortex');
num_out_vertices = size(output_mesh.vertices, 2);

if(num_out_vertices > num_vertices)
   error('Does not handle larger sampling output'); 
end

if(num_out_vertices == num_vertices)
    out_data = input_data;
else
    % Map out Closest Vertex
    index = MARS_findNV_kdTree(input_mesh.vertices, output_mesh.vertices);

    % Downsample data
    out_data = zeros(size(input_data, 1),  num_out_vertices);
    for i = 1:num_out_vertices
        if(output_mesh.MARS_label(i) == 2) %if output vertex is cortex
            out_data(:, i) = mean(input_data(:, index == i & input_mesh.MARS_label == 2), 2);
        else
            out_data(:, i) = input_data(:, i);
        end
    end
end

% save
eval([var_name ' = out_data;']);
save(output_file, var_name);

if(exit_flag)
  exit;
end
