function flag = CBIG_CheckValidSurfaceData(surface_template, lh_data, rh_data, ref_lh_vertices, ref_rh_vertices)

% flag = CBIG_CheckValidSurfaceData(surface_template, lh_data, rh_data, ref_lh_vertices, ref_rh_vertices)
%
% Check if the number of vertices in the input annotation match that of
% the reference vertices and the surface template in both hemispheres
% 
% Input:
%     - lh_data              : column vector containing surface data of
%                              the left hemisphere.
%     - rh_data              : column vector containing surface data of
%                              the right hemisphere.
%     - ref_lh_vertices      : column vector containing reference vertices
%                              of the left hemisphere.
%     - ref_rh_vertices      : column vector containing reference vertices
%                              of the right hemisphere
%     - surface_template     : name of the surface space of the input
%                              data.`surface_space` needs to correspond to
%                              number of vertices of the annotation in
%                              `lh_data`, `rh_data` and the reference annotations.
%        surface_template    |      number of vertices
%          fsaverage5        |            10242
%          fsaverage6        |            40962
%          fsaverage         |            163842
%          fs_LR_32k         |            32492
%          fs_LR_164k        |            163842
% Example:
%  flag = CBIG_CheckValidSurfaceData('fsaverage5', lh_data, rh_data, ref_lh_vertices, ref_rh_vertices)
%
%  Check if the length of vectors lh_data, rh_data, ref_lh_vertices and
%  ref_rh_vertices equal to each other and equal to 10242 (number of vertices
%  in fsaverage5 template). flag equals to 1 if true, and 0 otherwise.
%
% Written by Gia H. Ngo and CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md
    
    if size(lh_data,2) ~= 1
        error('Input argument ''lh_data'' should be a column vector');
    end
    if size(rh_data,2) ~= 1
        error('Input argument ''rh_data'' should be a column vector');
    end
    if size(ref_lh_vertices,2) ~= 1
        error('Input argument ''ref_lh_vertices'' should be a column vector');
    end
    if size(ref_rh_vertices,2) ~= 1
        error('Input argument ''ref_rh_vertices'' should be a column vector');
    end
  
    flag = numel(lh_data) == numel(rh_data);
    flag = flag && (numel(lh_data) == numel(ref_lh_vertices));
    flag = flag && (numel(ref_lh_vertices) == numel(ref_rh_vertices));
    
    switch(surface_template)
    case('fsaverage5')
        flag = flag && (numel(lh_data) == 10242);
    case ('fsaverage6')
        flag = flag && (numel(lh_data) == 40962);
    case('fsaverage')
        flag = flag && (numel(lh_data) == 163842);
    case('fs_LR_32k')
        flag = flag && (numel(lh_data) == 32492);
    case('fs_LR_164k')
        flag = flag && (numel(lh_data) == 163842);
    otherwise
        disp(['Invalid surface template: ' surface_template]);
        return;
    end
    
