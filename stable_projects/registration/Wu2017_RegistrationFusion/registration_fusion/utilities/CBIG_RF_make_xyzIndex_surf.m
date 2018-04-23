function CBIG_RF_make_xyzIndex_surf(reg_dir, output_dir, output_prefix)
% CBIG_make_RF_xyzIndex_surf(reg_dir, output_dir, output_prefix, faces_num)
%
% This function creates index files using a subject's lh.sphere.reg and 
% rh.sphere.reg from reg_dir
% The index files created are curvature files containing a vector of x/y/z 
% coordinates at each vertex in the subject's surface space, as mapped from 
% fsaverage spherical surface
%
% Input:
%     - reg_dir      :
%                      absolute/relative path to the input sphere.reg file
%     - output_dir   :
%                      absolute/relative path to directory where output should be stored
%     - output_prefix:
%                      desired prefix for the outputs
%
% Output:
%     - There is no function output.
%     - 6 index files are created in output_dir:
%           lh.xIndex_[output_prefix].index
%           rh.xIndex_[output_prefix].index
%           lh.yIndex_[output_prefix].index
%           rh.yIndex_[output_prefix].index
%           lh.zIndex_[output_prefix].index
%           rh.zIndex_[output_prefix].index
%
% Example:
% CBIG_make_RF_xyzIndex_surf('~/data/Sub0001_Ses1_FS/surf', 'results/index_surf/', 
%                       'fsaverage_to_Sub0001_Ses1', )
% This command generates index files containing fsaverage spherical coordinates 
% mapped to Sub0001's surface
%
% Written by Wu Jianxiao and CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md


if nargin < 3
    disp('usage: CBIG_RF_make_xyzIndex_surf(reg, output_dir, output_prefix)');
    return
end

for hemis = {'lh', 'rh'}
    hemi = hemis{1};
    reg = [reg_dir '/' hemi '.sphere.reg'];
    
    %Load sphere.reg file and extract the x/y/z coordinates
    coord = read_surf(reg);
    x = coord(:, 1);
    y = coord(:, 2);
    z = coord(:, 3);
    
    if size(find(coord==0), 1) > 0
        warning(['The sphere.reg contains ' size(find(coord==0,1)) ...
            'x/y/z coordinates which are zero.']);
    end
    
    %Save output index files (do not care about faces numbers since we
    %do not need it when reading the file
    write_curv([output_dir '/' hemi '.xIndex_' output_prefix '.index'], x, 1);
    write_curv([output_dir '/' hemi '.yIndex_' output_prefix '.index'], y, 1);
    write_curv([output_dir '/' hemi '.zIndex_' output_prefix '.index'], z, 1);
end

end
