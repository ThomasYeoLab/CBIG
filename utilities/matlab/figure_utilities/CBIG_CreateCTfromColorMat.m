function ct = CBIG_CreateCTfromColorMat(color_mat, struct_names)

% ct = CBIG_CreateCTfromColorMat(color_mat, struct_names)
%
% This function is used to generate 'MARS_ct' struct field for the mesh
% structure which can be obtained by CBIG_ReadNCAvgMesh.m
% Input:
%      -color_mat: 
%       Nx3 matrix where each coloum corresponds to r,g and b. N is the 
%       number of parcels, each parcel has a specific RGB value.
%
%      -struct_names:
%       Nx1 cell, each row corresponds to the structure name for each
%       parcel. (e.g. central sulcus and so on)
%
% Output:
%      -ct:
%       ct is a struct
%       ct.numEntries = number of entries
%       ct.orig_tab = name of original ct
%       ct.struct_names = list of structure names (e.g. central sulcus and so on)
%       ct.table = n x 5 matrix. 1st column is r, 2nd column is g, 3rd column
%       is b, 4th column is flag, 5th column is resultant integer values
%       calculated from r + g*2^8 + b*2^16 + flag*2^24. flag expected to be all 0
%
% Example:
% ct = CBIG_CreateCTfromColorMat(color_mat, struct_names_mat)
%
% Written by CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md


ct.numEntries = size(color_mat, 1);
ct.orig_tab = 'MyColorLUT';

if(nargin == 2)
   ct.struct_names = struct_names; 
else
   for i = 1:ct.numEntries 
       ct.struct_names{i} = ['NONAME' num2str(i-1)]; 
   end
end

ct.table = zeros(ct.numEntries, 5);
for i = 1:ct.numEntries
  r = color_mat(i, 1);
  g = color_mat(i, 2);
  b = color_mat(i, 3);
  
  ct.table(i, :) = [r g b 0 r + g*2^8 + b*2^16];
end



