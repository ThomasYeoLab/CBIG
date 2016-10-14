function CBIG_RemoveNeighborhoodFromProfile(hemi, profile_file, output_profile, radius)

% CBIG_RemoveNeighborhoodFromProfile(hemi, profile_file, output_profile, radius)
% This function is used to remove neighborhood vertices within radius for
% each vertex in correlation profile.
% hemi: 'lh' or 'rh'
% profile_file: correlation profile
% radius: radius used to find neighborhood vertices.
% Example:
% CBIG_RemoveNeighborhoodFromProfile('lh', profile_file, radius)
%
% Written by CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md


% load profile
x = MRIread(profile_file);
profile = reshape(x.vol, 10242, size(x.vol, 4));

% load neighborhood cell
neighborhood_cell=CBIG_ComputeNeighborhood(hemi, 'fsaverage5', radius, 'white');

% load mask
mesh_name = 'fsaverage5';
mask = 'cortex';
lh_avg_mesh = CBIG_ReadNCAvgMesh('lh', mesh_name, 'inflated', mask);
l1 = find(lh_avg_mesh.MARS_label(1:642) == 2); 
rh_avg_mesh = CBIG_ReadNCAvgMesh('rh', mesh_name, 'inflated', mask);
l2 = find(rh_avg_mesh.MARS_label(1:642) == 2); 

for i = 1:length(neighborhood_cell)
   
   index = zeros(1, length(l1)+length(l2));
   
   neighbors = neighborhood_cell{i};
   
   if(strcmp(hemi, 'lh'))
       [c, IA] = intersect(l1, neighbors);
       index(IA) = 1; 
   else
       [c, IA] = intersect(l2, neighbors);
       index(IA+length(l1)) = 1;  
   end
       
   profile(i, logical(index)) = nan;
end

x.vol = reshape(profile, size(x.vol));
MRIwrite(x, output_profile);
exit
