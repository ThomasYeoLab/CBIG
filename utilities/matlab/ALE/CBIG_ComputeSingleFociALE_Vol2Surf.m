function p = CBIG_ComputeSingleFociALE_Vol2Surf(surface_ras, foci_coor, res, sig)

% p = CBIG_ComputeSingleFociALE_Vol2Surf(surface_ras, foci_coor, res, sig)
%
% Compute the ALE scores of a given set of surface locations with respect to of a foci
%
% surface_ras  = 3 x N matrix, which is the N 3-D locations on the brain surface
% foci_coord   = 3 x 1 vector, 3-D location of an activation foci
% res          = volume of a brain voxel, set as 8mm^3 if no argument is supplied
% sig          = standard deviation of the Gaussian spatial probabilitiy distribution assumed at an
%                activation foci
% 
% p            = 1 x N vector of ALE scores of N surface locations
% __________________________________________________________________________________________________
%
% CBIG_ComputeSingleFociALE_Vol2Surf returns the ALE scores of a given set of surface locations given
% an activation foci, as described in Eickhoff et. al, 2009.
% ____________________________________________________________________________________________________
%
% Refs:
%
% Eickhoff et al., Human Brain Mapping, 2009. Coordinate-based activation likelihood estimation meta-analysis
% of neuroimaging data: A random-effects approach based on empirical estimates of spatial uncertainty
%
% Written by B.T.Thomas Yeo and CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md

  if size(foci_coor,2) ~= 1
      error('Input argument ''foci_coor'' should be a column vector');
  end

  if(~exist('res', 'var'))
     res = 8; % assume 2 x 2 x 2 mm
  end 

  if(~exist('sig', 'var'))
    sig = 12/(2*sqrt(2*log(2))); % 12mm fwhm
  end


  d = bsxfun(@minus, surface_ras, foci_coor);
  d2 = sum(d.^2, 1);

  p = exp(-d2/(2*(sig^2)))/(((2*pi)^(3/2))*(sig^3))*res;
