function [subcort_index, non_subcort_index] = CBIG_GetSubcortIndices

% Get Freesurfer subcortical indices.
% Reference: $FREESURFER_HOME/FreeSurferColorLUT.txt
%
%   [subcort_index, non_subcort_index] = CBIG_GetSubcortIndices
%   Output:
%       subcort_index       : column vector, subcortical index
%       non_subcort_index   : column vector, non subcortical index
%
% Written by CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md



% gray
subcort_index = [8; 9; 10; 11; 12; 13; 17; 18; 26; 27; 28; ...
                47; 48; 49; 50; 51; 52; 53; 54; 58; 59; 60];

non_subcort_index = [0; 2; 3; 4; 5; 7; 14; 15; 16; 24; ...
                    30; 31; 41; 42; 43; 44; 46; 62; 63; 72; 77; 80; 85; ...
                    251; 252; 253; 254; 255; transpose(1000:1035); transpose(2000:2035)];
