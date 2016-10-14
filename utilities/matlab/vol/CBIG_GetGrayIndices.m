function [gray_index, non_gray_index] = CBIG_GetGrayIndices

% Get Freesurfer gray matter indices.
% Reference: $FREESURFER_HOME/FreeSurferColorLUT.txt
%
%   [gray_index, non_gray_index] = CBIG_GetGrayIndices
%   Output:
%       gray_index : column vector, gray matter index
%       non_gray_index : column vector, non gray matter index
%
% Written by CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md



% gray
gray_index = [3; 8; 9; 10; 11; 12; 13; 17; 18; 26; 27; 42; ...
    47; 48; 49; 50; 51; 52; 53; 54; 58; 59; 77; 80; ...
    251; 252; 253; 254; 255; transpose(1000:1035); transpose(2000:2035)];

non_gray_index = [0; 2; 4; 5; 7; 14; 15; 16; 24; 28; 30; 31; 41; 43; 44; 46; 60; 62; 63; 72; 85];