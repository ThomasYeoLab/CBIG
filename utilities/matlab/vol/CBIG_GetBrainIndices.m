function [brain_index, non_brain_index] = CBIG_GetBrainIndices

% Get Freesurfer brain indices.
% Reference: $FREESURFER_HOME/FreeSurferColorLUT.txt
%
%   [brain_index, non_brain_index] = CBIG_GetBrainIndices
%   Output:
%       brain_index : column vector, brain index
%       non_brain_index : column vector, non brain index
%
% Written by CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md



% gray
brain_index = [2; 3; 7; 8; 9; 10; 11; 12; 13; 17; 18; 26; 27; 28; 41; 42; ...
    46; 47; 48; 49; 50; 51; 52; 53; 54; 58; 59; 60; 77; 80; ...
    251; 252; 253; 254; 255; transpose(1000:1035); transpose(2000:2035)];

non_brain_index = [0; 4; 5; 14; 15; 16; 24; 30; 31; 43; 44; 62; 63; 72; 85];
