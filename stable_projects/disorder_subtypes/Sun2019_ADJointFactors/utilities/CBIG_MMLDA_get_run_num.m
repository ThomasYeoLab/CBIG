function r = CBIG_MMLDA_get_run_num(in_dir, k)
% r = CBIG_MMLDA_get_run_num(in_dir, k)
%
% Get best run number based on visualization folder name.
%
% Input:
%   - in_dir: input directory which is the parent direcotry of "k*"
%   - k     : number of fators, e.g. 2, 3, 4
%
% Output:
%   - r     : number of run. e.g. 11, 20
%
% Example:
%   r = CBIG_MMLDA_get_run_num('~/visualizeFactors/', 3)
%
% Written by Nanbo Sun and CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md

rstruct = dir([in_dir sprintf('/k%s/r*', num2str(k))]);
run = rstruct.name;
r = str2num(run(2:end));
