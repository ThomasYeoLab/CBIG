function CBIG_ASDf_gamma2table(gamma_file, output_name)
% CBIG_ASDf_gamma2table(gamma_file, output_name)
%
% This function normalizes estimated gamma to get Pr(Factor|Participant)
%  (i.e., factor compositions of participants),
% and write to a text file specified by output_name.
%
% Input:
%     - gamma_file:
%           Gamma file name including the full path
%     - output_name:
%           Output factor composition file name including the full path
%
% Example:
%       CBIG_ASDf_gamma2table('~/example_output/estimate/k3/r10/final.gamma',
%       '~/example_output/visualizeFactors/k3/r10/factorComp.txt');
%
% Written by CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md

gamma = load(gamma_file);

% Normalize to a probability distribution
gamma_norm = bsxfun(@times, gamma, 1./(sum(gamma, 2)));

dlmwrite(output_name, gamma_norm, 'delimiter', ' ');
