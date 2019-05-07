% Written by Nanbo Sun and CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md
% normalize gamma file to probabilities

gamma = load(gamma_file);
gamma_norm = bsxfun(@times, gamma, 1./(sum(gamma, 2)));
gamma_norm = gamma_norm(:, [1, 3, 2]);
dlmwrite(prob_file, gamma_norm, 'delimiter', ',');
