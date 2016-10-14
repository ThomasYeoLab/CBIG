function prob = CBIG_get_prob(gamma_file)

% prob = CBIG_get_prob(gamma_file)
% % Subtracting alpha (optional)
% % From Blei:
% % "we (sometimes) subtract alpha to find a less smooth representation
% % of the topics within each document.  if you are looking for an
% % estimate of \theta_{doc} then you don't have to subtract it.  simply
% % normalizing the gamma gives the posterior expectation of the topic
% % proportions."
% other_file = strrep(gamma_file, '.gamma', '.other');
% fileID = fopen(other_file);
% C = textscan(fileID, '%s %f');
% alpha = C{1, 2}(3);
%
% Written by CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md

gamma = load(gamma_file);


% fclose(fileID);
% gamma = gamma-alpha;

gamma_norm = bsxfun(@times, gamma, 1./(sum(gamma, 2)));

prob = gamma_norm;
