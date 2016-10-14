function PI = CBIG_ComputeSoftParticipationIndex(prob_mat, corr_mat)

% PI = CBIG_ComputeSoftParticipationIndex(prob_mat, corr_mat)
%
% prob_mat = T x N  (probability that location n appears in topic t)
% corr_mat = N x N  (correlation between location n and location m, it is like weighted adjacency matrix)
% The self correlation and negative correlation are neglected inside the code.
%
% "Soft" participation coefficient. 
%
% "Hard" participation coefficient is defined by
% 1 - sum_{t=1}^{T} ( (k_{n,t}/k_{n}) ^2 ), where k_{n,t} is the degrees to
% cluster (or topic in this scenario) of node n, k_{n} is the total degree
% of node n. Please see:
% Guimera, Roger, and Luis A. Nunes Amaral. "Functional cartography of
% complex metabolic networks." Nature 433, no. 7028 (2005): 895-900.
%
% In "soft" participation coefficient, k_{n,t} is defined as 
% sum_{m!=n} (prob_{m,t} * corr_{n,m}), which is to sum the correlation
% between this node and any other node weighted by the probability that the
% other node belongs to cluster (topic) t. k_{n} is still defined by
% summing k_{n,t} over all clusters (topics).
%
% Written by CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md



T = size(prob_mat, 1);
N = size(prob_mat, 2);
[a, b] = meshgrid(1:N, 1:N);
corr_mat(a == b) = 0;
clear a; clear b;
corr_mat(corr_mat <= 0) = 0;

prob_mat = bsxfun(@times, prob_mat, 1./(squeeze(sum(prob_mat, 1))+eps));

K  = sum(corr_mat, 1);
KM = zeros(T, N);
for i = 1:T
   KM(i, :) = sum(bsxfun(@times, corr_mat, prob_mat(i, :)), 2)';
end

PI = 1 - sum(bsxfun(@times, KM, 1./K).^2, 1);







