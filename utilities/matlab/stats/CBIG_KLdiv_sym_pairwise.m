function div = CBIG_KLdiv_sym_pairwise(prob_mat)

% Calculate pairwise symmetric KL divergence for each row of input matrix.
%
%     div = CBIG_KLdiv_sym_pairwise(prob_mat)
%     Input:
%         prob_mat: N x T matrix, N is # of variables, T is # of outcomes
%     Output:
%         div     : N x N matrix, N is # of variables
%
% Written by CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md


N = size(prob_mat, 1);
div = zeros(N, N);
for i = 1:N
    for j = 1:N
        div(i, j) = CBIG_KLdiv_sym(prob_mat(i, :), prob_mat(j, :));
    end
end



