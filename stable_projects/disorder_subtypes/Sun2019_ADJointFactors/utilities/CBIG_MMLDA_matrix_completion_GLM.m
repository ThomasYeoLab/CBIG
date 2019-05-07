function M_recover = CBIG_MMLDA_matrix_completion_GLM(M, M_train)

% M_recover = CBIG_MMLDA_matrix_completion_GLM(M, M_train)
%
% Use GLM to do matrix completion for missing entries by using training matrix. 
% 
% Input:
%   - M         : N x B matrix with missing entries (NaN). N is the # of subjects. B is # of features.
%   - M_train   : N_train x B matrix without missing entries. N_train is # of training subjects.
%
% Output:
%   - M_recover : N x B recovered matrix without missing entries (NaN).
%
% Written by Nanbo Sun and CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md

M_recover = M;
row_miss_ind = find(isnan(sum(M, 2)) == 1);
for i = 1:length(row_miss_ind)
    row_score = M(row_miss_ind(i), :);
    col_miss_ind = find(isnan(row_score) == 1);
    for j = 1:length(col_miss_ind)
        X = [ones(size(M_train, 1), 1) M_train(:, ~isnan(row_score))];
        Y = M_train(:, col_miss_ind(j));
        b = (X'*X)\(X'*Y);
        M_recover(row_miss_ind(i), col_miss_ind(j)) = [1 row_score(:, ~isnan(row_score))] * b;
    end
end

end