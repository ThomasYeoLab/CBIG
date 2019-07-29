function R = myPLS_cov(X,Y,num_groups,subj_grouping)

% Compute PLS cross-covariance matrix (stacked)

for iter_group = 1:num_groups
    Ysel = Y(find(subj_grouping==iter_group),:);
    Xsel = X(find(subj_grouping==iter_group),:);
    
    R0 = Ysel.'*Xsel;

    if ~exist('R')
        R = R0;
    else
        R = [R; R0];
    end
end

