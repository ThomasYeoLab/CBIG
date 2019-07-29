function X = myPLS_norm(X,num_groups,subj_grouping,mode)

% Normalization

switch mode
    case 1
        X = zscore(X);
        
    case 2
        for iter_group = 1:num_groups
            idx = find(subj_grouping == iter_group);
            X(idx,:) = zscore(X(idx,:));
        end
        
    case 3
        X2 = sqrt(mean(X.^2,1));
        X = X./ repmat(X2,[size(X,1) 1]);
        
    case 4
        for iter_group = 1:num_groups
            idx = find(subj_grouping == iter_group);
            X2 = sqrt(mean(X(idx,:).^2,1));
            X(idx,:) = X(idx,:)./ repmat(X2,[size(X(idx,:),1) 1]);
        end
end
