function p = CBIG_lr_test(X_r, y, n, ll_u, dof)

if isempty(X_r) % only constant term
    mu = sum(y==1)/numel(y); % ML estimate
    ll_r = sum(y==1)*log(mu)+sum(y==0)*log(1-mu);
else
    [b_r, ~, ~] = glmfit(X_r, [y n], 'binomial', 'link', 'logit'); % log(µ/(1-µ)) = Xb
    yFit = glmval(b_r, X_r, 'logit', 'size', n);
    ll_r = sum(log(binopdf(y, n, yFit./n)));
end

% LR test
[~, p] = lratiotest(ll_u, ll_r, dof);