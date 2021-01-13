function y = wtmean(mat,weights)
%weighted mean(with NaN values removed)
    wmat = bsxfun(@times,~isnan(mat),weights);
    mat(isnan(mat)) = 0;
    swmat = mat .* wmat;
    y = sum(swmat(weights~=0,:),1)./sum(wmat,1);
end