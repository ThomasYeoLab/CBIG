function [cs matchInds scores] = performMatching(grpRes,indRes,matchingFlag)

% [cs matchInds scores] = performMatching(grpRes,indRes,matchingFlag)
%
% NOTE: This function is not exactly the same with Danial Lashikari's
%       function performMatching.m, which can be found here:
%       http://people.csail.mit.edu/danial/Site/Code.html
%
%       Our modification includes:
%       - Removing multiple cases of matchingFlag. Now only the case of
%       "matchingFlag == 0" is included.
%       - Removing two output arguments "mus" and "sps"
%
% Performs the matching of the group and individual clusters based on
% solving the bipartite graph matching problem between each individual and
% the group. grpRes and indRes are the group and individual results as
% created by discover.m.
%
% cs:           consistency scores for clustesrs
% matchInds:    matching indices of individual data to group clusters
% scores:       matching scores for each individual data

% matchingFlag specifies what kind of matching should be performed. The
% default is 0 which means matching based on the vectors of selectivity
% profiles. 

noSubjs = length(indRes);

mus = grpRes.mtc;
musnor = datanor(mus,'norm',2);

for i = 1:noSubjs

    sps(:,:,i) = indRes{i}.mtc;
    spsnor = datanor(sps(:,:,i),'norm',2);

    cors = musnor*spsnor';

    matching(:,:,i) = Hungarian(-cors);
    [temp matchInds(:,i)] = max(matching(:,:,i),[],2);

    sps_match(:,:,i) = sps(matchInds(:,i),:,i);

    scores(:,i) = diag(musnor*spsnor(matchInds(:,i),:)');

end

cs = mean(scores,2);


