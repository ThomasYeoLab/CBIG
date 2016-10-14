function [grpRes indRes] = discover(dat,noClus,noReps,dim)

% [grpRes indRes] = discover(dat,noClus,noReps,dim)
% 
% NOTE: This function is a modified version of Danial Lashikari's code
%       discover.m, which can be found here:
%       http://people.csail.mit.edu/danial/Site/Code.html
%
%       Our modification includes:
%       - Removing two input arguments meanFlag and epsilon. Now meanFlag
%       is always non-zero.
%       - Removing re-running direcClus and choosing the best likelihood
%
% Performs the analysis on the data (variable dat) with number of clusters
% indicated by variable noClus. The analysis is repeated noReps number of
% times with random initializations.
%
% Format of variable dat: 
% The variable has a cell structure. For subject i,
% dat{i} is an N X D dimensional matrix representing the fMRI response of
% the N voxels within the mask for subject i, to D different stimuli or
% time points.
%
% grpRes: clustering results for the concatenated group data
%
% indRes: clustering results for individual subject data. indRes has a cell
%         structure where indRes{i} presents the result for subject i (see the
%         structure of variable dat below).
%
% dim:    the actual dimensionality of data. If ny linear constraint that we
%         add to the data reduces dim from the value D. For instance, we
%         the vectors are all zero mean, then dim = D-1. If we remove the
%         mean 
%
% discover also has a number of other default values introduced running the
% clustering code (direcClus.m). For more information look a the help for
% direcClus.


d = size(dat{1},2);

if ~exist('dim')
    dim = d;
end

yG = [];
for s = 1:length(dat)
	yG = [yG; datanor(dat{s},'rmmn')];
end

grpRes = direcClus(yG,noClus,dim,noReps,0,0,0,1e-3,1);

for s = 1:length(dat)
	
	%indRes{s} = direcClus(datanor(dat{s},'rmmn'),noClus,dim,1,resG.lambda,resG.p,resG.mtc,1e-3,1);
    indRes{s} = direcClus(datanor(dat{s},'rmmn'),noClus,dim,noReps,0,0,0,1e-3,1);

end

