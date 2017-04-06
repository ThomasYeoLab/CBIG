function bestModel = CBIG_visualizeFactors(initFolders, GMMask, outDir)

% bestModel = CBIG_visualizeFactors(initFolders, GMMask, outDir)
%
% This function selects the best random initilization (that gives the
% highest likelihood) and visualizes the estimated factors as probabilistic
% atrophy maps.
%
% Input:
%     - initFolders:
%       Folders of all the random initializations; e.g., ~/R*
%     - GMMask:
%       Gray matter mask used in VBM
%     - outDir:
%       Output directory; e.g., ~/
%
% Output:
%     - bestModel:
%       Path to the best model estimated; for future inferences
%     - Factors and factor compositions from the best initialization:
%       Factors are in the form of probabilistic atrophy maps 
%   
% Example:
% See ../replicatePNAS/CBIG_LDA_wrapper.m
%
% Written by Xiuming Zhang and CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md

mkdir(outDir);
initDir = initFolders(1:find(initFolders=='/', 1, 'last'));

% Pick the initialization with maximum likelihood
runName = findMaxLikeRun(initFolders, outDir);
bestModel = [initDir runName '/final'];

% Plot the climbing of log-likelihood as this random initialization runs
% The log-likelihood should "go flat" if LDA has converged
plotLoglike(initDir, runName, outDir);

% Check if the current number of initializations are enough
% If many initializations yield similar (high correlations) results, enough
plotCorrWithBest(initFolders, runName, outDir);

% Write Pr(Voxel | Factor), i.e., exp(beta) in LDA, to 3D brain volumes
beta2brain(GMMask, initDir, runName, outDir);

% Write Pr(Factor | Subject), i.e., normalized gamma in LDA, to txt, where each row is a subject
gamma2table(initDir, runName, outDir);


function maxRunName = findMaxLikeRun(runFolders, outDir)
folderList = dir(runFolders);
runDir = runFolders(1:find(runFolders=='/', 1, 'last'));
maxLogLike = -inf;
figure;
hold on;
for idx = 1:numel(folderList)
    runName = folderList(idx).name;
    try
        logLikes = load([runDir runName '/likelihood.dat']);
        logLike = logLikes(end, 1);
        scatter(idx, logLike, 'b.');
        if logLike > maxLogLike
            maxLogLike = logLike;
            maxRunName = runName;
            maxRunID = idx;
        end
    catch
        warning([runName ': no likelihood.dat']);
    end
end
scatter(maxRunID, maxLogLike, 'ro');
hold off;
ylabel('Log-Likelihood');
xlabel('Random Initialization');
saveas(gcf, [outDir 'bestInit.png']);


function plotLoglike(runDir, maxRunName, outDir)
log_likelihood = load([runDir maxRunName '/likelihood.dat']);
figure;
plot(log_likelihood(:, 1), 'o-');
xlabel('Iteration');
ylabel('Log-Likelihood');
grid on;
saveas(gcf, [outDir 'convergence.png']);


function plotCorrWithBest(runFolders, maxRunName, outDir)
runDir = runFolders(1:find(runFolders=='/', 1, 'last'));
folderList = dir(runFolders);
noRuns = numel(folderList);
% Rank by likelihood
logLike = zeros(noRuns, 1);
for idx = 1:noRuns
    runName = folderList(idx).name;
    logLikes = load([runDir runName '/likelihood.dat']);
    logLike(idx) = logLikes(end, 1);
end
[~, ind_h2l] = sort(logLike, 'descend');
avgCorr = zeros(noRuns, 1);
for idx = 1:noRuns
    runName = folderList(ind_h2l(idx)).name;
    avgCorr(idx) = hunMatch_sameK(runDir, maxRunName, runName);
end
% Plot
figure;
plot(1:noRuns, avgCorr, '-o');
box off;
xlabel('Random Initializations');
ylabel('Correlation with the Best Run');
xlim([1, noRuns]);
ylim([0, 1]);
set(gca, 'YTick', 0:0.05:1);
grid on;
saveas(gcf, [outDir '/corrWithBest.png']);


function avgCorr = hunMatch_sameK(runDir, maxRunName, runName)
beta_best = exp(load([runDir maxRunName '/final.beta']));
beta = exp(load([runDir runName '/final.beta']));
%%% Reorder subtypes to obtain the maximal correlation coefficients
% Construct the COST matrix
%         pos1 pos2 ...
% topic1
% topic2
% ...
K = size(beta, 1);
costMat = zeros(K, K);
for rowIdx = 1:K
    for colIdx = 1:K
        % Assign beta (jobs, column) to beta_best (workers, row)
        corrMat = corrcoef(beta_best(rowIdx, :)', beta(colIdx, :)');
        costMat(rowIdx, colIdx) = 1-corrMat(1, 2);
    end
end
% Run the Hungarian matching algorithm
% order: each row (worker)'s matched column (job)
[order, ~] = munkres(costMat);
% Recompute the avergae correlation with sorted topics
corr = zeros(K, 1);
for idx = 1:K  
	corrMat = corrcoef(beta_best(idx, :)', beta(order(idx), :)');
    corr(idx) = corrMat(1, 2); 
end
avgCorr = mean(corr);


function beta2brain(GMMask, runDir, maxRunName, outDir)
mask = MRIread(GMMask);
mask = mask.vol;
mask_1d = reshape(mask, [1 numel(mask)]);
beta = load([runDir maxRunName '/final.beta']);
% For each topic
for topicIdx = 1:size(beta, 1)
    betaRow = exp(beta(topicIdx, :));
    % Convert back to 3D
    betaRow_1d = mask_1d; % insert 0's
    betaRow_1d(betaRow_1d==1) = betaRow; % 1's to real values
    betaRow_3d = reshape(betaRow_1d, size(mask));
    % Write values to MRI
    dummyMRI = MRIread(GMMask);
    dummyMRI.vol = betaRow_3d;
    MRIFilename = [outDir '/factor' num2str(topicIdx) '.nii.gz'];
    MRIwrite(dummyMRI, MRIFilename);
end


function gamma2table(runDir, maxRunName, outDir)
gamma = load([runDir maxRunName '/final.gamma']);
% Normalize to a probability distribution
gamma_norm = bsxfun(@times, gamma, 1./(sum(gamma, 2)));
dlmwrite([outDir 'factorCompositions.txt'], gamma_norm, 'delimiter', ' ');
