function structName_avgWinProb = CBIG_sortByAvgProbWinningFactor(label_structName_avgWinProb)

avgWinProb = cell2mat(label_structName_avgWinProb(:, 3));

[~, ind] = sort(avgWinProb, 'descend');

structName_avgWinProb = label_structName_avgWinProb(ind, [2 3]);
