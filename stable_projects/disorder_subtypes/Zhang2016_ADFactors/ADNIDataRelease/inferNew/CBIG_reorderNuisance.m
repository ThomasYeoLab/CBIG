function CBIG_reorderNuisance(brainListFile, nuisanceVarsFile, concatOrderFile)

% CBIG_reorderNuisance(brainListFile, nuisanceVarsFile, concatOrderFile)
%
% Helper function that reorders the lines in nuisanceVars so that now they
% match the line order of concatOrder. Before this function, they match the
% line order of brainList instead.
%
% Written by Xiuming Zhang and CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md

% Load brainList
fID = fopen(brainListFile);
brainList = textscan(fID, '%s');
fclose(fID);
brainList = brainList{1};

nuisanceVars = csvread(nuisanceVarsFile);

% Load concatOrder
fID = fopen(concatOrderFile);
concatOrder = textscan(fID, '%s');
fclose(fID);
concatOrder = concatOrder{1};

% Fetch nuisance variables for each image
% so that nuisanceVars and concatOrder correspond
nuisanceVars_reordered = nuisanceVars;
for idx = 1:numel(concatOrder)
    fName = concatOrder{idx};
    logIdx = ~cellfun(@isempty, strfind(brainList, fName));
    assert(sum(logIdx)==1);
    nuisanceVars_reordered(idx, :) = nuisanceVars(logIdx, :);
end
csvwrite([nuisanceVarsFile(1:end-4) '_matchConcatOrder.csv'], nuisanceVars_reordered);
