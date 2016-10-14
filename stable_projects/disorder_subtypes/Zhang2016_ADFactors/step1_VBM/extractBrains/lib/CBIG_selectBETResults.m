function CBIG_selectBETResults(BETOutFolderPattern)

% CBIG_selectBETResults(BETOutFolderPattern)
%
% This function selects the best (equivalently, the latest) BET result for
% each image. For example, if Image 6 is present in BET rounds 1, 2 and 3,
% then we should take the round 3 result as its final BET result.
%
% Input:
%     - BETOutFolderPattern:
%	Path to your BET output folders, for example, '~/outputs/VBM/brains*'
%
% Output:
%     - A text file listing out the selected BET result for each scan
%
% Written by Xiuming Zhang and CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md

% Parse directory out
tmp = strsplit(BETOutFolderPattern, '/');
noBits = numel(tmp{end});
BETOutDir = BETOutFolderPattern(1:end-noBits);
[~, BETOutDir_abs] = system(['readlink ' BETOutDir ' -f']); % relative to absolute

dirs = dir(BETOutFolderPattern);
dirs = dirs([dirs.isdir]); % directories only
[~, sortInd] = sort([dirs.datenum], 'descend'); % sort by date: new folders first

% Starting from the latest brain folder
brainList = {};
for idx = 1:numel(sortInd)
    searchPattern = [BETOutDir dirs(sortInd(idx)).name '/*_brain.nii.gz'];
    files = dir(searchPattern);
    for file = files'
        % If brainList doesn't yet contain this image
        if isempty(strfind([brainList{:}], file.name))
            brainList = [brainList; strcat(BETOutDir_abs, '/', dirs(sortInd(idx)).name, '/', file.name)];
        end
    end
end

fileID = fopen([BETOutDir 'brainList.txt'], 'w');
for idx = 1:numel(brainList)
    fprintf(fileID, '%s\n', brainList{idx});
end
fclose(fileID);
