function rid_icv_vol = CBIG_getVol(ridList, structNames)

% rid_icv_vol = CBIG_getVol(ridList, structNames)
% Get baseline structure volumes from FS stats
%
% Written by CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md


imgDir = '//share/users/imganalysis/yeolab/data/ADNI_mert/';

rid_icv_vol = zeros(numel(ridList), 3);

for idx = 1:numel(ridList)
    rid = ridList(idx);
    
    %% Read in stats from files
    
    % ICV
    asegFile = [imgDir sprintf('%04d', rid) '/stats/aseg.stats'];
    fileID = fopen(asegFile);
    tbl = textscan(fileID, '%s', 'delimiter', '\n');
    fclose(fileID);
    idxICV = find(not(cellfun('isempty', strfind(tbl{1, 1}, 'Intracranial Volume'))));
    line = tbl{1, 1}{idxICV};
    strs = strsplit(line, ',');
    icv = str2double(strs(4));
    
    % aseg
    fileID = fopen(asegFile);
    tbl = textscan(fileID, '%d %d %d %d %s %f %f %f %f %f %f', 'CommentStyle' ,'#');
    fclose(fileID);
    structName_vol_aseg = [tbl{5} num2cell(tbl{4})];
    
    % lh.aparc
    aparcFile = [imgDir sprintf('%04d', rid) '/stats/lh.aparc.stats'];
    fileID = fopen(aparcFile);
    tbl = textscan(fileID, '%s %d %d %d %f %f %f %f %d %f', 'CommentStyle' ,'#');
    fclose(fileID);
    structName_gmVol_lhAparc = [cellfun(@(x) ['ctx-lh-' x], tbl{1}, 'UniformOutput', false) num2cell(tbl{4})];
    
    % rh.aparc
    aparcFile = [imgDir sprintf('%04d', rid) '/stats/rh.aparc.stats'];
    fileID = fopen(aparcFile);
    tbl = textscan(fileID, '%s %d %d %d %f %f %f %f %d %f', 'CommentStyle' ,'#');
    fclose(fileID);
    structName_gmVol_rhAparc = [cellfun(@(x) ['ctx-rh-' x], tbl{1}, 'UniformOutput', false) num2cell(tbl{4})];
    
    % Concat all together
    structName_vol = [structName_vol_aseg; structName_gmVol_lhAparc; structName_gmVol_rhAparc];
    
    %% Fetch GM volume by structure name
    
    volSum = 0;
    for idx2 = 1:numel(structNames)
        structName = structNames{idx2};
        ind = strcmp(structName_vol(:, 1), structName);
        if sum(ind) == 1
            volSum = volSum+structName_vol{ind, 2};
        elseif sum(ind) == 0
            warning(sprintf('No stats found for %s', structName));
        else
            error(sprintf('More than one stats found for %s', structName));
        end
    end
    
    rid_icv_vol(idx, 1) = rid;
    rid_icv_vol(idx, 2) = icv;
    rid_icv_vol(idx, 3) = volSum;
end
