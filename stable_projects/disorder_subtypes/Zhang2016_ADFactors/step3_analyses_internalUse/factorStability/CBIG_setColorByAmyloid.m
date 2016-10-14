function rid_colors = CBIG_setColorByAmyloid(ridList)

% rid_colors = CBIG_setColorByAmyloid(ridList)
% Red a+; green a-; blue unknown

%
% Written by CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md

rid_viscode_dummyCol_amyloid = CBIG_get_amyloid(ridList);
rid_amyloid = cell2mat(rid_viscode_dummyCol_amyloid(:, [1 4]));

noSubjets = numel(ridList);

rid_colors = cell(noSubjets, 2);
for idx = 1:noSubjets
    rid = ridList(idx);
    rid_colors{idx, 1} = rid;
    ind = rid_amyloid(:, 1)==rid;
    if sum(ind) == 0 % no amyloid status
        rid_colors{idx, 2} = 'b';
    elseif sum(ind) == 1
        a = rid_amyloid(ind, 2);
        if a < 192 % a+
            rid_colors{idx, 2} = 'r';
        else % a-
            rid_colors{idx, 2} = 'g';
        end
    else
        error('More than two amyloid values found!');
    end
end