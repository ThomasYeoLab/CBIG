function rid_time_score = CBIG_compute_time(q, q_full)

% rid_time_score = CBIG_compute_time(q, q_full)
%
% Written by CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md

ridList_q_full = cell2mat(q_full(:, 1));

ind_sc_bl = ismember(q_full(:, 2), 'sc');
if sum(ind_sc_bl) == 0 % bl convention
    ind_sc_bl = ismember(q_full(:, 2), 'bl');
end

rid_time_score = zeros(size(q, 1), 3);

for idx = 1:size(rid_time_score, 1)
    rid = q{idx, 1};
    date_now = q{idx, 3};
    score = q{idx, 4};
    % Find 'sc' date
    idx_sc_bl = (ridList_q_full==rid)&ind_sc_bl;
    if sum(idx_sc_bl) == 1
        date_sc = q_full{idx_sc_bl, 3};
    else % Some rare subjects don't have baseline
        disp(q(cell2mat(q(:, 1))==rid, :));
        date_sc = input('Please enter estimated baseline date (yyyy-mm-dd): ', 's');
    end
    % Compute time
    noDays = datenum(date_now, 'yyyy-mm-dd')-datenum(date_sc, 'yyyy-mm-dd');
    noYears = noDays/365;
    % Pack
    rid_time_score(idx, :) = [rid noYears score];
end
