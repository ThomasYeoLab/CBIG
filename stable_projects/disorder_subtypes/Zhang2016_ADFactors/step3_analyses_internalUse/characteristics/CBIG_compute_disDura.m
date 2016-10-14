function rid_durations = CBIG_compute_disDura(subInfoFile, GMICVFile, K)

% rid_durations = CBIG_compute_disDura(subInfoFile, GMICVFile, K)
%%% Fetch onset year
%
% Written by CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md

fileID = fopen('//share/users/imganalysis/yeolab/data/ADNI_mert/scripts/PTDEMOG.csv', 'r');
rawData = textscan(fileID, '%*q%*q%q%*q%*q%*q%*q%*q%*q%*q%*q%*q%*q%*q%*q%*q%*q%*q%*q%*q%*q%*q%*q%*q%*q%*q%*q%q%*[^\n]',...
    'delimiter', {',', '\n'}, 'HeaderLines', 1);
fclose(fileID);
% RID string to integer
rid = rawData{1};
yyyyOnset = rawData{2};
for idx = 1:1:size(rid, 1)
    rid{idx} = str2double(rid{idx});
    yyyyOnset{idx} = str2double(yyyyOnset{idx});
end
rid_yyyyOnset = cell2mat([rid yyyyOnset]);


%%% Fetch baseline year
fileID = fopen('//share/users/imganalysis/yeolab/data/ADNI_mert/scripts/DXSUM_PDXCONV_ADNIALL.csv', 'r');
rawData = textscan(fileID, '%*q%*q%q%*q%q%*q%*q%*q%q%*[^\n]',...
    'delimiter', ',', 'HeaderLines', 1);
fclose(fileID);
% RID string to integer
rid = rawData{1};
for idx = 1:1:size(rid, 1)
    rid{idx} = str2double(rid{idx});
end
% VISCODE
viscode = rawData{2};
% EXAMDATE
examdate = rawData{3};
rid_viscode_examdate = [rid viscode examdate];
% 'bl' only
indBL = find(ismember(viscode, 'bl'));
rid_viscode_examdate = rid_viscode_examdate(indBL, :);
% Sort by RID
rid_viscode_examdate = sortrows(rid_viscode_examdate, 1);
% Only keep years
rid_yyyyExam = zeros(size(rid_viscode_examdate, 1), 2);
for idx = 1:size(rid_yyyyExam, 1)
    rid_yyyyExam(idx, 1) = rid_viscode_examdate{idx, 1};
    yyyymmdd = rid_viscode_examdate{idx, 3};
    rid_yyyyExam(idx, 2) = str2double(yyyymmdd(1:4));
end

%%% bl-AD only
[rid_diagnosis, ~, rid_prob, ~, ~, ~] = get_data(subInfoFile, GMICVFile, K);
rid_prob_AD = rid_prob(rid_diagnosis(:, 2)==3, :);
ridAD = rid_diagnosis(rid_diagnosis(:, 2)==3, 1);
indAD = ismember(rid_yyyyOnset(:, 1), ridAD);
rid_yyyyOnset_AD = rid_yyyyOnset(indAD, :);
indAD = ismember(rid_yyyyExam(:, 1), ridAD);
rid_yyyyExam_AD = rid_yyyyExam(indAD, :);
% Make unique and sort
[~, ia, ~] = unique(rid_yyyyOnset_AD(:, 1));
rid_yyyyOnset_AD = rid_yyyyOnset_AD(ia, :);
rid_yyyyExam_AD = sortrows(rid_yyyyExam_AD, 1);

%%% Compute disease duration
if isequal(rid_yyyyOnset_AD(:, 1), rid_yyyyExam_AD(:, 1))
    indMiss = rid_yyyyOnset_AD(:, 2)==-4;
    durations = rid_yyyyExam_AD(~indMiss, 2)-rid_yyyyOnset_AD(~indMiss, 2);
    rid_durations = [rid_yyyyOnset_AD(~indMiss, 1) durations];
else
    error('RID mismatched!');
end
