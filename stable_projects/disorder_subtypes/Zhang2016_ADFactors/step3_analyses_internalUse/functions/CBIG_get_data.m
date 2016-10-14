function [rid_diagnosis, rid_ICVoverGM, rid_prob, rid_age, rid_sex, rid_edu] = CBIG_get_data(subInfoFile, GMICVFile, K)

% [rid_diagnosis, rid_ICVoverGM, rid_prob, rid_age, rid_sex, rid_edu] = CBIG_get_data(subInfoFile, GMICVFile, K)
% Subject information
%
% Written by CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md

subInfo = csvread(subInfoFile);

% Diagnosis
rid_diagnosis = subInfo(:, [1 8]);
% Age
rid_age = subInfo(:, [1 7]);
% Gender
rid_sex = subInfo(:, [1 2]);
% Education
edu = CBIG_get_edu(subInfo(:, 1));
rid_edu = cell2mat(edu(:, [1 4]));

% ICV/GM_vol
rid_gm_icv = csvread(GMICVFile);
rid_ICVoverGM = [rid_gm_icv(:, 1) rid_gm_icv(:, 3)./rid_gm_icv(:, 2)];

r = dir(['../../lda/postprocess/188blAD/k' num2str(K) '/r*']);
% Probabilities of subtypes
% AD
prob = CBIG_get_prob(['../../lda/postprocess/188blAD/k' num2str(K) '/' r.name '/final.gamma']);
rid_prob_ad = [subInfo(subInfo(:, 8)==3, 1) prob];
% MCI
prob = CBIG_get_prob(['../../lda/outputs/k' num2str(K) r.name '_inf394MCI-gamma.dat']);
rid_prob_mci = [subInfo(subInfo(:, 8)==2, 1) prob];
% CN
prob = CBIG_get_prob(['../../lda/outputs/k' num2str(K) r.name '_inf228CN-gamma.dat']);
rid_prob_cn = [subInfo(subInfo(:, 8)==1, 1) prob];
% Consolidate and sort
rid_prob = [rid_prob_ad; rid_prob_mci; rid_prob_cn];
rid_prob = sortrows(rid_prob, 1);
