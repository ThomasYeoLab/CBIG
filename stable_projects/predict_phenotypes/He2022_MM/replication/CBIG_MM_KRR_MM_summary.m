function CBIG_MM_KRR_MM_summary(CBIG_CODE_DIR)

% CBIG_MM_KRR_MM_summary(CBIG_CODE_DIR)
% 
% This function runs the meta-matching with kernel ridge regression
% algorithm to all testing meta-set non-brain-imaging phenotypes with
% all K shot and random numbers.
%
% Inputs:
%   - CBIG_CODE_DIR
%     Full path of the ThomasYeoLab/CBIG repository in your local place
%     (https://github.com/ThomasYeoLab/CBIG).
% 
% Written by Tong He and CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md

%% add kernel regression code to path
addpath(genpath('../KRR_MM/'));

split = 'test';
base_dir = fullfile(CBIG_CODE_DIR, 'stable_projects', 'predict_phenotypes', 'He2022_MM');
input_dir = fullfile(getenv('CBIG_REPDATA_DIR'), 'stable_projects', 'predict_phenotypes', 'He2022_MM');
data_dir = fullfile(base_dir, 'replication', ['output_KRR_mm']);

ks = [10, 20, 50, 100, 200];
time = 100;
phe_tes = load_phe_from_txt(fullfile(input_dir, ['ukbb_' split '_final_phe_list.txt']));
phe_tes_mat = strrep(phe_tes, '.', '_');
phe_tes_mat = strrep(phe_tes_mat, '-', '_');
phe_tes_mat = strcat('x', phe_tes_mat);
if any(strcmp(phe_tes_mat, 'xage_2_0'))
    phe_tes_mat{strcmp(phe_tes_mat, 'xage_2_0')} = 'age_2_0';
end
phe_tra = load_phe_from_txt(fullfile(input_dir, 'ukbb_train_final_phe_list.txt'));
test_csv = fullfile(input_dir, ['ukbb_' split '_final.csv']);
test_csv = readtable(test_csv);
result_file = fullfile(data_dir, 'meta_result', ['meta_res_' split '.mat']);

meta_cor = zeros(length(phe_tes), time, length(ks));
meta_cod = zeros(length(phe_tes), time, length(ks));
meta_cod_k = zeros(length(phe_tes), time, length(ks));

y_pred = load_krr_base_res(data_dir, phe_tra);

for j = 1:length(ks)
    k = ks(j);
    for i = 1:time
        [opt_entry, opt_k, opt_test_corr, opt_test_cod] = ...
            CBIG_MM_KRR_MM(base_dir, data_dir, y_pred, test_csv,...
            phe_tra, phe_tes, phe_tes_mat, i, k, 'cod');
        meta_cor(:, i, j) = opt_test_corr;
        meta_cod(:, i, j) = opt_test_cod;
        meta_cod_k(:, i, j) = opt_k;
    end
end
save(result_file, 'meta_cor', 'meta_cod', 'meta_cod_k', 'phe_tes')
rmpath(genpath('../KRR_MM/'));
end


function phe_list = load_phe_from_txt(txt)
% phe_list = load_phe_from_txt(txt)
% 
% This function load list of phenotypes from txt file
%
% Inputs:
%   - txt
%     path of phenotypes list txt file, It should be a txt file that contains
%     #phenotypes of line, while each line is name of one phenotype.
%
% Outputs:
%	- phe_list = phenotypes list

    fileID = fopen(txt);
    temp = textscan(fileID, '%s');
    phe_list = temp{1}';
    fclose(fileID);
end

function y_pred = load_krr_base_res(data_dir, phe_list)
% y_pred = load_krr_base_res(data_dir, phe_list)
% 
% This function load base kernel ridge regression result from training
% meta-set
%
% Inputs:
%   - data_dir
%     name of result data directory
%   - phe_list
%     list of training meta-set phenotpyes
%
% Outputs:
%	- y_pred = predicted y from kernel ridge regression

    mat = load(fullfile(data_dir, 'ukbb_rng_num_1', 'final_result.mat'));
    y_pred = mat.y_predict_concat;
    tmp_ind = find(strcmp(phe_list, '31-0.0'));
    if logical(length(tmp_ind))
        mat = load(fullfile([data_dir '_sex'], 'ukbb_rng_num_1', 'final_result.mat'));
        y_pred_sex = mat.y_predict_concat;
        y_pred(:, tmp_ind + 1:end + 1) = y_pred(:, tmp_ind:end);
        y_pred(:, tmp_ind) = y_pred_sex;
    end
end