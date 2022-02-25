function CBIG_MM_KRR_classical_summary(CBIG_CODE_DIR, prefix)

% CBIG_MM_KRR_classical_summary(CBIG_CODE_DIR)
% 
% This function summarize the kernel ridge regression classical results
%
% Inputs:
%   - CBIG_CODE_DIR
%     Full path of the ThomasYeoLab/CBIG repository in your local place
%     (https://github.com/ThomasYeoLab/CBIG).
%
% Written by Tong He and CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md

% initialization
if ~exist('prefix','var')
    split = 'test';
    prefix = 'ukbb';
else
    split = 'diff_roi';
end

base_dir = fullfile(CBIG_CODE_DIR, ...
    'stable_projects/predict_phenotypes/He2022_MM');
input_dir = [getenv('CBIG_REPDATA_DIR'), '/stable_projects/predict_phenotypes/He2022_MM/'];

n_rng = 100;
ks = [10, 20, 50, 100, 200];
data_dir = fullfile(base_dir, 'replication', ['output_KRR_classical_' prefix]);
phe_list = fullfile(input_dir, [prefix '_' split '_final_phe_list.txt']);
fileID = fopen(phe_list);
temp = textscan(fileID, '%s');
entry_name = temp{1};

% output file
warning('off', 'MATLAB:MKDIR:DirectoryExists');
result_dir = fullfile(data_dir, 'final_result');
mkdir(result_dir)
result_file = fullfile(result_dir, ['krr_classical_res_' split '.mat']);

% load KRR result
meta_cor = zeros(length(entry_name), n_rng, length(ks));
meta_cod = zeros(length(entry_name), n_rng, length(ks));
for i_entry = 1:length(entry_name)
    entry = entry_name{i_entry};
    for i = 1:length(ks)
        k = ks(i);
        for rng_num = 1:n_rng
            result_mat = fullfile(data_dir, ['output_phe_' entry], ...
                [prefix, '_', entry, '_k_', num2str(k), '_rng_num_', ...
                num2str(rng_num)], 'final_result.mat');

            if ~exist(result_mat, 'file')
                error([result_mat, ' not exist!'])
            end

            tmp = load(result_mat);

            if isnan(tmp.optimal_acc)
                disp([result_mat ' corr is nan'])
            end

            if isnan(tmp.optimal_stats.predictive_COD)
                disp([result_mat ' cod is nan'])
            end

            meta_cod(i_entry, rng_num, i) = tmp.optimal_stats.predictive_COD;
            meta_cor(i_entry, rng_num, i) = tmp.optimal_acc;
        end
    end
end
phe_tes = entry_name';
save(result_file, 'meta_cor', 'meta_cod', 'phe_tes')
