function CBIG_MMM_KRR_classical_summary(CBIG_CODE_DIR, data_dir, phe_list, n_rng, prefix)

% CBIG_MMM_KRR_classical_summary(CBIG_CODE_DIR)
% 
% This function summarize the kernel ridge regression classical results
%
% Inputs:
%   - CBIG_CODE_DIR
%     Full path of the ThomasYeoLab/CBIG repository in your local place
%     (https://github.com/ThomasYeoLab/CBIG).
%   - data_dir
%     Full path of the input/output directory for the KRR classical.
%   - phe_list
%     Full path of phenotypes (non-brain-imaging phenotypes) list in dataset
%     that performed the KRR classical.
%   - n_rng
%     Number (integer) of random number generator repeats of KRR classical. It can be number or
%     string.
%   - prefix
%     str of the prefix for the dataset.
%
% Written by Pansheng Chen and CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md

% initialization
if ischar(n_rng)
    n_rng = str2num(n_rng);
end
ks = [10, 20, 50, 100, 200];

fileID = fopen(phe_list);
temp = textscan(fileID, '%s');
entry_name = temp{1};

% output file
warning('off', 'MATLAB:MKDIR:DirectoryExists');
result_dir = fullfile(data_dir, 'final_result');
mkdir(result_dir)
result_file = fullfile(result_dir, 'krr_classical_res.mat');

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
