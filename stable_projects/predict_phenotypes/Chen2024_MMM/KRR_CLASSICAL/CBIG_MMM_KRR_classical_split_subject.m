function CBIG_MMM_KRR_classical_split_subject(data_dir, subject_list_file, rng_num, k, y, prefix)

% CBIG_MMM_KRR_classical_split_subject(data_dir, subject_list_file, rng_num, k, y, prefix)
% 
% This function splits subjects into k subjects and remaining subjects.
%
% Inputs:
%   - data_dir
%     Path of the your output and intermediate values. You can
%     also change this to any place you want.
%   - subject_list_file
%     Path of subject list of setup_file dataset. It should be a txt file
%     that contains #subject of line, while each line is the subject id of
%   - rng_num
%     Number (integer) of random number generator
%   - k
%     Number (integer) of the K for K shot (participants) learning
%   - y
%     array of y, phenotype value
%   - prefix
%     dataset name prefix for the split file saving
% 
% Written by Pansheng Chen and CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md

rng(rng_num);

% read out subject file
fileID = fopen(subject_list_file);
temp = textscan(fileID, '%s');
subject_list = temp{1};
% Train/Validation/Test split
num_total = size(subject_list, 1);
num_train = k;
num_test = num_total - num_train;

% handle missing value
ind_real = ~isnan(y);
num_real = sum(ind_real);
if num_total ~= k
    % dont check if num_total is K for intentional usage
    if num_real < 1.5 * k
        error([num2str(num_real) ' is less than 1.5 * ' num2str(k)])
    end
end

subject_list_real = subject_list(ind_real);
subject_list_real = subject_list_real(randperm(num_real));
subject_list_nan = subject_list(isnan(y));

% generate list
subject_list_rand = cat(1, subject_list_real, subject_list_nan);
subject_list_train = sort(subject_list_rand(1:num_train));
subject_list_test = sort(subject_list_rand(end-num_test+1:end));

% generate index list
fold_index = zeros(num_total,1);

if ~isempty(subject_list_test)
    for i = 1:size(subject_list_test,1)
        fold_index = fold_index + (strcmp(subject_list_test(i), subject_list));
    end
end

% save them in strucuture
sub_fold = struct('fold_index', fold_index);

% save out
save(fullfile(data_dir, [prefix '_subject_split.mat']), 'sub_fold');

end
