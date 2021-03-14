function CBIG_pMFM_step2_generate_FCD_schaefer()

% This function is the wrapper to generate parcellated time serises and FCD
% for Desikan parcellation
%
% There is no input for this function as it can automatically get the
% output file from previous step.
% There is no output for this function as it will generate the output files
%
% Written by Kong Xiaolu and CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md

FCD_result_dir = '../output/FCD';
if ~exist(FCD_result_dir,'dir')
    mkdir(FCD_result_dir)
end

generate_training_FCD()
generate_validation_FCD()
generate_test_FCD()

end


function generate_training_FCD()

% This function is the wrapper to generate parcellated FCD matrices
% for Desikan parcellation for training set
%
% There is no input for this function as it can automatically get the
% output file from previous step.
% There is no output for this function as it will generate the output files
%
% Written by Kong Xiaolu and CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md

TC_dir = '../output/TC/train';
file_dir = dir(TC_dir);

FCD_train_dir = '../output/FCD/train';
if ~exist(FCD_train_dir,'dir')
    mkdir(FCD_train_dir)
end

for i = 3:size(file_dir,1)
    load([TC_dir file_dir(i).name], 'TC');
    TC_100 = TC;
    FCD_mat = FCD_plot(TC_100, 83);
    save([FCD_train_dir file_dir(i).name], 'FCD_mat')
end

end

function generate_validation_FCD()

% This function is the wrapper to generate parcellated FCD matrices
% for Desikan parcellation for validation set
%
% There is no input for this function as it can automatically get the
% output file from previous step.
% There is no output for this function as it will generate the output files
%
% Written by Kong Xiaolu and CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md

TC_dir = '../output/TC/validation';
file_dir = dir(TC_dir);

FCD_vali_dir = '../output/FCD/validation';
if ~exist(FCD_vali_dir,'dir')
    mkdir(FCD_vali_dir)
end

for i = 3:size(file_dir,1)
    load([TC_dir file_dir(i).name], 'TC');
    TC_100 = TC;
    FCD_mat = FCD_plot(TC_100, 83);
    save([FCD_vali_dir file_dir(i).name], 'FCD_mat')
end

end

function generate_test_FCD()

% This function is the wrapper to generate parcellated FCD matrices
% for Desikan parcellation for training set
%
% There is no input for this function as it can automatically get the
% output file from previous step.
% There is no output for this function as it will generate the output files
%
% Written by Kong Xiaolu and CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md

TC_dir = '../output/TC/test';
file_dir = dir(TC_dir);

FCD_test_dir = '../output/FCD/test';
if ~exist(FCD_test_dir,'dir')
    mkdir(FCD_test_dir)
end

for i = 3:size(file_dir,1)
    load([TC_dir file_dir(i).name], 'TC');
    TC_100 = TC;
    FCD_mat = FCD_plot(TC_100, 83);
    save([FCD_test_dir file_dir(i).name], 'FCD_mat')
end

end


function FCD_mat = FCD_plot(data_TC, window)

% This function is used to generate the parcellated FCD
% Input:
%   data_TC: parcellated TC
%   window:  window length of the sliding windows
% Output:
%   FCD_mat: parcellated FCD matrix

[dimension,timeline] = size(data_TC);
mask_tril = ~tril(ones(dimension,dimension));

for i = 1:timeline-window+1
    corr_swc = corr(data_TC(:,i:i+window-1)');
    corr_vect = corr_swc(mask_tril);
    corr_mat(:,i) = corr_vect;
end

FCD_mat = corr(corr_mat);

end




