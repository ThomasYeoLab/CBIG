function CBIG_pFIC_extrapolate_best_parameter(test_file_path, training_folder_path, training_file_name, ...
    output_path, myelin_full, rsfc_gradient_full)

% CBIG_pFIC_extrapolate_best_parameter(test_all_file_path, training_folder_path, output_path, ...
%   myelin_full, rsfc_gradient_full)
% 
% This function takes the linear coeffcients having the lowest validation
% cost and generates extrapolated parameters by using whole-cortex myelin
% and RSFC gradient.
% Input:
%   -test_file_path: absolute path of the .csv containing parameters and
%   associated test costs
%   -training_folder_path: absolute path of the training folder
%   -training_file_path: name of the training file (for example, 'training')
%   -output_path: asbolute path to save out the extrapolated parameters
%   -myelin_full: whole-cortex myelin
%   -rsfc_gradient_full: whole-cortex RSFC gradient
% Output:
%   -'extrapolated_parameter.csv' saved under <output_path>
% Written by Shaoshi Zhang and CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md


test_file = csvread(test_file_path);
ini_number = test_file(1, 1);
run_number = test_file(2, 1) + 1; %python counts from 0, matlab counts from 1

training_param_file = csvread([training_folder_path '/' training_file_name '_' num2str(ini_number) '.csv']);
training_param = training_param_file(1:10, run_number);
wee_coeff = training_param(1:3);
wei_coeff = training_param(4:6);
G = training_param(7);
sigma_coeff = training_param(8:10);

myelin = csvread(myelin_full);
rsfc_gradient = csvread(rsfc_gradient_full); 
num_roi = size(myelin, 1);
cmatrix = [ones(num_roi, 1) myelin rsfc_gradient];

wee = cmatrix * wee_coeff;
wei = cmatrix * wei_coeff;
sigma = cmatrix * sigma_coeff;

param = [wee; wei; G; sigma];
if ~exist(output_path, 'dir')
   mkdir(output_path)
end
dlmwrite([output_path '/extrapolated_parameter.csv'], param, 'delimiter', ',', 'precision', 15);

end
