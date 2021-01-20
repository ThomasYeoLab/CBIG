function [] = CBIG_gwMRF_generate_example_results(output_folder)

% This function generates the final results for the example provided in
% Alex's gwMRF stable project. Note that There will be some warnings from 
% the optimizer: "Warning: Neighbours array should be upper-triangular; 
% entries below the diagnonal will be ignored." You can ignore them.
% It should take less than 2 hours with 1 CPU for the example code to finish running.
%
% Input:
%   - output_folder: specify the path of the folder in which the example
%  results will be saved.
%
% Written by Xiaoxuan Yan and CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md

% generate example input file
cmd = ['sh ${CBIG_CODE_DIR}/stable_projects/brain_parcellation/Schaefer2018_LocalGlobal'...
    '/examples/example_input/CBIG_gwMRF_create_example_input_fullpaths.sh  ' output_folder];
system(cmd);

example_input_fullpaths = fullfile(output_folder,'example_input_fullpaths.csv');

code_dir = fullfile(getenv('CBIG_CODE_DIR'),'stable_projects','brain_parcellation',...
'Schaefer2018_LocalGlobal','Code');

% we shall run for 1 seed, 2 iterations, using the data of 2 subjects
cd(code_dir);
CBIG_gwMRF_build_data_and_perform_clustering(example_input_fullpaths,output_folder,1,2,50,50,5000,2,1,50000000,15);

end
