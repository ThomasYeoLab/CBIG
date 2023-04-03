classdef CBIG_hMRF_unit_test < matlab.unittest.TestCase
% Written by Xiaoxuan Yan and CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md
methods (Test)
    function test_example(testCase)
        % create output folder
        CBIG_CODE_DIR = getenv('CBIG_CODE_DIR');
        load(fullfile(CBIG_CODE_DIR, 'unit_tests', 'replace_unittest_flag'));
        cur_dir = fullfile(CBIG_CODE_DIR, 'stable_projects', 'brain_parcellation', 'Yan2023_homotopic', 'unit_tests');
        out_dir = fullfile(cur_dir, 'output');
        
        % run the example
        addpath(fullfile(CBIG_CODE_DIR, 'stable_projects', ...
            'brain_parcellation', 'Yan2023_homotopic', 'examples'));
        CBIG_hMRF_example_wrapper(out_dir);
        
        % check example results
        CBIG_hMRF_check_example_results(out_dir);
        
        % replace reference outputs if flag is 1
        ref_dir = fullfile(getenv('CBIG_CODE_DIR'), 'stable_projects', 'brain_parcellation',...
            'Yan2023_homotopic', 'examples', 'ref_results');
        
        % left hemi time matrix
        lh_time_mat = fullfile(out_dir, 'time_data', 'lh_time_matrix.mat');
        lh_time_mat_ref = fullfile(ref_dir, 'time_data', 'lh_time_matrix.mat');
        if(replace_unittest_flag)
            disp('Replacing unit test reference results for CBIG_hMRF_unit_test...');
            copyfile(lh_time_mat, lh_time_mat_ref);
        end
        
        % right hemi time matrix
        rh_time_mat = fullfile(out_dir, 'time_data', 'rh_time_matrix.mat');
        rh_time_mat_ref = fullfile(ref_dir, 'time_data', 'rh_time_matrix.mat');
        if(replace_unittest_flag)
            disp('Replacing unit test reference results for CBIG_hMRF_unit_test...');
            copyfile(rh_time_mat, rh_time_mat_ref);
        end

        % parcellation
        parcellation = fullfile(out_dir, 'parcellation', 'results',...
         '100parcels_C1.0e+02_K15_Wxyz1.5e+03_D10_A1_iterations_3_seed_835.mat');
        parcellation_ref = fullfile(ref_dir, 'parcellation_seed_835',...
         '100parcels_C1.0e+02_K15_Wxyz1.5e+03_D10_A1_iterations_3_seed_835.mat');
        if(replace_unittest_flag)
            disp('Replacing unit test reference results for CBIG_hMRF_unit_test...');
            copyfile(parcellation, parcellation_ref);
        end

        % remove the output directory
        rmdir(out_dir, 's')

        % remove path
        rmpath(fullfile(CBIG_CODE_DIR, 'stable_projects', ...
            'brain_parcellation', 'Yan2023_homotopic', 'examples'));
    end
end
end