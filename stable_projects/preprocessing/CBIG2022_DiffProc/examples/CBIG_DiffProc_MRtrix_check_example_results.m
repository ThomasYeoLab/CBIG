function [test_array] = CBIG_DiffProc_MRtrix_check_example_results(subj_output_dir, test_case)
% This function checks whether the example TBSS results are consistent for the
% unit test.
% Written by Leon Ooi and CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md

CBIG_TEST_DIR = getenv('CBIG_TESTDATA_DIR');
ref_dir = fullfile(CBIG_TEST_DIR,'stable_projects','preprocessing', ...
    'CBIG2022_DiffProc','MRtrix', 'ref_output', test_case);
% compare results
test_array = zeros(1,4);
% parcellations for comparison
parcellations = {'Schaefer2018_100Parcels_17Networks' 'Schaefer2018_400Parcels_17Networks'};
metrics = {'SIFT2' 'unfiltered'};
for i = 1:length(parcellations)
    for j = 1:length(metrics)
        test_connectome_dir = fullfile(subj_output_dir,strcat('connectome_', ...
            parcellations{i}, '_', metrics{j}, '.csv'));
        test_connectome = csvread(test_connectome_dir);
        ref_connectome_dir = fullfile(ref_dir, ...
            strcat('connectome_', parcellations{i}, '_', metrics{j}, '.csv'));
        ref_connectome = csvread(ref_connectome_dir);
        
        % check values for each subject
        if (corr2(test_connectome, ref_connectome) > 0.99)
                test_array((i-1)*length(metrics) + j) = 1;
            else
                fprintf('Correlation between test and reference (%s) %s is %f.\n', ...
                    parcellations{i}, metrics{j}, corr2(test_connectome, ref_connectome));
        end
    end
end
        
end