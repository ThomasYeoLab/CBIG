
function CBIG_preproc_FCmatrices_UnitTestCmp(preproc_out_dir, test_file_stem)

% CBIG_preproc_FCmatrices_UnitTestCmp(preproc_out_dir)
%
% Description of the function:
% This function is primarily used after running the
% CBIG_preproc_unit_tests_call_fMRI_preproc.csh (to generate the FC
% matrices of a single subject) to compare the FC matrices against the
% ground truth for the purpose of a Unit Test.
%
% Input:
%
% - preproc_out_dir:
%           a string. It contains the path of the directory that
%           stores the preprocessed data of the single
%           subject required for the Unit Test. This preproc_out_dir is the
%           same string passed to the
%           CBIG_preproc_unit_tests_call_fMRI_preproc.csh as an argument.
%
% - test_file_stem (optional):
%           a string appended to the generated FC_matrices filename.It is
%           used to specify which preprocessing pipelines were used
%           to generate the FC matrices. For example, a typical FC matrices
%           filename could be
%           "Sub1116_Ses1_rest_skip4_stc_mc_residc_interp_FDRMS0.2_
%           DVARS50_bp_0.009_0.08_fs6_sm6_all2all.mat". In this case, the
%           test_file_stem would be "rest_skip4_stc_mc_residc_interp_
%           FDRMS0.2_DVARS50_bp_0.009_0.08_fs6_sm6". This is an optional
%           input argument. If the user does not pass in this argument, the
%           test_file_stem will be set to the default "rest_skip4_stc_mc_
%           residc_interp_FDRMS0.2_DVARS50_bp_0.009_0.08_fs6_sm6".
%
% Output:
%
% - inequal_corr_log.txt:
%           A text file is generated and stored in the preproc_out_dir.
%           The text file states the FC matrices (eg. all2all) that
%           do not tally with the ground truth data. For each FC matrix,
%           the maximum difference will also be displayed in the text file.
%           If there are no differences at all between the generated
%           FC matrices and the ground truth, the maximum difference
%           displayed will be 0.
%
%
% - Displayed message:
%           After running the MATLAB function, a message will be displayed
%           to inform the user if there are differences or that no
%           differences are identified.
%
% Example:
%
% CBIG_preproc_FCmatrices_UnitTestCmp(fullfile(getenv('HOME'), 'storage', ...
% 'Pre_Proc_Unit_Test')
%
% The function takes in a string which contains the directory where the
% preprocessed data for the single subject are stored after running the
% CBIG_preproc_unit_tests_call_fMRI_preproc.csh.
%
%
%
% Written by CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md


if(~exist('test_file_stem', 'var'))
    test_file_stem = ['rest_skip8_stc_mc_sdc_me_residc_interp_FDRMS0.3_DVARS60'...
        '_bp_0.009_0.08_fs6_sm6'];
end

% Beginning of the main code
% Storing the ID of the 4 subjects chosen for the Unit Test and
% defining the 7 types of correlations within the brain.
subject_array = 'sub005';
corr_label_array = {'all2all', 'lh2lh', 'lh2rh', 'lh2subcort', 'rh2rh',...
    'rh2subcort', 'subcort2subcort'};

% num_ineq will be used to count the number of times an FC matrices is
% compared and showed different results from the ground truth FC matrices
num_ineq = 0;

% file handle to write in the specific subject and the correlation type
% which are different from the ground truth
fid = fopen(fullfile(preproc_out_dir, 'inequal_corr_log.txt'),'wt');

% true_path is the path to the FC matrices directory of the ground truth
% test_path is the path to the FC matrices directory of the test data
true_path = fullfile(getenv('CBIG_TESTDATA_DIR'), 'stable_projects', ...
    'preprocessing', 'CBIG_fMRI_Preproc2016', 'single_subject','data', ...
    subject_array, 'FC_metrics', 'Pearson_r');

test_path = fullfile(preproc_out_dir, subject_array, 'FC_metrics', 'Pearson_r');

% Loop through each FC matrix and compare
for j = 1: length(corr_label_array)
    file_path_true = [subject_array '_rest_skip8_stc_mc_sdc_me_residc_'...
        'interp_FDRMS0.3_DVARS60_bp_0.009_0.08_fs6_sm6_'...
        corr_label_array{1,j} '_Yan_400Parcels_Kong17Networks_order_with_19Subcortical' '.mat'];
    file_path_test = [subject_array '_' test_file_stem '_' ...
        corr_label_array{1,j} '_Yan_400Parcels_Kong17Networks_order_with_19Subcortical' '.mat' ];
    true_full_path = fullfile(true_path, file_path_true);
    test_full_path = fullfile(test_path, file_path_test);
    
    % Loading the FC matrices of the ground truth and the test data
    true_corr_mat_struct = load(true_full_path);
    test_corr_mat_struct = load(test_full_path);
    true_corr_mat = true_corr_mat_struct.corr_mat;
    test_corr_mat = test_corr_mat_struct.corr_mat;
    
    % check each element's difference between the 2 matrices
    test_equality = abs(true_corr_mat - test_corr_mat);
    
    % if the max difference between the 2 matrices is more than 1e-6, then
    % increment num_ineq by 1 and print the corresponding subject and the
    % specific correlation type into the textfile
    if max(max(test_equality)) > 1e-6
        num_ineq = num_ineq + 1;
        fprintf(fid, [corr_label_array{1,j} ': '  ...
            num2str(max(max(test_equality)))  '\n']);
    else
        fprintf(fid, [corr_label_array{1,j} ': '  ...
            num2str(0)  '\n']);
    end
end


% Display on MATLAB whether there are any differences between the ground
% truth matrices and the test data
if num_ineq ~= 0
    display(['Differences noted in one or more correlation matrices.' ...
        'Check the inequal_corr_log.txt file in your '...
        'preproc_out_dir for more details'])
    fprintf(fid, [ 'Num_differences: ' num2str(num_ineq)  '\n']);
else
    display('All correlation matrices are the same.')
    fprintf(fid, ['Num_differences: ' num2str(0)  '\n']);
end

fclose(fid);

end
