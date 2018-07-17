
function CBIG_preproc_FCmetrics_UnitTestComparison(testdir)

% CBIG_preproc_FCmetrics_UnitTestComparison(testdir)
%
% Description of the function:
% This function is primarily used after running the
% CBIG_preproc_unit_tests_FCmetrics_fMRI_preproc.csh (to generate the FC
% metrices of 4 subjects) to compare the FC metrices against the ground
% truth for the purpose of a Unit Test.
%
% Input:
%
% - testdir:
%           testdir is a string which contains the path of the directory
%           that stores the preprocessed data of all subjects required for
%           the Unit Test. This testdir is the same string passed to the
%           CBIG_preproc_unit_tests_FCmetrics_fMRI_preproc.csh as an arg.
%
% Output:
%
% - inequal_corr_log.txt:
%           This text file is generated and stored in the testdir. The
%           text file will state all the subjects and their corresponding
%           FC metrices (eg. all2all) which do not tally with the
%           ground truth data. If there are no differences at all between
%           the generated FC metrices and the ground truth, the text file
%           will state "All correlation matrices are the same".
%
% - Displayed message:
%           After running the MATLAB function, a message will be displayed
%           to inform the user if there are differences or that no 
%	    differences are identified.
%
% Example:
%   
% CBIG_preproc_FCmetrics_UnitTestComparison(['/data/users/usr/storage/'...
% 'Pre_Proc_Unit_Test'])
% 
% The function takes in a string which contains the directory where the
% preprocessed data for all 4 subjects are stored after running the 
% CBIG_preproc_unit_tests_FCmetrics_fMRI_preproc.csh.
% 
% 
%
% Written by CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md


% Beginning of the main code
% Storing the ID of the 4 subjects chosen for the Unit Test and 
% defining the 7 types of correlations within the brain. 
subject_array = ['Sub0017_Ses1'; 'Sub0735_Ses1'; 'Sub1155_Ses1';...
                'Sub1488_Ses1'];
corr_label_array = {'all2all', 'lh2lh', 'lh2rh', 'lh2subcort', 'rh2rh',...
                    'rh2subcort', 'subcort2subcort'};

% num_ineq will be used to count the number of times an FC metric is 
% compared and showed different results from the ground truth FC metric
num_ineq = 0;

% file handle to write in the specific subject and the correlation type
% which are different from the ground truth
fid = fopen([testdir '/inequal_corr_log.txt'],'wt');

% Loop through the 4 subjects
for i = 1:size(subject_array,1)

    % true_path is the path to the FC metrics directory of the ground truth
    % test_path is the path to the FC metrics directory of the test data
	true_path = ['/mnt/eql/yeo1/CBIG_private_unit_tests_data/' ...
                 'stable_projects/preprocessing/CBIG_fMRI_Preproc2016/' ...
                 'FCmetrics/' subject_array(i,:) '/FC_metrics/Pearson_r/'];
	test_path = [testdir '/' subject_array(i,:) '/FC_metrics/Pearson_r/'];

    % Count the number of FC metrices files in the true_path and test_path
    % directory
	true_num_files = size(dir([true_path '*.mat']));
	test_num_files = size(dir([test_path '*.mat']));
	
    % Loop through each FC metric and compare
	for j = 1: true_num_files
		file_path = [subject_array(i,:) '_rest_skip4_stc_mc_residc_'...
                    'interp_FDRMS0.2_DVARS50_bp_0.009_0.08_fs6_sm6_'...
                     corr_label_array{1,j} '.mat'];
		true_full_path = [true_path file_path];
		test_full_path = [test_path file_path];
        
        % Loading the FC metrics of the ground truth and the test data
		true_corr_mat_struct = load(true_full_path);
		test_corr_mat_struct = load(test_full_path);
		true_corr_mat = true_corr_mat_struct.corr_mat;
		test_corr_mat = test_corr_mat_struct.corr_mat;
        
        % if the ground truth and the test data are equal, test_equality
        % will be given a value of 1
		test_equality = isequal(true_corr_mat, test_corr_mat);
        
        % if ground truth and test data are not equal, then increment
        % num_ineq by 1 and print the corresponding subject and the
        % specific correlation type into the textfile
		if test_equality == 0
			num_ineq = num_ineq + 1;
			fprintf(fid, [subject_array(i,:) ' ' corr_label_array{1,j} ...
                    '\n']);
		end
	end

end


% Display on MATLAB whether there are any differences between the ground
% truth metrices and the test data after looping through all 4 subjects
if num_ineq ~= 0
	display(['Differences noted in one or more correlation matrices.' ...
            'Check the inequal_corr_log.txt file in your testdir for'...
            'more details'])
else
	display('All correlation matrices are the same.')
	fprintf(fid, 'No difference identified');
end

fclose(fid);

end

	
	
		
	



