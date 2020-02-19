classdef CBIG_preproc_single_subject_unit_test < matlab.unittest.TestCase
%
% Target project:
%                 CBIG_fMRI_Preproc2016
%
% Case design:
%                 CBIG_fMRI_Preproc2016 already have unit tests, here we
%                 call its "single_subject" unit test and automatically
%                 make judgement on whether the unit test is passed or
%                 failed
%
%                 For now, the stable projects' unit tests in our repo
%                 require manual check of the output txt files or images,
%                 making it incovenient for wrapper function to call them
%                 and automatically draw conclusions. As a result, we write
%                 some simple matlab test functions for these stable
%                 projects' unit tests.
%
% Written by Yang Qing and CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md

    methods (Test)
        function test_single_subject_Case(testCase)
            %% path setting
            addpath(fullfile(getenv('CBIG_CODE_DIR'), 'stable_projects', ...
                'preprocessing', 'CBIG_fMRI_Preproc2016', 'utilities')); 
            addpath(fullfile(getenv('CBIG_CODE_DIR'), 'stable_projects', ...
                'preprocessing', 'CBIG_fMRI_Preproc2016', 'unit_tests', 'single_subject'))
            UnitTestDir = fullfile(getenv('CBIG_CODE_DIR'), 'stable_projects',...
                'preprocessing', 'CBIG_fMRI_Preproc2016', 'unit_tests');
            OutputDir = fullfile(UnitTestDir, 'output', 'single_subject_Case'); 
            
            %create output dir (IMPORTANT)
            if(exist(OutputDir, 'dir'))
                rmdir(OutputDir, 's')
            end
            mkdir(OutputDir);
            
            
            %% call CBIG_fMRI_Preproc2016 single_subject unit test script to generate results
            cmd = [fullfile(UnitTestDir, 'single_subject', ...
                'CBIG_preproc_unit_tests_call_fMRI_preproc.csh'), ' ', OutputDir];
            system(cmd); % this will submit a job to HPC
            
            
            %% we need to periodically check whether the job has finished or not             
            cmdout = 1;
            while(cmdout ~= 0)
                cmd = 'qstat | grep preproc | grep `whoami` | wc -l';
                [~, cmdout] = system(cmd);
                % after job finishes, cmdout should be 0
                cmdout = str2num(cmdout(1: end-1));
                pause(60); % sleep for 1min and check again
            end

            
            %% check surf files
            pipe_dir1 = fullfile(getenv('CBIG_TESTDATA_DIR'), 'stable_projects', 'preprocessing', ...
                'CBIG_fMRI_Preproc2016', 'single_subject', 'data');
            pipe_name1 = 'gt';
            pipe_stem1 = '_rest_skip4_stc_mc_residc_interp_FDRMS0.2_DVARS50_bp_0.009_0.08_fs6_sm6_fs5';

            pipe_dir2 = OutputDir;
            pipe_name2 = 'user-test';
            pipe_stem2 = '_rest_skip4_stc_mc_residc_interp_FDRMS0.2_DVARS50_bp_0.009_0.08_fs6_sm6_fs5';

            subject_id = 'Sub1116_Ses1';
            runs = {'002', '003'};

            output_dir = fullfile(OutputDir, 'compare_output');
            
            for i = 1: length(runs)                
                CBIG_preproc_compare_two_pipelines(pipe_dir1, pipe_name1, ...
                    pipe_stem1, pipe_dir2, pipe_name2, pipe_stem2, subject_id, ...
                    runs{i}, output_dir, 'surf');
            end
            
            % check surf run 002
            corr_file = fullfile(OutputDir, 'compare_output', 'Sub1116_Ses1',...
                '002', 'gt_user-test_corr_surf_stat.txt');
            corr_result = importdata(corr_file);
            corr_result = regexp(corr_result{3}, ':', 'split'); % we look at the min corr
            corr_result = corr_result(2);
            corr_result = str2num(corr_result{1});
            assert(corr_result > 0.99, ...
                sprintf('surface min_corr of run002 is less than 0.99! The value is %f \n', ...
                corr_result))
            
            % check surf run 003
            corr_file = fullfile(OutputDir, 'compare_output', 'Sub1116_Ses1',...
                '003', 'gt_user-test_corr_surf_stat.txt');
            corr_result = importdata(corr_file);
            corr_result = regexp(corr_result{3}, ':', 'split'); % we look at the min corr
            corr_result = corr_result(2);
            corr_result = str2num(corr_result{1});
            assert(corr_result > 0.99, ...
                sprintf('surface min_corr of run003 is less than 0.99! The value is %f \n', ...
                corr_result))

            
            %% check volume files
            pipe_dir1 = fullfile(getenv('CBIG_TESTDATA_DIR'), 'stable_projects', 'preprocessing', ...
                'CBIG_fMRI_Preproc2016', 'single_subject', 'data');
            pipe_name1 = 'gt';
            pipe_stem1 = '_rest_skip4_stc_mc_residc_interp_FDRMS0.2_DVARS50_bp_0.009_0.08_MNI2mm_sm6_finalmask';

            pipe_dir2 = OutputDir;
            pipe_name2 = 'user-test';
            pipe_stem2 = '_rest_skip4_stc_mc_residc_interp_FDRMS0.2_DVARS50_bp_0.009_0.08_MNI2mm_sm6_finalmask';

            subject_id = 'Sub1116_Ses1';
            runs = {'002', '003'};

            output_dir = fullfile(OutputDir, 'compare_output');
            
            for i = 1:length(runs)                
                CBIG_preproc_compare_two_pipelines(pipe_dir1, pipe_name1, ...
                    pipe_stem1, pipe_dir2, pipe_name2, pipe_stem2, ...
                    subject_id, runs{i}, output_dir, 'vol');
            end
            
            % check vol run 002
            corr_file = fullfile(OutputDir, 'compare_output', 'Sub1116_Ses1',...
                '002', 'gt_user-test_corr_vol_stat.txt');
            corr_result = importdata(corr_file);
            corr_result = regexp(corr_result{3}, ':', 'split'); % we look at the min corr
            corr_result = corr_result(2);
            corr_result = str2num(corr_result{1});
            assert(corr_result > 0.99, ...
                sprintf('volume min_corr of run002 is less than 0.99! The value is %f \n', ...
                corr_result))
            
            % check vol run 003
            corr_file = fullfile(OutputDir, 'compare_output', ...
                'Sub1116_Ses1/003/gt_user-test_corr_vol_stat.txt');
            corr_result = importdata(corr_file);
            corr_result = regexp(corr_result{3}, ':', 'split'); % we look at the min corr
            corr_result = corr_result(2);
            corr_result = str2num(corr_result{1});
            assert(corr_result > 0.99, ...
                sprintf('volume min_corr of run003 is less than 0.99! The value is %f \n', ...
                corr_result))

            %% check FC matrices
            CBIG_preproc_FCmatrices_UnitTestCmp(OutputDir)
            
            % check FC matrices
            FC_cmp_file = fullfile(OutputDir, 'inequal_corr_log.txt');
            FC_cmp = importdata(FC_cmp_file);
            size_FC_cmp = length(FC_cmp.textdata);
            dif = FC_cmp.data(1:(size_FC_cmp - 1));
            labels = FC_cmp.textdata;
            [v,I] = max(dif);
            max_label = labels(I);
            message1 = sprintf('Differences identified in %i FC_matrices! \n', ...
                FC_cmp.data(size_FC_cmp));
            message2 = sprintf('The %s FC_matrices have the largest difference of %d. \n', ...
                max_label{1}, v);
            message3 = sprintf('For more details, refer to inequal_corr_log.txt in OutputDir. \n');
            assert(v ==0, [message1 message2 message3])
                           
            % remove intermediate output data (IMPORTANT)
            rmdir(OutputDir, 's');
            rmpath(fullfile(getenv('CBIG_CODE_DIR'), 'stable_projects', ...
            'preprocessing', 'CBIG_fMRI_Preproc2016', 'utilities')); 
            rmpath(fullfile(getenv('CBIG_CODE_DIR'), 'stable_projects', ...
                'preprocessing', 'CBIG_fMRI_Preproc2016', 'unit_tests', 'single_subject'))
            
        end
        
        
    end
end
