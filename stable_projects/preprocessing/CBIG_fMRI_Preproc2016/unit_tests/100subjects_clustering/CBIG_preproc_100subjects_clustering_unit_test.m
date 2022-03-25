classdef CBIG_preproc_100subjects_clustering_unit_test < matlab.unittest.TestCase
%
% Target project:
%                 CBIG_fMRI_Preproc2016
%
% Case design:
%                 CBIG_fMRI_Preproc2016 already have unit tests, here we
%                 call its "100subjects_clustering" unit test and automatically
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
% Written by Yang Qing, Shaoshi Zhang and CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md

    methods (Test)
        function test_100_subjects_Case(testCase)
            %% path setting
            addpath(fullfile(getenv('CBIG_CODE_DIR'), 'stable_projects', ...
            'preprocessing', 'CBIG_fMRI_Preproc2016', 'utilities')); 
            addpath(fullfile(getenv('CBIG_CODE_DIR'), 'stable_projects', ...
                'preprocessing', 'CBIG_fMRI_Preproc2016', 'unit_tests', '100subjects_clustering'))
            UnitTestDir = fullfile(getenv('CBIG_CODE_DIR'), 'stable_projects',...
                'preprocessing', 'CBIG_fMRI_Preproc2016', 'unit_tests');
            OutputDir = fullfile(UnitTestDir, 'output', '100_subjects_Case'); 
            load(fullfile(getenv('CBIG_CODE_DIR'), 'unit_tests', 'replace_unittest_flag'));
            
            %create output dir (IMPORTANT)
            if(exist(OutputDir, 'dir'))
                rmdir(OutputDir, 's')
            end
            mkdir(OutputDir);
			%% generate fmrinii list
			cmd = [fullfile(UnitTestDir, '100subjects_clustering', ...
				'CBIG_preproc_unit_tests_generate_fmrinii_list.sh'), ' ', OutputDir];
			system(cmd);
            
            %% call CBIG_fMRI_Preproc2016 100 subjects unit test script to preprocess subjects
            cmd = [fullfile(UnitTestDir, '100subjects_clustering', ...
                'CBIG_preproc_unit_tests_preprocess_100subjects.csh'), ' ', OutputDir];
            system(cmd); % this will submit a bunch of jobs to HPC
            
            %% periodically check whether the job has finished or not             
            cmdout = 1;
            while(cmdout ~= 0)
                cmd = 'ssh headnode "qstat | grep prep_100sub_ut | grep `whoami` | wc -l"';
                [~, cmdout] = system(cmd);
                % after job finishes, cmdout should be 0
                cmdout = str2num(cmdout(1: end-1));
                pause(60); % sleep for 1min and check again
            end
            
            %% call CBIG_fMRI_Preproc2016 100 subjects unit test script to obtain clustering results
            cmd = [fullfile(UnitTestDir, '100subjects_clustering', ...
                'CBIG_preproc_unit_tests_general_cluster_GSP_80_low_motion+20_w_censor.csh'), ...
                ' ', OutputDir, ' ', OutputDir];
            system(cmd);
            
            %% periodically check whether the job has finished or not             
            cmdout = 1;
            while(cmdout ~= 0)
                cmd = 'ssh headnode "qstat | grep clust_100sub_ut | grep `whoami` | wc -l"';
                [~, cmdout] = system(cmd);
                cmdout = str2num(cmdout(1: end-1));
                pause(60);
            end
            
            if(replace_unittest_flag)
                disp('Replacing 100 subjects clustering unit test results...')
                disp('Make sure that the refernece directory and its parent directory have write permission!')
                preproc_ref_dir = fullfile(getenv('CBIG_TESTDATA_DIR'), 'stable_projects', 'preprocessing', ...
                    'CBIG_fMRI_Preproc2016', '100subjects_clustering', 'preproc_out');
                clustering_ref_dir = fullfile(getenv('CBIG_TESTDATA_DIR'), 'stable_projects' ,'preprocessing', ...
                    'CBIG_fMRI_Preproc2016', '100subjects_clustering', 'clustering');
                clustering_output_dir = fullfile(getenv('CBIG_CODE_DIR'), 'stable_projects', 'preprocessing', ...
                    'CBIG_fMRI_Preproc2016/unit_tests', 'output', '100_subjects_Case', 'clustering');
                % replace clustering results
                movefile(clustering_output_dir, clustering_ref_dir);
               
                % load subject list and replace preprocessing results
                fid = fopen(fullfile(getenv('CBIG_TESTDATA_DIR'), 'stable_projects', 'preprocessing', ...
                    'CBIG_fMRI_Preproc2016', '100subjects_clustering', 'subject_list.txt'));
                subject_list = textscan(fid, '%s');
                subject_list = subject_list{1};
                for i = 1:length(subject_list)
                    subject_id = subject_list{i};
                    preproc_out_dir = fullfile(getenv('CBIG_CODE_DIR'), 'stable_projects', 'preprocessing', ...
                        'CBIG_fMRI_Preproc2016/unit_tests', 'output', '100_subjects_Case', subject_id);
                    movefile(preproc_out_dir, fullfile(preproc_ref_dir, subject_id))                    
                end
            else
                %% check clustering results
                CBIG_preproc_unit_tests_cmp_clusters(fullfile(OutputDir, ...
                    'clustering/GSP_80_low_mt_20_w_censor_clusters017_scrub.mat'), OutputDir);
            end
            
            % remove intermediate output data (IMPORTANT)
            rmdir(OutputDir, 's');
            rmpath(fullfile(getenv('CBIG_CODE_DIR'), 'stable_projects', ...
            'preprocessing', 'CBIG_fMRI_Preproc2016', 'utilities')); 
            rmpath(fullfile(getenv('CBIG_CODE_DIR'), 'stable_projects', ...
                'preprocessing', 'CBIG_fMRI_Preproc2016', 'unit_tests', '100subjects_clustering'))
            
        end
        
        
    end
end
