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
    % Written by Yang Qing, Zhang Shaoshi, Lyu Xingyu and CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md
    
    methods (Test)
        function test_single_subject_Case(testCase)
            %% path setting
            addpath(fullfile(getenv('CBIG_CODE_DIR'), 'stable_projects', ...
                'preprocessing', 'CBIG_fMRI_Preproc2016', 'utilities'));
            addpath(fullfile(getenv('CBIG_CODE_DIR'), 'stable_projects', ...
                'preprocessing', 'CBIG_fMRI_Preproc2016', 'unit_tests', 'single_subject'))
            UnitTestDir = fullfile(getenv('CBIG_CODE_DIR'), 'stable_projects',...
                'preprocessing', 'CBIG_fMRI_Preproc2016', 'unit_tests');
            OutputDir = fullfile(UnitTestDir, 'output', 'single_subject_multi_echo_Case');
            load(fullfile(getenv('CBIG_CODE_DIR'), 'unit_tests', 'replace_unittest_flag'));
            
            %create output dir (IMPORTANT)
            if(exist(OutputDir, 'dir'))
                rmdir(OutputDir, 's')
            end
            mkdir(OutputDir);
            
            
            %% call CBIG_fMRI_Preproc2016 single_subject unit test script to generate results
            cmd = [fullfile(UnitTestDir, 'single_subject', ...
                'CBIG_preproc_unit_tests_call_fMRI_preproc.sh'), ' ', OutputDir];
            system(cmd); % this will submit a job to HPC
            
            
            %% we need to periodically check whether the job has finished or not
            cmdout = 1;
            while(cmdout ~= 0)
                cmd = 'ssh headnode "qstat | grep preproc | grep `whoami` | wc -l"';
                [~, cmdout] = system(cmd);
                % after job finishes, cmdout should be 0
                cmdout = str2num(cmdout(1: end-1));
                pause(60); % sleep for 1min and check again
            end
            
            if(replace_unittest_flag)
                disp('Replacing single subject multi echo preprocessing unit test results...')
                disp('Make sure that the reference directory and its parent directory have write permission!')
                ref_dir = fullfile(getenv('CBIG_TESTDATA_DIR'), 'stable_projects', 'preprocessing', ...
                    'CBIG_fMRI_Preproc2016', 'single_subject', 'data');
                ref_dir_subject = fullfile(ref_dir,'sub005');
                output_dir = fullfile(getenv('CBIG_CODE_DIR'), 'stable_projects', 'preprocessing', ...
                    'CBIG_fMRI_Preproc2016', 'unit_tests', 'output', 'single_subject_multi_echo_Case', 'sub005');
                if(exist(ref_dir_subject, 'dir'))
                    rmdir(ref_dir_subject, 's')
                end
                movefile(output_dir, ref_dir)
            else
                %% check surf files
                pipe_dir1 = fullfile(getenv('CBIG_TESTDATA_DIR'), 'stable_projects', 'preprocessing', ...
                    'CBIG_fMRI_Preproc2016', 'single_subject', 'data');
                pipe_name1 = 'gt';
                pipe_stem1 = '_rest_skip8_stc_mc_sdc_me_residc_interp_FDRMS0.3_DVARS60_bp_0.009_0.08_fs6_sm6_fs5';
                
                pipe_dir2 = OutputDir;
                pipe_name2 = 'user-test';
                pipe_stem2 = '_rest_skip8_stc_mc_sdc_me_residc_interp_FDRMS0.3_DVARS60_bp_0.009_0.08_fs6_sm6_fs5';
                
                subject_id = 'sub005';
                runs = {'001'};
                
                
                output_dir = fullfile(OutputDir, 'compare_output');
                
                for i = 1: length(runs)
                    CBIG_preproc_compare_two_pipelines(pipe_dir1, pipe_name1, ...
                        pipe_stem1, pipe_dir2, pipe_name2, pipe_stem2, subject_id, ...
                        runs{i}, output_dir, 'surf');
                end
                
                % check surf run 001
                corr_file = fullfile(OutputDir, 'compare_output', 'sub005',...
                    '001', 'gt_user-test_corr_surf_stat.txt');
                corr_result = importdata(corr_file);
                corr_result = regexp(corr_result{3}, ':', 'split'); % we look at the min corr
                corr_result = corr_result(2);
                corr_result = str2num(corr_result{1});
                assert(corr_result > 0.99, ...
                    sprintf('surface min_corr of run001 is less than 0.99! The value is %f \n', ...
                    corr_result))                
                
                %% check volume files
                pipe_dir1 = fullfile(getenv('CBIG_TESTDATA_DIR'), 'stable_projects', 'preprocessing', ...
                    'CBIG_fMRI_Preproc2016', 'single_subject', 'data');
                pipe_name1 = 'gt';
                pipe_stem1 = '_rest_skip8_stc_mc_sdc_me_residc_interp_FDRMS0.3_DVARS60_bp_0.009_0.08_MNI2mm_sm6_finalmask';
                
                pipe_dir2 = OutputDir;
                pipe_name2 = 'user-test';
                pipe_stem2 = '_rest_skip8_stc_mc_sdc_me_residc_interp_FDRMS0.3_DVARS60_bp_0.009_0.08_MNI2mm_sm6_finalmask';
                
                subject_id = 'sub005';
                runs = {'001'};
                
                output_dir = fullfile(OutputDir, 'compare_output');
                
                for i = 1:length(runs)
                    CBIG_preproc_compare_two_pipelines(pipe_dir1, pipe_name1, ...
                        pipe_stem1, pipe_dir2, pipe_name2, pipe_stem2, ...
                        subject_id, runs{i}, output_dir, 'vol');
                end
                
                % check vol run 001
                corr_file = fullfile(OutputDir, 'compare_output', 'sub005',...
                    '001', 'gt_user-test_corr_vol_stat.txt');
                corr_result = importdata(corr_file);
                corr_result = regexp(corr_result{3}, ':', 'split'); % we look at the min corr
                corr_result = corr_result(2);
                corr_result = str2num(corr_result{1});
                assert(corr_result > 0.99, ...
                    sprintf('volume min_corr of run001 is less than 0.99! The value is %f \n', ...
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
            end
            
            % remove intermediate output data (IMPORTANT)
            rmdir(OutputDir, 's');
            rmpath(fullfile(getenv('CBIG_CODE_DIR'), 'stable_projects', ...
                'preprocessing', 'CBIG_fMRI_Preproc2016', 'utilities'));
            rmpath(fullfile(getenv('CBIG_CODE_DIR'), 'stable_projects', ...
                'preprocessing', 'CBIG_fMRI_Preproc2016', 'unit_tests', 'single_subject'))
            
        end
        
        function test_motion_filter(testCase)
            addpath(fullfile(getenv('CBIG_CODE_DIR'), 'stable_projects', ...
                'preprocessing', 'CBIG_fMRI_Preproc2016', 'utilities'));
            addpath(fullfile(getenv('CBIG_CODE_DIR'), 'stable_projects', ...
                'preprocessing', 'CBIG_fMRI_Preproc2016', 'unit_tests', 'single_subject'))
            UnitTestDir = fullfile(getenv('CBIG_CODE_DIR'), 'stable_projects',...
                'preprocessing', 'CBIG_fMRI_Preproc2016', 'unit_tests');
            OutputDir = fullfile(UnitTestDir, 'output', 'test_motion_filter');
            load(fullfile(getenv('CBIG_CODE_DIR'), 'unit_tests', 'replace_unittest_flag'));
            
            %create output dir (IMPORTANT)
            if(exist(OutputDir, 'dir'))
                rmdir(OutputDir, 's')
            end
            mkdir(OutputDir);
            
            
            %% call CBIG_fMRI_Preproc2016 single_subject unit test script to generate results
            cmd = [fullfile(UnitTestDir, 'single_subject', ...
                'CBIG_preproc_unit_tests_test_motion_filtering.csh'), ' ', OutputDir];
            system(cmd); % this will submit a job to HPC
            
            
            %% we need to periodically check whether the job has finished or not
            cmdout = 1;
            while(cmdout ~= 0)
                cmd = 'ssh headnode "qstat | grep preproc | grep `whoami` | wc -l"';
                [~, cmdout] = system(cmd);
                % after job finishes, cmdout should be 0
                cmdout = str2num(cmdout(1: end-1));
                pause(60); % sleep for 1min and check again
            end
            sub_id = 'sub-NDARBF851NH6';
            ref_dir = fullfile(getenv('CBIG_TESTDATA_DIR'), 'stable_projects', 'preprocessing', ...
                'CBIG_fMRI_Preproc2016', 'single_subject', 'test_motion_filter');
            output_dir = fullfile(getenv('CBIG_CODE_DIR'), 'stable_projects', 'preprocessing', ...
                'CBIG_fMRI_Preproc2016', 'unit_tests', 'output', 'test_motion_filter',sub_id);
            
            if(replace_unittest_flag)
                disp('Replacing single subject preprocessing unit test results...')
                disp('Make sure that the reference directory and its parent directory have write permission!')
                ref_dir_subject = fullfile(ref_dir,'sub-NDARBF851NH6');
                if(exist(ref_dir_subject, 'dir'))
                    rmdir(ref_dir_subject, 's')
                end
                movefile(output_dir, ref_dir)
            else
                ref_dir = fullfile(getenv('CBIG_TESTDATA_DIR'), 'stable_projects', 'preprocessing', ...
                    'CBIG_fMRI_Preproc2016', 'single_subject', 'test_motion_filter', sub_id);
                runs = {'001','002'};
                for i = 1:length(runs)
                    motion_ref_file = fullfile(ref_dir,'bold',runs{i},...
                        [sub_id '_bld' runs{i} '_rest_skip8_stc_mc.par']);
                    motion_ref = load(motion_ref_file);
                    
                    motion_out_file = fullfile(output_dir,'bold',runs{i},...
                        [sub_id '_bld' runs{i} '_rest_skip8_stc_mc.par']);
                    motion_out = load(motion_out_file);
                    
                    max_diff = max(abs(motion_ref(:)-motion_out(:)));
                    assert(max_diff < 1e-6, 'maximum difference of filtered motion parameters greater than 1e-6');
                end
            end
            
            % remove intermediate output data (IMPORTANT)
            
            rmdir(OutputDir, 's');
            rmpath(fullfile(getenv('CBIG_CODE_DIR'), 'stable_projects', ...
                'preprocessing', 'CBIG_fMRI_Preproc2016', 'utilities'));
            rmpath(fullfile(getenv('CBIG_CODE_DIR'), 'stable_projects', ...
                'preprocessing', 'CBIG_fMRI_Preproc2016', 'unit_tests', 'single_subject'))
        end
        
        function test_aCompCor_deoblique(testCase)
            %% path setting
            addpath(fullfile(getenv('CBIG_CODE_DIR'), 'stable_projects', ...
                'preprocessing', 'CBIG_fMRI_Preproc2016', 'utilities'));
            addpath(fullfile(getenv('CBIG_CODE_DIR'), 'stable_projects', ...
                'preprocessing', 'CBIG_fMRI_Preproc2016', 'unit_tests', 'single_subject'))
            UnitTestDir = fullfile(getenv('CBIG_CODE_DIR'), 'stable_projects',...
                'preprocessing', 'CBIG_fMRI_Preproc2016', 'unit_tests');
            OutputDir = fullfile(UnitTestDir, 'output', 'aCompCor_deoblique');
            load(fullfile(getenv('CBIG_CODE_DIR'), 'unit_tests', 'replace_unittest_flag'));
            
            %create output dir (IMPORTANT)
            if(exist(OutputDir, 'dir'))
                rmdir(OutputDir, 's')
            end
            mkdir(OutputDir);
            
            
            %% call CBIG_fMRI_Preproc2016 single_subject unit test script to generate results
            cmd = [fullfile(UnitTestDir, 'single_subject', ...
                'CBIG_preproc_unit_tests_call_fMRI_preproc_aCompCor_deoblique.csh'), ' ', OutputDir];
            system(cmd); % this will submit a job to HPC
            
            
            %% we need to periodically check whether the job has finished or not
            cmdout = 1;
            while(cmdout ~= 0)
                cmd = 'ssh headnode "qstat | grep preproc | grep `whoami` | wc -l"';
                [~, cmdout] = system(cmd);
                % after job finishes, cmdout should be 0
                cmdout = str2num(cmdout(1: end-1));
                pause(60); % sleep for 1min and check again
            end
            
            if(replace_unittest_flag)
                disp('Replacing aCompCor and deoblique preprocessing unit test results...')
                disp('Make sure that the reference directory and its parent directory have write permission!')
                ref_dir = fullfile(getenv('CBIG_TESTDATA_DIR'), 'stable_projects', 'preprocessing', ...
                    'CBIG_fMRI_Preproc2016', 'single_subject', 'aCompCor_deoblique');
                output_dir = fullfile(getenv('CBIG_CODE_DIR'), 'stable_projects', 'preprocessing', ...
                    'CBIG_fMRI_Preproc2016', 'unit_tests', 'output', 'aCompCor_deoblique', 'Sub1116_Ses1');
                ref_dir_subject = fullfile(ref_dir,'Sub1116_Ses1');
                if(exist(ref_dir_subject, 'dir'))
                    rmdir(ref_dir_subject, 's')
                end
                movefile(output_dir, ref_dir)
            else
                %% check volume files
                pipe_dir1 = fullfile(getenv('CBIG_TESTDATA_DIR'), 'stable_projects', 'preprocessing', ...
                    'CBIG_fMRI_Preproc2016', 'single_subject', 'aCompCor_deoblique');
                pipe_stem1 = '_rest_skip4_stc_mc_resid';
                
                pipe_dir2 = OutputDir;
                pipe_stem2 = '_rest_skip4_stc_mc_resid';
                
                subject_id = 'Sub1116_Ses1';
                runs = {'002', '003'};
                
                % check deoblique
                vol_test_file = fullfile(pipe_dir2, subject_id, 'bold', runs{1}, ...
                    [subject_id '_bld' runs{1} '_rest.nii.gz']);
                
                cmd = ['3dinfo -is_oblique ' vol_test_file];
                [~, obl_test] = system(cmd);
                if(contains(obl_test, '1'))
                    assert('Output data is oblique.');
                end
                
                cmd = ['3dinfo -orient ' vol_test_file];
                [~, ori_test] = system(cmd);
                if(~contains(ori_test, 'RPI'))
                    assert(['Orientation is ' ori_test]);
                end
                
                % check final volume
                for i = 1:length(runs)
                    vol_ref_file = fullfile(pipe_dir1,subject_id,'bold',runs{i},...
                        [subject_id '_bld' runs{i} pipe_stem1 '.nii.gz']);
                    vol_ref = MRIread(vol_ref_file);
                    vol_ref = vol_ref.vol;
                    
                    vol_test_file = fullfile(pipe_dir2,subject_id,'bold',runs{i},...
                        [subject_id '_bld' runs{i} pipe_stem2 '.nii.gz']);
                    vol_test = MRIread(vol_test_file);
                    vol_test = vol_test.vol;
                    
                    max_diff = max(abs(vol_ref(:)-vol_test(:)));
                    assert(max_diff < 1e-6, 'maximum difference greater than 1e-6');
                end
            end
            
            % remove intermediate output data (IMPORTANT)
            rmdir(OutputDir, 's');
            rmpath(fullfile(getenv('CBIG_CODE_DIR'), 'stable_projects', ...
                'preprocessing', 'CBIG_fMRI_Preproc2016', 'utilities'));
            rmpath(fullfile(getenv('CBIG_CODE_DIR'), 'stable_projects', ...
                'preprocessing', 'CBIG_fMRI_Preproc2016', 'unit_tests', 'single_subject'))
            
        end
        
        function test_sdc_oppo_PED(testCase)
            %% path setting
            addpath(fullfile(getenv('CBIG_CODE_DIR'), 'stable_projects', ...
                'preprocessing', 'CBIG_fMRI_Preproc2016', 'utilities'));
            addpath(fullfile(getenv('CBIG_CODE_DIR'), 'stable_projects', ...
                'preprocessing', 'CBIG_fMRI_Preproc2016', 'unit_tests', 'single_subject'))
            UnitTestDir = fullfile(getenv('CBIG_CODE_DIR'), 'stable_projects',...
                'preprocessing', 'CBIG_fMRI_Preproc2016', 'unit_tests');
            OutputDir = fullfile(UnitTestDir, 'output', 'sdc_oppo_PED');
            load(fullfile(getenv('CBIG_CODE_DIR'), 'unit_tests', 'replace_unittest_flag'));
            
            %create output dir (IMPORTANT)
            if(exist(OutputDir, 'dir'))
                rmdir(OutputDir, 's')
            end
            mkdir(OutputDir);
            
            
            %% call CBIG_fMRI_Preproc2016 single_subject unit test script to generate results
            cmd = [fullfile(UnitTestDir, 'single_subject', ...
                'CBIG_preproc_unit_tests_call_fMRI_preproc_sdc_oppo_PED.csh'), ' ', OutputDir];
            system(cmd); % this will submit a job to HPC
            
            
            %% we need to periodically check whether the job has finished or not
            cmdout = 1;
            while(cmdout ~= 0)
                cmd = 'ssh headnode "qstat | grep preproc | grep `whoami` | wc -l"';
                [~, cmdout] = system(cmd);
                % after job finishes, cmdout should be 0
                cmdout = str2num(cmdout(1: end-1));
                pause(60); % sleep for 1min and check again
            end
            
            if(replace_unittest_flag)
                disp('Replacing spatial distortion correction preprocessing unit test results...')
                disp('Make sure that the reference directory and its parent directory have write permission!')
                ref_dir = fullfile(getenv('CBIG_TESTDATA_DIR'), 'stable_projects', 'preprocessing', ...
                    'CBIG_fMRI_Preproc2016', 'single_subject', 'data');
                output_dir = fullfile(getenv('CBIG_CODE_DIR'), 'stable_projects', 'preprocessing', ...
                    'CBIG_fMRI_Preproc2016', 'unit_tests', 'output', 'sdc_oppo_PED', 'sub-NDARBF851NH6');
                ref_dir_subject = fullfile(ref_dir,'sub-NDARBF851NH6');
                if(exist(ref_dir_subject, 'dir'))
                    rmdir(ref_dir_subject, 's')
                end
                movefile(output_dir, ref_dir)
            else
                %% check volume files
                pipe_dir1 = fullfile(getenv('CBIG_TESTDATA_DIR'), 'stable_projects', 'preprocessing', ...
                    'CBIG_fMRI_Preproc2016', 'single_subject', 'data');
                pipe_stem1 = '_rest_skip8_stc_mc_sdc';
                
                pipe_dir2 = OutputDir;
                pipe_stem2 = '_rest_skip8_stc_mc_sdc';
                
                subject_id = 'sub-NDARBF851NH6';
                runs = {'001','002'};
                
                for i = 1:length(runs)
                    vol_ref_file = fullfile(pipe_dir1,subject_id,'bold',runs{i},...
                        [subject_id '_bld' runs{i} pipe_stem1 '.nii.gz']);
                    vol_ref = MRIread(vol_ref_file);
                    vol_ref = vol_ref.vol;
                    
                    vol_test_file = fullfile(pipe_dir2,subject_id,'bold',runs{i},...
                        [subject_id '_bld' runs{i} pipe_stem2 '.nii.gz']);
                    vol_test = MRIread(vol_test_file);
                    vol_test = vol_test.vol;
                    
                    max_diff = max(abs(vol_ref(:)-vol_test(:)));
                    assert(max_diff < 1e-6, 'maximum difference greater than 1e-6 in run %s', runs{i});
                end
            end
            
            % remove intermediate output data (IMPORTANT)
            rmdir(OutputDir, 's');
            rmpath(fullfile(getenv('CBIG_CODE_DIR'), 'stable_projects', ...
                'preprocessing', 'CBIG_fMRI_Preproc2016', 'utilities'));
            rmpath(fullfile(getenv('CBIG_CODE_DIR'), 'stable_projects', ...
                'preprocessing', 'CBIG_fMRI_Preproc2016', 'unit_tests', 'single_subject'))
            
        end
    end
end
