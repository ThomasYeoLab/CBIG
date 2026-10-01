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
    % Written by CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md
    
    methods (Test)        
        function test_motion_filter(testCase)
            addpath(fullfile(getenv('CBIG_CODE_DIR'), 'stable_projects', ...
                'preprocessing', 'CBIG_fMRI_Preproc2016', 'utilities'));
            addpath(fullfile(getenv('CBIG_CODE_DIR'), 'stable_projects', ...
                'preprocessing', 'CBIG_fMRI_Preproc2016', 'unit_tests', 'single_subject'))
            UnitTestDir = fullfile(getenv('CBIG_CODE_DIR'), 'stable_projects',...
                'preprocessing', 'CBIG_fMRI_Preproc2016', 'unit_tests');
            OutputDir = fullfile(UnitTestDir, 'output', 'test_motion_filter');
            load(fullfile(getenv('CBIG_CODE_DIR'), 'unit_tests', 'replace_unittest_flag.txt'));
            
            %create output dir (IMPORTANT)
            if(exist(OutputDir, 'dir'))
                rmdir(OutputDir, 's')
            end
            mkdir(OutputDir);
            
            
            %% call CBIG_fMRI_Preproc2016 single_subject unit test script to generate results
            cmd = [fullfile(UnitTestDir, 'single_subject', ...
                'CBIG_preproc_unit_tests_test_motion_filtering.csh'), ' ', OutputDir];
            system(cmd); % this will submit a job to HPC
            
            
            sub_id = 'sub-NDARBF851NH6';
            out_files = { ...
                fullfile(OutputDir, sub_id, 'bold', '001', [sub_id '_bld001_rest_skip8_stc_mc.par']), ...
                fullfile(OutputDir, sub_id, 'bold', '002', [sub_id '_bld002_rest_skip8_stc_mc.par'])};
            CBIG_check_job_status('CBIG_prep', 'out', out_files);
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
            load(fullfile(getenv('CBIG_CODE_DIR'), 'unit_tests', 'replace_unittest_flag.txt'));
            
            %create output dir (IMPORTANT)
            if(exist(OutputDir, 'dir'))
                rmdir(OutputDir, 's')
            end
            mkdir(OutputDir);
            
            
            %% call CBIG_fMRI_Preproc2016 single_subject unit test script to generate results
            cmd = [fullfile(UnitTestDir, 'single_subject', ...
                'CBIG_preproc_unit_tests_call_fMRI_preproc_aCompCor_deoblique.csh'), ' ', OutputDir];
            system(cmd); % this will submit a job to HPC
            
            
            subject_id = 'Sub1116_Ses1';
            out_files = { ...
                fullfile(OutputDir, subject_id, 'bold', '002', ...
                [subject_id '_bld002_rest_skip4_stc_mc_resid.nii.gz']), ...
                fullfile(OutputDir, subject_id, 'bold', '003', ...
                [subject_id '_bld003_rest_skip4_stc_mc_resid.nii.gz'])};
            CBIG_check_job_status('CBIG_prep', 'out', out_files);

            if(replace_unittest_flag)
                disp('Replacing aCompCor and deoblique preprocessing unit test results...')
                disp('Make sure that the reference directory and its parent directory have write permission!')
                ref_dir = fullfile(getenv('CBIG_TESTDATA_DIR'), 'stable_projects', 'preprocessing', ...
                    'CBIG_fMRI_Preproc2016', 'single_subject', 'data');
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
                    'CBIG_fMRI_Preproc2016', 'single_subject', 'data');
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
            load(fullfile(getenv('CBIG_CODE_DIR'), 'unit_tests', 'replace_unittest_flag.txt'));
            
            %create output dir (IMPORTANT)
            if(exist(OutputDir, 'dir'))
                rmdir(OutputDir, 's')
            end
            mkdir(OutputDir);
            
            
            %% call CBIG_fMRI_Preproc2016 single_subject unit test script to generate results
            cmd = [fullfile(UnitTestDir, 'single_subject', ...
                'CBIG_preproc_unit_tests_call_fMRI_preproc_sdc_oppo_PED.csh'), ' ', OutputDir];
            system(cmd); % this will submit a job to HPC
            
            
            subject_id = 'sub-NDARBF851NH6';
            out_files = { ...
                fullfile(OutputDir, subject_id, 'bold', '001', ...
                [subject_id '_bld001_rest_skip8_stc_mc_sdc.nii.gz']), ...
                fullfile(OutputDir, subject_id, 'bold', '002', ...
                [subject_id '_bld002_rest_skip8_stc_mc_sdc.nii.gz'])};
            CBIG_check_job_status('CBIG_prep', 'out', out_files);

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
        function test_simple_ME_case(testCase)
            %% path setting
            addpath(fullfile(getenv('CBIG_CODE_DIR'), 'stable_projects', ...
                'preprocessing', 'CBIG_fMRI_Preproc2016', 'utilities'));
            addpath(fullfile(getenv('CBIG_CODE_DIR'), 'stable_projects', ...
                'preprocessing', 'CBIG_fMRI_Preproc2016', 'unit_tests', 'single_subject'))
            UnitTestDir = fullfile(getenv('CBIG_CODE_DIR'), 'stable_projects',...
                'preprocessing', 'CBIG_fMRI_Preproc2016', 'unit_tests');
            OutputDir = fullfile(UnitTestDir, 'output', 'single_subject_simple_ME_case');
            load(fullfile(getenv('CBIG_CODE_DIR'), 'unit_tests', 'replace_unittest_flag.txt'));

            %create output dir (IMPORTANT)
            if(exist(OutputDir, 'dir'))
                rmdir(OutputDir, 's')
            end
            mkdir(OutputDir);


            %% call CBIG_fMRI_Preproc2016 single_subject unit test script to generate results
            cmd = [fullfile(UnitTestDir, 'single_subject', ...
                'CBIG_preproc_unit_tests_call_fMRI_preproc_motion_correction.sh'), ' ', OutputDir];
            system(cmd); % this will submit a job to HPC


            subject_id = 'TMR00032_TMSD';
            pipe_stem_final = '_rest_skip4_mc_me_residc_interp_FDRMS0.2_DVARS75_bp_0.009_0.08';
            out_files = { ...
                fullfile(OutputDir, subject_id, 'bold', '001', ...
                [subject_id '_bld001' pipe_stem_final '.nii.gz']), ...
                fullfile(OutputDir, subject_id, 'bold', '002', ...
                [subject_id '_bld002' pipe_stem_final '.nii.gz'])};
            CBIG_check_job_status('CBIG_prep', 'out', out_files);

            runs = {'001', '002'};
            echos = {'e1', 'e2', 'e3'};

            if(replace_unittest_flag)
                disp('Replacing single subject multi echo preprocessing unit test results...')
                disp('Make sure that the reference directory and its parent directory have write permission!')
                ref_dir = fullfile(getenv('CBIG_TESTDATA_DIR'), 'stable_projects', 'preprocessing', ...
                    'CBIG_fMRI_Preproc2016', 'single_subject', 'data');
                ref_dir_subject = fullfile(ref_dir,'TMR00032_TMSD');
                output_dir = fullfile(getenv('CBIG_CODE_DIR'), 'stable_projects', 'preprocessing', ...
                    'CBIG_fMRI_Preproc2016', 'unit_tests', 'output', 'single_subject_simple_ME_case', subject_id);
                if(exist(ref_dir_subject, 'dir'))
                    rmdir(ref_dir_subject, 's')
                end
                movefile(output_dir, ref_dir)
            else
                %% check surf files
                pipe_dir1 = fullfile(getenv('CBIG_TESTDATA_DIR'), 'stable_projects', 'preprocessing', ...
                    'CBIG_fMRI_Preproc2016', 'single_subject', 'data');
                pipe_name1 = 'gt';
                pipe_stem1 = '_rest_skip4_mc_me_residc_interp_FDRMS0.2_DVARS75_bp_0.009_0.08_fs6_sm6_fs5';

                pipe_dir2 = OutputDir;
                pipe_name2 = 'user-test';
                pipe_stem2 = '_rest_skip4_mc_me_residc_interp_FDRMS0.2_DVARS75_bp_0.009_0.08_fs6_sm6_fs5';

                output_dir = fullfile(OutputDir, 'compare_output');

                for i = 1: length(runs)
                    CBIG_preproc_compare_two_pipelines(pipe_dir1, pipe_name1, ...
                        pipe_stem1, pipe_dir2, pipe_name2, pipe_stem2, subject_id, ...
                        runs{i}, output_dir, 'surf');
                    corr_file = fullfile(OutputDir, 'compare_output', subject_id,...
                        runs{i}, 'gt_user-test_corr_surf_stat.txt');
                    corr_result = importdata(corr_file);
                    corr_result = regexp(corr_result{3}, ':', 'split'); % we look at the min corr
                    corr_result = corr_result(2);
                    corr_result = str2num(corr_result{1});
                    assert(corr_result > 0.99, ...
                        sprintf(['surface min_corr of run ' runs{i} ' is less than 0.99! The value is %f \n'], ...
                        corr_result))
                end


                %% check volume files after motion correction
                pipe_stem1 = '_rest_skip4_mc';
                pipe_stem2 = '_rest_skip4_mc';


                for i = 1:length(runs)
                    for j = 1:length(echos)
                        disp(['Comparing results for run ' runs{i} ', ' echos{j} '...'])
                        vol_file1 = fullfile(pipe_dir1, subject_id, 'bold', runs{i}, ...
                            [subject_id '_bld' runs{i} '_' echos{j} pipe_stem1 '.nii.gz']);
                        vol_file2 = fullfile(pipe_dir2, subject_id, 'bold', runs{i}, ...
                            [subject_id '_bld' runs{i} '_' echos{j} pipe_stem2 '.nii.gz']);
                        [corr_vol, vol_size] = CBIG_preproc_compare_two_vols(vol_file1, vol_file2);
                        corr_vol_ex = corr_vol(~isnan(corr_vol));
                        corr_result = min(corr_vol_ex);
                        assert(corr_result > 0.99, ...
                            sprintf(['volume min_corr of run ' runs{i} ', ' echos{j} ...
                            ' ...is less than 0.99! ...The value is %f \n'], corr_result))
                    end
                end


                %% check final volume files
                pipe_stem1 = '_rest_skip4_mc_me_residc_interp_FDRMS0.2_DVARS75_bp_0.009_0.08';
                pipe_stem2 = '_rest_skip4_mc_me_residc_interp_FDRMS0.2_DVARS75_bp_0.009_0.08';

                for i = 1:length(runs)
                    disp(['Comparing results for run ' runs{i} '...'])
                    vol_file1 = fullfile(pipe_dir1, subject_id, 'bold', runs{i}, ...
                        [subject_id '_bld' runs{i} pipe_stem1 '.nii.gz']);
                    vol_file2 = fullfile(pipe_dir2, subject_id, 'bold', runs{i}, ...
                        [subject_id '_bld' runs{i} pipe_stem2 '.nii.gz']);
                    [corr_vol, vol_size] = CBIG_preproc_compare_two_vols(vol_file1, vol_file2);
                    corr_vol_ex = corr_vol(~isnan(corr_vol));
                    corr_result = min(corr_vol_ex);
                    assert(corr_result > 0.99, ...
                        sprintf(['volume min_corr of run ' runs{i} ' is less than 0.99! The value is %f \n'], ...
                        corr_result))
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
