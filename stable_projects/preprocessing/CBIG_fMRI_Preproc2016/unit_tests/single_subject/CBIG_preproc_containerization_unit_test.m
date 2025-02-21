classdef CBIG_preproc_containerization_unit_test < matlab.unittest.TestCase
    %
    % Target project:
    %                 CBIG_fMRI_Preproc2016
    %
    % Case design:
    %                 This unit test is similar to CBIG_preproc_single_subject_unit_test.m,
    %                 except that it is used to test the containerizated version of
    %                 of the pipeline.
    %
    %                 To be able to run this test, you need to first ensure that
    %                 the containerized version of the pipeline has been pushed to
    %                 the Docker Hub repository thomasyeolab/cbig_fmri_preproc2016.
    %
    % Written by Fang Tian and CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md
    
    methods (Test)
        function test_docker_single_subject_Case(testCase)
            %% path setting
            addpath(fullfile(getenv('CBIG_CODE_DIR'), 'stable_projects', ...
                'preprocessing', 'CBIG_fMRI_Preproc2016', 'utilities'));
            addpath(fullfile(getenv('CBIG_CODE_DIR'), 'stable_projects', ...
                'preprocessing', 'CBIG_fMRI_Preproc2016', 'unit_tests', 'single_subject'))
            UnitTestDir = fullfile(getenv('CBIG_CODE_DIR'), 'stable_projects',...
                'preprocessing', 'CBIG_fMRI_Preproc2016', 'unit_tests');
            OutputDir = fullfile(UnitTestDir, 'output', 'single_subject_multi_echo_Case_docker');

            %create output dir (IMPORTANT)
            if(exist(OutputDir, 'dir'))
                rmdir(OutputDir, 's')
            end
            mkdir(OutputDir);
            
            
            %% call Singularity version CBIG_fMRI_Preproc2016 single_subject unit test script to generate results
            cmd = [fullfile(UnitTestDir, 'single_subject', ...
                'CBIG_preproc_unit_tests_call_fMRI_preproc_docker.sh'), ' ', OutputDir];
            system(cmd); % this will submit a job to HPC

            %% check surf files
            check_files('surf', OutputDir, ...
                '_rest_skip8_stc_mc_sdc_me_residc_interp_FDRMS0.3_DVARS60_bp_0.009_0.08_fs6_sm6_fs5');

            %% check volume files
            check_files('vol', OutputDir, ...
                '_rest_skip8_stc_mc_sdc_me_residc_interp_FDRMS0.3_DVARS60_bp_0.009_0.08_MNI2mm_sm6_finalmask');

            %% Check FC matrices
            check_FC_matrices(OutputDir)

            % remove intermediate output data (IMPORTANT)
            rmdir(OutputDir, 's');
            rmpath(fullfile(getenv('CBIG_CODE_DIR'), 'stable_projects', ...
                'preprocessing', 'CBIG_fMRI_Preproc2016', 'utilities'));
            rmpath(fullfile(getenv('CBIG_CODE_DIR'), 'stable_projects', ...
                'preprocessing', 'CBIG_fMRI_Preproc2016', 'unit_tests', 'single_subject'))

        end

        function test_singularity_single_subject_Case(testCase)
            %% path setting
            addpath(fullfile(getenv('CBIG_CODE_DIR'), 'stable_projects', ...
                'preprocessing', 'CBIG_fMRI_Preproc2016', 'utilities'));
            addpath(fullfile(getenv('CBIG_CODE_DIR'), 'stable_projects', ...
                'preprocessing', 'CBIG_fMRI_Preproc2016', 'unit_tests', 'single_subject'))
            UnitTestDir = fullfile(getenv('CBIG_CODE_DIR'), 'stable_projects',...
                'preprocessing', 'CBIG_fMRI_Preproc2016', 'unit_tests');
            OutputDir = fullfile(UnitTestDir, 'output', 'single_subject_multi_echo_Case_singularity');

            % create output dir (IMPORTANT)
            if(exist(OutputDir, 'dir'))
                rmdir(OutputDir, 's')
            end
            mkdir(OutputDir);
            
            
            %% call Singularity version CBIG_fMRI_Preproc2016 single_subject unit test script to generate results
            cmd = [fullfile(UnitTestDir, 'single_subject', ...
                'CBIG_preproc_unit_tests_call_fMRI_preproc_singularity.sh'), ' ', OutputDir];
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

            %% check surf files
            check_files('surf', OutputDir, ...
                '_rest_skip8_stc_mc_sdc_me_residc_interp_FDRMS0.3_DVARS60_bp_0.009_0.08_fs6_sm6_fs5');

            %% check volume files
            check_files('vol', OutputDir, ...
                '_rest_skip8_stc_mc_sdc_me_residc_interp_FDRMS0.3_DVARS60_bp_0.009_0.08_MNI2mm_sm6_finalmask');

            %% Check FC matrices
            check_FC_matrices(OutputDir)

            % remove intermediate output data (IMPORTANT)
            rmdir(OutputDir, 's');
            rmpath(fullfile(getenv('CBIG_CODE_DIR'), 'stable_projects', ...
                'preprocessing', 'CBIG_fMRI_Preproc2016', 'utilities'));
            rmpath(fullfile(getenv('CBIG_CODE_DIR'), 'stable_projects', ...
                'preprocessing', 'CBIG_fMRI_Preproc2016', 'unit_tests', 'single_subject'))

        end
    end
end


function check_files(file_type, OutputDir, pipe_stem)
    pipe_dir1 = fullfile(getenv('CBIG_TESTDATA_DIR'), 'stable_projects', 'preprocessing', ...
        'CBIG_fMRI_Preproc2016', 'single_subject', 'data');
    pipe_name1 = 'gt';
    
    pipe_dir2 = OutputDir;
    pipe_name2 = 'user-test';
    
    subject_id = 'sub005';
    runs = {'001'};
    
    output_dir = fullfile(OutputDir, 'compare_output');
    
    for i = 1:length(runs)
        CBIG_preproc_compare_two_pipelines(pipe_dir1, pipe_name1, ...
            pipe_stem, pipe_dir2, pipe_name2, pipe_stem, subject_id, ...
            runs{i}, output_dir, file_type);
    end
    
    % check correlation result
    corr_file = fullfile(OutputDir, 'compare_output', 'sub005', ...
        '001', sprintf('gt_user-test_corr_%s_stat.txt', file_type));
    corr_result = importdata(corr_file);
    corr_result = regexp(corr_result{3}, ':', 'split'); % we look at the min corr
    corr_result = corr_result(2);
    corr_result = str2num(corr_result{1});
    assert(corr_result > 0.99, ...
        sprintf('%s min_corr of run001 is less than 0.99! The value is %f \n', ...
        file_type, corr_result))
end

function check_FC_matrices(OutputDir)
    % Run the FC matrices unit test comparison
    CBIG_preproc_FCmatrices_UnitTestCmp(OutputDir)

    % Check FC matrices
    FC_cmp_file = fullfile(OutputDir, 'inequal_corr_log.txt');
    FC_cmp = importdata(FC_cmp_file);
    size_FC_cmp = length(FC_cmp.textdata);
    dif = FC_cmp.data(1:(size_FC_cmp - 1));
    labels = FC_cmp.textdata;
    [v,I] = max(dif);
    max_label = labels(I);
    message1 = sprintf('Differences identified in %i FC_matrices! \n', FC_cmp.data(size_FC_cmp));
    message2 = sprintf('The %s FC_matrices have the largest difference of %d. \n', max_label{1}, v);
    message3 = sprintf('For more details, refer to inequal_corr_log.txt in OutputDir. \n');
    assert(v == 0, [message1 message2 message3])
end
