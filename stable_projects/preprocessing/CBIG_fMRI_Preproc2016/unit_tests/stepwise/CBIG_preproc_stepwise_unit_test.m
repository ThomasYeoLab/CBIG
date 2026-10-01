classdef CBIG_preproc_stepwise_unit_test < matlab.unittest.TestCase
    %
    % Target project:
    %                 CBIG_fMRI_Preproc2016
    %
    % Case design:
    %                 Tests each preprocessing step independently. For each step:
    %                   - Input is taken from the reference data of the previous
    %                     step (or raw input_truncate data for the first step)
    %                   - The step script is run locally on a temp copy
    %                   - Output is compared against the reference data for
    %                     that step (max absolute difference < 1e-6)
    %
    %                 Subject: sub005 (multi-echo, run 001) following
    %                 the prepro.config pipeline stages.
    %
    % Written by XUE Aihuiping and CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md

    properties
        TestDataDir  % .../single_subject reference data root
        CodeDir      % CBIG_fMRI_Preproc2016 project root
        UnitTestDir
        SubjectID    % sub005
        Run          % '001'
        Echoes       % {'e1', 'e2', 'e3'}
    end

    methods (TestMethodSetup)
        function setup(testCase)
            ss_dir = fullfile(getenv('CBIG_TESTDATA_DIR'), 'stable_projects', ...
                'preprocessing', 'CBIG_fMRI_Preproc2016', 'single_subject');

            testCase.TestDataDir = ss_dir;
            testCase.CodeDir = fullfile(getenv('CBIG_CODE_DIR'), ...
                'stable_projects', 'preprocessing', 'CBIG_fMRI_Preproc2016');
            testCase.UnitTestDir = fullfile(testCase.CodeDir, 'unit_tests');
            testCase.SubjectID = 'sub005';
            testCase.Run = '001';
            testCase.Echoes = {'e1', 'e2', 'e3'};

            addpath(fullfile(testCase.CodeDir, 'utilities'));
        end
    end

    methods (TestMethodTeardown)
        function teardown(testCase)
            rmpath(fullfile(testCase.CodeDir, 'utilities'));
        end
    end

    methods (Test)

        function test_deoblique(testCase)
            % Input:  raw echo files from input_truncate/sub005/
            % Script: CBIG_preproc_deoblique.sh -i <raw> -o <tmp>
            % Ref:    data/sub005/bold/001/<subject>_bld001_e{1,2,3}_rest.nii.gz

            input_dir = fullfile(testCase.TestDataDir, 'data', 'input_truncate', testCase.SubjectID);
            ref_dir = fullfile(testCase.TestDataDir, 'data', testCase.SubjectID, 'bold', testCase.Run);
            output_dir = fullfile(testCase.UnitTestDir, 'output', 'stepwise_deoblique');
            script = fullfile(testCase.CodeDir, 'utilities', 'CBIG_preproc_deoblique.sh');

            if(exist(output_dir, 'dir'))
                rmdir(output_dir, 's')
            end
            mkdir(output_dir);

            fprintf('\n');
            disp('========== CBIG_preproc_deoblique ==========');
            echo = testCase.Echoes{1};
            file_name = [testCase.SubjectID '_bld' testCase.Run '_' echo '_rest.nii.gz'];
            input_file = fullfile(input_dir, file_name);
            output_file = fullfile(output_dir, file_name);
            ref_file = fullfile(ref_dir, file_name);

            cmd = [script ' -i ' input_file ' -o ' output_file];
            disp(cmd);
            [status, ~] = system(cmd);
            assert(status == 0, 'Deoblique script failed for echo %s', echo);

            [~, out_hd] = system(['fslhd ' output_file]);
            [~, ref_hd] = system(['fslhd ' ref_file]);
            out_lines = strsplit(out_hd, newline);
            ref_lines = strsplit(ref_hd, newline);
            out_xfm = out_lines(startsWith(out_lines, 'qto_xyz') | startsWith(out_lines, 'sto_xyz'));
            ref_xfm = ref_lines(startsWith(ref_lines, 'qto_xyz') | startsWith(ref_lines, 'sto_xyz'));
            disp('Reference transform:'); disp(strjoin(ref_xfm, newline));
            disp('Output transform:');    disp(strjoin(out_xfm, newline));
            assert(isequal(out_xfm, ref_xfm), 'Transform header mismatch for echo %s', echo);

            rmdir(output_dir, 's');
        end

        function test_skip(testCase)
            % Input:  deobliqued echo files from data/sub005/bold/001/
            % Script: CBIG_preproc_skip.csh -s sub005 -d <tmp> -bld '001'
            %         -BOLD_stem _rest -skip 8 -echo_number 3
            % Ref:    data/sub005/bold/001/<subject>_bld001_e{1,2,3}_rest_skip8.nii.gz

            ref_bold_dir = fullfile(testCase.TestDataDir, 'data', testCase.SubjectID, ...
                'bold', testCase.Run);
            output_dir = fullfile(testCase.UnitTestDir, 'output', 'stepwise_skip');
            script = fullfile(testCase.CodeDir, 'CBIG_preproc_skip.csh');

            if(exist(output_dir, 'dir'))
                rmdir(output_dir, 's')
            end
            mkdir(output_dir);

            fprintf('\n');
            disp('========== CBIG_preproc_skip ==========');

            bold_dir = fullfile(output_dir, testCase.SubjectID, 'bold', testCase.Run);
            log_dir = fullfile(output_dir, testCase.SubjectID, 'logs');
            mkdir(bold_dir);
            mkdir(log_dir);

            for i = 1:length(testCase.Echoes)
                file_name = [testCase.SubjectID '_bld' testCase.Run '_' testCase.Echoes{i} '_rest.nii.gz'];
                copyfile(fullfile(ref_bold_dir, file_name), fullfile(bold_dir, file_name));
            end

            cmd = [script ' -s ' testCase.SubjectID ' -d ' output_dir ...
                ' -bld ''' testCase.Run ''' -BOLD_stem _rest -skip 8 -echo_number 3'];
            disp(cmd);
            [status, ~] = system(cmd);
            assert(status == 0, 'Skip script failed');

            for i = 1:length(testCase.Echoes)
                echo  = testCase.Echoes{i};
                file_name = [testCase.SubjectID '_bld' testCase.Run '_' echo '_rest_skip8.nii.gz'];
                out = MRIread(fullfile(bold_dir, file_name));
                ref = MRIread(fullfile(ref_bold_dir, file_name));
                disp(['Skip ' echo ]);
                n_out = size(out.vol, 4); n_ref = size(ref.vol, 4);
                disp(['  Output: ' num2str(n_out) ' frames, reference: ' num2str(n_ref) ' frames.']);
                assert(isequal(size(out.vol), size(ref.vol)), ...
                    'Skip volume size mismatch for %s', echo);
            end

            rmdir(output_dir, 's');
        end

        function test_stc(testCase)
            % Input:  skip8 echo files from data/sub005/bold/001/
            % Script: CBIG_preproc_fslslicetimer.csh -s sub005 -d <output_dir> -bld '001'
            %         -BOLD_stem _rest_skip8 -echo_number 3 -slice_timing e1.txt,e2.txt,e3.txt
            % Ref:    data/sub005/bold/001/<subject>_bld001_e{1,2,3}_rest_skip8_stc.nii.gz

            ref_bold_dir = fullfile(testCase.TestDataDir, 'data', testCase.SubjectID, ...
                'bold', testCase.Run);
            output_dir = fullfile(testCase.UnitTestDir, 'output', 'stepwise_stc');
            script = fullfile(testCase.CodeDir, 'CBIG_preproc_fslslicetimer.csh');

            if exist(output_dir, 'dir')
                rmdir(output_dir, 's')
            end
            mkdir(output_dir);

            fprintf('\n');
            disp('========== CBIG_preproc_fslslicetimer ==========');

            bold_dir = fullfile(output_dir, testCase.SubjectID, 'bold', testCase.Run);
            log_dir  = fullfile(output_dir, testCase.SubjectID, 'logs');
            mkdir(bold_dir);
            mkdir(log_dir);

            for i = 1:length(testCase.Echoes)
                file_name = [testCase.SubjectID '_bld' testCase.Run '_' testCase.Echoes{i} '_rest_skip8.nii.gz'];
                copyfile(fullfile(ref_bold_dir, file_name), fullfile(bold_dir, file_name));
            end

            slice_timing_dir = fullfile(testCase.TestDataDir, 'scripts', 'slice_timing');
            st_files = cellfun(@(e) fullfile(slice_timing_dir, ...
                [testCase.SubjectID '_slice_timing_' e '.txt']), testCase.Echoes, ...
                'UniformOutput', false);

            cmd = [script ' -s ' testCase.SubjectID ' -d ' output_dir ...
                ' -bld ''' testCase.Run ''' -BOLD_stem _rest_skip8 -echo_number 3' ...
                ' -slice_timing ' strjoin(st_files, ',')];
            disp(cmd);
            [status, ~] = system(cmd);
            assert(status == 0, 'STC script failed');

            mask_idx = loadBrainMask(testCase);
            for i = 1:length(testCase.Echoes)
                echo  = testCase.Echoes{i};
                file_name = [testCase.SubjectID '_bld' testCase.Run '_' echo '_rest_skip8_stc.nii.gz'];
                disp(['STC ' echo]);
                stats = CBIG_compare_volumes(fullfile(bold_dir, file_name), ...
                    fullfile(ref_bold_dir, file_name), fullfile(output_dir, 'compare_output'), mask_idx);
                assert(stats.max_diff < 1e-6, ...
                    'STC max diff for %s = %e (threshold 1e-6)', echo, stats.max_diff);
            end

            rmdir(output_dir, 's');
        end

        function test_mc(testCase)
            % Input:  stc echo files from data/sub005/bold/001/
            % Script: CBIG_preproc_fslmcflirt_outliers.csh -s sub005 -d <output_dir> -bld '001'
            %         -BOLD_stem _rest_skip8_stc -echo_number 3
            %         -FD_th 0.3 -DV_th 60 -discard-run 50 -rm-seg 5
            % Ref:    data/sub005/bold/001/<subject>_bld001_e{1,2,3}_rest_skip8_stc_mc.nii.gz
            %         data/sub005/bold/001/<subject>_bld001_e1_rest_skip8_stc_mc.par

            ref_bold_dir = fullfile(testCase.TestDataDir, 'data', testCase.SubjectID, ...
                'bold', testCase.Run);
            output_dir = fullfile(testCase.UnitTestDir, 'output', 'stepwise_mc');
            script = fullfile(testCase.CodeDir, 'CBIG_preproc_fslmcflirt_outliers.csh');

            if exist(output_dir, 'dir')
                rmdir(output_dir, 's')
            end
            mkdir(output_dir);

            fprintf('\n');
            disp('========== CBIG_preproc_fslmcflirt_outliers ==========');

            bold_dir = fullfile(output_dir, testCase.SubjectID, 'bold', testCase.Run);
            log_dir  = fullfile(output_dir, testCase.SubjectID, 'logs');
            mkdir(bold_dir);
            mkdir(log_dir);

            % The MC script reads sub005.bold to get the run list for applying
            % transforms to e2/e3; the main pipeline creates it but we must do it here
            fid = fopen(fullfile(log_dir, [testCase.SubjectID '.bold']), 'w');
            fprintf(fid, '%s\n', testCase.Run);
            fclose(fid);

            for i = 1:length(testCase.Echoes)
                file_name = [testCase.SubjectID '_bld' testCase.Run '_' testCase.Echoes{i} '_rest_skip8_stc.nii.gz'];
                copyfile(fullfile(ref_bold_dir, file_name), fullfile(bold_dir, file_name));
            end

            cmd = [script ' -s ' testCase.SubjectID ' -d ' output_dir ...
                ' -bld ''' testCase.Run ''' -BOLD_stem _rest_skip8_stc -echo_number 3' ...
                ' -FD_th 0.3 -DV_th 60 -discard-run 50 -rm-seg 5'];
            disp(cmd);
            [status, ~] = system(cmd);
            assert(status == 0, 'MC script failed');

            mask_idx = loadBrainMask(testCase);
            for i = 1:length(testCase.Echoes)
                echo  = testCase.Echoes{i};
                file_name = [testCase.SubjectID '_bld' testCase.Run '_' echo '_rest_skip8_stc_mc.nii.gz'];
                disp(['MC ' echo]);
                stats = CBIG_compare_volumes(fullfile(bold_dir, file_name), ...
                    fullfile(ref_bold_dir, file_name), fullfile(output_dir, 'compare_output'), mask_idx);
                assert(stats.max_diff < 1e-6, 'MC max diff for %s = %e (threshold 1e-6)', echo, stats.max_diff);
            end

            % Check motion parameters (e1 only — MC is run on e1 and applied to all echoes)
            par_file_name = [testCase.SubjectID '_bld' testCase.Run '_e1_rest_skip8_stc_mc.par'];
            out_par  = load(fullfile(bold_dir, par_file_name));
            ref_par  = load(fullfile(ref_bold_dir, par_file_name));
            par_diff = max(abs(out_par(:) - ref_par(:)));
            disp(['MC par max diff: ' num2str(par_diff)]);
            assert(par_diff < 1e-6, 'MC .par max diff = %e (threshold 1e-6)', par_diff);

            rmdir(output_dir, 's');
        end

        function test_sdc(testCase)
            % Input:  pre-MC echo files + e1 MC .cat (motion correction matrices)
            % Script: CBIG_preproc_spatial_distortion_correction.csh -s sub005 -d <output_dir>
            %         -bld '001' -echo_number 3 -BOLD_stem _rest_skip8_stc
            %         -fpm mag+phasediff -m magnitude1.nii.gz -p phase_diff.nii.gz
            %         -delta 2.46 -ees 0.25 -te 12
            % Ref:    data/sub005/bold/001/<subject>_bld001_e{1,2,3}_rest_skip8_stc_mc_sdc.nii.gz

            ref_bold_dir = fullfile(testCase.TestDataDir, 'data', testCase.SubjectID, ...
                'bold', testCase.Run);
            input_truncate_dir = fullfile(testCase.TestDataDir, 'data', 'input_truncate', testCase.SubjectID);
            output_dir = fullfile(testCase.UnitTestDir, 'output', 'stepwise_sdc');
            script = fullfile(testCase.CodeDir, 'CBIG_preproc_spatial_distortion_correction.csh');

            if exist(output_dir, 'dir'); rmdir(output_dir, 's'); end
            mkdir(output_dir);

            fprintf('\n');
            disp('========== CBIG_preproc_spatial_distortion_correction ==========');

            bold_dir = fullfile(output_dir, testCase.SubjectID, 'bold', testCase.Run);
            log_dir  = fullfile(output_dir, testCase.SubjectID, 'logs');
            mkdir(bold_dir); mkdir(log_dir);

            % Copy pre-MC volumes (SDC applies MC+SDC simultaneously to reduce interpolations)
            for i = 1:length(testCase.Echoes)
                file_name = [testCase.SubjectID '_bld' testCase.Run '_' testCase.Echoes{i} '_rest_skip8_stc.nii.gz'];
                copyfile(fullfile(ref_bold_dir, file_name), fullfile(bold_dir, file_name));
            end

            % Copy e1 MC concatenated transformation matrices (required by applywarp --premat)
            cat_file = [testCase.SubjectID '_bld' testCase.Run '_e1_rest_skip8_stc_mc.cat'];
            copyfile(fullfile(ref_bold_dir, cat_file), fullfile(bold_dir, cat_file));

            mag_f   = fullfile(input_truncate_dir, 'magnitude1.nii.gz');
            phase_f = fullfile(input_truncate_dir, 'phase_diff.nii.gz');

            cmd = [script ' -s ' testCase.SubjectID ' -d ' output_dir ...
                ' -bld ''' testCase.Run ''' -echo_number 3 -BOLD_stem _rest_skip8_stc' ...
                ' -fpm mag+phasediff -m ' mag_f ' -p ' phase_f ...
                ' -delta 2.46 -ees 0.25 -te 12'];
            disp(cmd);
            [status, ~] = system(cmd);
            assert(status == 0, 'SDC script failed');

            mask_idx = loadBrainMask(testCase);

            % Compare field map: phase_difference_rad_unwrap_sec.nii.gz (3D, rad/s)
            sdc_dir     = fullfile(output_dir, testCase.SubjectID, 'bold', 'sdc');
            ref_sdc_dir = fullfile(testCase.TestDataDir, 'data', testCase.SubjectID, 'bold', 'sdc');
            fmap_name   = 'phase_difference_rad_unwrap_sec.nii.gz';
            disp('Field map');
            fmap_stats = CBIG_compare_volumes(fullfile(sdc_dir, fmap_name), ...
                fullfile(ref_sdc_dir, fmap_name), fullfile(output_dir, 'compare_output'), mask_idx);
            assert(fmap_stats.max_diff < 1e-6, ...
                'Field map max diff = %e (threshold 1e-6)', fmap_stats.max_diff);

            for i = 1:length(testCase.Echoes)
                echo  = testCase.Echoes{i};
                file_name = [testCase.SubjectID '_bld' testCase.Run '_' echo '_rest_skip8_stc_mc_sdc.nii.gz'];
                disp(['SDC ' echo]);
                stats = CBIG_compare_volumes(fullfile(bold_dir, file_name), ...
                    fullfile(ref_bold_dir, file_name), fullfile(output_dir, 'compare_output'), mask_idx);
                assert(stats.max_diff < 1e-6, 'SDC max diff for %s = %e (threshold 1e-6)', echo, stats.max_diff);
            end

            rmdir(output_dir, 's');
        end

        function test_me_denoise(testCase)
            % Input:  SDC-corrected echo files from data/sub005/bold/001/
            % Script: CBIG_preproc_multiecho_denoise.csh -s sub005 -d <output_dir> -bld '001'
            %         -BOLD_stem _rest_skip8_stc_mc_sdc -echo_number 3 -echo_time 12,30.11,48.22
            % Ref:    data/sub005/bold/001/sub005_bld001_rest_skip8_stc_mc_sdc_me.nii.gz

            ref_bold_dir = fullfile(testCase.TestDataDir, 'data', testCase.SubjectID, ...
                'bold', testCase.Run);
            output_dir = fullfile(testCase.UnitTestDir, 'output', 'stepwise_me_denoise');
            script = fullfile(testCase.CodeDir, 'CBIG_preproc_multiecho_denoise.csh');

            if exist(output_dir, 'dir'); rmdir(output_dir, 's'); end
            mkdir(output_dir);

            fprintf('\n');
            disp('========== CBIG_preproc_multiecho_denoise ==========');

            bold_dir = fullfile(output_dir, testCase.SubjectID, 'bold', testCase.Run);
            log_dir  = fullfile(output_dir, testCase.SubjectID, 'logs');
            qc_dir   = fullfile(output_dir, testCase.SubjectID, 'qc');
            mkdir(bold_dir); mkdir(log_dir); mkdir(qc_dir);

            for i = 1:length(testCase.Echoes)
                echo = testCase.Echoes{i};
                file_name = [testCase.SubjectID '_bld' testCase.Run '_' echo '_rest_skip8_stc_mc_sdc.nii.gz'];
                copyfile(fullfile(ref_bold_dir, file_name), fullfile(bold_dir, file_name));
            end

            % tedana requires the CBIG_py3 conda env; use conda run to activate it
            conda_run = [getenv('HOME') '/storage/miniconda/bin/conda run -n CBIG_py3 '];
            cmd = [conda_run script ' -s ' testCase.SubjectID ' -d ' output_dir ...
                ' -bld ''' testCase.Run ''' -BOLD_stem _rest_skip8_stc_mc_sdc' ...
                ' -echo_number 3 -echo_time 12,30.11,48.22'];
            disp(cmd);
            [status, ~] = system(cmd);
            assert(status == 0, 'ME denoise script failed');

            file_name = [testCase.SubjectID '_bld' testCase.Run '_rest_skip8_stc_mc_sdc_me.nii.gz'];
            mask_idx = loadBrainMask(testCase);
            disp('ME denoise');
            stats = CBIG_compare_volumes(fullfile(bold_dir, file_name), ...
                fullfile(ref_bold_dir, file_name), fullfile(output_dir, 'compare_output'), mask_idx);
            assert(stats.max_diff < 1e-6, 'ME denoise max diff = %e (threshold 1e-6)', stats.max_diff);

            rmdir(output_dir, 's');
        end

        function test_bbregister(testCase)
            % Input:  ME combined volume from data/sub005/bold/001/ (already masked by reference run)
            % Script: CBIG_preproc_bbregister.csh -s sub005 -d <output_dir>
            %         -anat_s sub005_FS -anat_d <testdata>/data
            %         -bld '001' -BOLD_stem _rest_skip8_stc_mc_sdc_me
            % Ref:    data/sub005/bold/001/sub005_bld001_rest_skip8_stc_mc_sdc_me.nii.gz (masked)
            %         data/sub005/bold/001/sub005_bld001_rest_skip8_stc_mc_sdc_me_reg.dat.mincost

            ref_bold_dir = fullfile(testCase.TestDataDir, 'data', testCase.SubjectID, ...
                'bold', testCase.Run);
            output_dir = fullfile(testCase.UnitTestDir, 'output', 'stepwise_bbregister');
            script = fullfile(testCase.CodeDir, 'CBIG_preproc_bbregister.csh');

            if exist(output_dir, 'dir'); rmdir(output_dir, 's'); end
            mkdir(output_dir);

            fprintf('\n');
            disp('========== CBIG_preproc_bbregister ==========');

            bold_dir = fullfile(output_dir, testCase.SubjectID, 'bold', testCase.Run);
            log_dir  = fullfile(output_dir, testCase.SubjectID, 'logs');
            qc_dir   = fullfile(output_dir, testCase.SubjectID, 'qc');
            mkdir(bold_dir); mkdir(log_dir); mkdir(qc_dir);

            file_name = [testCase.SubjectID '_bld' testCase.Run '_rest_skip8_stc_mc_sdc_me.nii.gz'];
            copyfile(fullfile(ref_bold_dir, file_name), fullfile(bold_dir, file_name));

            % FreeSurfer recon-all data lives in the test data dir (read-only, no copy needed)
            anat_dir = fullfile(testCase.TestDataDir, 'data');
            anat_s   = 'sub005_FS';

            cmd = [script ' -s ' testCase.SubjectID ' -d ' output_dir ...
                ' -anat_s ' anat_s ' -anat_d ' anat_dir ...
                ' -bld ''' testCase.Run ''' -BOLD_stem _rest_skip8_stc_mc_sdc_me'];
            disp(cmd);
            [status, ~] = system(cmd);
            assert(status == 0, 'BBregister script failed');

            % Compare registration matrix (.dat, 4x4 matrix starting at line 5)
            dat_file = [testCase.SubjectID '_bld' testCase.Run '_rest_skip8_stc_mc_sdc_me_reg.dat'];
            out_M = testCase.readRegDat(fullfile(bold_dir, dat_file));
            ref_M = testCase.readRegDat(fullfile(ref_bold_dir, dat_file));
            dat_diff = max(abs(out_M(:) - ref_M(:)));
            disp(['BBregister reg.dat max diff: ' num2str(dat_diff)]);
            assert(dat_diff < 1e-4, 'BBregister reg.dat max diff = %e (threshold 1e-4)', dat_diff);

            % Compare registration cost (single float in .dat.mincost)
            cost_file = [testCase.SubjectID '_bld' testCase.Run '_rest_skip8_stc_mc_sdc_me_reg.dat.mincost'];
            out_cost = load(fullfile(bold_dir, cost_file));
            ref_cost = load(fullfile(ref_bold_dir, cost_file));
            cost_diff = abs(out_cost(1) - ref_cost(1));
            disp(['BBregister mincost diff: ' num2str(cost_diff)]);
            assert(cost_diff < 1e-4, 'BBregister mincost diff = %e (threshold 1e-4)', cost_diff);

            rmdir(output_dir, 's');
        end

        function test_regress(testCase)
            % Input:  brain-masked ME volume + reg.dat + mc.par + motion outlier file
            % Script: CBIG_preproc_regression.csh -s sub005 -d <output_dir>
            %         -anat_s sub005_FS -anat_d <testdata>/data -bld '001'
            %         -BOLD_stem _rest_skip8_stc_mc_sdc_me
            %         -REG_stem _rest_skip8_stc_mc_sdc_me_reg
            %         -MASK_stem _rest_skip8_stc_mc_sdc_me
            %         -OUTLIER_stem _FDRMS0.3_DVARS60_motion_outliers.txt
            %         -whole_brain -wm -csf -motion12_itamar -detrend_method detrend
            %         -per_run -censor -polynomial_fit 1
            % Ref:    data/sub005/bold/001/sub005_bld001_rest_skip8_stc_mc_sdc_me_residc.nii.gz

            ref_bold_dir = fullfile(testCase.TestDataDir, 'data', testCase.SubjectID, ...
                'bold', testCase.Run);
            ref_qc_dir = fullfile(testCase.TestDataDir, 'data', testCase.SubjectID, 'qc');
            output_dir = fullfile(testCase.UnitTestDir, 'output', 'stepwise_regress');
            script = fullfile(testCase.CodeDir, 'CBIG_preproc_regression.csh');

            if exist(output_dir, 'dir'); rmdir(output_dir, 's'); end
            mkdir(output_dir);

            fprintf('\n');
            disp('========== CBIG_preproc_regression ==========');

            bold_dir = fullfile(output_dir, testCase.SubjectID, 'bold', testCase.Run);
            log_dir  = fullfile(output_dir, testCase.SubjectID, 'logs');
            qc_dir   = fullfile(output_dir, testCase.SubjectID, 'qc');
            mkdir(bold_dir); mkdir(log_dir); mkdir(qc_dir);

            % fMRI volume (brain-masked by bbregister step)
            me_file = [testCase.SubjectID '_bld' testCase.Run '_rest_skip8_stc_mc_sdc_me.nii.gz'];
            copyfile(fullfile(ref_bold_dir, me_file), fullfile(bold_dir, me_file));

            % Registration file (needed by create_mask for WM/CSF/WB masks)
            reg_file = [testCase.SubjectID '_bld' testCase.Run '_rest_skip8_stc_mc_sdc_me_reg.dat'];
            copyfile(fullfile(ref_bold_dir, reg_file), fullfile(bold_dir, reg_file));

            % Motion parameters (script finds via glob *_mc.par)
            par_file = [testCase.SubjectID '_bld' testCase.Run '_e1_rest_skip8_stc_mc.par'];
            copyfile(fullfile(ref_bold_dir, par_file), fullfile(bold_dir, par_file));

            % Motion outlier file (for -censor flag)
            outlier_file = [testCase.SubjectID '_bld' testCase.Run '_FDRMS0.3_DVARS60_motion_outliers.txt'];
            copyfile(fullfile(ref_qc_dir, outlier_file), fullfile(qc_dir, outlier_file));

            anat_dir = fullfile(testCase.TestDataDir, 'data');
            anat_s   = 'sub005_FS';
            outlier_stem = '_FDRMS0.3_DVARS60_motion_outliers.txt';

            cmd = [script ' -s ' testCase.SubjectID ' -d ' output_dir ...
                ' -anat_s ' anat_s ' -anat_d ' anat_dir ...
                ' -bld ''' testCase.Run ''' -BOLD_stem _rest_skip8_stc_mc_sdc_me' ...
                ' -REG_stem _rest_skip8_stc_mc_sdc_me_reg' ...
                ' -MASK_stem _rest_skip8_stc_mc_sdc_me' ...
                ' -OUTLIER_stem ' outlier_stem ...
                ' -whole_brain -wm -csf -motion12_itamar -detrend_method detrend' ...
                ' -per_run -censor -polynomial_fit 1'];
            disp(cmd);
            [status, ~] = system(cmd);
            assert(status == 0, 'Regress script failed');

            file_name = [testCase.SubjectID '_bld' testCase.Run '_rest_skip8_stc_mc_sdc_me_residc.nii.gz'];
            mask_idx = loadBrainMask(testCase);
            disp('Regress');
            stats = CBIG_compare_volumes(fullfile(bold_dir, file_name), ...
                fullfile(ref_bold_dir, file_name), fullfile(output_dir, 'compare_output'), mask_idx);
            assert(stats.max_diff < 1e-6, 'Regress max diff = %e (threshold 1e-6)', stats.max_diff);

            % Compare regressor files
            ref_reg_dir = fullfile(testCase.TestDataDir, 'data', testCase.SubjectID, 'bold', 'regression');
            out_reg_dir = fullfile(output_dir, testCase.SubjectID, 'bold', 'regression');
            reg_files = { ...
                [testCase.SubjectID '_bld' testCase.Run '_mc_regressor.txt'], ...
                [testCase.SubjectID '_bld' testCase.Run '_WB_WM_CSF_regressors.txt'], ...
                [testCase.SubjectID '_bld' testCase.Run '_all_regressors.txt'] };
            for i = 1:length(reg_files)
                out_r = load(fullfile(out_reg_dir, reg_files{i}));
                ref_r = load(fullfile(ref_reg_dir, reg_files{i}));
                reg_diff = max(abs(out_r(:) - ref_r(:)));
                disp([reg_files{i} ' max diff: ' num2str(reg_diff)]);
                assert(reg_diff < 1e-6, '%s max diff = %e (threshold 1e-6)', reg_files{i}, reg_diff);
            end

            rmdir(output_dir, 's');
        end

        function test_censor(testCase)
            % Input:  residualised volume + reg.dat + motion outlier file
            % Script: CBIG_preproc_censor.csh -s sub005 -d <output_dir>
            %         -anat_s sub005_FS -anat_d <testdata>/data -bld '001'
            %         -BOLD_stem _rest_skip8_stc_mc_sdc_me_residc
            %         -REG_stem _rest_skip8_stc_mc_sdc_me_reg
            %         -OUTLIER_stem _FDRMS0.3_DVARS60_motion_outliers.txt -FD_metric FDRMS
            % Ref:    data/sub005/bold/001/<subject>_bld001_rest_..._residc_interp_FDRMS0.3_DVARS60.nii.gz
            %         data/sub005/bold/001/<subject>_bld001_rest_..._residc_interp_inter_FDRMS0.3_DVARS60.nii.gz

            ref_bold_dir = fullfile(testCase.TestDataDir, 'data', testCase.SubjectID, ...
                'bold', testCase.Run);
            ref_qc_dir = fullfile(testCase.TestDataDir, 'data', testCase.SubjectID, 'qc');
            output_dir = fullfile(testCase.UnitTestDir, 'output', 'stepwise_censor');
            script = fullfile(testCase.CodeDir, 'CBIG_preproc_censor.csh');

            if exist(output_dir, 'dir'); rmdir(output_dir, 's'); end
            mkdir(output_dir);

            fprintf('\n');
            disp('========== CBIG_preproc_censor ==========');

            bold_dir = fullfile(output_dir, testCase.SubjectID, 'bold', testCase.Run);
            log_dir  = fullfile(output_dir, testCase.SubjectID, 'logs');
            qc_dir   = fullfile(output_dir, testCase.SubjectID, 'qc');
            mkdir(bold_dir); mkdir(log_dir); mkdir(qc_dir);

            % Residualised fMRI volume
            residc_file = [testCase.SubjectID '_bld' testCase.Run '_rest_skip8_stc_mc_sdc_me_residc.nii.gz'];
            copyfile(fullfile(ref_bold_dir, residc_file), fullfile(bold_dir, residc_file));

            % Registration file (for brain/GM mask creation inside censor script)
            reg_file = [testCase.SubjectID '_bld' testCase.Run '_rest_skip8_stc_mc_sdc_me_reg.dat'];
            copyfile(fullfile(ref_bold_dir, reg_file), fullfile(bold_dir, reg_file));

            % Motion outlier file (accessed as ../../qc/<file> from bold/001/)
            outlier_file = [testCase.SubjectID '_bld' testCase.Run '_FDRMS0.3_DVARS60_motion_outliers.txt'];
            copyfile(fullfile(ref_qc_dir, outlier_file), fullfile(qc_dir, outlier_file));

            anat_dir     = fullfile(testCase.TestDataDir, 'data');
            anat_s       = 'sub005_FS';
            outlier_stem = '_FDRMS0.3_DVARS60_motion_outliers.txt';

            cmd = [script ' -s ' testCase.SubjectID ' -d ' output_dir ...
                ' -anat_s ' anat_s ' -anat_d ' anat_dir ...
                ' -bld ''' testCase.Run ''' -BOLD_stem _rest_skip8_stc_mc_sdc_me_residc' ...
                ' -REG_stem _rest_skip8_stc_mc_sdc_me_reg' ...
                ' -OUTLIER_stem ' outlier_stem ' -FD_metric FDRMS -nocleanup'];
            disp(cmd);
            [status, ~] = system(cmd);
            assert(status == 0, 'Censor script failed');

            interp_stem = [testCase.SubjectID '_bld' testCase.Run ...
                '_rest_skip8_stc_mc_sdc_me_residc_interp'];
            out_files = { ...
                [interp_stem '_FDRMS0.3_DVARS60.nii.gz'], ...
                [interp_stem '_inter_FDRMS0.3_DVARS60.nii.gz'] };
            labels = {'interp', 'interp_inter'};

            mask_idx = loadBrainMask(testCase);
            for i = 1:length(out_files)
                disp(['Censor ' labels{i}]);
                stats = CBIG_compare_volumes(fullfile(bold_dir, out_files{i}), ...
                    fullfile(ref_bold_dir, out_files{i}), fullfile(output_dir, 'compare_output'), mask_idx);
                assert(stats.max_diff < 1e-6, 'Censor %s max diff = %e (threshold 1e-6)', labels{i}, stats.max_diff);
            end

            rmdir(output_dir, 's');
        end

        function test_bandpass(testCase)
            % Input:  censor-interpolated volume + motion outlier file
            % Script: CBIG_preproc_bandpass_fft.csh -s sub005 -d <output_dir> -bld '001'
            %         -BOLD_stem _rest_skip8_stc_mc_sdc_me_residc_interp_FDRMS0.3_DVARS60
            %         -OUTLIER_stem _FDRMS0.3_DVARS60_motion_outliers.txt
            %         -low_f 0.009 -high_f 0.08 -detrend
            % Ref:    data/sub005/bold/001/<subject>_bld001_..._interp_FDRMS0.3_DVARS60_bp_0.009_0.08.nii.gz

            ref_bold_dir = fullfile(testCase.TestDataDir, 'data', testCase.SubjectID, ...
                'bold', testCase.Run);
            ref_qc_dir = fullfile(testCase.TestDataDir, 'data', testCase.SubjectID, 'qc');
            output_dir = fullfile(testCase.UnitTestDir, 'output', 'stepwise_bandpass');
            script = fullfile(testCase.CodeDir, 'CBIG_preproc_bandpass_fft.csh');

            if exist(output_dir, 'dir'); rmdir(output_dir, 's'); end
            mkdir(output_dir);

            fprintf('\n');
            disp('========== CBIG_preproc_bandpass_fft ==========');

            bold_dir = fullfile(output_dir, testCase.SubjectID, 'bold', testCase.Run);
            log_dir  = fullfile(output_dir, testCase.SubjectID, 'logs');
            qc_dir   = fullfile(output_dir, testCase.SubjectID, 'qc');
            mkdir(bold_dir); mkdir(log_dir); mkdir(qc_dir);

            input_file = [testCase.SubjectID '_bld' testCase.Run ...
                '_rest_skip8_stc_mc_sdc_me_residc_interp_FDRMS0.3_DVARS60.nii.gz'];
            copyfile(fullfile(ref_bold_dir, input_file), fullfile(bold_dir, input_file));

            outlier_file = [testCase.SubjectID '_bld' testCase.Run '_FDRMS0.3_DVARS60_motion_outliers.txt'];
            copyfile(fullfile(ref_qc_dir, outlier_file), fullfile(qc_dir, outlier_file));

            bold_stem    = '_rest_skip8_stc_mc_sdc_me_residc_interp_FDRMS0.3_DVARS60';
            outlier_stem = '_FDRMS0.3_DVARS60_motion_outliers.txt';

            cmd = [script ' -s ' testCase.SubjectID ' -d ' output_dir ...
                ' -bld ''' testCase.Run ''' -BOLD_stem ' bold_stem ...
                ' -OUTLIER_stem ' outlier_stem ' -low_f 0.009 -high_f 0.08 -detrend'];
            disp(cmd);
            [status, ~] = system(cmd);
            assert(status == 0, 'Bandpass script failed');

            file_name = [testCase.SubjectID '_bld' testCase.Run ...
                '_rest_skip8_stc_mc_sdc_me_residc_interp_FDRMS0.3_DVARS60_bp_0.009_0.08.nii.gz'];
            mask_idx = loadBrainMask(testCase);
            disp('Bandpass');
            stats = CBIG_compare_volumes(fullfile(bold_dir, file_name), ...
                fullfile(ref_bold_dir, file_name), fullfile(output_dir, 'compare_output'), mask_idx);
            assert(stats.max_diff < 1e-6, 'Bandpass max diff = %e (threshold 1e-6)', stats.max_diff);

            rmdir(output_dir, 's');
        end

        function test_native2fsaverage(testCase)
            % Input:  bandpass-filtered volume + reg.dat
            % Script: CBIG_preproc_native2fsaverage.csh -s sub005 -d <output_dir>
            %         -anat_s sub005_FS -anat_d <testdata>/data -bld '001'
            %         -BOLD_stem _rest_skip8_stc_mc_sdc_me_residc_interp_FDRMS0.3_DVARS60_bp_0.009_0.08
            %         -REG_stem _rest_skip8_stc_mc_sdc_me_reg
            %         -proj fsaverage6 -sm 6 -down fsaverage5
            % Ref:    data/sub005/surf/{lh,rh}.<subject>_bld001_..._fs6{,_sm6,_sm6_fs5}.nii.gz

            ref_bold_dir = fullfile(testCase.TestDataDir, 'data', testCase.SubjectID, ...
                'bold', testCase.Run);
            ref_surf_dir = fullfile(testCase.TestDataDir, 'data', testCase.SubjectID, 'surf');
            output_dir   = fullfile(testCase.UnitTestDir, 'output', 'stepwise_native2fsaverage');
            script = fullfile(testCase.CodeDir, 'CBIG_preproc_native2fsaverage.csh');

            if exist(output_dir, 'dir'); rmdir(output_dir, 's'); end
            mkdir(output_dir);

            fprintf('\n');
            disp('========== CBIG_preproc_native2fsaverage ==========');

            bold_dir = fullfile(output_dir, testCase.SubjectID, 'bold', testCase.Run);
            log_dir  = fullfile(output_dir, testCase.SubjectID, 'logs');
            mkdir(bold_dir); mkdir(log_dir);

            bold_stem = '_rest_skip8_stc_mc_sdc_me_residc_interp_FDRMS0.3_DVARS60_bp_0.009_0.08';
            bp_file  = [testCase.SubjectID '_bld' testCase.Run bold_stem '.nii.gz'];
            copyfile(fullfile(ref_bold_dir, bp_file), fullfile(bold_dir, bp_file));

            reg_file = [testCase.SubjectID '_bld' testCase.Run '_rest_skip8_stc_mc_sdc_me_reg.dat'];
            copyfile(fullfile(ref_bold_dir, reg_file), fullfile(bold_dir, reg_file));

            anat_dir = fullfile(testCase.TestDataDir, 'data');
            anat_s   = 'sub005_FS';

            cmd = [script ' -s ' testCase.SubjectID ' -d ' output_dir ...
                ' -anat_s ' anat_s ' -anat_d ' anat_dir ...
                ' -bld ''' testCase.Run ''' -BOLD_stem ' bold_stem ...
                ' -REG_stem _rest_skip8_stc_mc_sdc_me_reg' ...
                ' -proj fsaverage6 -sm 6 -down fsaverage5'];
            disp(cmd);
            [status, ~] = system(cmd);
            assert(status == 0, 'native2fsaverage script failed');

            out_surf_dir = fullfile(output_dir, testCase.SubjectID, 'surf');
            base = [testCase.SubjectID '_bld' testCase.Run bold_stem];
            disp('_fs6_sm6_fs5');
            lh_out = fullfile(out_surf_dir, ['lh.' base '_fs6_sm6_fs5.nii.gz']);
            rh_out = fullfile(out_surf_dir, ['rh.' base '_fs6_sm6_fs5.nii.gz']);
            lh_ref = fullfile(ref_surf_dir,  ['lh.' base '_fs6_sm6_fs5.nii.gz']);
            rh_ref = fullfile(ref_surf_dir,  ['rh.' base '_fs6_sm6_fs5.nii.gz']);
            stats = CBIG_compare_surfaces(lh_out, rh_out, lh_ref, rh_ref, ...
                fullfile(output_dir, 'compare_output'));
            assert(stats.max_diff < 1e-6, 'native2fsaverage max diff = %e (threshold 1e-6)', stats.max_diff);

            rmdir(output_dir, 's');
        end

        function test_fc_metrics(testCase)
            % Input:  surface sm6 files (lh/rh) + bp volume + func.aseg mask + outlier file
            % Script: CBIG_preproc_FCmetrics_wrapper.csh -s sub005 -d <output_dir> -bld '001'
            %         -BOLD_stem _rest_skip8_stc_mc_sdc_me_residc_interp_FDRMS0.3_DVARS60_bp_0.009_0.08
            %         -SURF_stem _rest_skip8_stc_mc_sdc_me_residc_interp_FDRMS0.3_DVARS60_bp_0.009_0.08_fs6_sm6
            %         -OUTLIER_stem _FDRMS0.3_DVARS60_motion_outliers.txt -Pearson_r -censor -nocleanup
            % Ref:    data/sub005/FC_metrics/Pearson_r/sub005_rest_..._fs6_sm6_<pair>_Yan_400Parcels_....mat

            ref_bold_dir = fullfile(testCase.TestDataDir, 'data', testCase.SubjectID, ...
                'bold', testCase.Run);
            ref_surf_dir = fullfile(testCase.TestDataDir, 'data', testCase.SubjectID, 'surf');
            ref_mask_dir = fullfile(testCase.TestDataDir, 'data', testCase.SubjectID, 'bold', 'mask');
            ref_qc_dir   = fullfile(testCase.TestDataDir, 'data', testCase.SubjectID, 'qc');
            ref_fc_dir   = fullfile(testCase.TestDataDir, 'data', testCase.SubjectID, ...
                'FC_metrics', 'Pearson_r');
            output_dir   = fullfile(testCase.UnitTestDir, 'output', 'stepwise_fc_metrics');
            script = fullfile(testCase.CodeDir, 'CBIG_preproc_FCmetrics_wrapper.csh');

            if exist(output_dir, 'dir'); rmdir(output_dir, 's'); end
            mkdir(output_dir);

            fprintf('\n');
            disp('========== CBIG_preproc_FCmetrics_wrapper ==========');

            bold_dir = fullfile(output_dir, testCase.SubjectID, 'bold', testCase.Run);
            mask_dir = fullfile(output_dir, testCase.SubjectID, 'bold', 'mask');
            surf_dir = fullfile(output_dir, testCase.SubjectID, 'surf');
            log_dir  = fullfile(output_dir, testCase.SubjectID, 'logs');
            qc_dir   = fullfile(output_dir, testCase.SubjectID, 'qc');
            mkdir(bold_dir); mkdir(mask_dir); mkdir(surf_dir); mkdir(log_dir); mkdir(qc_dir);

            bold_stem = '_rest_skip8_stc_mc_sdc_me_residc_interp_FDRMS0.3_DVARS60_bp_0.009_0.08';
            surf_stem = [bold_stem '_fs6_sm6'];

            % Bandpass volume
            bp_fname = [testCase.SubjectID '_bld' testCase.Run bold_stem '.nii.gz'];
            copyfile(fullfile(ref_bold_dir, bp_fname), fullfile(bold_dir, bp_fname));

            % Surface files (lh and rh, fs6_sm6)
            hemis = {'lh', 'rh'};
            for h = 1:length(hemis)
                sf = [hemis{h} '.' testCase.SubjectID '_bld' testCase.Run surf_stem '.nii.gz'];
                copyfile(fullfile(ref_surf_dir, sf), fullfile(surf_dir, sf));
            end

            % Subcortical aseg mask (for 19 subcortical label extraction)
            copyfile(fullfile(ref_mask_dir, [testCase.SubjectID '.func.aseg.nii']), ...
                fullfile(mask_dir, [testCase.SubjectID '.func.aseg.nii']));

            % Motion outlier file (for -censor)
            outlier_stem  = '_FDRMS0.3_DVARS60_motion_outliers.txt';
            outlier_file = [testCase.SubjectID '_bld' testCase.Run outlier_stem];
            copyfile(fullfile(ref_qc_dir, outlier_file), fullfile(qc_dir, outlier_file));

            cmd = [script ' -s ' testCase.SubjectID ' -d ' output_dir ...
                ' -bld ''' testCase.Run ''' -BOLD_stem ' bold_stem ...
                ' -SURF_stem ' surf_stem ...
                ' -OUTLIER_stem ' outlier_stem ' -Pearson_r -censor -nocleanup'];
            disp(cmd);
            [status, ~] = system(cmd);
            assert(status == 0, 'FC metrics script failed');

            out_fc_dir = fullfile(output_dir, testCase.SubjectID, 'FC_metrics', 'Pearson_r');
            parcel_suffix = '_Yan_400Parcels_Kong17Networks_order_with_19Subcortical.mat';
            fc_stem = [testCase.SubjectID surf_stem];

            file_name = [fc_stem '_all2all' parcel_suffix];
            out = load(fullfile(out_fc_dir, file_name));
            ref = load(fullfile(ref_fc_dir, file_name));
            diff_mat = abs(out.corr_mat(:) - ref.corr_mat(:));
            disp('FC all2all');
            disp(['  max diff: '    num2str(max(diff_mat)) ...
                  '  min diff: '    num2str(min(diff_mat)) ...
                  '  median diff: ' num2str(median(diff_mat)) ...
                  '  mean diff: '   num2str(mean(diff_mat))]);
            assert(max(diff_mat) < 1e-6, 'FC all2all max diff = %e (threshold 1e-6)', max(diff_mat));

            rmdir(output_dir, 's');
        end

        function test_native2mni(testCase)
            % Input:  bp volume + reg.dat (bbregister) + pre-computed ANTs vol files
            % Script: CBIG_preproc_native2mni_ants.csh -s sub005 -d <output_dir>
            %         -anat_s sub005_FS -anat_d <testdata>/data -bld '001'
            %         -BOLD_stem _rest_skip8_stc_mc_sdc_me_residc_interp_FDRMS0.3_DVARS60_bp_0.009_0.08
            %         -REG_stem _rest_skip8_stc_mc_sdc_me_reg
            %         -final_mask ${FSL_DIR}/data/standard/MNI152_T1_2mm_brain_mask_dil.nii.gz
            %         -iter_affine 5x5x2 -iter_SyN 5x5x2
            % Note:   Providing pre-computed norm_MNI1mm warp+affine in vol/ skips ANTs registration.
            %         Default sm=6, so output is _MNI2mm_sm6.nii.gz
            % Ref:    data/sub005/vol/sub005_bld001_..._MNI2mm_sm6{,_finalmask}.nii.gz

            bold_stem = '_rest_skip8_stc_mc_sdc_me_residc_interp_FDRMS0.3_DVARS60_bp_0.009_0.08';

            ref_bold_dir = fullfile(testCase.TestDataDir, 'data', testCase.SubjectID, ...
                'bold', testCase.Run);
            ref_vol_dir  = fullfile(testCase.TestDataDir, 'data', testCase.SubjectID, 'vol');
            output_dir   = fullfile(testCase.UnitTestDir, 'output', 'stepwise_native2mni');
            script = fullfile(testCase.CodeDir, 'CBIG_preproc_native2mni_ants.csh');

            if exist(output_dir, 'dir'); rmdir(output_dir, 's'); end
            mkdir(output_dir);

            fprintf('\n');
            disp('========== CBIG_preproc_native2mni_ants ==========');

            bold_dir = fullfile(output_dir, testCase.SubjectID, 'bold', testCase.Run);
            vol_dir  = fullfile(output_dir, testCase.SubjectID, 'vol');
            log_dir  = fullfile(output_dir, testCase.SubjectID, 'logs');
            mkdir(bold_dir); mkdir(vol_dir); mkdir(log_dir);

            % Bandpass volume (input BOLD)
            bp_file = [testCase.SubjectID '_bld' testCase.Run bold_stem '.nii.gz'];
            copyfile(fullfile(ref_bold_dir, bp_file), fullfile(bold_dir, bp_file));

            % BBregister .dat file (used by tkregister2 to convert to ITK format)
            reg_file = [testCase.SubjectID '_bld' testCase.Run '_rest_skip8_stc_mc_sdc_me_reg.dat'];
            copyfile(fullfile(ref_bold_dir, reg_file), fullfile(bold_dir, reg_file));

            % Pre-computed ANTs vol files (warp+affine presence skips re-registration)
            for file_name = {'norm.nii.gz', 'norm_MNI1mm0GenericAffine.mat', ...
                         'norm_MNI1mm1Warp.nii.gz', 'norm_MNI1mm1InverseWarp.nii.gz'}
                copyfile(fullfile(ref_vol_dir, file_name{1}), fullfile(vol_dir, file_name{1}));
            end

            anat_dir = fullfile(testCase.TestDataDir, 'data');
            anat_s   = 'sub005_FS';
            fsl_dir  = getenv('FSL_DIR');
            final_mask = fullfile(fsl_dir, 'data', 'standard', ...
                'MNI152_T1_2mm_brain_mask_dil.nii.gz');

            cmd = [script ' -s ' testCase.SubjectID ' -d ' output_dir ...
                ' -anat_s ' anat_s ' -anat_d ' anat_dir ...
                ' -bld ''' testCase.Run ''' -BOLD_stem ' bold_stem ...
                ' -REG_stem _rest_skip8_stc_mc_sdc_me_reg' ...
                ' -final_mask ' final_mask ...
                ' -iter_affine 5x5x2 -iter_SyN 5x5x2'];
            disp(cmd);
            [status, ~] = system(cmd);
            assert(status == 0, 'native2mni script failed');

            % When -final_mask is used with nocleanup=0 (default), the script deletes
            % _MNI2mm_sm6.nii.gz after masking; only _MNI2mm_sm6_finalmask.nii.gz survives.
            out_vol_dir = fullfile(output_dir, testCase.SubjectID, 'vol');
            out_stem = [testCase.SubjectID '_bld' testCase.Run bold_stem];
            file_name = [out_stem '_MNI2mm_sm6_finalmask.nii.gz'];

            disp('native2mni MNI2mm_sm6_finalmask');
            stats = CBIG_compare_volumes(fullfile(out_vol_dir, file_name), ...
                fullfile(ref_vol_dir, file_name), fullfile(output_dir, 'compare_output'));
            assert(stats.max_diff < 1e-6, '%s max diff = %e (threshold 1e-6)', file_name, stats.max_diff);

            rmdir(output_dir, 's');
        end

        function test_mc_notch_filter(testCase)
            % Test notch-filter MC for sub-NDARBF851NH6 (single-echo, 2 runs).
            % Input:  _rest_skip8_stc.nii.gz per run
            % Script: CBIG_preproc_fslmcflirt_outliers.csh with -low_f 0.31 -high_f 0.43 -spline_final
            % Ref:    test_motion_filter/sub-NDARBF851NH6/bold/{001,002}/<sub>_bld<run>_rest_skip8_stc_mc.par

            sub_id = 'sub-NDARBF851NH6';
            runs   = {'001', '002'};
            ref_base_dir = fullfile(testCase.TestDataDir, 'test_motion_filter', sub_id);
            output_dir = fullfile(testCase.UnitTestDir, 'output', 'stepwise_mc_notch');
            script = fullfile(testCase.CodeDir, 'CBIG_preproc_fslmcflirt_outliers.csh');

            if exist(output_dir, 'dir'); rmdir(output_dir, 's'); end
            mkdir(output_dir);

            fprintf('\n');
            disp('========== CBIG_preproc_fslmcflirt_outliers (notch filter) ==========');

            log_dir = fullfile(output_dir, sub_id, 'logs');
            mkdir(log_dir);

            % Write .bold run-list file (needed by script to apply transforms across runs)
            fid = fopen(fullfile(log_dir, [sub_id '.bold']), 'w');
            for i = 1:length(runs)
                fprintf(fid, '%s\n', runs{i});
            end
            fclose(fid);

            % Copy input STC volumes for each run
            for i = 1:length(runs)
                run = runs{i};
                bold_dir = fullfile(output_dir, sub_id, 'bold', run);
                mkdir(bold_dir);
                in_file = [sub_id '_bld' run '_rest_skip8_stc.nii.gz'];
                copyfile(fullfile(ref_base_dir, 'bold', run, in_file), ...
                         fullfile(bold_dir, in_file));
            end

            run_str = strjoin(runs, ' ');
            cmd = [script ' -s ' sub_id ' -d ' output_dir ...
                ' -bld ''' run_str '''' ...
                ' -BOLD_stem _rest_skip8_stc' ...
                ' -FD_th 0.3 -DV_th 60 -discard-run 50 -rm-seg 5 -spline_final' ...
                ' -low_f 0.31 -high_f 0.43'];
            disp(cmd);
            [status, ~] = system(cmd);
            assert(status == 0, 'MC notch-filter script failed');

            for i = 1:length(runs)
                run = runs{i};
                par_name = [sub_id '_bld' run '_rest_skip8_stc_mc.par'];
                out_par = load(fullfile(output_dir, sub_id, 'bold', run, par_name));
                ref_par = load(fullfile(ref_base_dir, 'bold', run, par_name));
                par_diff = max(abs(out_par(:) - ref_par(:)));
                disp(['MC notch par run ' run ' max diff: ' num2str(par_diff)]);
                assert(par_diff < 1e-6, ...
                    'MC notch .par run %s max diff = %e (threshold 1e-6)', run, par_diff);
            end

            rmdir(output_dir, 's');
        end

        function test_regress_aCompCor(testCase)
            % Test aCompCor regression for Sub1116_Ses1 (single-echo, 2 runs).
            % Input per run: _rest_skip4_stc_mc.nii.gz, .par, _reg.dat; FDRMS/DVARS outliers from qc/
            % Script: CBIG_preproc_regression.csh with -aCompCor -aCompCor_nPCs 5 -wm -csf -motion12_itamar -per_run
            % Ref: data/Sub1116_Ses1/bold/{002,003}/<sub>_bld<run>_rest_skip4_stc_mc_resid.nii.gz

            sub_id   = 'Sub1116_Ses1';
            runs     = {'002', '003'};
            anat_dir = fullfile(testCase.TestDataDir, 'data');
            ref_base_dir = fullfile(testCase.TestDataDir, 'data', sub_id);
            output_dir = fullfile(testCase.UnitTestDir, 'output', 'stepwise_regress_aCompCor');
            script = fullfile(testCase.CodeDir, 'CBIG_preproc_regression.csh');

            if exist(output_dir, 'dir'); rmdir(output_dir, 's'); end
            mkdir(output_dir);

            fprintf('\n');
            disp('========== CBIG_preproc_regression (aCompCor) ==========');

            log_dir = fullfile(output_dir, sub_id, 'logs');
            qc_dir  = fullfile(output_dir, sub_id, 'qc');
            mkdir(log_dir); mkdir(qc_dir);

            bold_stem = '_rest_skip4_stc_mc';
            for i = 1:length(runs)
                run = runs{i};
                bold_dir    = fullfile(output_dir, sub_id, 'bold', run);
                ref_run_dir = fullfile(ref_base_dir, 'bold', run);
                mkdir(bold_dir);

                vol_file = [sub_id '_bld' run bold_stem '.nii.gz'];
                copyfile(fullfile(ref_run_dir, vol_file), fullfile(bold_dir, vol_file));

                par_file = [sub_id '_bld' run bold_stem '.par'];
                copyfile(fullfile(ref_run_dir, par_file), fullfile(bold_dir, par_file));

                reg_file = [sub_id '_bld' run bold_stem '_reg.dat'];
                copyfile(fullfile(ref_run_dir, reg_file), fullfile(bold_dir, reg_file));

                outlier_file = [sub_id '_bld' run '_FDRMS0.2_DVARS50_motion_outliers.txt'];
                copyfile(fullfile(ref_base_dir, 'qc', outlier_file), fullfile(qc_dir, outlier_file));
            end

            run_str = strjoin(runs, ' ');
            anat_s  = [sub_id '_FS'];

            cmd = [script ' -s ' sub_id ' -d ' output_dir ...
                ' -anat_s ' anat_s ' -anat_d ' anat_dir ...
                ' -bld ''' run_str '''' ...
                ' -BOLD_stem ' bold_stem ...
                ' -REG_stem ' bold_stem '_reg' ...
                ' -MASK_stem ' bold_stem ...
                ' -OUTLIER_stem _FDRMS0.2_DVARS50_motion_outliers.txt' ...
                ' -aCompCor -aCompCor_nPCs 5 -wm -wm_max_erode 3' ...
                ' -csf -csf_max_erode 1 -motion12_itamar -per_run'];
            disp(cmd);
            [status, ~] = system(cmd);
            assert(status == 0, 'Regress aCompCor script failed');

            % Load brainmask from reference data (this subject uses .loosebrainmask naming)
            mask_vol = MRIread(fullfile(ref_base_dir, 'bold', 'mask', ...
                [sub_id '.loosebrainmask.bin.nii.gz']));
            mask_idx = logical(mask_vol.vol(:));

            ref_reg_dir = fullfile(ref_base_dir, 'bold', 'regression');
            out_reg_dir = fullfile(output_dir, sub_id, 'bold', 'regression');
            reg_files_tmpl = { ...
                '_mc_regressor.txt', ...
                '_WM_CSF_regressors.txt', ...
                '_aCompCor_regressor.txt', ...
                '_all_regressors.txt' };

            for i = 1:length(runs)
                run = runs{i};
                out_run_dir = fullfile(output_dir, sub_id, 'bold', run);
                ref_run_dir = fullfile(ref_base_dir, 'bold', run);
                file_name = [sub_id '_bld' run bold_stem '_resid.nii.gz'];
                disp(['Regress aCompCor run ' run]);
                stats = CBIG_compare_volumes(fullfile(out_run_dir, file_name), ...
                    fullfile(ref_run_dir, file_name), fullfile(output_dir, 'compare_output'), mask_idx);
                assert(stats.max_diff < 1e-5, ...
                    'Regress aCompCor run %s max diff = %e (threshold 1e-5)', run, stats.max_diff);

                for j = 1:length(reg_files_tmpl)
                    reg_name = [sub_id '_bld' run reg_files_tmpl{j}];
                    out_r = load(fullfile(out_reg_dir, reg_name));
                    ref_r = load(fullfile(ref_reg_dir, reg_name));
                    reg_diff = max(abs(out_r(:) - ref_r(:)));
                    disp([reg_name ' max diff: ' num2str(reg_diff)]);
                    assert(reg_diff < 1e-5, '%s max diff = %e (threshold 1e-5)', reg_name, reg_diff);
                end
            end

            rmdir(output_dir, 's');
        end

        function test_sdc_oppo_ped(testCase)
            % Test opposite phase-encode SDC for sub-NDARBF851NH6 (single-echo, 2 runs).
            % Input per run: _rest_skip8_stc.nii.gz, _rest_skip8_stc_mc.cat
            %               AP/PA EPI from data/input_truncate/
            % Script: CBIG_preproc_spatial_distortion_correction.csh with -fpm oppo_PED
            % Ref: data/sub-NDARBF851NH6/bold/{001,002}/<sub>_bld<run>_rest_skip8_stc_mc_sdc.nii.gz

            sub_id = 'sub-NDARBF851NH6';
            runs   = {'001', '002'};
            ref_base_dir = fullfile(testCase.TestDataDir, 'data', sub_id);
            output_dir = fullfile(testCase.UnitTestDir, 'output', 'stepwise_sdc_oppo_ped');
            script = fullfile(testCase.CodeDir, 'CBIG_preproc_spatial_distortion_correction.csh');

            if exist(output_dir, 'dir'); rmdir(output_dir, 's'); end
            mkdir(output_dir);

            fprintf('\n');
            disp('========== CBIG_preproc_spatial_distortion_correction (oppo-PED) ==========');

            log_dir = fullfile(output_dir, sub_id, 'logs');
            mkdir(log_dir);

            for i = 1:length(runs)
                run = runs{i};
                bold_dir    = fullfile(output_dir, sub_id, 'bold', run);
                ref_run_dir = fullfile(ref_base_dir, 'bold', run);
                mkdir(bold_dir);

                stc_file = [sub_id '_bld' run '_rest_skip8_stc.nii.gz'];
                copyfile(fullfile(ref_run_dir, stc_file), fullfile(bold_dir, stc_file));

                cat_file = [sub_id '_bld' run '_rest_skip8_stc_mc.cat'];
                copyfile(fullfile(ref_run_dir, cat_file), fullfile(bold_dir, cat_file));
            end

            input_dir = fullfile(testCase.TestDataDir, 'data', 'input_truncate');
            ap_file = fullfile(input_dir, [sub_id '_dir-AP_acq-fMRI_epi.nii.gz']);
            pa_file = fullfile(input_dir, [sub_id '_dir-PA_acq-fMRI_epi.nii.gz']);
            run_str = strjoin(runs, ' ');

            cmd = [script ' -s ' sub_id ' -d ' output_dir ...
                ' -bld ''' run_str '''' ...
                ' -BOLD_stem _rest_skip8_stc' ...
                ' -fpm oppo_PED' ...
                ' -j_minus ' ap_file ' -j_plus ' pa_file ...
                ' -j_minus_trt 0.0414992 -j_plus_trt 0.0414992' ...
                ' -ees .589987000 -te 30.00'];
            disp(cmd);
            [status, ~] = system(cmd);
            assert(status == 0, 'SDC oppo-PED script failed');

            for i = 1:length(runs)
                run = runs{i};
                out_run_dir = fullfile(output_dir, sub_id, 'bold', run);
                ref_run_dir = fullfile(ref_base_dir, 'bold', run);
                file_name = [sub_id '_bld' run '_rest_skip8_stc_mc_sdc.nii.gz'];
                disp(['SDC oppo-PED run ' run]);
                stats = CBIG_compare_volumes(fullfile(out_run_dir, file_name), ...
                    fullfile(ref_run_dir, file_name), fullfile(output_dir, 'compare_output'));
                assert(stats.max_diff < 1e-6, ...
                    'SDC oppo-PED run %s max diff = %e (threshold 1e-6)', run, stats.max_diff);
            end

            rmdir(output_dir, 's');
        end

    end

    methods (Access = private)

        function mask_idx = loadBrainMask(testCase)
            % Return logical column vector of brain voxels (functional space).
            mask_dir = fullfile(testCase.TestDataDir, 'data', testCase.SubjectID, 'bold', 'mask');
            mask     = MRIread(fullfile(mask_dir, [testCase.SubjectID '.brainmask.bin.nii.gz']));
            mask_idx = logical(mask.vol(:));
        end

    end

    methods (Static, Access = private)

        function M = readRegDat(dat_file)
            % Read the 4x4 registration matrix from a FreeSurfer .dat file.
            % Format: 4 header lines, then 4 rows of the matrix, then 'round'.
            fid = fopen(dat_file, 'r');
            for k = 1:4; fgetl(fid); end  % skip subject, voxel sizes, intensity
            M = zeros(4, 4);
            for k = 1:4
                M(k, :) = sscanf(fgetl(fid), '%f', 4)';
            end
            fclose(fid);
        end

    end
end
