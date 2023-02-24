classdef CBIG_LiGSR_unit_test < matlab.unittest.TestCase
%
% Target project:
%                 Li2019_GSR
%
% Case design:
%                 Li2019_GSR already has unit tests, here we call its
%                 "intelligence_score" unit test and automatically make
%                 judgement on whether the unit test is passed or failed
%
%                 For now, the stable projects' unit tests in our repo
%                 require manual check of the output txt files or images,
%                 making it incovenient for wrapper function to call them
%                 and automatically draw conclusions. As a result, we write
%                 some simple matlab test functions for these stable
%                 projects' unit tests.
%
% Written by Jingwei Li and CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md

methods (Test)
    function test_examples(testCase)
        addpath(fullfile(getenv('CBIG_CODE_DIR'), 'stable_projects', 'preprocessing', ...
            'Li2019_GSR', 'examples', 'scripts'));
        replace_unit_test = load(fullfile(getenv('CBIG_CODE_DIR'), 'unit_tests', 'replace_unittest_flag'));
        
        % generate example results
        CBIG_LiGSR_generate_example_results
        
        % replace unit test if flag is 1
        if replace_unit_test
            % display differences
            disp("Replacing unit test results for Li2019_GSR, test_examples")
            % differences in variance component model
            exm_dir = fullfile(getenv('CBIG_CODE_DIR'), 'stable_projects', 'preprocessing', ...
            'Li2019_GSR', 'examples');
            outdir = fullfile(exm_dir, 'output');
            refdir = fullfile(exm_dir, 'ref_output');
            nbehaviors = 2;
            jk_seeds = 10;
            d = 2;
            fprintf('1. Variance Component model:\n')
            for seed = 1:jk_seeds
                fprintf('\tJackknife #%d ...\n', seed)
                for b = 1:nbehaviors
                    fprintf('\t\tBehavior %d ...\n', b)

                    outfile = fullfile(outdir, 'VarianceComponentModel', ['del' num2str(d) '_set' ...
                        num2str(seed)], ['m2_QuantileNorm_Behavior_' num2str(b) '.mat']);
                    reffile = fullfile(refdir, 'VarianceComponentModel', ['del' num2str(d) '_set' ...
                        num2str(seed)], ['m2_QuantileNorm_Behavior_' num2str(b) '.mat']);
                    out = load(outfile);
                    ref = load(reffile);

                    reffields = fieldnames(ref.morpho);
                    outfields = fieldnames(out.morpho);
                    if length(reffields) ~= length(outfields)
                        fprintf('\t\t\tOutput ''morpho'' structure has changed.\n')
                    end
        
                    for i = 1:length(reffields)
                        curr_outfield = getfield(out.morpho, outfields{i});
                        curr_reffield = getfield(ref.morpho, reffields{i});
                        if length(reffields) ~= length(outfields)
                            fprintf('\t\t\tstructure field %s size has changed.\n',reffields{i})
                        end
                        if all(isnan(curr_outfield)) 
                            if ~all(isnan(curr_outfield))
                                fprintf('\t\t\tstructure field %s is now NaN, but was not in old reference.\n', ...
                                    reffields{i})
                            end
                        else
                            maxdif = max(abs(curr_reffield(:) - curr_outfield(:)));
                            fprintf('\t\t\tstructure field %s differed by (max abs diff) %f\n.', ...
                                reffields{i}, maxdif)
                        end
                    end
                end
            end

            fprintf('\tFull set ...\n')
            for b = 1:nbehaviors
                fprintf('\t\tBehavior %d ...\n', b)

                outfile = fullfile(outdir, 'VarianceComponentModel', 'fullset', ...
                    ['m2_QuantileNorm_Behavior_' num2str(b) '.mat']);
                reffile = fullfile(refdir, 'VarianceComponentModel', 'fullset', ...
                    ['m2_QuantileNorm_Behavior_' num2str(b) '.mat']);
                out = load(outfile);
                ref = load(reffile);
    
                reffields = fieldnames(ref.morpho);
                outfields = fieldnames(out.morpho);
                if length(reffields) ~= length(outfields)
                        fprintf('\t\t\tOutput ''morpho'' structure has changed.\n')
                end

                for i = 1:length(reffields)
                        curr_outfield = getfield(out.morpho, outfields{i});
                        curr_reffield = getfield(ref.morpho, reffields{i});
                        if length(reffields) ~= length(outfields)
                            fprintf('\t\t\tstructure field %s size has changed.\n',reffields{i})
                        end
                        if all(isnan(curr_outfield)) 
                            if ~all(isnan(curr_outfield))
                                fprintf('\t\t\tstructure field %s is now NaN, but was not in old reference.\n', ...
                                    reffields{i})
                            end
                        else
                            maxdif = max(abs(curr_reffield(:) - curr_outfield(:)));
                            fprintf('\t\t\tstructure field %s differed by (max abs diff) %f\n.', ...
                                reffields{i}, maxdif)
                        end
                    end
            end

            % differences in kernel regression
            fprintf('2. Kernel regression:\n')
            outfile = fullfile(outdir, 'KernelRidgeRegression', 'final_result.mat');
            reffile = fullfile(refdir, 'KernelRidgeRegression', 'final_result.mat');
            out = load(outfile);
            ref = load(reffile);

            out = orderfields(out,ref);
            reffields = fieldnames(ref);
            outfields = fieldnames(out);
            if length(reffields) ~= length(outfields)
                        fprintf('Output structure has changed.\n')
            end

            for i = 1:length(reffields)
                curr_outfield = getfield(out, outfields{i});
                curr_reffield = getfield(ref, reffields{i});
                if ~isequal(size(curr_reffield), size(curr_outfield))
                    fprintf('\t\t\tvariable %s size has changed.\n',reffields{i})
                end
                if(strcmp(reffields{i}, 'optimal_stats') || strcmp(reffields{i}, 'optimal_kernel'))
                    subfields = fieldnames(curr_reffield);
                    for j = 1:length(subfields)
                        sub_reffield = getfield(curr_reffield, subfields{j});
                        sub_outfield = getfield(curr_outfield, subfields{j});

                        if all(reshape(isnan(sub_outfield), [numel(sub_outfield) 1]))
                            if ~all(reshape(isnan(sub_reffield), [numel(sub_reffield) 1]))
                                fprintf('\t\t\tstructure field %s %s is now NaN, but was not in old reference.\n', ...
                                    reffields{i}, subfields{j})
                            end
                        else
                            maxdif = max(abs(sub_reffield(:) - sub_outfield(:)));
                            fprintf('\t\t\tstructure field %s %s differed by (max abs diff) %f\n.', ...
                                reffields{i}, subfields{j}, maxdif)
                        end
                    end
                else
                    if all(reshape(isnan(curr_outfield), [numel(curr_outfield) 1]))
                            if ~all(reshape(isnan(curr_reffield), [numel(curr_reffield) 1]))
                                fprintf('\t\t\tstructure field %s is now NaN, but was not in old reference.\n', ...
                                    reffields{i})
                            end
                        else
                            maxdif = max(abs(curr_reffield(:) - curr_outfield(:)));
                            fprintf('\t\t\tstructure field %s differed by (max abs diff) %f\n.', ...
                                reffields{i}, maxdif)
                    end
                end
            end
            
            % copy and replace reference
            copyfile(fullfile(exm_dir, 'output'), fullfile(exm_dir, 'ref_output') )
        end
        
        % check results
        CBIG_LiGSR_check_example_results
        
        rmpath(fullfile(getenv('CBIG_CODE_DIR'), 'stable_projects', 'preprocessing', ...
            'Li2019_GSR', 'examples', 'scripts'));
    end
    
    function test_intelligence_score_LME_GSP_Case(testCase)
        %% path setting
        addpath(fullfile(getenv('CBIG_CODE_DIR'), 'stable_projects', 'preprocessing', ...
            'Li2019_GSR', 'unit_tests', 'intelligence_score', 'scripts'));
        OutputDir = fullfile(getenv('CBIG_CODE_DIR'), 'stable_projects', 'preprocessing', ...
            'Li2019_GSR', 'unit_tests', 'output', 'intelligence_score_LME_GSP_Case');
        replace_unit_test = load(fullfile(getenv('CBIG_CODE_DIR'), 'unit_tests', 'replace_unittest_flag'));
        
        % create output dir (IMPORTANT)
        if(exist(OutputDir, 'dir'))
            rmdir(OutputDir, 's')
        end
        mkdir(OutputDir);
        
        %% call Li2019_GSR intelligence_score shell script (variance component model & GSP data)
        cmdfun = fullfile(getenv('CBIG_CODE_DIR'), 'stable_projects', 'preprocessing', 'Li2019_GSR', ...
            'unit_tests', 'intelligence_score', 'scripts', 'CBIG_LiGSR_LME_unittest_intelligence_score_GSP.sh');
        cmd = [cmdfun, ' ', OutputDir];
        system(cmd) % this will submit 2 jobs to HPC
        
        %% we need to periodically check whether the job has finished or not
        cmd='sh $CBIG_CODE_DIR/utilities/scripts/CBIG_check_job_status -n LiGSRUT_ME';
        system(cmd)
        
        %% compare two preprocessing pipelines
        CBIG_LiGSR_LME_unittest_intelligence_score_cmp2pipe_GSP( OutputDir );
        
        %% compare result with ground truth
        % no more assert command needed (they are already written inside 
        % CBIG_LiGSR_LME_unittest_intelligence_score_cmp_w_reference_GSP)
        % replace unit test if flag is 1
        if replace_unit_test
            disp("Replacing unit test results for Li2019_GSR, test_intelligence_score_LME_GSP_Case")
            % display differences
            ref_dir = fullfile(getenv('CBIG_TESTDATA_DIR'), ...
                'stable_projects', 'preprocessing', 'Li2019_GSR');
            ref_result = fullfile(ref_dir, 'intelligence_score', 'VarianceComponentModel', ...
                'GSP', 'ref_output', 'compare_2pipe', 'allstats_cmp2pipelines.mat');            
            ref = load(ref_result);
            load(fullfile(OutputDir, 'compare_2pipe', 'allstats_cmp2pipelines.mat') );
            disp(['Difference in percentage improvement is ' num2str(abs(ref.perc_improv - perc_improv))])
            disp(['Difference in jackknife mean is ' num2str(abs(ref.m_jack - m_jack))])
            disp(['Difference in jackknife variance is ' num2str(abs(ref.v_jack - v_jack))])
            disp(['Difference in "#traits whose IQR > 0" is ' num2str(abs(ref.IQR_pos - IQR_pos))])
            disp(['Difference in "#traits whose IQR < 0" is ' num2str(abs(ref.IQR_neg - IQR_neg))])
            disp(['Difference in "#traits whose median > 0" is ' num2str(abs(ref.med_pos - med_pos))])
            disp(['Difference in "#traits whose median < 0" is ' num2str(abs(ref.med_neg - med_neg))])
            
            % save new reference
            save(ref_result, 'perc_improv', 'm_jack', 'v_jack', 'IQR_pos', 'IQR_neg', 'med_pos', 'med_neg');          
        end
        CBIG_LiGSR_LME_unittest_intelligence_score_cmp_w_reference_GSP( fullfile(OutputDir, ...
            'compare_2pipe', 'allstats_cmp2pipelines.mat') );
        
        rmpath(fullfile(getenv('CBIG_CODE_DIR'), 'stable_projects', 'preprocessing', ...
            'Li2019_GSR', 'unit_tests', 'intelligence_score', 'scripts'));
    end
    
    function test_intelligence_score_LME_HCP_Case(testCase)
        %% path setting
        addpath(fullfile(getenv('CBIG_CODE_DIR'), 'stable_projects', 'preprocessing', ...
            'Li2019_GSR', 'unit_tests', 'intelligence_score', 'scripts'));
        OutputDir = fullfile(getenv('CBIG_CODE_DIR'), 'stable_projects', 'preprocessing', ...
            'Li2019_GSR', 'unit_tests', 'output', 'intelligence_score_LME_HCP_Case');
        replace_unit_test = load(fullfile(getenv('CBIG_CODE_DIR'), 'unit_tests', 'replace_unittest_flag'));
        
        % create output dir (IMPORTANT)
        if(exist(OutputDir, 'dir'))
            rmdir(OutputDir, 's')
        end
        mkdir(OutputDir);
        
        %% call Li2019_GSR intelligence_score shell script (variance component model & HCP data)
        cmdfun = fullfile(getenv('CBIG_CODE_DIR'), 'stable_projects', 'preprocessing', 'Li2019_GSR', ...
            'unit_tests', 'intelligence_score', 'scripts', 'CBIG_LiGSR_LME_unittest_PMAT_HCP.sh');
        cmd = [cmdfun, ' ', OutputDir];
        system(cmd);
        
        %% we need to periodically check whether the job has finished or not
        cmd='sh $CBIG_CODE_DIR/utilities/scripts/CBIG_check_job_status -n LiGSRUT_ME';
        system(cmd)
        
        %% compare two preprocessing pipelines
        CBIG_LiGSR_LME_unittest_PMAT_cmp2pipe_HCP( OutputDir );
        
        %% compare result with ground truth
        % no more assert command needed (they are already written inside 
        % CBIG_LiGSR_LME_unittest_PMAT_cmp_w_reference_HCP)
        % replace unit test if flag is 1
        if replace_unit_test
            disp("Replacing unit test results for Li2019_GSR, test_intelligence_score_LME_HCP_Case")
            % display differences
            ref_dir = fullfile(getenv('CBIG_TESTDATA_DIR'), ...
                'stable_projects', 'preprocessing', 'Li2019_GSR');
            ref_result = fullfile(ref_dir, 'intelligence_score', 'VarianceComponentModel', 'HCP', ...
                'ref_output', 'compare_2pipe', 'allstats_cmp2pipelines.mat');            
            ref = load(ref_result);
            load(fullfile(OutputDir, 'compare_2pipe', 'allstats_cmp2pipelines.mat') );
            disp(['Difference in percentage improvement is ' num2str(abs(ref.perc_improv - perc_improv))])
            disp(['Difference in jackknife mean is ' num2str(abs(ref.m_jack - m_jack))])
            disp(['Difference in jackknife variance is ' num2str(abs(ref.v_jack - v_jack))])
            disp(['Difference in "#traits whose IQR > 0" is ' num2str(abs(ref.IQR_pos - IQR_pos))])
            disp(['Difference in "#traits whose IQR < 0" is ' num2str(abs(ref.IQR_neg - IQR_neg))])
            disp(['Difference in "#traits whose median > 0" is ' num2str(abs(ref.med_pos - med_pos))])
            disp(['Difference in "#traits whose median < 0" is ' num2str(abs(ref.med_neg - med_neg))])
            
            % save new reference
            save(ref_result, 'perc_improv', 'm_jack', 'v_jack', 'IQR_pos', 'IQR_neg', 'med_pos', 'med_neg');          
        end
        CBIG_LiGSR_LME_unittest_PMAT_cmp_w_reference_HCP( fullfile(OutputDir, ...
            'compare_2pipe', 'allstats_cmp2pipelines.mat') );
        
        rmpath(fullfile(getenv('CBIG_CODE_DIR'), 'stable_projects', 'preprocessing', ...
            'Li2019_GSR', 'unit_tests', 'intelligence_score', 'scripts'));
    end
    
    function test_intelligence_score_KRR_GSP_Case(testCase)
        %% path setting
        addpath(fullfile(getenv('CBIG_CODE_DIR'), 'stable_projects', 'preprocessing', ...
            'Li2019_GSR', 'unit_tests', 'intelligence_score', 'scripts'));
        OutputDir = fullfile(getenv('CBIG_CODE_DIR'), 'stable_projects', 'preprocessing', ...
            'Li2019_GSR', 'unit_tests', 'output', 'intelligence_score_KRR_GSP_Case');
        replace_unit_test = load(fullfile(getenv('CBIG_CODE_DIR'), 'unit_tests', 'replace_unittest_flag'));
        
        % create output dir (IMPORTANT)
        if(exist(OutputDir, 'dir'))
            rmdir(OutputDir, 's')
        end
        mkdir(OutputDir);
        
        %% call Li2019_GSR intelligence_score shell script (kernel regression method & GSP data)
        cmdfun = fullfile(getenv('CBIG_CODE_DIR'), 'stable_projects', 'preprocessing', 'Li2019_GSR', ...
            'unit_tests', 'intelligence_score', 'scripts', 'CBIG_LiGSR_KRR_unittest_intelligence_score_GSP.sh');
        cmd = [cmdfun, ' ', OutputDir];
        system(cmd); % this will submit 6 jobs to HPC
        
        %% we need to periodically check whether the job has finished or not
        cmd='sh $CBIG_CODE_DIR/utilities/scripts/CBIG_check_job_status -n LiGSRUT_KR';
        system(cmd)
        
        %% compare two preprocessing pipelines
        CBIG_LiGSR_KRR_unittest_intelligence_score_cmp2pipe_GSP( OutputDir );
        
        %% compare result with ground truth
        % no more assert command needed (they are already written inside 
        % CBIG_LiGSR_KRR_unittest_intelligence_score_cmp_w_reference_GSP)
        % replace unit test if flag is 1
        if replace_unit_test
            disp("Replacing unit test results for Li2019_GSR, test_intelligence_score_KRR_GSP_Case")
            % display differences
            ref_dir = fullfile(getenv('CBIG_TESTDATA_DIR'), ...
                'stable_projects', 'preprocessing', 'Li2019_GSR');
            ref_result = fullfile(ref_dir, 'intelligence_score', 'KernelRidgeRegression', ...
                'GSP', 'ref_output', 'compare_2pipe', 'final_result.mat');          
            ref = load(ref_result);
            load(fullfile(OutputDir, 'compare_2pipe', 'final_result.mat') );
            disp(['Difference in accuracy with GSR is ' num2str(max(max(abs(ref.acc_GSR_mean - acc_GSR_mean))))])
            disp(['Difference in baseline accuracies is ' ...
                num2str(max(max(abs(ref.acc_Baseline_mean - acc_Baseline_mean))))])
            disp(['Difference in accuracy between pipelines is ' ...
                num2str(max(max(abs(ref.mean_acc_dif - mean_acc_dif))))])
                       
            % save new reference
            save(ref_result, 'acc_GSR_mean', 'acc_Baseline_mean', 'mean_acc_dif');          
        end
        CBIG_LiGSR_KRR_unittest_intelligence_score_cmp_w_reference_GSP( fullfile(OutputDir, ...
            'compare_2pipe', 'final_result.mat') );
        
        rmpath(fullfile(getenv('CBIG_CODE_DIR'), 'stable_projects', 'preprocessing', ...
            'Li2019_GSR', 'unit_tests', 'intelligence_score', 'scripts'));
    end
    
    function test_intelligence_score_KRR_HCP_Case(testCase)
        %% path setting
        addpath(fullfile(getenv('CBIG_CODE_DIR'), 'stable_projects', 'preprocessing', ...
            'Li2019_GSR', 'unit_tests', 'intelligence_score', 'scripts'));
        OutputDir = fullfile(getenv('CBIG_CODE_DIR'), 'stable_projects', 'preprocessing', ...
            'Li2019_GSR', 'unit_tests', 'output', 'intelligence_score_KRR_HCP_Case');
        replace_unit_test = load(fullfile(getenv('CBIG_CODE_DIR'), 'unit_tests', 'replace_unittest_flag'));
        
        % create output dir (IMPORTANT)
        if(exist(OutputDir, 'dir'))
            rmdir(OutputDir, 's')
        end
        mkdir(OutputDir);
        
        %% call Li2019_GSR intelligence_score shell script (kernel regression method & HCP data)
        cmdfun = fullfile(getenv('CBIG_CODE_DIR'), 'stable_projects', 'preprocessing', 'Li2019_GSR', ...
            'unit_tests', 'intelligence_score', 'scripts', 'CBIG_LiGSR_KRR_unittest_PMAT_HCP.sh');
        cmd = [cmdfun, ' ', OutputDir];
        system(cmd); % this will submit 6 jobs to HPC
        
        %% we need to periodically check whether the job has finished or not
        cmd='sh $CBIG_CODE_DIR/utilities/scripts/CBIG_check_job_status -n LiGSRUT_KR';
        system(cmd)
        
        %% compare two preprocessing pipelines
        CBIG_LiGSR_KRR_unittest_PMAT_cmp2pipe_HCP( OutputDir );
        
        %% compare result with ground truth
        % no more assert command needed (they are already written inside 
        % CBIG_LiGSR_KRR_unittest_PMAT_cmp_w_reference_HCP)
        % replace unit test if flag is 1
        if replace_unit_test
            disp("Replacing unit test results for Li2019_GSR, test_intelligence_score_KRR_HCP_Case")
            % display differences
            ref_dir = fullfile(getenv('CBIG_TESTDATA_DIR'), ...
                'stable_projects', 'preprocessing', 'Li2019_GSR');
            ref_result = fullfile(ref_dir, 'intelligence_score', 'KernelRidgeRegression', ...
                'HCP', 'ref_output', 'compare_2pipe', 'final_result.mat');         
            ref = load(ref_result);
            load(fullfile(OutputDir, 'compare_2pipe', 'final_result.mat') );
            disp(['Difference in accuracy with GSR is ' num2str(max(max(abs(ref.acc_GSR_mean - acc_GSR_mean))))])
            disp(['Difference in baseline accuracies is ' ...
                num2str(max(max(abs(ref.acc_Baseline_mean - acc_Baseline_mean))))])
            disp(['Difference in accuracy between pipelines is ' ...
                num2str(max(max(abs(ref.mean_acc_dif - mean_acc_dif))))])
                       
            % save new reference
            save(ref_result, 'acc_GSR_mean', 'acc_Baseline_mean', 'mean_acc_dif');          
        end
        CBIG_LiGSR_KRR_unittest_PMAT_cmp_w_reference_HCP( fullfile(OutputDir, ...
            'compare_2pipe', 'final_result.mat') );
        
        rmpath(fullfile(getenv('CBIG_CODE_DIR'), 'stable_projects', 'preprocessing', ...
            'Li2019_GSR', 'unit_tests', 'intelligence_score', 'scripts'));
    end
    
    function test_intelligence_score_LRR_GSP_Case(testCase)
        %% path setting
        addpath(fullfile(getenv('CBIG_CODE_DIR'), 'stable_projects', 'preprocessing', ...
            'Li2019_GSR', 'unit_tests', 'intelligence_score', 'scripts'));
        gpso_dir = fullfile(getenv('CBIG_CODE_DIR'), 'external_packages', 'matlab', ...
            'non_default_packages', 'Gaussian_Process');
        OutputDir = fullfile(getenv('CBIG_CODE_DIR'), 'stable_projects', 'preprocessing', ...
            'Li2019_GSR', 'unit_tests', 'output', 'intelligence_score_LRR_GSP_Case');
        replace_unit_test = load(fullfile(getenv('CBIG_CODE_DIR'), 'unit_tests', 'replace_unittest_flag'));
        
        % create output dir (IMPORTANT)
        if(exist(OutputDir, 'dir'))
            rmdir(OutputDir, 's')
        end
        mkdir(OutputDir);
        
        
        %% call Li2019_GSR intelligence_score shell script (linear ridge regression method & GSP data)
        cmdfun = fullfile(getenv('CBIG_CODE_DIR'), 'stable_projects', 'preprocessing', 'Li2019_GSR', ...
            'unit_tests', 'intelligence_score', 'scripts', 'CBIG_LiGSR_LRR_unittest_intelligence_score_GSP.sh');
        cmd = [cmdfun, ' ', gpso_dir, ' ', OutputDir];
        system(cmd); % this will submit 6 jobs to HPC
        
        %% we need to periodically check whether the job has finished or not
        cmd='sh $CBIG_CODE_DIR/utilities/scripts/CBIG_check_job_status -n LiGSRUT_LR';
        system(cmd)
        
        %% compare two preprocessing pipelines
        CBIG_LiGSR_LRR_unittest_intelligence_score_cmp2pipe_GSP( OutputDir );
        
        %% compare result with ground truth
        % no more assert command needed (they are already written inside 
        % CBIG_LiGSR_LRR_unittest_intelligence_score_cmp_w_reference_GSP)
        % replace unit test if flag is 1
        if replace_unit_test
            disp("Replacing unit test results for Li2019_GSR, test_intelligence_score_LRR_GSP_Case")
            % display differences
            ref_dir = fullfile(getenv('CBIG_TESTDATA_DIR'), ...
                'stable_projects', 'preprocessing', 'Li2019_GSR');
            ref_result = fullfile(ref_dir, 'intelligence_score', 'LinearRidgeRegression', ...
                'GSP', 'ref_output', 'compare_2pipe', 'final_result.mat');        
            ref = load(ref_result);
            load(fullfile(OutputDir, 'compare_2pipe', 'final_result.mat') );
            disp(['Difference in accuracy with GSR is ' num2str(max(max(abs(ref.acc_GSR_mean - acc_GSR_mean))))])
            disp(['Difference in baseline accuracies is ' ...
                num2str(max(max(abs(ref.acc_Baseline_mean - acc_Baseline_mean))))])
            disp(['Difference in accuracy between pipelines is ' ...
                num2str(max(max(abs(ref.mean_acc_dif - mean_acc_dif))))])
                       
            % save new reference
            save(ref_result, 'acc_GSR_mean', 'acc_Baseline_mean', 'mean_acc_dif');          
        end
        CBIG_LiGSR_LRR_unittest_intelligence_score_cmp_w_reference_GSP( fullfile(OutputDir, ...
            'compare_2pipe', 'final_result.mat') );
        
        rmpath(fullfile(getenv('CBIG_CODE_DIR'), 'stable_projects', 'preprocessing', ...
            'Li2019_GSR', 'unit_tests', 'intelligence_score', 'scripts'));
    end
    
    function test_intelligence_score_LRR_HCP_Case(testCase)
        %% path setting
        addpath(fullfile(getenv('CBIG_CODE_DIR'), 'stable_projects', 'preprocessing', ...
            'Li2019_GSR', 'unit_tests', 'intelligence_score', 'scripts'));
        gpso_dir = fullfile(getenv('CBIG_CODE_DIR'), 'external_packages', 'matlab', ...
            'non_default_packages', 'Gaussian_Process');
        OutputDir = fullfile(getenv('CBIG_CODE_DIR'), 'stable_projects', 'preprocessing', ...
            'Li2019_GSR', 'unit_tests', 'output', 'intelligence_score_LRR_HCP_Case');
        replace_unit_test = load(fullfile(getenv('CBIG_CODE_DIR'), 'unit_tests', 'replace_unittest_flag'));
        
        % create output dir (IMPORTANT)
        if(exist(OutputDir, 'dir'))
            rmdir(OutputDir, 's')
        end
        mkdir(OutputDir);
        
        %% call Li2019_GSR intelligence_score shell script (linear ridge regression method & HCP data)
        cmdfun = fullfile(getenv('CBIG_CODE_DIR'), 'stable_projects', 'preprocessing', 'Li2019_GSR', ...
            'unit_tests', 'intelligence_score', 'scripts/CBIG_LiGSR_LRR_unittest_PMAT_HCP.sh');
        cmd = [cmdfun, ' ', gpso_dir, ' ', OutputDir];
        system(cmd); % this will submit 6 jobs to HPC
        
        %% we need to periodically check whether the job has finished or not
        cmd='sh $CBIG_CODE_DIR/utilities/scripts/CBIG_check_job_status -n LiGSRUT_LR';
        system(cmd)
        
        %% compare two preprocessing pipelines
        CBIG_LiGSR_LRR_unittest_PMAT_cmp2pipe_HCP( OutputDir );
        
        %% compare result with ground truth
        % no more assert command needed (they are already written inside 
        % CBIG_LiGSR_LRR_unittest_PMAT_cmp_w_reference_HCP)
        % replace unit test if flag is 1
        if replace_unit_test
            disp("Replacing unit test results for Li2019_GSR, test_intelligence_score_LRR_HCP_Case")
            % display differences
            ref_dir = fullfile(getenv('CBIG_TESTDATA_DIR'), ...
                'stable_projects', 'preprocessing', 'Li2019_GSR');
            ref_result = fullfile(ref_dir, 'intelligence_score', 'LinearRidgeRegression', ...
                'HCP', 'ref_output', 'compare_2pipe', 'final_result.mat');      
            ref = load(ref_result);
            load(fullfile(OutputDir, 'compare_2pipe', 'final_result.mat') );
            disp(['Difference in accuracy with GSR is ' num2str(max(max(abs(ref.acc_GSR_mean - acc_GSR_mean))))])
            disp(['Difference in baseline accuracies is ' ...
                num2str(max(max(abs(ref.acc_Baseline_mean - acc_Baseline_mean))))])
            disp(['Difference in accuracy between pipelines is ' ...
                num2str(max(max(abs(ref.mean_acc_dif - mean_acc_dif))))])
                       
            % save new reference
            save(ref_result, 'acc_GSR_mean', 'acc_Baseline_mean', 'mean_acc_dif');          
        end
        CBIG_LiGSR_LRR_unittest_PMAT_cmp_w_reference_HCP( fullfile(OutputDir, ...
            'compare_2pipe', 'final_result.mat') );
        
        rmpath(fullfile(getenv('CBIG_CODE_DIR'), 'stable_projects', 'preprocessing', ...
            'Li2019_GSR', 'unit_tests', 'intelligence_score', 'scripts'));
    end
end

end