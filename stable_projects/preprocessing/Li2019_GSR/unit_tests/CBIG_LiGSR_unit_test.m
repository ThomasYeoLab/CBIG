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
    function test_intelligence_score_LME_GSP_Case(testCase)
        %% path setting
        addpath(fullfile(getenv('CBIG_CODE_DIR'), 'stable_projects', 'preprocessing', ...
            'Li2019_GSR', 'unit_tests', 'intelligence_score', 'scripts'));
        OutputDir = fullfile(getenv('CBIG_CODE_DIR'), 'stable_projects', 'preprocessing', ...
            'Li2019_GSR', 'unit_tests', 'output', 'intelligence_score_LME_GSP_Case');
        
        % create output dir (IMPORTANT)
        if(exist(OutputDir, 'dir'))
            rmdir(OutputDir, 's')
        end
        mkdir(OutputDir);
        
        [~, whoami] = system('whoami');
        whoami = whoami(1:end-1);
        
        %% call Li2019_GSR intelligence_score shell script (variance component model & GSP data)
        cmd = ['$CBIG_CODE_DIR/stable_projects/preprocessing/Li2019_GSR/unit_tests/intelligence_score/', ...
            'scripts/CBIG_LiGSR_LME_unittest_intelligence_score_GSP.sh ', OutputDir];
        system(cmd) % this will submit 2 jobs to HPC
        
        %% we need to periodically check whether the job has finished or not
        cmdout = 1;
        while(cmdout~=0)
            % the job name is  specified as CBIG_LiGSR_LME_unittest_intelligence_score_GSP
            cmd = ['qstat -a | grep CBIG_LiGSR_LME | grep ' whoami ' | wc -l']; 
            [~, cmdout] = system(cmd);
            cmdout = str2num(cmdout(1: end-1)); % after the job is finished, cmdout should be 0
            
            pause(20); % sleep for 20s and check again
        end
        
        %% compare two preprocessing pipelines
        CBIG_LiGSR_LME_unittest_intelligence_score_cmp2pipe_GSP( OutputDir );
        
        %% compare result with ground truth
        % no more assert command needed (they are already written inside 
        % CBIG_LiGSR_LME_unittest_intelligence_score_cmp_w_reference_GSP)
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
        
        % create output dir (IMPORTANT)
        if(exist(OutputDir, 'dir'))
            rmdir(OutputDir, 's')
        end
        mkdir(OutputDir);
        
        [~, whoami] = system('whoami');
        whoami = whoami(1:end-1);
        
        %% call Li2019_GSR intelligence_score shell script (variance component model & HCP data)
        cmd = ['$CBIG_CODE_DIR/stable_projects/preprocessing/Li2019_GSR/unit_tests/intelligence_score/'...
            'scripts/CBIG_LiGSR_LME_unittest_PMAT_HCP.sh ', OutputDir];
        system(cmd);
        
        %% we need to periodically check whether the job has finished or not
        cmdout = 1;
        while(cmdout~=0)
            % the job name is  specified as CBIG_LiGSR_LME_unittest_PMAT_HCP
            cmd = ['qstat -a | grep CBIG_LiGSR_LME | grep ' whoami ' | wc -l']; 
            [~, cmdout] = system(cmd);
            cmdout = str2num(cmdout(1: end-1)); % after the job is finished, cmdout should be 0
            
            pause(20); % sleep for 20s and check again
        end
        
        %% compare two preprocessing pipelines
        CBIG_LiGSR_LME_unittest_PMAT_cmp2pipe_HCP( OutputDir );
        
        %% compare result with ground truth
        % no more assert command needed (they are already written inside 
        % CBIG_LiGSR_LME_unittest_PMAT_cmp_w_reference_HCP)
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
        
        % create output dir (IMPORTANT)
        if(exist(OutputDir, 'dir'))
            rmdir(OutputDir, 's')
        end
        mkdir(OutputDir);
        
        [~, whoami] = system('whoami');
        whoami = whoami(1:end-1);
        
        %% call Li2019_GSR intelligence_score shell script (kernel regression method & GSP data)
        cmd = ['$CBIG_CODE_DIR/stable_projects/preprocessing/Li2019_GSR/unit_tests/intelligence_score/',...
            'scripts/CBIG_LiGSR_KRR_unittest_intelligence_score_GSP.sh ', OutputDir];
        system(cmd); % this will submit 6 jobs to HPC
        
        %% we need to periodically check whether the job has finished or not
        cmdout = 1;
        while(cmdout~=0)
            % the job name is  specified as CBIG_LiGSR_KRR_unittest_intelligence_score_GSP
            cmd = ['qstat -a | grep CBIG_LiGSR_KRR | grep ' whoami ' | wc -l']; 
            [~, cmdout] = system(cmd);
            cmdout = str2num(cmdout(1: end-1)); % after the job is finished, cmdout should be 0
            
            pause(20); % sleep for 20s and check again
        end
        
        %% compare two preprocessing pipelines
        CBIG_LiGSR_KRR_unittest_intelligence_score_cmp2pipe_GSP( OutputDir );
        
        %% compare result with ground truth
        % no more assert command needed (they are already written inside 
        % CBIG_LiGSR_KRR_unittest_intelligence_score_cmp_w_reference_GSP)
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
        
        % create output dir (IMPORTANT)
        if(exist(OutputDir, 'dir'))
            rmdir(OutputDir, 's')
        end
        mkdir(OutputDir);
        
        [~, whoami] = system('whoami');
        whoami = whoami(1:end-1);
        
        %% call Li2019_GSR intelligence_score shell script (kernel regression method & HCP data)
        cmd = ['$CBIG_CODE_DIR/stable_projects/preprocessing/Li2019_GSR/unit_tests/intelligence_score/',...
            'scripts/CBIG_LiGSR_KRR_unittest_PMAT_HCP.sh ', OutputDir];
        system(cmd); % this will submit 6 jobs to HPC
        
        %% we need to periodically check whether the job has finished or not
        cmdout = 1;
        while(cmdout~=0)
            % the job name is  specified as CBIG_LiGSR_KRR_unittest_PMAT_HCP
            cmd = ['qstat -a | grep CBIG_LiGSR_KRR | grep ' whoami ' | wc -l']; 
            [~, cmdout] = system(cmd);
            cmdout = str2num(cmdout(1: end-1)); % after the job is finished, cmdout should be 0
            
            pause(20); % sleep for 20s and check again
        end
        
        %% compare two preprocessing pipelines
        CBIG_LiGSR_KRR_unittest_PMAT_cmp2pipe_HCP( OutputDir );
        
        %% compare result with ground truth
        % no more assert command needed (they are already written inside 
        % CBIG_LiGSR_KRR_unittest_PMAT_cmp_w_reference_HCP)
        CBIG_LiGSR_KRR_unittest_PMAT_cmp_w_reference_HCP( fullfile(OutputDir, ...
            'compare_2pipe', 'final_result.mat') );
        
        rmpath(fullfile(getenv('CBIG_CODE_DIR'), 'stable_projects', 'preprocessing', ...
            'Li2019_GSR', 'unit_tests', 'intelligence_score', 'scripts'));
    end
end

end