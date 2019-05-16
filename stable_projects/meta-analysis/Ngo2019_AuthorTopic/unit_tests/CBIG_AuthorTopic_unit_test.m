classdef CBIG_AuthorTopic_unit_test < matlab.unittest.TestCase
  %
  % Wrapper function to unit-test a coordinate-based meta-analysis with the author-topic model.
  % The inference was performed with the following hard-coded hyperparameters
  %  - K = 2:3 - number of cognitive components to be estimated is 2 and 3
  %  - alpha = 100, eta = 0.01 - hyperparameters of the Dirichlet's distribution
  %  - seeds = 1:2 - 2 reinitializations per K
  %  - Input data is a sampled set of activation foci of self-generated thought dataset
  %    in ./unit_test_sample_coordinates.txt
  %
  %  Output:
  %  - Formatted input data is saved at ./output/data/SelfGeneratedThought_CVBData.mat
  %  - Estimates of the author-topic model's parameters are saved at
  %    ./output/outputs/K<K>/alpha100_eta0.01/params_K<K>_SEED<seed>.mat
  %  - Visualization of the 2-component solution is saved under ./output/figures directory
  %  - BIC files and figures are saved under ./output/BIC directory
  %  - Note that if the unit test pass, the ./output folder will be
  %  removed, otherwise it will be kept for people to the check the results.
  %
  % Example:
  %   runtests('CBIG_AuthorTopic_unit_test.m')
  %   Perfoming unit-testing of a coordinate-based meta-analysis with the author-topic
  %   model using self-generated thought data
  % Written by Gia H. Ngo and CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md
  
  
    methods (Test)
        function test_Case(testCase)
            %% path setting            
            UnitTestDir = fullfile(getenv('CBIG_CODE_DIR'), 'stable_projects', 'meta-analysis', ...
                'Ngo2019_AuthorTopic', 'unit_tests');
            OutputDir = fullfile(UnitTestDir, 'output'); % this output dir is case specific
            
            % create output dir (IMPORTANT)
            if(exist(OutputDir, 'dir'))
                rmdir(OutputDir, 's')
            end
            mkdir(OutputDir);
            
            %% generate results
            
            tic
            
            % add paths to functions specific to author-topic model
            utilitiesDir = fullfile(getenv('CBIG_CODE_DIR'), 'stable_projects', 'meta-analysis', ...
                'Ngo2019_AuthorTopic', 'utilities');
            addpath(fullfile(utilitiesDir, 'preprocessing'));
            addpath(fullfile(utilitiesDir, 'inference'));
            addpath(fullfile(utilitiesDir, 'BIC'));
 
            % prepare dataset
            textFilePath = fullfile(UnitTestDir, 'unit_test_sample_coordinates.txt'); % sample data
            dataDirPath = fullfile(OutputDir, 'data');
            dataFileName = 'SelfGeneratedThought_CVBData.mat';
            system(['mkdir -p ' dataDirPath]);
            CBIG_AuthorTopic_GenerateCVBDataFromText(textFilePath, dataDirPath, dataFileName);
            data = load(fullfile(dataDirPath, dataFileName));
            assert(size(data.act, 1) == 30, 'Wrong number of experiments');
            assert(size(data.act, 2) == 284100, 'Wrong number of voxels in the brain');
            assert(size(data.taskByExp, 1) == 7, 'Wrong number of tasks');
            assert(size(data.taskByExp, 2) == size(data.act, 1), ...
                'Mismatched number of experiments');
            
            % inference
            allKs = 2:3; % set allKs to wider range number of components to discover the full range of solutions
            alpha = 100;
            eta = 0.01;
            doSmoothing = 1;
            cvbData = fullfile(dataDirPath, dataFileName);
            seeds = 1:2; % set seeds = 1:1000 to obtain stable estimates

            for K = allKs
                for seed = seeds
                    CBIG_AuthorTopic_RunInference(seed, K, alpha, eta, doSmoothing, OutputDir, cvbData);
                end
            end
            
            % find best solutions from all reinitializations
            outputsDir = fullfile(OutputDir, 'outputs');
            for K = allKs
                CBIG_AuthorTopic_ComputeBestSolution(outputsDir, K, seeds, alpha, eta);

                bestSolution = load(fullfile(OutputDir, 'outputs', 'bestSolution', ...
                    ['alpha' num2str(alpha) '_eta' num2str(eta)], ...
                    ['BestSolution_K' num2str(K, '%03d') '.mat']));
 
                assert(size(bestSolution.params.beta, 1) == K, ...
                    'Wrong number of componets in beta');
                assert(size(bestSolution.params.beta, 2) == size(data.act, 2), ...
                    'Mismatched number of voxels in beta');
                assert(size(bestSolution.params.theta, 2) == K, ...
                    'Wrong number of componets in theta');
                assert(size(bestSolution.params.theta, 1) == size(data.taskByExp, 1), ...
                    'Mismatched number of tasks in theta');
            end
            
            
            bestSolutionDir = fullfile(OutputDir, 'outputs', 'bestSolution',...
                ['alpha' num2str(alpha) '_eta' num2str(eta)]);
            inputFile = fullfile(bestSolutionDir, 'BestSolution_K002.mat');     
            refTheta = [0.3712, 0.5404, 0.4980, 0.7596, 0.1926, 0.7443, 0.0840; ...
                0.6288, 0.4596, 0.5020, 0.2404, 0.8074, 0.2557, 0.9160]';
            bestK2Solution = load(inputFile);
            thetaDiff = bestK2Solution.params.theta - refTheta;
            maxThetaDiff = max(thetaDiff(:));
            assert(maxThetaDiff < 1e-3, ...
                sprintf('Maximum difference in estimateed theta at K = 2: %f', maxThetaDiff));
 
            % compute BIC measures
            maskPath = fullfile(dataDirPath, 'mask', 'expMask.nii.gz');
            bicDir = fullfile(OutputDir, 'BIC');
            CBIG_AuthorTopic_ComputeBIC(allKs, maskPath, bestSolutionDir, bicDir);
            assert(exist(fullfile(bicDir, 'BIC_3dFWHMx.eps'), 'file') == 2, 'Missing BIC plot');
            
            % clean up paths
            rmpath(fullfile(utilitiesDir, 'preprocessing'));
            rmpath(fullfile(utilitiesDir, 'inference'));
            rmpath(fullfile(utilitiesDir, 'BIC'));
            
            toc
            
            %% remove intermediate output data (IMPORTANT)
            rmdir(OutputDir, 's');
            close all;
        end
              
    end
end


