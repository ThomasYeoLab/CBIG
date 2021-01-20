classdef CBIG_MMLDA_unit_test < matlab.unittest.TestCase
% Written by Nanbo Sun and CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md

    methods (Test)
        function test_example(testCase)
            % create output folder
            CBIG_CODE_DIR = getenv('CBIG_CODE_DIR');
            cur_dir = fullfile(CBIG_CODE_DIR, 'stable_projects/disorder_subtypes/Sun2019_ADJointFactors/unit_tests');
            out_dir = fullfile(cur_dir, 'output');
            mkdir(out_dir)
            
            % run the example
            addpath(fullfile(CBIG_CODE_DIR, 'stable_projects/disorder_subtypes/Sun2019_ADJointFactors/examples'))
            queue = 'circ-spool'
            CBIG_MMLDA_example_wrapper(out_dir, queue)
            
            % compare the results
            params = {'beta1' 'beta2' 'gamma' 'likelihood' 'inferred gamma'};
            flags = CBIG_MMLDA_check_example_results(out_dir);
            for n = 1:length(flags)
                assert(flags(n), sprintf('Result differ for %s', params{n}))
            end

            % remove the output directory
            rmdir(out_dir, 's')
            rmpath(fullfile(CBIG_CODE_DIR, 'stable_projects/disorder_subtypes/Sun2019_ADJointFactors/examples'))
        end
    end
end