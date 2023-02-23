classdef CBIG_Yeo2011_unit_test < matlab.unittest.TestCase
% Written by XUE Aihuiping and CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md

    methods (Test)
        function test_example(testCase)
            CBIG_CODE_DIR = getenv('CBIG_CODE_DIR');
            load(fullfile(CBIG_CODE_DIR, 'unit_tests', 'replace_unittest_flag'));
            % Create output folder
            script_dir = fileparts(mfilename('fullpath'));
            pos_v = strfind(script_dir, filesep);
            proj_dir = fullfile(script_dir(1:pos_v(length(pos_v)) - 1));
            output_dir = fullfile(script_dir, 'output');
            if(exist(output_dir, 'dir'))
                rmdir(output_dir, 's');
            end
            mkdir(output_dir);
            
            % Run the example
            example_dir = fullfile(proj_dir, 'examples');
            addpath(fullfile(example_dir, 'scripts'));
            CBIG_Yeo2011_generate_example_results(output_dir);
            
            if(replace_unittest_flag)
                disp('Replacing unit test reference results for CBIG_Yeo2011_unit_test...');
                CBIG_Yeo2011_check_example_results(output_dir);
                out_file = fullfile(output_dir, 'clustering', 'HNU_example_clusters017_scrub.mat');
                ref_file = fullfile(example_dir, 'results', 'HNU_example_clusters017_scrub.mat');
                copyfile(out_file, ref_file);
            else
                % Compare results
                assert(CBIG_Yeo2011_check_example_results(output_dir), sprintf('Result check failed.'));
            end
            rmdir(output_dir, 's');
            rmpath(fullfile(example_dir, 'scripts'));
        end
    end
end
