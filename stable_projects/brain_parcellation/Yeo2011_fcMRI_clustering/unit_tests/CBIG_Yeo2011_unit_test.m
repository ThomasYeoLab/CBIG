classdef CBIG_Yeo2011_unit_test < matlab.unittest.TestCase
% Written by XUE Aihuiping and CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md

    methods (Test)
        function test_example(testCase)
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
            
            % Compare results
            assert(CBIG_Yeo2011_check_example_results(output_dir), sprintf('Result check failed.'));
            rmdir(output_dir, 's');
            rmpath(fullfile(example_dir, 'scripts'));
        end
    end
end
