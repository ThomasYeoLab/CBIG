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
                rmdir(output_dir, 's')
            end
            mkdir(output_dir)
            % Run the example
            example_dir = fullfile(proj_dir, 'examples');
            cmd = [fullfile(example_dir, 'scripts', 'CBIG_Yeo2011_example.csh'), ' ', output_dir];
            system(cmd);
            % Compare results
            out_file = fullfile(script_dir, 'output', 'clustering', 'HNU_example_clusters017_scrub.mat');
            out = load(out_file);
            ref_file = fullfile(example_dir, 'results', 'HNU_example_clusters017_scrub.mat');
            ref = load(ref_file);
            assert(isequal(out.lh_labels, ref.lh_labels), sprintf('LH Labels are different.'));
            assert(isequal(out.rh_labels, ref.rh_labels), sprintf('RH Labels are different.'));
            rmdir(output_dir, 's');
        end
    end
end
