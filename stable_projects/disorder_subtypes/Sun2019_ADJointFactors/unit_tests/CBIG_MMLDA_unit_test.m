classdef CBIG_MMLDA_unit_test < matlab.unittest.TestCase
% Written by Nanbo Sun and CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md

    methods (Test)
        function test_example(testCase)
            CBIG_CODE_DIR = getenv('CBIG_CODE_DIR');
            load(fullfile(CBIG_CODE_DIR, 'unit_tests','replace_unittest_flag'));
            cur_dir = fullfile(CBIG_CODE_DIR, 'stable_projects/disorder_subtypes/', ...
                'Sun2019_ADJointFactors/unit_tests');
            addpath(fullfile(CBIG_CODE_DIR, 'stable_projects/disorder_subtypes/', ...
                'Sun2019_ADJointFactors/examples'))
            if replace_unittest_flag
                % create output folder
                out_dir = fullfile(cur_dir, 'output');
                mkdir(out_dir)

                % run the example
                queue = 'circ-spool'
                CBIG_MMLDA_example_wrapper(out_dir, queue)

                out_res_dir = fullfile(out_dir,'estimation/ADNI2_bl_AD_meanCNstdALL_plus1_2sub/k2/r3');
                
                % copy results to reference directory
                ref_dir = fullfile(CBIG_CODE_DIR,'stable_projects/disorder_subtypes/', ...
                    'Sun2019_ADJointFactors/examples/correct_output');
                
                ref_res_dir = fullfile(ref_dir,'estimation/ADNI2_bl_AD_meanCNstdALL_plus1_2sub/k2/r3');
                file_list = {'final.beta1' 'final.beta2' 'final.gamma' 'final.other' 'likelihood.dat'};
                for ii = 1:length(file_list)
                    filename = file_list{ii};
                    source = fullfile(out_res_dir,filename);
                    destination = fullfile(ref_res_dir,filename);
                    copyfile(source,destination)
                end
                
                % remove the output directory
                rmdir(out_dir, 's')
            else
                % create output folder
                out_dir = fullfile(cur_dir, 'output');
                mkdir(out_dir)

                % run the example
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
            end
            rmpath(fullfile(CBIG_CODE_DIR, 'stable_projects/disorder_subtypes/', ...
                'Sun2019_ADJointFactors/examples'))    
        end
    end
end