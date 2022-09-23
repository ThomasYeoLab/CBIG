function [test_array] = CBIG_DiffProc_diffusionQC_check_example_result(test_struct,test_case)
% This function checks whether the example diffusion QC results are consistent for the
% unit test.
% Written by Leon Ooi and CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md

CBIG_TEST_DIR = getenv('CBIG_TESTDATA_DIR');
ref_dir = fullfile(CBIG_TEST_DIR,'stable_projects', 'preprocessing', ...
        'CBIG2022_DiffProc', 'diffusionQC');
            
switch test_case
    case "basic"
        test_array = zeros(1,5);
        
        % test whether text file has same values
        load(fullfile(ref_dir, 'ref_output', 'diffusionQC_sub-1_basic.mat'));
        ref = qc_vals;
        
        % test whether text file has same values
        abs_diff = abs(sum(ref) - sum(test_struct.qc_vals));
        if(abs_diff < 1e-6)
            test_array(1) = 1;
        else
            fprintf('QC results are different from reference file by %f.\n', abs_diff);
        end
        
        % read fdt images and compare to reference
        test = MRIread(test_struct.FA_dir);
        ref = MRIread(fullfile(ref_dir, 'ref_output', 'dti_FA.nii.gz'));
        % check size
        if (isequal(size(test.vol), size(ref.vol))) 
            test_array(2) = 1;
        else
            fprintf('Size of test and reference FA images are different.\n');
        end
        % check values
        if (abs(sum(test.vol,'all') - sum(ref.vol,'all')) < 1e-6) 
            test_array(3) = 1;
        else
            fprintf('Test and reference FA values are different.\n');
        end
       
        clear test ref

        test = MRIread(test_struct.SSE_dir);
        ref = MRIread(fullfile(ref_dir, 'ref_output', 'dti_sse.nii.gz'));
        % check size
        if (isequal(size(test.vol), size(ref.vol))) 
            test_array(4) = 1;
        else
            fprintf('Size of test and reference SSE images are different.\n');
        end
        % check values
        if (abs(sum(test.vol,'all') - sum(ref.vol,'all')) < 1e-6) 
            test_array(5) = 1;
        else
            fprintf('Test and reference SSE values are different.\n');
        end
        
        case "rounding"
            test_array = zeros(1,1);
        
            % test whether text file has same values
            load(fullfile(ref_dir, 'ref_output', 'diffusionQC_sub-1_round.mat'));
            ref = qc_vals;

            % test whether text file has same values
            abs_diff = abs(sum(ref) - sum(test_struct.qc_vals));
            if(abs_diff < 1e-6)
               test_array(1) = 1;
            else
               fprintf('QC results are different from reference file by %f.\n', abs_diff);
            end
            
        case "single_shell"
            test_array = zeros(1,9);
        
            % test whether text file has same values
            load(fullfile(ref_dir, 'ref_output', 'diffusionQC_sub-1_single_shell.mat'));
            ref = qc_vals;

            % test whether text file has same values
            abs_diff = abs(sum(ref) - sum(test_struct.qc_vals));
            if(abs_diff < 1e-6)
               test_array(1) = 1;
            else
               fprintf('QC results are different from reference file by %f.\n', abs_diff);
            end
            
            % read single shell images
            test = dlmread(test_struct.bval_dir);
            ref = dlmread(fullfile(ref_dir, 'ref_output', 'sub-1_single_shell_1000.bval'));
            if(abs(sum(ref) - sum(test)) < 1e-6)
                test_array(2) = 1;
            else
                fprintf('bvals are different from reference file.\n');
            end
            test = dlmread(test_struct.bvec_dir);
            ref = dlmread(fullfile(ref_dir, 'ref_output', 'sub-1_single_shell_1000.bvec'));
            if(abs(sum(sum(ref) - sum(test))) < 1e-6)
                test_array(3) = 1;
            else
                fprintf('bvecs results are different from reference file.\n');
            end      
            test = MRIread(test_struct.single_shell_dir);
            ref = MRIread(fullfile(ref_dir, 'ref_output', 'sub-1_single_shell_1000.nii.gz'));
            % check size
            if(isequal(size(test.vol), size(ref.vol)))
                test_array(4) = 1;
            else
                fprintf('Size of test and reference single shell images are different.\n');
            end  
            % check values
            if(abs(sum(test.vol,'all') - sum(ref.vol,'all')) < 1e-6)
                test_array(5) = 1;
            else
                fprintf('Test and reference single shell image values are different.\n');
            end  

        
            % read fdt images and compare to reference
            test = MRIread(test_struct.FA_dir);
            ref = MRIread(fullfile(ref_dir, 'ref_output', 'dti_FA_ss.nii.gz'));
            % check size
            if (isequal(size(test.vol), size(ref.vol))) 
                test_array(6) = 1;
            else
                fprintf('Size of test and reference FA images are different.\n');
            end
            % check values
            if (abs(sum(test.vol,'all') - sum(ref.vol,'all')) < 1e-6) 
                test_array(7) = 1;
            else
                fprintf('Test and reference FA values are different.\n');
            end

            clear test ref

            test = MRIread(test_struct.SSE_dir);
            ref = MRIread(fullfile(ref_dir, 'ref_output', 'dti_sse_ss.nii.gz'));
            % check size
            if (isequal(size(test.vol), size(ref.vol))) 
                test_array(8) = 1;
            else
                fprintf('Size of test and reference SSE images are different.\n');
            end
            % check values
            if (abs(sum(test.vol,'all') - sum(ref.vol,'all')) < 1e-6) 
                test_array(9) = 1;
            else
                fprintf('Test and reference SSE values are different.\n');
            end
end
        
end