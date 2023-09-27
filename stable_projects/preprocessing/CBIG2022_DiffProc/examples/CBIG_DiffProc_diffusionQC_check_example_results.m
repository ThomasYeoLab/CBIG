function [test_array] = CBIG_DiffProc_diffusionQC_check_example_results(test_struct)
% This function checks whether the example diffusion QC results are consistent
% with the generated example.
%
% Input:
% test_struct - A struct containing qc_txt, FA_dir and SSE_dir, which refer
%               to the paths to the corresponding qc text file and output
%               files from fdt.
%               E.g. s.qc_txt = fullfile('QC_output', 'sub1_QC.txt');
%                    s.FA_dir = fullfile('fdt', 'dti_FA.nii.gz');
%                    s.SSE_dir = fullfile('fdt', 'dti_sse.nii.gz');
%
% Output:
% test_array - A array with 5 entries. 1 indicates a pass and 0 is a fail.
%              test_array(1): QC text file values are similar
%              test_array(2): FA image volume sizes are equal
%              test_array(3): FA image values are similar
%              test_array(4): SSE image volume sizes are equal
%              test_array(5): SSE image values are similar
%
% Written by Leon Ooi and CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md

% set directories
CBIG_CODE_DIR = getenv('CBIG_CODE_DIR');
ref_dir = fullfile(CBIG_CODE_DIR,'stable_projects', 'preprocessing', ...
        'CBIG2022_DiffProc', 'examples', 'diffusionQC_example_data', ...
        'ref_output');

% compare results
test_array = zeros(1,5);

% test whether text file has same values
test_qc = read_QC_file(test_struct.qc_txt);
ref_qc = read_QC_file(fullfile(ref_dir, 'QC_output', 'sub-A00059845_QC.txt'));
        
% test whether text file has same values
abs_diff = abs(sum(ref_qc.qc_vals) - sum(test_qc.qc_vals));
if(abs_diff < 1e-6)
    test_array(1) = 1;
else
    fprintf('QC results are different from reference file by %f.\n', abs_diff);
end

% read fdt images and compare to reference
test = MRIread(test_struct.FA_dir);
ref = MRIread(fullfile(ref_dir, 'fdt', 'dti_FA.nii.gz'));
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
ref = MRIread(fullfile(ref_dir, 'fdt', 'dti_sse.nii.gz'));
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
               
end

function s = read_QC_file(txt_file)
% read QC file and compare to reference
fid = fopen(txt_file,'r'); 
tline = fgetl(fid);
i = 1;
while ischar(tline)
    tline = fgetl(fid);
    try
        s_qc_all{:,i} = split(tline, ": ");
        s.qc_vals(:,i) = str2num(s_qc_all{:,i}{2});
    end
    i = i + 1;
end
fclose(fid);
end