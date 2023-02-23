function tf = CBIG_IndCBM_check_example_results(output_dir)

% CBIG_IndCBM_check_example_results(output_dir)
%
% This function checks if the generated example results are identical with
% the reference files.
%
% Input:
%   - output_dir: The output directory with generated example results.
%
% Written by XUE Aihuiping and CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md
    
CBIG_CODE_DIR = getenv('CBIG_CODE_DIR');

ref_dir = fullfile(CBIG_CODE_DIR, 'stable_projects', ...
    'brain_parcellation', 'Xue2021_IndCerebellum', 'examples', 'ref_output');

ref_sub1 = ft_read_cifti(fullfile(ref_dir, 'sub1', 'IndCBM_parcellation_top100.dlabel.nii'));
ref_sub2 = ft_read_cifti(fullfile(ref_dir, 'sub2', 'IndCBM_parcellation_top100.dlabel.nii'));

new_sub1 = ft_read_cifti(fullfile(output_dir, 'parcellation', 'sub1', 'IndCBM_parcellation_top100.dlabel.nii'));
new_sub2 = ft_read_cifti(fullfile(output_dir, 'parcellation', 'sub2', 'IndCBM_parcellation_top100.dlabel.nii'));

diff_sub1 = sum(ref_sub1.dscalar ~= new_sub1.dscalar);
diff_sub2 = sum(ref_sub2.dscalar ~= new_sub2.dscalar);

tf = true;
if(diff_sub1 ~= 0)
    sprintf('%d labels are different for subject 1 ', diff_sub1);
    tf = false;
end
if(diff_sub2 ~= 0)
    sprintf('%d labels are different for subject 2 ', diff_sub2);
    tf = false;
end
if(tf)
    fprintf('Your example results are correct!\n');
end
fclose('all');

end
