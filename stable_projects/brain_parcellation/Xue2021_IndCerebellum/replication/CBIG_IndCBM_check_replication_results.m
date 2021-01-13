function  CBIG_IndCBM_check_replication_results(output_dir)

% CBIG_IndCBM_check_replication_results(output_dir)
%
% This function checks if the replication results are identical with the 
% reference files.
%
% Input:
%   - output_dir: The output directory with generated replication results.
%
% Written by XUE Aihuiping and CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md
    
CBIG_IndCBM_REP_DIR = getenv('CBIG_IndCBM_REP_DIR');
if(isempty(CBIG_IndCBM_REP_DIR))
    error('CBIG_IndCBM_REP_DIR does not exist. Please check your config file.')
end

fid = fopen(fullfile(CBIG_IndCBM_REP_DIR, 'rep_sub_list.txt'), 'r');
sub_list = textscan(fid,'%s');
sub_list = sub_list{1};
fclose(fid);
set_list = {'all'; 'discovery'; 'replication'};
new_dir = fullfile(output_dir, 'parcellation');
ref_dir = fullfile(CBIG_IndCBM_REP_DIR, 'ref_replication_results');

for i = 1:length(sub_list)
    for j = 1:length(set_list)
        ref_sub = ft_read_cifti(fullfile(ref_dir, [sub_list{i} '_' set_list{j} '_IndCBM_parcellation.dlabel.nii']));
        new_sub = ft_read_cifti(fullfile(new_dir, [sub_list{i} '_' set_list{j} '_IndCBM_parcellation.dlabel.nii']));
        diff_sub = sum(ref_sub.dscalar ~= new_sub.dscalar);
        assert(diff_sub == 0, sprintf(['%d labels are different for ' sub_list{i} ', ' set_list{j} '.'], diff_sub));
    end
end
fprintf('Your replication results are correct!\n')
fclose('all');

end
