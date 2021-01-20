function tf = CBIG_Yeo2011_check_example_results(output_dir)

% result = CBIG_Yeo2011_check_example_results(output_dir)
%
% This function compares the example results saved in <output_dir> with the
% reference results.
%
% Input:
%     - output_dir:
%       The output directory where example results are saved.
%
% Output:
%     - tf:
%       The flag indicates whether the example results saved in 
%       <output_dir> match the reference results. If example results match 
%       reference results within tolerance, tf = 1, otherwise tf = 0.
%
% Written by XUE Aihuiping and CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md

script_dir = fileparts(mfilename('fullpath'));
pos_v = strfind(script_dir, filesep);
example_dir = fullfile(script_dir(1:pos_v(length(pos_v)) - 1));
if(~exist(output_dir, 'dir'))
    error('Result folder does not exist')
end

% Compare results
out_file = fullfile(output_dir, 'clustering', 'HNU_example_clusters017_scrub.mat');
if(~exist(out_file, 'file'))
    error('Result file does not exist')
end
out = load(out_file);
ref_file = fullfile(example_dir, 'results', 'HNU_example_clusters017_scrub.mat');
ref = load(ref_file);
tf = true;
lh_diff = sum(out.lh_labels ~= ref.lh_labels);
if(lh_diff == 0)
    disp('LH Labels are the same');
elseif(lh_diff/length(ref.lh_labels) < 0.0005) % Less then 0.05% of total voxels are different.
    disp([num2str(lh_diff) ' LH Labels are different.']);
    warning('The small difference might be caused by different running environments.');
else % Too many voxels are different.
    disp([num2str(lh_diff) ' LH Labels are different.']);
    tf = false;
end
rh_diff = sum(out.rh_labels ~= ref.rh_labels);
if(rh_diff == 0)
    disp('RH Labels are the same');
elseif(rh_diff/length(ref.rh_labels) < 0.0005) % Less then 0.05% of total voxels are different.
    disp([num2str(rh_diff) ' RH Labels are different.']);
    warning('The small difference might be caused by different running environments.');
else % Too many voxels are different.
    disp([num2str(rh_diff) ' RH Labels are different.']);
    tf = false;
end

end
