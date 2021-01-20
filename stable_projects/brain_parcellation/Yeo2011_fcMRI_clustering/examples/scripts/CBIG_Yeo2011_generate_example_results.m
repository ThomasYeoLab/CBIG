function CBIG_Yeo2011_generate_example_results(output_dir)

% CBIG_Yeo2011_generate_example_results(output_dir)
%
% This function is a matlab wrapper to run the csh script to generate the
% example results for Yeo2011. 
% See README file for more details about the example: 
% $CBIG_CODE_DIR/stable_projects/brain_parcellation/Yeo2011_fcMRI_clustering/examples/README.md
%
% Input:
%     - output_dir:
%       The output directory where example results are saved.
%
% Written by XUE Aihuiping and CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md

script_dir = fileparts(mfilename('fullpath'));
if(~exist(output_dir, 'dir'))
    mkdir(output_dir)
end

% Run the example
cmd = [fullfile(script_dir, 'CBIG_Yeo2011_example.csh'), ' ', output_dir];
system(cmd);

out_file = fullfile(output_dir, 'clustering', 'HNU_example_clusters017_scrub.mat');
if(~exist(out_file, 'file'))
    error('Fail to generate final output')
end

end


