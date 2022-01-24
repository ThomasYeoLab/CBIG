function CBIG_TRBPC_example_multiKRR_get_inputs(outdir)

% This function generate some input files to run teh multikKRR example
% Written by Jianzhong Chen and CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md

example_dir = fullfile(getenv('CBIG_CODE_DIR'),'stable_projects','predict_phenotypes','ChenTam2022_TRBPC','examples');
feature_files = cell(2,1);
feature_files{1} = fullfile(example_dir,'input','FC_rs.mat');
feature_files{2} = fullfile(example_dir,'input','FC_nback.mat');
save(fullfile(outdir, 'FC_all.mat'), 'feature_files');

end
