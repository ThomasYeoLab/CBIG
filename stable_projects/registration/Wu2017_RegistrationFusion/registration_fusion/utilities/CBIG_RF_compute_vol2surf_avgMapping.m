function CBIG_RF_compute_vol2surf_avgMapping(sub_list, input_dir, input_prefix, output_dir, output_prefix)
% CBIG_RF_compute_vol2surf_avgMapping(sub_list, input_dir, input_prefix, output_dir, output_prefix)
%
% This function computes volume-to-surface avgerage mapping of x/y/z index 
% across a given list of subjects.
%
% Input:
%     - sub_list     :
%                      absolute/relative path to subject list file, 
%                        which contains a subject name each line
%     - input_dir    :
%                      absolute/relative path to directory where input can be found
%                      inputs are index files projected to fsaverage, 
%                        which can be read by MRIread()
%     - input_prefix :
%                      a fixed prefix for each subject's mapping, 
%                        which may specify the approach and template used
%                      For example, an index file projected to left hemisphere of 
%                        fsaverage through Subject 0001 is named: 
%                        lh.xIndex_input_prefix_to_Sub0001_Ses1_FS_to_fsaverage.nii.gz
%     - output_dir   :
%                      absolute/relative path to directory where output should be stored
%     - output_prefix:
%                      desired prefix for the outputs
%
% Output:
%     - There is no function output.
%     - 2 mapping files are created in output_dir:
%           lh.avgMapping_[output_prefix].mat
%           rh.avgMapping_[output_prefix].mat
%
% Example:
% CBIG_RF_compute_vol2surf_avgMapping('~/data/GSP_subjectnames.csv', '../results/index_fsaverage/', 
%               'RF_M3Z_MNI', '../results/mappings/', 'RF_M3Z_MNI2fsaverage_avg1490Sub')
% This command reads in subject mappings from  ../results/index_fsaverage directory with 
% prefix of 'RF_M3Z_MNI'. An average mapping is generated across the subjects provided in 
% GSP_subjectnames.csv and put into ../results/mappings directory.
%
% Written by Wu Jianxiao and CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md

if nargin < 5
    disp('usage: CBIG_RF_compute_vol2surf_avgMaping(sub_list, input_dir, input_prefix, output_dir, output_prefix)');
    return
end

%Read subject names
fid = fopen(sub_list, 'r');
i = 1;
while ~feof(fid)
    sub_names{i} = fgetl(fid);
    i = i + 1;
end
fclose(fid);

%Loop through each subject in each hemisphere
for hemis = {'lh', 'rh'}
    hemi = hemis{1};
    
    %Read and combine x, y, and z index volume results for each subject
    for i = 1:length(sub_names)
        disp([num2str(i) ': ' sub_names{i} ' in ' hemi]);
        x = MRIread([input_dir '/' hemi '.xIndex_' input_prefix '_to_' sub_names{i} '_to_fsaverage.nii.gz']);
        y = MRIread([input_dir '/' hemi '.yIndex_' input_prefix '_to_' sub_names{i} '_to_fsaverage.nii.gz']);
        z = MRIread([input_dir '/' hemi '.zIndex_' input_prefix '_to_' sub_names{i} '_to_fsaverage.nii.gz']);
        if i == 1
            ras = [x.vol; y.vol; z.vol];
        else
            ras = ras + [x.vol; y.vol; z.vol];
        end
    end
    
    %Get average of the results and save to mat file
    ras = ras ./ length(sub_names);
    save([output_dir '/' hemi '.avgMapping_' output_prefix '.mat'], 'ras');
end

end
