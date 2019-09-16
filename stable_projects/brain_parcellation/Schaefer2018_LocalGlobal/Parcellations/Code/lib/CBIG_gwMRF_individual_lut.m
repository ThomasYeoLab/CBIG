function CBIG_gwMRF_individual_lut(output_dir)

% CBIG_gwMRF_individual_lut(output_dir)
%
% This function generates the lookup table for individual subjects.
%
% Input:
%      -output_dir: 
%       Path of the output folder, There should already be annot files in
%       fsaverage space '?h.Schaefer2018_?Parcels_?Networks_order.annot'
%       for all resolutions, both hemisphere saved under 
%       'output_dir/FreeSurfer5.3/fsaverage/label'
% 
% Example:
% CBIG_gwMRF_individual_lut(output_dir)
% 
% Written by XUE Aihuiping and CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md

CBIG_CODE_DIR = getenv('CBIG_CODE_DIR');
annot_dir = fullfile(output_dir, 'FreeSurfer5.3', 'fsaverage', 'label');
lut_dir = fullfile(output_dir, 'project_to_individual');
if(~exist(lut_dir,'dir'))
    mkdir(lut_dir);
end

for k = 7:10:17
    for p = 100:100:1000
        lut_file = fullfile(lut_dir, ['Schaefer2018_' num2str(p) 'Parcels_' num2str(k) 'Networks_order_LUT.txt']);
        subcortical = fullfile(CBIG_CODE_DIR, 'stable_projects', 'brain_parcellation', ...
            'Schaefer2018_LocalGlobal', 'Parcellations', 'Code', 'input', 'subcortical_lut.txt');

        if(exist(lut_file, 'file'))
            delete(lut_file);
        end
        copyfile(subcortical, lut_file);

        lh_annot = fullfile(annot_dir, ['lh.Schaefer2018_' num2str(p) 'Parcels_' num2str(k) 'Networks_order.annot']);
        [~, ~, colortable] = read_annotation(lh_annot);

        fid = fopen(lut_file,'a');
        line = '1000    LH_Background+FreeSurfer_Defined_Medial_Wall  1   1   1   0';
        fprintf(fid, [line '\n']);
        p_2 = p / 2;
        for i = 1:p_2
            line='                                                                  0';
            line(1:4) = sprintf('%4d', 1000+i);
            name = colortable.struct_names{i+1};
            line(9:8+length(name)) = name;
            line(55:66) = sprintf('%-4d%-4d%-4d', colortable.table(i+1,1), colortable.table(i+1,2), ...
                colortable.table(i+1,3));
            fprintf(fid, [line '\n']);
        end

        rh_annot = fullfile(annot_dir, ['rh.Schaefer2018_' num2str(p) 'Parcels_' num2str(k) 'Networks_order.annot']);
        [~, ~, colortable] = read_annotation(rh_annot);
        fprintf(fid, '\n');
        line='2000    RH_Background+FreeSurfer_Defined_Medial_Wall  1   1   1   0';
        fprintf(fid, [line '\n']);
        for i=1:p_2
            line = '                                                                  0';
            line(1:4) = sprintf('%4d', 2000 + i);
            name = colortable.struct_names{i + 1};
            line(9:8+length(name)) = name;
            line(55:66) = sprintf('%-4d%-4d%-4d', colortable.table(i+1,1), colortable.table(i+1,2), ...
                colortable.table(i+1,3));
            fprintf(fid, [line '\n']);
        end
        fclose(fid);
    end
end

end
