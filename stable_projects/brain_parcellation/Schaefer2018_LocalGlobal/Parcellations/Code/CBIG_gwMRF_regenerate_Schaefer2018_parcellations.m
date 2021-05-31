function CBIG_gwMRF_regenerate_Schaefer2018_parcellations(output_dir)

% CBIG_gwMRF_regenerate_Schaefer2018_parcellations(output_dir)
% 
% This function uses current Schaefer2018 parcellation to regenerate new
% version of Schaefer2018 parcellation under given folder. 
%
% In our previous versions (v0.8.1-Schaefer2018_LocalGlobal and earlier),
% there are some parcels incorrectly labeled. This function uses 
% the Schaefer2018 parcellation files in CBIG_Private
% (v0.8.1-Schaefer2018_LocalGlobal) to regenerate a corrected version 
% (v0.14.3-Update_Yeo2011_Schaefer2018_labelname).
% 
% You can modify this function when you want to make another update on 
% Schaefer2018 parcellation. This function refers to the files in the repo. 
% Please be careful about the version when you modify this function to do 
% further updates. We recommend you to first make a pull request to modify
% the code and then replace the parcellation files. 
% 
% Input:
%      -output_dir: 
%       Full path of the output folder. Sub folders "Code", "Freesurfer5.3",
%        "HCP", "project_to_individual" "MNI" would be created under the 
%       output_dir to store the respective output files.
%
% Example:
% CBIG_gwMRF_regenerate_Schaefer2018_parcellations(output_dir)
%       
% Written by XUE Aihuiping and CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md

CBIG_CODE_DIR = getenv('CBIG_CODE_DIR');
Schaefer2018_dir = fullfile(CBIG_CODE_DIR, 'stable_projects', 'brain_parcellation', ...
    'Schaefer2018_LocalGlobal', 'Parcellations');
addpath(fullfile(Schaefer2018_dir, 'Code', 'lib'));
grow_dir = fullfile(CBIG_CODE_DIR, 'stable_projects', 'brain_parcellation', ...
    'Yeo2011_fcMRI_clustering', '1000subjects_reference', ...
    'Yeo_JNeurophysiol11_SplitLabels', 'grow_boundary', 'annot');

mkdir(fullfile(output_dir, 'FreeSurfer5.3', 'fsaverage5', 'label'));
mkdir(fullfile(output_dir, 'FreeSurfer5.3', 'fsaverage6', 'label'));
mkdir(fullfile(output_dir, 'FreeSurfer5.3', 'fsaverage', 'label'));
mkdir(fullfile(output_dir, 'HCP', 'fslr32k', 'cifti'));
mkdir(fullfile(output_dir, 'MNI'));

% Adjust the order of the components for 17 networks
lh_reorder_idx = [1:28,30,29,31:57];
rh_reorder_idx =[1:30,32,31,33:57];

for k = 7:10:17    
    % Read the growed Yeo2011 split components parcellation
    lh_yeosplit_annot = fullfile(grow_dir, ['lh.Yeo2011_' num2str(k) 'Networks_growed.annot']);
    rh_yeosplit_annot = fullfile(grow_dir, ['rh.Yeo2011_' num2str(k) 'Networks_growed.annot']);
    
    % Define the name of the networks
    if(k==7)
        network_name = [{'Medial_Wall'}; {'Vis'}; {'SomMot'}; {'DorsAttn'}; ...
            {'SalVentAttn'}; {'Limbic'}; {'Cont'}; {'Default'}];
    elseif(k==17)
        network_name = [{'Medial_Wall'}; {'VisCent'}; {'VisPeri'}; {'SomMotA'}; ...
            {'SomMotB'}; {'DorsAttnA'}; {'DorsAttnB'}; {'SalVentAttnA'}; ...
            {'SalVentAttnB'}; {'LimbicA'}; {'LimbicB' }; {'ContC'}; {'ContA'}; ...
            {'ContB'}; {'TempPar'}; {'DefaultC'}; {'DefaultA'}; {'DefaultB'}];
    end
    
    for p = 100:100:1000
        disp(['Reordering Schaefer2018 parcellation on fsaverage, ' num2str(p) ' Parcels, ' num2str(k) ' Networks.']);
        % Original Schaefer parcellation
        lh_alex_annot = fullfile(CBIG_CODE_DIR, 'stable_projects', 'brain_parcellation', ...
            'Schaefer2018_LocalGlobal', 'Parcellations', 'FreeSurfer5.3', 'fsaverage', 'label', ...
            ['lh.Schaefer2018_' num2str(p) 'Parcels_' num2str(k) 'Networks_order.annot']);
        rh_alex_annot = fullfile(CBIG_CODE_DIR, 'stable_projects', 'brain_parcellation', ...
            'Schaefer2018_LocalGlobal', 'Parcellations', 'FreeSurfer5.3', 'fsaverage', 'label', ...
            ['rh.Schaefer2018_' num2str(p) 'Parcels_' num2str(k) 'Networks_order.annot']);
        
        % Match the parcels with Yeo2011 split components
        if(k==7)
            [lh_label, lh_table] = CBIG_gwMRF_match_yeo2011(lh_alex_annot, lh_yeosplit_annot, ...
                'lh', network_name, k, p);
            [rh_label, rh_table] = CBIG_gwMRF_match_yeo2011(rh_alex_annot, rh_yeosplit_annot, ...
                'rh', network_name, k, p);
        else
            [lh_label, lh_table] = CBIG_gwMRF_match_yeo2011(lh_alex_annot, lh_yeosplit_annot, ...
                'lh', network_name, k, p, lh_reorder_idx);
            [rh_label, rh_table] = CBIG_gwMRF_match_yeo2011(rh_alex_annot, rh_yeosplit_annot, ...
                'rh', network_name, k, p, rh_reorder_idx);
        end
        
        % Add information about the parcellation
        lh_table.orig_tab=['Schaefer2018_' num2str(p) 'Parcels_' num2str(k) 'Networks'];
        rh_table.orig_tab=['Schaefer2018_' num2str(p) 'Parcels_' num2str(k) 'Networks'];
        
        % Write new Schaefer parcellation
        lh_output_annot = fullfile(output_dir, 'FreeSurfer5.3', 'fsaverage', 'label', ...
            ['lh.Schaefer2018_' num2str(p) 'Parcels_' num2str(k) 'Networks_order.annot']);
        rh_output_annot = fullfile(output_dir, 'FreeSurfer5.3', 'fsaverage', 'label', ...
            ['rh.Schaefer2018_' num2str(p) 'Parcels_' num2str(k) 'Networks_order.annot']);
        vertices = 0:(length(lh_label)-1);        
        write_annotation(lh_output_annot, vertices', lh_label', lh_table);
        write_annotation(rh_output_annot, vertices', rh_label', rh_table);
        
        % Do not change the parcel order for 100-400, 17 networks and
        % 100-800, 1000, 7 networks because after applying the new
        % algorithm, the network labels were not changed for these
        % resolutions. In the future updates, we might change the order of
        % other resolutions. Please modify the code below. 
        if((p <= 400 && k == 17) || (p ~= 900 && k == 7))
            [lh_label_old, lh_table_modify] = keep_original_order(lh_alex_annot, lh_output_annot);
            [rh_label_old, rh_table_modify] = keep_original_order(rh_alex_annot, rh_output_annot);
            write_annotation(lh_output_annot, vertices', lh_label_old, lh_table_modify);
            write_annotation(rh_output_annot, vertices', rh_label_old, rh_table_modify);
        end
        
        % Use labels on fsaverage to create files for fsaverage5/6, fsLR and MNI
        create_all_files(output_dir, k, p);
        
        % Initially, after projection to fsLR space, it was noticed that 
        % the boundaries of the parcels have shifted. This could be due to 
        % changes of the system and software version. Therefore, we read 
        % the old Schaefer Parcellation in fsLR space and reordered them 
        % using the new ordering. This will prevent the boundaries of the 
        % parcels from shifting in the new version as compared to the old 
        % version.
        disp(['Reordering cifti files, ' num2str(p) ' Parcels, ' num2str(k) ' Networks.']);
        reorder_hcp(output_dir, k, p);
    end
end

% Create lookup table for individual projection
CBIG_gwMRF_individual_lut(output_dir);

rmpath(fullfile(Schaefer2018_dir, 'Code', 'lib'));

end

function [label_old, table_modify] = keep_original_order(old_annot, new_annot)

% [label_old, table_modify] = keep_original_order(old_annot, new_annot)
% 
% This function uses the new annot file to rename the color table of the  
% old annot file while keeping the parcel order and label value exactly the
% same as the old annot file.  
% 
% Input:
%      -old_annot: 
%       Old Schaefer2018 parcellation annot file
% 
%      -new_annot: 
%       New Schaefer2018 parcellation annot file
% 
% Output:
%      -label_old:
%       Labels read from Old Schaefer2018 parcellation annot file
% 
%      -table_modify:
%       Color table with old label values and new parcel names
%
% Example:
% [label_old, table_modify] = keep_original_order(old_annot, new_annot)

[~, label_old, table_old]=read_annotation(old_annot);
[~, label_new, table_new]=read_annotation(new_annot);

list = unique(label_old);
list(list == table_old.table(1,5)) = [];

table_modify = table_old;
for i = 1:length(list)
    index_old = find(table_old.table(:, 5) == list(i));
    index_new = find(table_new.table(:, 5) == mode(label_new(label_old == list(i))));
    table_modify.struct_names{index_old} = table_new.struct_names{index_new};
end

table_modify.orig_tab = table_new.orig_tab;

end

function reorder_hcp(output_dir, k, p)

% reorder_hcp(output_dir, k, p)
% 
% This fucntion compares the old annot and the new annot files to reorder 
% the old cifti files.
% Possibly due to changes of the system and software version in our server, 
% if we directly project the parcellations in fsaverage space to fs_LR
% space, the boundaries are slightly shifted compared to the ones in our 
% repo (v0.8.1-Schaefer2018_LocalGlobal).
% To keep the boundaries the same as the previous version, we read the old
% Schaefer Parcellation in fsLR space and reordered them using the new 
% ordering. 
% 
% Input:
%      -output_dir: 
%       Path of the output folder
% 
%      -k: 
%       A Scalar. Total number of networks.
% 
%      -p: 
%       A Scalar. Total number of parcels.
%
% Example:
% reorder_hcp(output_dir, 17, 400)

CBIG_CODE_DIR = getenv('CBIG_CODE_DIR');
old_dlabel_file = fullfile(CBIG_CODE_DIR, 'stable_projects', 'brain_parcellation', ...
    'Schaefer2018_LocalGlobal', 'Parcellations', 'HCP', 'fslr32k', 'cifti', ...
    ['Schaefer2018_' num2str(p) 'Parcels_' num2str(k) 'Networks_order.dlabel.nii']);
old_dlabel = ft_read_cifti(old_dlabel_file);
output_name = fullfile(output_dir, 'HCP', 'fslr32k', 'cifti', ...
    ['Schaefer2018_' num2str(p) 'Parcels_' num2str(k) 'Networks_order']);
lh_old_annot = fullfile(CBIG_CODE_DIR, 'stable_projects', 'brain_parcellation', ...
    'Schaefer2018_LocalGlobal', 'Parcellations', 'FreeSurfer5.3', 'fsaverage', 'label', ...
    ['lh.Schaefer2018_' num2str(p) 'Parcels_' num2str(k) 'Networks_order.annot']);
rh_old_annot = fullfile(CBIG_CODE_DIR, 'stable_projects', 'brain_parcellation', ...
    'Schaefer2018_LocalGlobal', 'Parcellations', 'FreeSurfer5.3', 'fsaverage', 'label', ...
    ['rh.Schaefer2018_' num2str(p) 'Parcels_' num2str(k) 'Networks_order.annot']);
lh_new_annot = fullfile(output_dir, 'FreeSurfer5.3', 'fsaverage', 'label', ...
    ['lh.Schaefer2018_' num2str(p) 'Parcels_' num2str(k) 'Networks_order.annot']);
rh_new_annot = fullfile(output_dir, 'FreeSurfer5.3', 'fsaverage', 'label', ...
    ['rh.Schaefer2018_' num2str(p) 'Parcels_' num2str(k) 'Networks_order.annot']);
index = CBIG_gwMRF_index_trans_btwn2versions(lh_old_annot, rh_old_annot, lh_new_annot, rh_new_annot);
load(fullfile(CBIG_CODE_DIR, 'stable_projects', 'brain_parcellation', ...
    'Schaefer2018_LocalGlobal', 'Parcellations', 'Code', 'input', 'common_cifti.mat'));
new_dlabel = old_dlabel;
for i=1:p
    new_dlabel.parcels(old_dlabel.parcels == index(i)) = i;
end
% Exclude medial wall vertices defined by HCP
lh_mesh_fslr_32k=CBIG_read_fslr_surface('lh','fs_LR_32k','inflated','medialwall.annot');
rh_mesh_fslr_32k=CBIG_read_fslr_surface('rh','fs_LR_32k','inflated','medialwall.annot');
new_dlabel.parcels([lh_mesh_fslr_32k.MARS_label == 1;rh_mesh_fslr_32k.MARS_label == 1]) = 0;


ft_write_cifti(output_name, new_dlabel, 'parameter', 'parcels');
input = [output_name, '.dscalar.nii'];
label_list_file = [output_name,'_info.txt'];
output = [output_name, '.dlabel.nii'];
command = ['wb_command -cifti-label-import ', input, ' ', label_list_file, ' ', output];
system(command);

end

function create_all_files(output_dir, k, p)

% create_all_files(output_dir, k, p)
% 
% This fucntion use the parcellation on fsaverage to generate files for 
% fsaverage5/6, fsLR and MNI. 
% 
% Input:
%      -output_dir: 
%       Path of the output folder
% 
%      -k: 
%       A Scalar. Total number of networks.
% 
%      -p: 
%       A Scalar. Total number of parcels.
%
% Example:
% create_all_files(output_dir, 17, 400)

lh_annot = fullfile(output_dir, 'FreeSurfer5.3', 'fsaverage', 'label', ...
    ['lh.Schaefer2018_' num2str(p) 'Parcels_' num2str(k) 'Networks_order.annot']);
rh_annot = fullfile(output_dir, 'FreeSurfer5.3', 'fsaverage', 'label', ...
    ['rh.Schaefer2018_' num2str(p) 'Parcels_' num2str(k) 'Networks_order.annot']);

[~, lh, lh_s] = read_annotation(lh_annot);
[~, rh, rh_s] = read_annotation(rh_annot);

% check for usage, we will ignore the first parcel so we start with 2
for i = 2:p/2+1
    if(length(find(lh == lh_s.table(i, 5))) < 1)
        warning('parcel in lh not used %i',i)
    end
    if(length(find(rh == rh_s.table(i, 5))) < 1)
        warning('parcel in rh not used %i',i)
    end
end

% Create parcellation on fsaverage5/6
disp(['Writing annot files for fsaverage5 and fsaverage6, ' num2str(p) ' Parcels, ' num2str(k) ' Networks.']);
write_annotation(fullfile(output_dir, 'FreeSurfer5.3', 'fsaverage5', 'label', ...
    ['lh.Schaefer2018_',num2str(p),'Parcels_',num2str(k) ,'Networks_order.annot']), 0:10241, lh(1:10242), lh_s);
write_annotation(fullfile(output_dir, 'FreeSurfer5.3', 'fsaverage5', 'label', ...
    ['rh.Schaefer2018_',num2str(p),'Parcels_',num2str(k) ,'Networks_order.annot']), 0:10241, rh(1:10242), rh_s);

write_annotation(fullfile(output_dir, 'FreeSurfer5.3', 'fsaverage6', 'label', ...
    ['lh.Schaefer2018_',num2str(p),'Parcels_',num2str(k) ,'Networks_order.annot']), 0:40961, lh(1:40962), lh_s);
write_annotation(fullfile(output_dir, 'FreeSurfer5.3', 'fsaverage6', 'label', ...
    ['rh.Schaefer2018_',num2str(p),'Parcels_',num2str(k) ,'Networks_order.annot']), 0:40961, rh(1:40962), rh_s);

nparcels = p / 2;
lh_new = zeros(size(lh));
rh_new = zeros(size(rh));
for l = 1:nparcels
    lh_new(lh == lh_s.table(l+1, 5)) = l;
    rh_new(rh == rh_s.table(l+1, 5)) = l;
end

% Create parcellation on fsLR
disp(['Writing cifti files, ' num2str(p) ' Parcels, ' num2str(k) ' Networks.']);
system('rm -rf ~/temp123/');
[lh_fslr32k, rh_fslr32k] = CBIG_project_fsaverage2fsLR(lh_new, rh_new, 'fsaverage6', ...
    'label', '~/temp123', '20160827');

% Exclude medial wall vertices defined by HCP
lh_mesh_fslr_32k=CBIG_read_fslr_surface('lh','fs_LR_32k','inflated','medialwall.annot');
rh_mesh_fslr_32k=CBIG_read_fslr_surface('rh','fs_LR_32k','inflated','medialwall.annot');
lh_fslr32k(lh_mesh_fslr_32k.MARS_label == 1) = 0;
rh_fslr32k(rh_mesh_fslr_32k.MARS_label == 1) = 0;

cifti_name = fullfile(output_dir, 'HCP', 'fslr32k', 'cifti', ...
    ['Schaefer2018_',num2str(p),'Parcels_',num2str(k) ,'Networks_order']);
CBIG_gwMRF_write_cifti_from_annot(lh_annot, rh_annot, cifti_name, p/2, lh_fslr32k, rh_fslr32k);

% Create parcellation in MNI
disp(['Writing MNI files, ' num2str(p) ' Parcels, ' num2str(k) ' Networks.']);
mni_name = fullfile(output_dir, 'MNI', ['Schaefer2018_',num2str(p),'Parcels_',num2str(k) ,'Networks_order']);
project_to_MNI(lh_new, rh_new, mni_name);
txt_name = fullfile(output_dir, 'MNI', ['Schaefer2018_',num2str(p),'Parcels_',num2str(k) ,'Networks_order.txt']);
lut_name = fullfile(output_dir, 'MNI', ['Schaefer2018_',num2str(p),'Parcels_',num2str(k) ,'Networks_order_lut.txt']);
cell2csv(txt_name, [num2cell(1:2*nparcels)', [lh_s.struct_names(2:end); rh_s.struct_names(2:end)], ...
    num2cell([lh_s.table(2:nparcels+1, 1:3); rh_s.table(2:nparcels+1, 1:3)]), num2cell(zeros(2*nparcels, 1))], '\t');
cell2csv(lut_name, [num2cell(1:2*nparcels)',...
    num2cell([lh_s.table(2:nparcels+1, 1:3)./255; rh_s.table(2:nparcels+1, 1:3)./255]),...
    [lh_s.struct_names(2:end); rh_s.struct_names(2:end)]], '\t');

%NOTE the command below does not work for fsleyes. We have deprecated it.
%CBIG_gwMRF_create_FSL_LUT(lut_name, [lh_s.table(2:end, 1:3); rh_s.table(2:end, 1:3)]);

end

function project_to_MNI(lh_label, rh_label, filename)

% project_to_MNI(lh_label, rh_label, filename)
% 
% This function projects from fsaverage to MNI using
% CBIG_Projectfsaverage2MNI.
% 
% Input:
%      -lh_label:
%       Vector, label of the left hemisphere on fsaverage or fsaverge6
%
%      -rh_label:
%       Vector, label of the right hemisphere on fsaverage or fsaverge6
%
%      -filename:
%       Path of the output filename
% 
% Example:
% project_to_MNI(lh_label, rh_label, filename)

rh_label(rh_label > 0) = rh_label(rh_label > 0) + max(lh_label);

if(size(lh_label, 1) ~= 1)
   lh_label = lh_label';
end

if(size(rh_label, 1) ~= 1)
    rh_label = rh_label';
end

if (max(size(lh_label)) == 40962)
    lh_mesh7 = CBIG_ReadNCAvgMesh('lh', 'fsaverage', 'sphere', 'cortex');
    rh_mesh7 = CBIG_ReadNCAvgMesh('rh', 'fsaverage', 'sphere', 'cortex');
    lh_mesh6 = CBIG_ReadNCAvgMesh('lh', 'fsaverage6', 'sphere', 'cortex');
    rh_mesh6 = CBIG_ReadNCAvgMesh('rh', 'fsaverage6', 'sphere', 'cortex');
    lh_labels7 = MARS_NNInterpolate_kdTree(lh_mesh7.vertices, lh_mesh6, lh_label);
    rh_labels7 = MARS_NNInterpolate_kdTree(rh_mesh7.vertices, rh_mesh6, rh_label);
else
    lh_labels7 = lh_label;
    rh_labels7 = rh_label;
end

if(nargin == 3)
    output = CBIG_Projectfsaverage2MNI(lh_labels7', rh_labels7');
else
    error('provide correct input arguments')
end

MRIwrite(output, [filename, '.nii.gz']);
MNItemplates_dir = fullfile(getenv('CBIG_CODE_DIR'), 'data', 'templates', 'volume', 'FSL5.0.8_MNI_templates');
system(['mri_vol2vol --targ ', fullfile(MNItemplates_dir, 'MNI152_T1_2mm_brain.nii.gz'), ...
    ' --regheader --mov ',[filename,'.nii.gz'],' --o ',filename,'_FSLMNI152_2mm.nii.gz']);
system(['mri_vol2vol --targ ', fullfile(MNItemplates_dir, 'MNI152_T1_1mm_brain.nii.gz'), ...
    ' --regheader --mov ',[filename,'.nii.gz'],' --o ',filename,'_FSLMNI152_1mm.nii.gz']); 
system(['flirt -applyisoxfm 3 -interp nearestneighbour -in ',filename,'_FSLMNI152_2mm.nii.gz -ref ', ...
    fullfile(MNItemplates_dir, 'MNI152_T1_3mm_brain.nii.gz -out '), filename, '_3mm.nii.gz'])
system(['3dcalc -overwrite -a ', filename,'_3mm.nii.gz -expr ''step(a)''  -prefix ',filename,'_3mm_mask.nii.gz']);
system(['rm -f ', filename,'*.lta']);
system(['rm -f ', filename,'*.reg']);
system(['rm -f ', filename,'*mask*']);
system(['mv ', [filename,'.nii.gz'], ' ', [filename,'_conform.nii.gz']]);

vol=MRIread([filename,'_FSLMNI152_2mm.nii.gz']);
a = hist(reshape(vol.vol, 1, 109*91*91), max(rh_label) + 1);
disp(['Size_of_smallest_parcel = ' num2str(min(a))]);

end
