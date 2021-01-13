function vol2surf_fc = CBIG_IndCBM_compute_vol2surf_fc(lh_surf_file, rh_surf_file, vol_file, template_file, output_file)

% vol2surf_fc = CBIG_IndCBM_compute_vol2surf_fc(lh_surf_file, rh_surf_file, vol_file, template_file, output_file)
%
% This function computes the functional connectivity between the cerebellum
% in the volume and the cerebral cortex on the surface. The resolution
% should match your created cifti template. For template details see
% CBIG_IndCBM_create_template.sh.
%
% Input:
%
%     - lh_surf_file:
%           One single run: path of the lh surface time series (on surface)
%           Example: CBIG_CODE_DIR/data/example_data/CoRR_HNU/subj01/ ...
%           subj01_sess1/surf/lh.subj01_sess1_bld002_rest_skip4_stc_mc_ ...
%           residc_interp_FDRMS0.2_DVARS50_bp_0.009_0.08_fs6_sm6_fs5.nii.gz
%           Multiple runs: text file including all runs. Each row is the 
%           file path of a run.
%
%     - rh_surf_file:
%           One single run: path of the rh surface time series (on surface)
%           Example: CBIG_CODE_DIR/data/example_data/CoRR_HNU/subj01/ ...
%           subj01_sess1/surf/rh.subj01_sess1_bld002_rest_skip4_stc_mc_ ...
%           residc_interp_FDRMS0.2_DVARS50_bp_0.009_0.08_fs6_sm6_fs5.nii.gz
%           Multiple runs: text file including all runs. Each row is the 
%           file path of a run.
%
%     - vol_file:
%           One single run: path of the volume time series. 
%           Example: CBIG_CODE_DIR/stable_projects/brain_parcellation/ ...
%           Xue2021_IndCerebellum/examples/input/vol/sub1/ ...
%           sub1_sess1_vol_4mm.nii.gz
%           Multiple runs: text file including all runs. Each row is the 
%           file path of a run.
%
%     - template_file: 
%           Cifti template dscalar file specifying the cerebellar voxels.
%           Please create this file before you run this function. See 
%           CBIG_IndCBM_create_template.sh to create this template file. 
%           Example: ./examples/example_files/Sub1_fsaverage5_cerebellum_template.dscalar.nii
%
% Optional input:
%
%     - output_file: 
%           Path of output file. If given, 'vol2surf_fc' will be saved.
%
% Output:
%
%     - vol2surf_fc:
%           M x N functional connectivity matrix. 
%           M: cerebellar voxels, same order with the cifti template.
%           N: cerebral cortical vertices, lh and rh concatenated.
%
% Example:
% vol2surf_fc = CBIG_IndCBM_compute_vol2surf_fc('proj/lh_run1.nii.gz', 'proj/rh_run1.nii.gz', ...
% 'proj/vol_run1.nii.gz', 'proj/template.dscalar.nii')
% CBIG_IndCBM_compute_vol2surf_fc('proj/lh_list.txt', 'proj/rh_list.txt', 'proj/vol_list.txt', ...
% 'proj/template.dscalar.nii', 'proj/sub1_fc')
%
% Written by XUE Aihuiping and CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md

if(contains(vol_file,'.nii')) % Single run
    vol_list{1} = vol_file;
    lh_surf_list{1} = lh_surf_file;
    rh_surf_list{1} = rh_surf_file;
else % Multiple runs
    fid = fopen(vol_file);
    vol_list = textscan(fid, '%s');
    vol_list = vol_list{1};
    fclose(fid);
    fid = fopen(lh_surf_file);
    lh_surf_list = textscan(fid, '%s');
    lh_surf_list = lh_surf_list{1};
    fclose(fid);
    fid = fopen(rh_surf_file);
    rh_surf_list = textscan(fid, '%s');
    rh_surf_list = rh_surf_list{1};
    fclose(fid);
end

N = length(vol_list);
if(~isequal(N, length(lh_surf_list)) || ~isequal(N, length(rh_surf_list)))
    error('List lengths for vol, lh and rh are not consistent.');
end
disp([num2str(N) ' runs in total.']);

% Find cerebellum structure
cifti = ft_read_cifti(template_file);
for i = 1:length(cifti.brainstructurelabel)
    if(strcmp(cifti.brainstructurelabel{i}, 'CEREBELLUM_LEFT'))
        cbm(1) = i;
    elseif(strcmp(cifti.brainstructurelabel{i}, 'CEREBELLUM_RIGHT'))
        cbm(2) = i;
    end
end

% Read ras information from template and convert to matlab index
vol = MRIread(vol_list{1});
lh_cbm_ras = cifti.pos(cifti.brainstructure == cbm(1), :);
lh_cbm_mask = zeros(size(vol.vol,  2), size(vol.vol,1), size(vol.vol, 3));
v_num = length(lh_cbm_ras);
for i = 1:v_num
    ras = lh_cbm_ras(i,:)';
    vox = CBIG_ConvertRas2Vox(ras, vol.vox2ras);
    vox = ceil(vox([2 1 3]));
    lh_cbm_mask(vox(2), vox(1), vox(3)) = cbm(1);
end
rh_cbm_ras = cifti.pos(cifti.brainstructure == cbm(2), :);
rh_cbm_mask = zeros(size(vol.vol, 2), size(vol.vol, 1), size(vol.vol, 3));
v_num = length(rh_cbm_ras);
for i = 1:v_num
    ras = rh_cbm_ras(i,:)';
    vox = CBIG_ConvertRas2Vox(ras, vol.vox2ras);
    vox = ceil(vox([2 1 3]));
    rh_cbm_mask(vox(2), vox(1), vox(3)) = cbm(2);
end
lh_vol_index = find(lh_cbm_mask == cbm(1));
rh_vol_index = find(rh_cbm_mask == cbm(2));
vol_index = [lh_vol_index; rh_vol_index];
disp(['Left hemisphere of cerebellum: ' num2str(length(lh_vol_index)) ' voxels.']);
disp(['Right hemisphere of cerebellum: ' num2str(length(rh_vol_index)) ' voxels.']);

disp(['Computing vol to surf... ']);
for i =1:N
    disp(['Run: ' num2str(i)]);
    disp(['Reading lh surface data: ' lh_surf_list{i}]);
    lh = MRIread(lh_surf_list{i});
    lh_data = reshape(lh.vol, lh.nvoxels, lh.nframes);
    disp(['Reading rh surface data: ' rh_surf_list{i}]);
    rh = MRIread(rh_surf_list{i});
    rh_data = reshape(rh.vol, rh.nvoxels, rh.nframes);
    surf_data = [lh_data; rh_data]';
    clear lh rh lh_data rh_data
    disp(['Reading volume data: ' vol_list{i}]);
    if(i ~= 1)
        vol = MRIread(vol_list{i});
    end
    vol_data = permute(vol.vol,  [2,1,3,4]);
    vol_data = reshape(vol_data, vol.nvoxels, vol.nframes);
    clear vol
    vol_data = vol_data(vol_index, :);
    vol_data = vol_data';
    disp('Computing correlation...');
    r_vol2surf = single(CBIG_corr(vol_data, surf_data));
    clear surf_data vol_data
    z_vol2surf = CBIG_StableAtanh(r_vol2surf);
    clear r_vol2surf
    if(i == 1)
        vol2surf_fc = z_vol2surf;
    else
        vol2surf_fc = vol2surf_fc + z_vol2surf;
    end
    clear z_vol2surf
end
vol2surf_fc = vol2surf_fc / N;

% Save 'vol2surf' if 'output_name' is given
if(nargin==5)
    disp(['Saving correlation: ' output_file]);
    save(output_file, 'vol2surf_fc', '-v7.3');
    disp(['Correlation saved at: ' output_file]);
end

end
