function CBIG_Yeo2011_ProjectSplitLabels2MNI1mm

% CBIG_Yeo2011_ProjectSplitLabels2MNI1mm
% Written by CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md

cortex_mask  = MRIread(fullfile(getenv('CBIG_CODE_DIR'), 'data', 'extras', 'surf2surf_gui_data', 'CorrespondenceFreeSurferVolSurfSpace', 'coord_vol2surf', 'MNI_cortex_estimate.150.nii.gz'));
vol2surf_map = MRIread(fullfile(getenv('CBIG_CODE_DIR'), 'data', 'extras', 'surf2surf_gui_data', 'CorrespondenceFreeSurferVolSurfSpace', 'coord_vol2surf', '1000sub.FSL_MNI152.1mm.full_vertex_map.500.fsaverage5.nii.gz'));

lh_parms.SUBJECTS_DIR = fullfile(pwd, '..');
lh_parms.hemi = 'lh';
lh_parms.SUBJECT = 'fsaverage5';
lh_parms.read_surface = @MARS2_readSbjMesh;
lh_parms.radius = 100;
lh_parms.unfoldBool = 0;
lh_parms.flipFacesBool = 1;
lh_parms.surf_filename = 'inflated';
lh_parms.metric_surf_filename = 'inflated';
lh_parms.data_filename_cell = {'sulc', 'curv'};

rh_parms.SUBJECTS_DIR = fullfile(pwd, '..');
rh_parms.hemi = 'rh';
rh_parms.SUBJECT = 'fsaverage5';
rh_parms.read_surface = @MARS2_readSbjMesh;
rh_parms.radius = 100;
rh_parms.unfoldBool = 0;
rh_parms.flipFacesBool = 1;
rh_parms.surf_filename = 'inflated';
rh_parms.metric_surf_filename = 'inflated';
rh_parms.data_filename_cell = {'sulc', 'curv'};

mask_index = find(cortex_mask.vol ~= 0);
for cluster = [7 17]
    
    lh_parms.annot_filename = ['Yeo2011_' num2str(cluster) 'Networks_N1000.split_components.annot'];
    lh_avg_mesh = MARS2_readSbjMesh(lh_parms);
    lh_avg_mesh.vertices = lh_avg_mesh.metricVerts;
    
    rh_parms.annot_filename = ['Yeo2011_' num2str(cluster) 'Networks_N1000.split_components.annot'];
    rh_avg_mesh = MARS2_readSbjMesh(rh_parms);
    rh_avg_mesh.vertices = rh_avg_mesh.metricVerts;
    
    
    output = cortex_mask;
    for i = mask_index'
        vertex_id = vol2surf_map.vol(i);
        
        if(vertex_id > 0)
            output.vol(i) = lh_avg_mesh.MARS_label(vertex_id) - 1;
        else
            output.vol(i) = rh_avg_mesh.MARS_label(abs(vertex_id)) - 1;
        end
    end
    
    MRIwrite(output, ['../MNI152/' 'Yeo2011_' num2str(cluster) 'Networks_N1000.split_components.FSL_MNI152_FreeSurferConformed_1mm.nii.gz']);
    
    for i = 2:length(lh_avg_mesh.MARS_ct.struct_names)
       disp([num2str(i) ': ' num2str(sum(output.vol(:) == (i-1)))]); 
    end
    
    % Write freesurfer colormap
    fid = fopen(['../MNI152/' num2str(cluster) 'Networks_ColorLUT_freeview.txt'], 'w');
    fprintf(fid, '%3d %40s %3d %3d %3d   0\n', 0, 'NONE', 0, 0, 0);
    for i = 2:length(lh_avg_mesh.MARS_ct.struct_names)
        fprintf(fid, '%3d %40s %3d %3d %3d   0\n', i - 1, lh_avg_mesh.MARS_ct.struct_names{i}, lh_avg_mesh.MARS_ct.table(i, 1), lh_avg_mesh.MARS_ct.table(i, 2), lh_avg_mesh.MARS_ct.table(i, 3));
    end
    fclose(fid);
    
    
    % Write fslview colormap
    fid = fopen(['../MNI152/' num2str(cluster) 'Networks_ColorLUT_fslview.lut'], 'w');
    fprintf(fid, '%%!VEST-LUT\n');
    fprintf(fid, '%%%%BeginInstance\n');
    fprintf(fid, '<<\n');
    fprintf(fid, '/SavedInstanceClassName /ClassLUT \n');
    fprintf(fid, '/PseudoColorMinimum 0.00 \n');
    fprintf(fid, '/PseudoColorMaximum 1.00 \n');
    fprintf(fid, '/PseudoColorMinControl /Low \n');
    fprintf(fid, '/PseudoColorMaxControl /High \n');
    fprintf(fid, '/PseudoColormap [\n');
    
    for i = 1:length(lh_avg_mesh.MARS_ct.struct_names)
        fprintf(fid, '<-color{%.6f,%.6f,%.6f}->\n', double(lh_avg_mesh.MARS_ct.table(i, 1))/255, double(lh_avg_mesh.MARS_ct.table(i, 2))/255, double(lh_avg_mesh.MARS_ct.table(i, 3))/255);
    end
    
    fprintf(fid, ']\n');
    fprintf(fid, '>>\n');
    fprintf(fid, '\n');
    fprintf(fid, '%%%%EndInstance\n');
    fprintf(fid, '%%%%EOF\n');
    fclose(fid);
    
end


