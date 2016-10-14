function dice = ComputeDiceBetweenLabels(surface, orig_surface, label1, label2)

% e.g. ComputeDiceBetweenLabels('<SUBJECTS_DIR>/<SUBJECT>/surf/lh.sphere', '<SUBJECTS_DIR>/<SUBJECT>/surf/lh.white', 'lateraloccipital_gt', 'lateraloccipital')
% e.g. ComputeDiceBetweenLabels('<SUBJECTS_DIR>/<SUBJECT>/surf/lh.sphere', '<SUBJECTS_DIR>/<SUBJECT>/surf/lh.white', '<dir_name>/lh.lateraloccipital_gt.label', '<dir_name>/lh.lateraloccipital.label')

indices = strfind(surface, filesep);
hemi = surface(indices(end)+1:indices(end)+2);
SUBJECTS_DIR = surface(1:indices(end-2)-1);
subject = surface(indices(end-2)+1:indices(end-1)-1);
surf_dir = surface(indices(end-1)+1:indices(end)-1);
if(strcmp(surf_dir, 'surf') == 0)
   error(['ComputeDiceBetweenLabels assumes surfaces are in the "surf" directory of the subject, not ' surf_dir]);
end

sbj_parms.hemi = hemi;
sbj_parms.SUBJECTS_DIR = SUBJECTS_DIR;
sbj_parms.SUBJECT = subject;
sbj_parms.surf_filename = surface(indices(end)+4:end);

indices = strfind(orig_surface, filesep);
sbj_parms.metric_surf_filename = orig_surface(indices(end)+4:end);


if(isempty(strfind(label1, filesep)))
    sbj_parms.label_filename_cell = {label1};
    sbjMesh1 = MARS2_readSbjMesh(sbj_parms);
else
    sbjMesh1 = MARS2_readSbjMesh(sbj_parms);
    
    MARS_label = ones(1, size(sbjMesh1.vertices, 2));
    l = read_label([], label1);
    MARS_label(l(:, 1)+1) = 2;
    MARS_ct.numEntries = int32(2);
    
    sbjMesh1.MARS_ct = MARS_ct;
    sbjMesh1.MARS_label = int32(MARS_label);
end

    
if(isempty(strfind(label2, filesep)))    
    sbj_parms.label_filename_cell = {label2};
    sbjMesh2 = MARS2_readSbjMesh(sbj_parms);
else
    sbjMesh2 = MARS2_readSbjMesh(sbj_parms);
    
    MARS_label = ones(1, size(sbjMesh2.vertices, 2));
    l = read_label([], label2);
    MARS_label(l(:, 1)+1) = 2;
    MARS_ct.numEntries = int32(2);
    
    sbjMesh2.MARS_ct = MARS_ct;
    sbjMesh2.MARS_label = int32(MARS_label);
end


sbjMesh1.metricVerts = sbjMesh1.metricVerts/sbjMesh1.surface_scaling_factor;
weights = MARS_computeWeightsOfMeshVertices(sbjMesh1, 'metric');
[dice_vec] = MARS_computeDice(sbjMesh1, sbjMesh2.MARS_label, 1, weights);
dice = dice_vec(2);

