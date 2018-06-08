function CBIG_VisualizeSurfaceAnnotationInFreeview(lh_annot_file, rh_annot_file, surface_template, label, output_dir)

% CBIG_VisualizeSurfaceAnnotationInFreeview(lh_annot_file, rh_annot_file, surface_template, label, output_dir)
%
% Automatically view the surface annotation of the left and right hemisphere
% in Freeview, take screenshot of the four views (lateral, medial, dorsal, ventral)
% of each hemisphere, and stitch them togther into either a grid or
% serial layout.
% Note that only the inflated brain is used for visualization
%
% Input:
%   - lh_annot_file       : .annot file containing the surface annotation
%                           of the left hemisphere
%   - rh_annot_file       : .annot file containing the surface annotation
%                           of the right hemisphere
%   - surface_template    : name of the surface space of the input
%                           data.`surface_template` needs to correspond to
%                           number of vertices of the annotation
%     surface_template    |      number of vertices
%       fsaverage5        |            10242
%       fsaverage6        |            40962
%       fsaverage         |            163842
%       fslr_32k          |            32492
%       fslr_164k         |            163842
%   - label               : common label of the output images
%   - output_dir          : directory to store the output images
% Output:
%   - Images of the different views of the projected annotation are stored
%     under output_dir. The images are named as '[label]_[hemi]_[view].png',
%     in which [label] is the input argument, [hemi] is either 'lh' or 'rh',
%     [view] is one of the following 'lateral', 'medial', 'dorsal', 'ventral'.
%     A subfolder 'tmp' under output_dir also stores the intermediate files.
%
% Example:
%   CBIG_VisualizeSurfaceAnnotationInFreeview('/data/Work/area8_visualization/area8_lh.annot', ...
%          '/data/Work/area8_visualization/area8_rh.annot', 'fsaverage', 'area8', '/data/Work/area8_visualization/output')
%   The annotation 'area8_lh.annot' and 'area8_rh.annot' in fsaverage surface space will be
%   viewed in FreeView. The output images are stored with the format
%   'area8_[hemi]_[view].png' under /data/Work/area8_visualization/output'
%
% Written by Gia H. Ngo and CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md

    curr_dir = pwd;
    tmp_dir = fullfile(output_dir, 'tmp');
    system(['mkdir -p ' tmp_dir]);
    
    cd(output_dir);

    % colorscale of the figure 
    colorscale_file = fullfile(output_dir, 'colorscale.png');
    system(CBIG_TrimFreeviewScreenshotCommand(colorscale_file, colorscale_file));
    
    % retrieve path to file containing the underlaying brain surface
    if (strcmp(surface_template, 'fsaverage5') || strcmp(surface_template, 'fsaverage6') || strcmp(surface_template, 'fsaverage'))
        lh_underlay_file = fullfile(getenv('FREESURFER_HOME'), 'subjects', surface_template, 'surf', 'lh.inflated');
        rh_underlay_file = fullfile(getenv('FREESURFER_HOME'), 'subjects', surface_template, 'surf', 'rh.inflated');
    elseif (strcmp(surface_template, 'fs_LR_32k') || strcmp(surface_template, 'fs_LR_164k'))
        lh_underlay_file = fullfile(getenv('CBIG_CODE_DIR'), 'data', 'templates', 'surface', surface_template, 'surf', 'lh.inflated');
        rh_underlay_file = fullfile(getenv('CBIG_CODE_DIR'), 'data', 'templates', 'surface', surface_template, 'surf', 'rh.inflated');
    else
        disp(['Invalid surface template: ' surface_template]);
        return;
    end

    % LEFT VIEWS 
    % left lateral view 
    raw_lh_lateral_file = fullfile(tmp_dir, 'lh_lateral.raw.png');
    lh_lateral_file = fullfile(output_dir, [label '_lh_lateral.png']);
    tmp_lh_lateral_file = fullfile(tmp_dir, 'lh_lateral.tmp.png');
    system(CBIG_ViewLHAnnotationInFreeview(lh_annot_file, lh_underlay_file, 0, 0, raw_lh_lateral_file));
    system(CBIG_ReplaceBlackByWhiteBkgCommand(raw_lh_lateral_file, tmp_lh_lateral_file));

    % left medial view 
    raw_lh_medial_file = fullfile(tmp_dir, 'lh_medial.raw.png');
    lh_medial_file = fullfile(output_dir, [label '_lh_medial.png']);
    tmp_lh_medial_file = fullfile(tmp_dir, 'lh_medial.tmp.png');
    system(CBIG_ViewLHAnnotationInFreeview(lh_annot_file, lh_underlay_file, 180, 0, raw_lh_medial_file));
    system(CBIG_ReplaceBlackByWhiteBkgCommand(raw_lh_medial_file, tmp_lh_medial_file));

    % left dorsal view 
    raw_lh_dorsal_file = fullfile(tmp_dir, 'lh_dorsal.raw.png');
    lh_dorsal_file = fullfile(output_dir, [label '_lh_dorsal.png']);
    tmp_lh_dorsal_file = fullfile(tmp_dir, 'lh_dorsal.tmp.png');
    system(CBIG_ViewLHAnnotationInFreeview(lh_annot_file, lh_underlay_file, 180, 90, raw_lh_dorsal_file));
    system(CBIG_ReplaceBlackByWhiteBkgCommand(raw_lh_dorsal_file, tmp_lh_dorsal_file));

    % left ventral view
    raw_lh_ventral_file = fullfile(tmp_dir, 'lh_ventral.raw.png');
    lh_ventral_file = fullfile(output_dir, [label '_lh_ventral.png']);
    tmp_lh_ventral_file = fullfile(tmp_dir, 'lh_ventral.tmp.png');
    system(CBIG_ViewLHAnnotationInFreeview(lh_annot_file, lh_underlay_file, 180, -90, raw_lh_ventral_file));
    system(CBIG_ReplaceBlackByWhiteBkgCommand(raw_lh_ventral_file, tmp_lh_ventral_file));

   
    % RIGHT VIEWS 
    % right lateral view
    raw_rh_lateral_file = fullfile(tmp_dir, 'rh_lateral.raw.png');
    rh_lateral_file = fullfile(output_dir, [label '_rh_lateral.png']);
    tmp_rh_lateral_file = fullfile(tmp_dir, 'rh_lateral.tmp.png');
    system(CBIG_ViewRHAnnotationInFreeview(rh_annot_file, rh_underlay_file, 180, 0, raw_rh_lateral_file));
    system(CBIG_ReplaceBlackByWhiteBkgCommand(raw_rh_lateral_file, tmp_rh_lateral_file));

    % right medial view
    raw_rh_medial_file = fullfile(tmp_dir, 'rh_medial.raw.png');
    rh_medial_file = fullfile(output_dir, [label '_rh_medial.png']);
    tmp_rh_medial_file = fullfile(tmp_dir, 'rh_medial.tmp.png');
    system(CBIG_ViewRHAnnotationInFreeview(rh_annot_file, rh_underlay_file, 0, 0, raw_rh_medial_file));
    system(CBIG_ReplaceBlackByWhiteBkgCommand(raw_rh_medial_file, tmp_rh_medial_file));

    % right dorsal view 
    raw_rh_dorsal_file = fullfile(tmp_dir, 'rh_dorsal.raw.png');
    rh_dorsal_file = fullfile(output_dir, [label '_rh_dorsal.png']);
    tmp_rh_dorsal_file = fullfile(tmp_dir, 'rh_dorsal.tmp.png');
    system(CBIG_ViewRHAnnotationInFreeview(rh_annot_file, rh_underlay_file, 180, 90, raw_rh_dorsal_file));
    system(CBIG_ReplaceBlackByWhiteBkgCommand(raw_rh_dorsal_file, tmp_rh_dorsal_file));

    % right ventral view
    raw_rh_ventral_file = fullfile(tmp_dir, 'rh_ventral.raw.png');
    rh_ventral_file = fullfile(output_dir, [label '_rh_ventral.png']);
    tmp_rh_ventral_file = fullfile(tmp_dir, 'rh_ventral.tmp.png');
    system(CBIG_ViewRHAnnotationInFreeview(rh_annot_file, rh_underlay_file, 180, -90, raw_rh_ventral_file));
    system(CBIG_ReplaceBlackByWhiteBkgCommand(raw_rh_ventral_file, tmp_rh_ventral_file));


    % trim the borders
    system(CBIG_TrimFreeviewScreenshotCommand(tmp_lh_lateral_file, lh_lateral_file));
    system(CBIG_TrimFreeviewScreenshotCommand(tmp_lh_medial_file, lh_medial_file));
    system(CBIG_TrimFreeviewScreenshotCommand(tmp_lh_dorsal_file, lh_dorsal_file));
    system(CBIG_TrimFreeviewScreenshotCommand(tmp_lh_ventral_file, lh_ventral_file));
    
    system(CBIG_TrimFreeviewScreenshotCommand(tmp_rh_lateral_file, rh_lateral_file));
    system(CBIG_TrimFreeviewScreenshotCommand(tmp_rh_medial_file, rh_medial_file));
    system(CBIG_TrimFreeviewScreenshotCommand(tmp_rh_dorsal_file, rh_dorsal_file));
    system(CBIG_TrimFreeviewScreenshotCommand(tmp_rh_ventral_file, rh_ventral_file));

    % pad the views with borders to align the views
    dorsal_file = fullfile(output_dir, [label '_dorsal.png']);
    ventral_file = fullfile(output_dir, [label '_ventral.png']);
    tmp_dorsal_file = fullfile(tmp_dir, 'dorsal.tmp.png');
    tmp_ventral_file = fullfile(tmp_dir, 'ventral.tmp.png');
    
    system(CBIG_PadImage(lh_dorsal_file, lh_dorsal_file, 0, 5));
    system(CBIG_PadImage(lh_ventral_file, lh_ventral_file, 0, 5));
    system(CBIG_PadImage(rh_dorsal_file, rh_dorsal_file, 0, 5));
    system(CBIG_PadImage(rh_ventral_file, rh_ventral_file, 0, 5));
    system(CBIG_CombineImagesVertically(lh_dorsal_file, rh_dorsal_file, dorsal_file));
    system(CBIG_CombineImagesVertically(rh_ventral_file, lh_ventral_file, ventral_file));

    % GRID LAYOUT 
    % stich the views to a grid layout
    tmp_grid_output_path = fullfile(tmp_dir, [label '.grid.tmp.png']);
    grid_output_path = fullfile(output_dir, [label '.grid.png']);
    
    [lm_w, ~] = CBIG_GetMinimumDimensionsOfImages({lh_lateral_file, lh_medial_file, rh_lateral_file, rh_medial_file});
    
    CBIG_ResizeImage(lh_lateral_file, lm_w, 0, tmp_lh_lateral_file);
    CBIG_ResizeImage(lh_medial_file, lm_w, 0, tmp_lh_medial_file);
    
    CBIG_ResizeImage(rh_lateral_file, lm_w, 0, tmp_rh_lateral_file);
    CBIG_ResizeImage(rh_medial_file, lm_w, 0, tmp_rh_medial_file);

    [~, lm_h] = CBIG_GetMinimumDimensionsOfImages({lh_lateral_file, lh_medial_file, rh_lateral_file, rh_medial_file});
    CBIG_ResizeImage(dorsal_file, 0, lm_h, tmp_dorsal_file);
    CBIG_ResizeImage(ventral_file, 0, lm_h, tmp_ventral_file);
    
    [dv_w, dv_h] = CBIG_GetMinimumDimensionsOfImages({tmp_dorsal_file, tmp_ventral_file});
    CBIG_ResizeImage(tmp_dorsal_file, dv_w, 0, tmp_dorsal_file);
    CBIG_ResizeImage(tmp_ventral_file, dv_w, 0, tmp_ventral_file);
    
    
    system(CBIG_PadImage(tmp_lh_lateral_file, tmp_lh_lateral_file));
    system(CBIG_PadImage(tmp_lh_medial_file, tmp_lh_medial_file));
    system(CBIG_PadImage(tmp_rh_lateral_file, tmp_rh_lateral_file));
    system(CBIG_PadImage(tmp_rh_medial_file, tmp_rh_medial_file));
    system(CBIG_PadImage(tmp_dorsal_file, tmp_dorsal_file));
    system(CBIG_PadImage(tmp_ventral_file, tmp_ventral_file));
    
    % Set width of the color scale to 40% of one of the views
    tmp_colorscale_file = fullfile(tmp_dir, 'colorscale.tmp.png');
    colorscale_width = int16(lm_w * 0.4);
    CBIG_ResizeImage(colorscale_file, colorscale_width, 0, tmp_colorscale_file);
    
    system(CBIG_StitchBrainScreenshotsIntoGrid(tmp_lh_lateral_file, tmp_lh_medial_file, ...
        tmp_rh_lateral_file, tmp_rh_medial_file, ...
        tmp_dorsal_file, tmp_ventral_file, tmp_grid_output_path));
    system(CBIG_InsertColorscaleToImage(tmp_colorscale_file, tmp_grid_output_path, grid_output_path));

    % SERIAL LAYOUT 
    % stich the views to a serial layout
    tmp_series_output_path = fullfile(tmp_dir, [label 'series.tmp.png']);
    series_output_path = fullfile(output_dir, [label '.series.png']);

    [~, lm_h] = CBIG_GetMinimumDimensionsOfImages({lh_lateral_file, lh_medial_file, ...
        rh_lateral_file, rh_medial_file, ...
        dorsal_file, ventral_file});
    CBIG_ResizeImage(lh_lateral_file, 0, lm_h, tmp_lh_lateral_file);
    CBIG_ResizeImage(lh_medial_file, 0, lm_h, tmp_lh_medial_file);
    CBIG_ResizeImage(rh_lateral_file, 0, lm_h, tmp_rh_lateral_file);
    CBIG_ResizeImage(rh_medial_file, 0, lm_h, tmp_rh_medial_file);
    CBIG_ResizeImage(dorsal_file, 0, lm_h, tmp_dorsal_file);
    CBIG_ResizeImage(ventral_file, 0, lm_h, tmp_ventral_file);
    
    system(CBIG_PadImage(tmp_lh_lateral_file, [tmp_lh_lateral_file '.pad'], 5, 0));
    system(CBIG_PadImage(tmp_lh_medial_file, tmp_lh_medial_file, 5, 0));
    system(CBIG_PadImage(tmp_rh_lateral_file, tmp_rh_lateral_file, 5, 0));
    system(CBIG_PadImage(tmp_rh_medial_file, tmp_rh_medial_file, 5, 0));
    system(CBIG_PadImage(tmp_dorsal_file, tmp_dorsal_file, 5, 0));
    system(CBIG_PadImage(tmp_ventral_file, tmp_ventral_file, 5, 0));
    
     % Set width of the color scale to 40% of one of the views
    colorscale_width = int16(lm_w * 0.4);
    CBIG_ResizeImage(colorscale_file, colorscale_width, 0, tmp_colorscale_file);
    
    system(CBIG_StitchBrainScreenshotsIntoSeries(tmp_lh_lateral_file, tmp_lh_medial_file, ...
        tmp_rh_lateral_file, tmp_rh_medial_file, ...
        tmp_dorsal_file, tmp_ventral_file, tmp_series_output_path));
    system(CBIG_InsertColorscaleToImage(tmp_colorscale_file, tmp_series_output_path, series_output_path, 0.58, 0.49, 0.9));
    
    cd(curr_dir);
