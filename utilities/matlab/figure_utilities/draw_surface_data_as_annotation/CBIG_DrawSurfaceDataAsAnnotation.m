function CBIG_DrawSurfaceDataAsAnnotation(lh_data, rh_data, ...
    abs_path_to_lh_ref_annot, abs_path_to_rh_ref_annot, ref_medialwall_label, ...
    surface_template, abs_path_to_output_dir, label, ...
    colorscheme, discretization_res, min_thresh, max_thresh)

% CBIG_DrawSurfaceDataAsAnnotation(lh_data, rh_data, ...
%   abs_path_to_lh_ref_annot, abs_path_to_rh_ref_annot, ref_medialwall_label, ...
%   surface_template, abs_path_to_output_dir, label, ...
%   colorscheme, discretization_res, min_thresh, max_thresh)
%
% Wrapper function to draw surface data as surface annotation and capture screenshots of their visualization in FreeView
%
% Input:
%   Compulsory:
%     - lh_data              : column or row vector containing surface data of
%                              the left hemisphere.
%     - rh_data              : column or row vector containing surface data of
%                              the right hemisphere.
%     - abs_path_to_lh_ref_annot: absolute path to the reference
%                                 annotation file for the left hemisphere.
%     - abs_path_to_rh_ref_annot: absolute path to the reference
%                                 annotation file for the right hemisphere
%                                 The reference annotation is used to
%                                 mark the medial wall in the
%                                 visualization. The reference
%                                 annotations needs to be in the same surface
%                                 space as the input data.
%     - ref_medialwall_label : label assigned to medial wall in the reference
%                                 annotation files
%     - surface_template     : name of the surface space of the input
%                              data.`surface_template` needs to correspond to
%                              number of vertices of the annotation in
%                              `lh_data`, `rh_data` and the reference annotations.
%        surface_template    |      number of vertices
%          fsaverage5        |            10242
%          fsaverage6        |            40962
%          fsaverage         |            163842
%          fs_LR_32k         |            32492
%          fs_LR_164k        |            163842
%     - abs_path_to_output_dir: absolute path to the output folder.
%     - label                : common label of the output files.
%   Optional:
%     - colorscheme          : string specifying the colorscheme.
%                              'clear_brain' (default): Human Connectome
%                               Workbench's clear_brain color scheme
%                              'parula': Matlab's parula color scheme
%                              'hsv': Matlab's HSV color scheme
%     - discretization_res   : number of discrete intensity levels.
%                              Default: 28
%     - min_thresh           : minimum threshold of the visualized data.
%                              Values below min_thresh are excluded.
%                              Default: 1e-5
%     - max_thresh           : maximum threshold of the visualized data.
%                              Values above max_thresh are excluded.
%                              Default: 1
%
% Example:
%   abs_path_to_lh_annot_file = fullfile(getenv('FREESURFER_HOME'), ...
%                  'subjects', 'fsaverage5', 'label', 'lh.aparc.annot');
%   abs_path_to_rh_annot_file = fullfile(getenv('FREESURFER_HOME'), ...
%                  'subjects', 'fsaverage5', 'label', 'rh.aparc.annot');
%   ref_medial_wall_label = 0;
%   CBIG_DrawSurfaceDataAsAnnotation(lh_data, rh_data, 'fsaverage5', ...
%           abs_path_to_lh_annot_file, abs_path_to_rh_annot_file,
%           ref_medial_wall_label, '/data/Work/area8_visualization', 'area8')
%   Draw fsaverage5 surface data saved in lh_data and rh_data and save
%   the resulting images under /data/Work/area8_visualization with the 
%   label area8 prefixed in all file names. The reference annotation
%   files are the default fsaverage5 aparc.annot parcellation, with the
%   medialwall labelled as 0. The values in the resulting images are
%   discretized into 28 levels, with the minimum at 1e-5 and maximum
%   thresholded at 1
%
% Written by Gia H. Ngo and CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md

% This function does not need vector check because the function itself
% contains checking statement.

    if nargin < 12
        max_thresh = 1;
    else
        if ischar(max_thresh)
            max_thresh = str2num(max_thresh);
        end
    end

    if nargin < 11
        min_thresh = 1e-5;
    else
        if ischar(min_thresh)
            min_thresh = str2num(min_thresh);
        end
    end
    
    if nargin < 10
        discretization_res = 28;
    else
        if ischar(discretization_res)
            discretization_res = str2num(discretization_res);
        end
    end
    
    if nargin < 9
        colorscale = CBIG_GenerateClearbrainColorscale(discretization_res, min_thresh, max_thresh, abs_path_to_output_dir);
    elseif ~isnumeric(colorscheme)
        if strcmp(colorscheme, 'truncated')
            colorscale = CBIG_GenerateTruncatedClearbrainColorscale(discretization_res, min_thresh, max_thresh, abs_path_to_output_dir);
        elseif strcmp(colorscheme, 'hsv')
            colorscale = CBIG_GenerateHSVColorscale(discretization_res, min_thresh, max_thresh, abs_path_to_output_dir);
        elseif strcmp(colorscheme, 'parula')
            colorscale = CBIG_GenerateParulaColorscale(discretization_res, min_thresh, max_thresh, abs_path_to_output_dir);
        elseif strcmp(colorscheme, 'jet')
            colorscale = CBIG_GenerateJetColorscale(discretization_res, min_thresh, max_thresh, abs_path_to_output_dir);
        elseif strcmp(colorscheme, 'clear_brain') || strcmp(colorscheme, 'default')
            colorscale = CBIG_GenerateClearbrainColorscale(discretization_res, min_thresh, max_thresh, abs_path_to_output_dir);
        end
    end

    tmp_dir = fullfile(abs_path_to_output_dir, 'tmp');
    system(['mkdir -p ' tmp_dir]);
    
    if size(lh_data, 1) == 1
        lh_data = lh_data';
    end
    
    if size(rh_data, 1) == 1
        rh_data = rh_data';
    end
    
    % read reference annotations
    [lh_vertices, ref_lh_labels, ~] = read_annotation(abs_path_to_lh_ref_annot);
    [rh_vertices, ref_rh_labels, ~] = read_annotation(abs_path_to_rh_ref_annot);
    
    % check inputs
    if ~CBIG_CheckValidSurfaceData(surface_template, lh_data, rh_data, lh_vertices, rh_vertices)
        disp('Mismatched number of vertices in the input data, reference annotation and surface template!');
        return;
    end

    % convert the input surface data to discretized annotation
    [lh_labels, lh_colortable] = CBIG_ConvertSingleHemiSurfaceDataToDiscretizedAnnotation(lh_data, ...
        discretization_res, colorscale, min_thresh, max_thresh);
    [rh_labels, rh_colortable] = CBIG_ConvertSingleHemiSurfaceDataToDiscretizedAnnotation(rh_data, ...
        discretization_res, colorscale, min_thresh, max_thresh);
    % mark the medial wall on the final annotation. The medial wall is delineated
    % in the reference annotation
    lh_final_annot_file = fullfile(tmp_dir, ['lh.' label '.annot']);
    rh_final_annot_file = fullfile(tmp_dir, ['rh.' label '.annot']);
    
    CBIG_AnnotateSingleHemiMedialWall(lh_vertices, lh_labels, lh_colortable, ref_lh_labels, ref_medialwall_label, lh_final_annot_file);
    CBIG_AnnotateSingleHemiMedialWall(rh_vertices, rh_labels, rh_colortable, ref_rh_labels, ref_medialwall_label, rh_final_annot_file);
    
    % automatically view the annotation in FreeView and capture the
    % screenshots
    CBIG_VisualizeSurfaceAnnotationInFreeview(lh_final_annot_file, rh_final_annot_file, ...
        surface_template, label, abs_path_to_output_dir);
    
    close all;
