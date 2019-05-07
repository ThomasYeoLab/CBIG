function [labels, colortable] = CBIG_ConvertSingleHemiSurfaceDataToDiscretizedAnnotation(data, ...
    discretization_res, colorscale, min_thresh, max_thresh)  

% [labels, colortable] = CBIG_ConvertSingleHemiSurfaceDataToDiscretizedAnnotation(data, ..
%   discretization_res, colorscale, min_thresh, max_thresh)  
%
% Convert surface data of one  hemisphere to discretized surface annotation
%
% Input:
%   - data              : original surface data of one hemisphere represented
%                         as a columnvector
%   - discretization_res: number of distinct discrete annotation labels
%   - colorscale        : a colorscale in Matlab's format that would be
%                         discretized to produce the annotation's colortable
%   - min_thresh        : minimum threshold of the output values
%   - max_thresh        : maximum threshold of the output values
% Output:
%   - labels            : annotation labels assigned to the vertices
%   - colortable        : colortable for the annotated labels
%
% Example:
%   [labels, colortable] = CBIG_ConvertSingleHemiSurfaceDataToDiscretizedAnnotation(lh_data, ..
%        28, 'clear_brain', 1e-5, 5e-5);
%   Discretize a single hemisphere's surface data saved in
%   lh_data. The values in the surface data are discretized into
%   28 levels, with the minimum at 1e-5 and maximum
%   thresholded at 5e-5
%
% Written by Gia H. Ngo and CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md
  
    % cap the maximum value in the data at the threshold's upper limit
    if max_thresh < max(data)
        data(data >= max_thresh) = max_thresh;
    end
    if min(data) < min_thresh
        data(data <= min_thresh) = min_thresh;
    end
    
    % discretize non-zero values
    nonzero_values = data(data >= min_thresh);
    value_step = (max_thresh - min_thresh) / (discretization_res-1);
    binranges = min_thresh:value_step:max_thresh;
    [~, bin_ind] = histc(nonzero_values, binranges);

    % generate corresponding colortable for the new labels
    colortable = CBIG_GenerateAnnotationColortable(discretization_res, colorscale);
    
    % underlay the new annotation with gray color saved at the end of the new
    % colortable
    underlay_rgb_sum = colortable.table(end, 5);
    
    labels = underlay_rgb_sum * ones(length(data), 1);
    labels(data >= min_thresh) = colortable.table(bin_ind, 5);
