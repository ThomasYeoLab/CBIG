function CBIG_TRBPC_subcortical_heatmap(path_data, lim_cog, lim_pers, lim_mh)
% CBIG_TRBPC_subcortical_heatmap(path_data, lim_cog, lim_pers, lim_mh)
%
% This function will make heatmaps for the subcortical ROIs
%
% Required inputs:
% - path_data: a string to a directory that contains a .mat file called 
%         'vec_relevance_sum_conjunction_subcortical.mat'; This directory
%         will also serve as the output directory
% - lim_cog: a vector of two numbers specifying the minimum and maximum values
%         for the heatmap for the cognitive cluster
% - lim_pers: a vector of two numbers specifying the minimum and maximum values
%         for the heatmap for the personality cluster
% - lim_mh: a vector of two numbers specifying the minimum and maximum values
%         for the heatmap for the mental health cluster
%
% Written by CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md

%% load data
path_mat = [path_data filesep 'vec_relevance_sum_conjunction_subcortical.mat'];
load(path_mat)

%% output directory
path_out = [path_data filesep 'subcortical_heatmap'];
if ~exist(path_out), mkdir(path_out); end

%% set variables
clus = {'cog','pers','ment_health'};
vals = {'pos','neg'};

% labels in original order
roi_labels_o = {'Cerebellum_Left';'Thalamus_Left';'Caudate_Left';'Putamen_Left';...
    'Pallidum_Left';'BrainStem';'Hippocampus_Left';'Amygdala_Left';'Accumbens_Left'; ...
    'DiencephalonVentral_Left';'Cerebellum_Right';'Thalamus_Right';'Caudate_Right';...
    'Putamen_Right';'Pallidum_Right';'Hippocampus_Right';'Amygdala_Right';...
    'Accumbens_Right';'DiencephalonVentral_Right'};

% labels in new order
roi_labels_n = {'Accumbens_Left';'Accumbens_Right';'Amygdala_Left';'Amygdala_Right';...
    'BrainStem';'Caudate_Left';'Caudate_Right';'Cerebellum_Left';'Cerebellum_Right';...
    'DiencephalonVentral_Left';'DiencephalonVentral_Right';'Hippocampus_Left';...
    'Hippocampus_Right';'Pallidum_Left';'Pallidum_Right';'Putamen_Left';...
    'Putamen_Right';'Thalamus_Left';'Thalamus_Right'};

% get new index for labels
Index = zeros(numel(roi_labels_o),1);

for kk = 1:numel(roi_labels_o)
    Index(kk) = find(strcmp(roi_labels_n(kk),roi_labels_o));
end

%% plot
for vv = 1:length(vals)
    val = vals{vv};
    for cc = 1:length(clus)
        clu = clus{cc};
        % get the vector
        vec_tmp = vec_subcortical.(val).(clu);
        % plot the heatmap
        h = heatmap(vec_tmp(Index)');
        h.YDisplayLabels = roi_labels_n;
        % set color scale
        if isequal(val,'pos')
            col_scale = CBIG_GenerateColorscale('red_yell');
        else
            col_scale = CBIG_GenerateColorscale('blue_cyan');
        end
        colormap(col_scale)
        % set the right limits for each behavior
        if isequal(clu,'cog')
            caxis(lim_cog)
        elseif isequal(clu,'pers')
            caxis(lim_pers)
        else
            caxis(lim_mh)
        end
        % save the figure
        fname = [path_out filesep 'subcort_heatmap_' clu '_' val '.svg'];
        saveas(h,fname)
        close all
    end
end
end

%% custom color scale function
function colorscale = CBIG_GenerateColorscale(cmap)
    res = 140;
    if strcmp(cmap, 'red_yell')
        RGB_new = [255 255 255;
            255 236 94;
            232 208 28;
            181 148 0;
            181 97  0;
            181 51  0;
            181 0   0;
            117 0   0;
            0   0   0];
    end
    if strcmp(cmap, 'blue_cyan')
        RGB_new = [255 255 255;
            179 255 251;
            43  224 215;
            43  182 224;
            43  155 224;
            43  118 224;
            13  75  184;
            0   43  117;
            0   0   0]; 
    end
        
    RGB_COLORSCALE = RGB_new;
    orig_colorscale_length = size(RGB_COLORSCALE, 1);
    colorscale = [];
    step = round(res / orig_colorscale_length);
    
    count = 1;
    for i = 1:orig_colorscale_length-1
        for j = 0:step-1
            colorscale(count, :) = RGB_COLORSCALE(i, :) * (step - j) + RGB_COLORSCALE(i+1, :) * j;
            colorscale(count, :) = colorscale(count, :) * 1.0 / step;
            count = count + 1;
        end
    end
    colorscale = [colorscale; RGB_COLORSCALE(orig_colorscale_length, :)];
    
    pad = (res - size(colorscale, 1)) / 2;
    for i = 1:pad
        colorscale = [RGB_COLORSCALE(1, :); colorscale];
    end
    
    pad = res - size(colorscale, 1);
    for i = 1:pad
        colorscale = [colorscale; RGB_COLORSCALE(orig_colorscale_length, :)];
    end
    
    colorscale = flipud(colorscale) / 255;
end
