function CBIG_ASDf_visualizeFactorComp(sub_info_file, factorComp, K, output_name)
% CBIG_ASDf_visualizeFactorComp(sub_info_file, factorComp, K, output_name)
%
% This function plots participants' factor compositions for K = 2 and 3.
% For K = 2, the visualization is created using histogram; for K = 3, the
% visualization is created on a triangle using Barycentric coordinates.
% If output_name is specified, the plot will be saved.
%
% Input:
%     - sub_info_file:
%           absolute path to .csv file containing subject demographics info
%     - factorComp:
%           NxK matrix, where N is the number of subjets and K is the
%           number of factors. The i-th row represents the i-th subject's
%           factor composition which sums to 1.
%     - K:
%           Integer. Number of factors
%     - output_name (optional):
%           Output file name (including full path) that will be used to save the plot
%
% Example:
%     CBIG_ASDf_visualizeFactorComp('~/subInfo.csv', factorComp, 3);
%     Draws factorComp on a triangle using Barycentric coordinates. Each
%     dot on the plot represents one subject, and the face color of the dot
%     is different for different sites. Plot will not be saved.
%
% Written by CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md

%% get site info
[id_sites, id_dx] = CBIG_ASDf_getSubData(sub_info_file); % get sites and diagnoses of all subjects
sites_all = id_sites(:,2);
dx_all = cell2mat(id_dx(:,2));

sites_asd = sites_all(dx_all == 1);
sites_name = unique(sites_asd);

%% define marker type and color, 14 sites in total
C_EDGE = 'k';
marker_c = [18, 41, 255; 
            0, 3, 120; 
            217, 250, 0; 
            94, 136, 181; 
            203, 111, 0; 
            181, 240, 92; 
            237, 102, 245; 
            107, 23, 255; 
            21, 126, 0; 
            66, 255, 99; 
            61, 146, 202; 
            191, 154, 89; 
            61, 198, 0; 
            230, 0, 0];
marker_c = marker_c ./ 255;

%% plot
figure;
if K == 2
    %--------------- Scatter
    hist(factorComp(:,1), 0.025:0.05:0.975); % first topic
    h = findobj(gca, 'Type', 'patch');
    set(h(1), 'FaceColor', C_AP, 'EdgeColor', 'w');
    ax = gca;
    set(ax, 'XTick', [0.5]);
    set(ax, 'XLim', [0 1]);
    set(ax, 'TickDir', 'out');
    %     set(ax, 'YLim', [0 10]); %hard-coded here
    box off;
    %!!! Annotation
    text('String', {'F2'}, 'FontSize', 20, 'Position', [0, 0]);
    text('String', {'F1'}, 'FontSize', 20, 'Position', [1, 0]);
elseif K == 3
    %--------------- Create triangulation
    P = [0 1/sqrt(3); 0.5 -0.5/sqrt(3); -0.5 -0.5/sqrt(3)];
    T = [1 2 3]; % connectivity
    TR = triangulation(T, P); % create the triangulation
    %--------------- Draw triangle
    hold on;
    lw = 2;
    plot([P(1, 1) P(2, 1)], [P(1, 2) P(2, 2)], 'Color', C_EDGE, 'LineWidth', lw);
    plot([P(2, 1) P(3, 1)], [P(2, 2) P(3, 2)], 'Color', C_EDGE, 'LineWidth', lw);
    plot([P(3, 1) P(1, 1)], [P(3, 2) P(1, 2)], 'Color', C_EDGE, 'LineWidth', lw);
    %!!! Annotation
    text('String', {'F1'}, 'FontSize', 20, 'Position', P(1, :));
    text('String', {'F2'}, 'FontSize', 20, 'Position', P(2, :));
    text('String', {'F3'}, 'FontSize', 20, 'Position', P(3, :));
    for j = 1:length(sites_name)
        indices = strcmp(sites_asd, sites_name(j));
        factorComp_currSite = factorComp(indices,:);
        %--------------- Scatter
        ti = ones(size(factorComp_currSite, 1), 1);
        cartCoor = barycentricToCartesian(TR, ti, factorComp_currSite);
        scatter(cartCoor(:, 1), cartCoor(:, 2), 100, marker_c(j,:), 'filled');
%         legend(sites_name(j));
        hold on;
    end
    axis off;
    legends = ['axis1'; 'axis2'; 'axis3'; sites_name];
    legend(legends);
    
else
    error('Not configured!');
end

%% Save the plot
if nargin > 3 && ~isempty(output_name)
    hgexport(gcf, output_name);
    eps2xxx([output_name '.eps'], {'png'});
end
