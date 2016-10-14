function CBIG_visualize_factorComp(subInfoFile, GMICVFile, K, GRP)

% CBIG_visualize_factorComp(subInfoFile, GMICVFile, K, GRP)
%
% Written by CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md

C_EDGE = 'k';
C_AP = 'r';


%--------------- Get barycentric coordinates (probability vector)
% Subtypes estimated at ADNI 1 baseline
[rid_dx, ~, rid_prob, ~, ~, ~] = get_data(subInfoFile, GMICVFile, K);
rid = CBIG_select_predementedGrp(rid_dx, GRP);
baryCoor = rid_prob(ismember(rid_prob(:, 1), rid), 2:end);

if K == 2
    figure;
    %--------------- Scatter
    hist(baryCoor(:, 1), 0.025:0.05:0.975); % first topic
    h = findobj(gca, 'Type', 'patch');
    set(h(1), 'FaceColor', C_AP, 'EdgeColor', 'w');
    ax = gca;
    set(ax, 'XTick', [0.5]);
    set(ax, 'XLim', [0 1]);
    set(ax, 'YTick', [0 5 10 15 20], 'FontSize', 20);
    set(ax, 'YLim', [0 20]);
    box off;
    %!!! Annotation
    fprintf('!!! Subtype order assumed to be: cortical, subcortical+temporal.\n');
    fprintf('!!! Make sure this is true in the gamma file.\n');
    text('String', {'Temporal+Subcortical'}, 'FontSize', 20, 'Position', [0, 0]);
    text('String', {'Cortical'}, 'FontSize', 20, 'Position', [1, 0]);
elseif K == 3
    %--------------- Create triangulation
    P = [0 1/sqrt(3); 0.5 -0.5/sqrt(3); -0.5 -0.5/sqrt(3)];
    T = [1 2 3]; % connectivity
    TR = triangulation(T, P); % create the triangulation
    %--------------- Draw triangle
    figure;
    hold on;
    lw = 2;
    plot([P(1, 1) P(2, 1)], [P(1, 2) P(2, 2)], 'Color', C_EDGE, 'LineWidth', lw);
    plot([P(2, 1) P(3, 1)], [P(2, 2) P(3, 2)], 'Color', C_EDGE, 'LineWidth', lw);
    plot([P(3, 1) P(1, 1)], [P(3, 2) P(1, 2)], 'Color', C_EDGE, 'LineWidth', lw);
    %!!! Annotation
    fprintf('!!! Subtype order assumed to be: cortical, temporal, subcortical.\n');
    fprintf('!!! Make sure this is true in the gamma file.\n');
    text('String', {'Cortical'}, 'FontSize', 20, 'Position', P(1, :));
    text('String', {'Temporal'}, 'FontSize', 20, 'Position', P(2, :));
    text('String', {'Subcortical'}, 'FontSize', 20, 'Position', P(3, :));
    %--------------- Scatter
    ti = ones(size(baryCoor, 1), 1);
    cartCoor = barycentricToCartesian(TR, ti, baryCoor);
    scatter(cartCoor(:, 1), cartCoor(:, 2), 100, C_AP, 'filled');
    hold off;
    axis off;
elseif K == 4
    P = [0 1/sqrt(3) 0; 0.5 -0.5/sqrt(3) 0; -0.5 -0.5/sqrt(3) 0; 0 0 sqrt(2/3)];
    %--------------- Draw tetrahedron
    figure;
    hold on;
    lw = 2;
    plot3([P(1, 1) P(2, 1)], [P(1, 2) P(2, 2)], [P(1, 3) P(2, 3)], 'Color', C_EDGE, 'LineWidth', lw);
    plot3([P(2, 1) P(3, 1)], [P(2, 2) P(3, 2)], [P(2, 3) P(3, 3)], 'Color', C_EDGE, 'LineWidth', lw);
    plot3([P(3, 1) P(1, 1)], [P(3, 2) P(1, 2)], [P(3, 3) P(1, 3)], 'Color', C_EDGE, 'LineWidth', lw);
    plot3([P(1, 1) P(4, 1)], [P(1, 2) P(4, 2)], [P(1, 3) P(4, 3)], 'Color', C_EDGE, 'LineWidth', lw);
    plot3([P(2, 1) P(4, 1)], [P(2, 2) P(4, 2)], [P(2, 3) P(4, 3)], 'Color', C_EDGE, 'LineWidth', lw);
    plot3([P(3, 1) P(4, 1)], [P(3, 2) P(4, 2)], [P(3, 3) P(4, 3)], 'Color', C_EDGE, 'LineWidth', lw);
    %!!! Annotation
    fprintf('!!! Subtype order assumed to be: parietal, subcortical, frontal, temporal.\n');
    fprintf('!!! Make sure this is true in the gamma file.\n');
    text('String', {'Parietal'}, 'FontSize', 20, 'Position', P(1, :));
    text('String', {'Subcortical'}, 'FontSize', 20, 'Position', P(2, :));
    text('String', {'Frontal'}, 'FontSize', 20, 'Position', P(3, :));
    text('String', {'Temporal'}, 'FontSize', 20, 'Position', P(4, :));
    %--------------- Scatter
    % Barycentric to cartesian
    cartCoor = zeros(size(baryCoor, 1), 3);
    for idx = 1:size(cartCoor, 1)
        w = baryCoor(idx, :);
        cartCoor(idx, :) = w(1)*P(1, :)+w(2)*P(2, :)+w(3)*P(3, :)+w(4)*P(4, :);
    end
    scatter3(cartCoor(:, 1), cartCoor(:, 2), cartCoor(:, 3), 100, C_AP, 'filled');
    hold off;
    axis off;
    view(153, 34);
else
    error('Not configured!');
end
