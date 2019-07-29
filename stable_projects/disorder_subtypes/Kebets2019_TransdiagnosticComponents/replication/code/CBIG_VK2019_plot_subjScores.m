function CBIG_VK2019_plot_subjScores(Lx,Ly,names_groups,grouping,signif_LC)
%
% This function plots the RSFC and behavior composite scores of all 
% significant latent components (LCs). The diagnostic group of the subjects
% is indicated. Correlations between RSFC & behavior composite scores is
% displayed.
%
% Inputs:
% - Lx                  : N x L matrix, N is #subjects, L#components, 
% RSFC composite scores
% - Ly                  : N x L matrix, behavior composite scores
% - names_groups        : string, names of diagnostic groups
% - grouping            : N x 1 vector, subject (diagnostic) grouping 
%                         e.g. grouping = [1,2,3]
% - signif_LC           : significant latent components to plot
%
% Outputs:
% - corr_LxLy           : correlation between RSFC and behavior composite
% scores
%
% Written by Valeria Kebets and CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md

nGroups = numel(names_groups);

colors = {'b','r','c','g','m','y','w','k'}; % Matlab colors
plot_colors = colors(1:nGroups); % select as many colors as groups

disp('Correlations between RSFC and behavioral composite scores');
for iter_lc = 1:size(signif_LC,1)
    this_lc = signif_LC(iter_lc);
    
    figure;
    for iter_group = 1:numel(names_groups)
        plot(Lx(find(grouping==iter_group),this_lc),...
            Ly(find(grouping==iter_group),this_lc),[plot_colors{iter_group} '.'],'MarkerSize',10);
        hold on
    end
    hold off
    title(['LC' num2str(this_lc) ' - Correlations between RSFC and behavioral composite scores']);
    
    legend(names_groups,'Location','southeast');
    xlabel('RSFC composite scores');
    ylabel('Behavioral composite scores');
    
    set(gcf,'Color','w');
    set(gca,'Box','off');
    
    [corr_LxLy(iter_lc),~] = corr(Lx(:,this_lc),Ly(:,this_lc));
    disp(['LC' num2str(this_lc) ': r = ' num2str(corr_LxLy(iter_lc),'%0.2f')]);
    
end