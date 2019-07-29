function CBIG_VK2019_plot_loadings(out_dir,LC_behav_loadings,behav_names,LC_RSFC_loadings,nRois,signif_LC) 
%  
% This function plots the RSFC and behavior loadings of all 
% significant latent components (LCs). Behavior loadings are plotted as
% barcharts and RSFC loadings are plotted as correlation matrices.
%
% Inputs:
% - out_dir             : output directory where plots are saved
% - LC_behav_loadings   : B x L matrix, B is #behavior, L is #components, behavior loadings
% - behav_names         : string, names of behavior variables
% - LC_RSFC_loadings    : M x L matrix, M is #RSFC, RSFC loadings
% - nRois               : number of ROIs in RSFC matrix (e.g. 419)
% - signif_LC           : significant latent components to plot (e.g. [1,2])
%
% Written by Valeria Kebets and CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md

nBehav = size(LC_behav_loadings,1);

cd(out_dir);

% Plot behavioral loadings
for iter_lc = 1:length(signif_LC)
    this_lc = signif_LC(iter_lc);
    
    figure;
    bar(LC_behav_loadings(:,this_lc));    
    xticks(1:nBehav);
    xticklabels(behav_names);
    set(gca,'TickLabelInterpreter','none','FontSize',6,'Box','off');
    set(gcf,'Color','w');       
    xtickangle(45);
    xlabel('Behavioral variables');
    ylabel('Correlation');
    title(['LC' num2str(this_lc) ' - Behavioral loadings']);
    name_fig = ['LC' num2str(this_lc) '_behav_loadings.jpg'];
    saveas(gcf,name_fig);
end


% Plot RSFC loadings
for iter_lc = 1:length(signif_LC)
    this_lc = signif_LC(iter_lc);
    
    clear CM myMin myMax myLim lh2lh lh2rh rh2rh lh2subcor rh2subcor subcor2subcor name_fig
        
    CM = jVecToSymmetricMat(LC_RSFC_loadings(:,this_lc),nRois,1);
    
    % Separate FC matrix into 3 subparts (lh2lh, lh2rh, rh2rh)
    lh2lh = CM(1:200,1:200);
    lh2rh = CM(1:200,201:400);
    rh2rh = CM(201:400,201:400);
    lh2subcor = CM(1:200,401:419);
    rh2subcor = CM(201:400,401:419);
    subcor2subcor = CM(401:419,401:419);
    
    myLim = 0.3;
    name_fig = ['LC' num2str(this_lc) '_RSFC_loadings'];
    
    CBIG_Plot_Schaefer400_17Networks19SubcorRearrCorrMat_WhiteGrid...
        (lh2lh, lh2rh, rh2rh, lh2subcor, rh2subcor, subcor2subcor, [-myLim myLim], name_fig);    
    
end