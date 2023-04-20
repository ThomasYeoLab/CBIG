function CBIG_PlotCorrMatNetOrder(res, corr_mat, netorder, scalelim)

% CBIG_PlotCorrMatNetOrder(res, corr_mat, netorder, scalelim)
%
% This function is a wrapper script which visualizes ROI2ROI correlation matrix of Schaefer2018 
% and Yeo2011 with several pre-defined ROI order. The correlation matrix <corr_mat> can contain
% subcortical ROIs.
% The ROIs will be orderred based on a network structure which is defined as a cell structure:
%  {{'Net1','Net2'},
%   {'Net3'},
%   {'Net4'}}
% A more general function CBIG_PlotCorrMat_general.m can be used to visualize any ROI order.
%
% Input:
%   - res: (scalar)
%       The number of CORTICAL ROIs. For example:
%       res = 400 (not 419!) for Schaefer400+19subcortical
%       res = 100 for Schaefer100 without subcortical
%       res = 114 for Yeo2011 17 networks with split components.
%
%   - corr_mat: (KxK matrix)
%       The ROI2ROI correlation matrix. K = res if there is no subcortical ROI, else K = res+M if
%       there are M subcortical ROIs. The corr_mat is assumes to be formated as follows:
%       If there is no subcortical ROI:
%       corr_mat = [ lh2lh lh2rh;
%                    rh2lh rh2rh; ]
%       If there are subcortical ROIs:
%       corr_mat = [ lh2lh       lh2rh       lh2subcor;
%                    rh2lh       rh2rh       rh2subcor;
%                    subcor2lh   subcor2rh   subcor2subcor;]                    
%
%   - netorder: string
%       The predefined order type: 
%       1) 'Schaefer_Yeo17': The Schaefer parcellations with Yeo2011 17network ordering 
%       2) 'Schaefer_Yeo7': The Schaefer parcellations with Yeo2011 7network ordering 
%       3) 'Schaefer_Kong17': The Schaefer parcellations with Kong2022 17network ordering 
%       4) 'Yeo17': The Yeo2011 17 network parcellations with split components
%       5) 'Yeo7': The Yeo2011 7 network parcellations with split components 
%       6) 'Yan_Yeo17': The Yan parcellations with Yeo2011 17network ordering 
%       7) 'Yan_Yeo7': The Yan parcellations with Yeo2011 7network ordering 
%       8) 'Yan_Kong17': The Yan parcellations with Kong2022 17network ordering 
%       Check <network_order_struct> variable in the current script for ordering details.
%       
%   - scalelim: ([lim_min limmax])
%       The range of the colormap. The <scalelim> should be a 1x2 vector where <lim_min> and <lim_max>
%       are the mininum and maximum value of the range. By default, <lim_min> is the -1*maximum
%       absolute value of corr_mat, <lim_max> is the maximum abosulute value of corr_mat.
%
% Examples:
% corr_mat = rand(419,419);
% CBIG_PlotCorrMatNetOrder(400, corr_mat, 'Schaefer_Yeo17', [-0.4 0.4]);
%
% Written by Ruby Kong and CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md

%% set color limit
if ((nargin < 4) || (isempty(scalelim)))
    collim = max(max(abs(corr_mat)));
    scalelim = [-1*collim, 1*collim];
end

switch netorder
case 'Schaefer_Yeo17'
    group_labels = ft_read_cifti(fullfile(getenv('CBIG_CODE_DIR'), 'stable_projects',...
        'brain_parcellation', 'Schaefer2018_LocalGlobal', 'Parcellations', 'HCP', 'fslr32k',...
        'cifti', ['Schaefer2018_' num2str(res) 'Parcels_17Networks_order.dlabel.nii']), 'mapname', 'array');
    parcelname = group_labels.dlabellabel';
    network_order_struct = {{'TempPar'},
                            {'DefaultC', 'DefaultB', 'DefaultA'},
                            {'ContC', 'ContB', 'ContA'},
                            {'LimbicA', 'LimbicB'},
                            {'SalVentAttnB', 'SalVentAttnA'},
                            {'DorsAttnB', 'DorsAttnA'},
                            {'SomMotB', 'SomMotA'},
                            {'VisPeri', 'VisCent'}};

case 'Schaefer_Yeo7'
    group_labels = ft_read_cifti(fullfile(getenv('CBIG_CODE_DIR'), 'stable_projects',...
        'brain_parcellation', 'Schaefer2018_LocalGlobal', 'Parcellations', 'HCP', 'fslr32k',...
        'cifti', ['Schaefer2018_' num2str(res) 'Parcels_7Networks_order.dlabel.nii']), 'mapname', 'array');
    parcelname = group_labels.dlabellabel';
    network_order_struct = {{'Default'},
                            {'Cont'},
                            {'Limbic'},
                            {'SalVentAttn'},
                            {'DorsAttn'},
                            {'SomMot'},
                            {'Vis'}};

case 'Schaefer_Kong17'
    group_labels = ft_read_cifti(fullfile(getenv('CBIG_CODE_DIR'), 'stable_projects',...
        'brain_parcellation', 'Schaefer2018_LocalGlobal', 'Parcellations', 'HCP', 'fslr32k', 'cifti',...
        ['Schaefer2018_' num2str(res) 'Parcels_Kong2022_17Networks_order.dlabel.nii']), 'mapname', 'array');
    parcelname = group_labels.dlabellabel';
    network_order_struct = {{'DefaultC', 'DefaultB', 'DefaultA'},
                            {'ContC', 'ContB', 'ContA'},
                            {'Language'},
                            {'SalVenAttnB', 'SalVenAttnA'},
                            {'DorsAttnB', 'DorsAttnA'},
                            {'Aud'},
                            {'SomMotB', 'SomMotA'}
                            {'VisualC', 'VisualB', 'VisualA'}};                           

case 'Yeo17'
    group_labels = ft_read_cifti(fullfile(getenv('CBIG_CODE_DIR'), 'stable_projects',...
        'brain_parcellation', 'Yeo2011_fcMRI_clustering', '1000subjects_reference',...
        'Yeo_JNeurophysiol11_SplitLabels', 'fs_LR32k',...
        'Yeo2011_17Networks.split_components.dlabel.nii'), 'mapname', 'array');
    parcelname = group_labels.split_componentslabel';
    network_order_struct = {{'TempPar'},
                            {'DefaultC', 'DefaultB', 'DefaultA'},
                            {'ContC', 'ContB', 'ContA'},
                            {'LimbicA', 'LimbicB'},
                            {'SalVentAttnB', 'SalVentAttnA'},
                            {'DorsAttnB', 'DorsAttnA'},
                            {'SomMotB', 'SomMotA'},
                            {'VisPeri', 'VisCent'}};
                            
case 'Yeo7'
    group_labels = ft_read_cifti(fullfile(getenv('CBIG_CODE_DIR'), 'stable_projects',...
        'brain_parcellation', 'Yeo2011_fcMRI_clustering', '1000subjects_reference',...
        'Yeo_JNeurophysiol11_SplitLabels', 'fs_LR32k',...
        'Yeo2011_7Networks.split_components.dlabel.nii'), 'mapname', 'array');
    parcelname = group_labels.split_componentslabel';
    network_order_struct = {{'Default'},
                            {'Cont'},
                            {'Limbic'},
                            {'SalVentAttn'},
                            {'DorsAttn'},
                            {'SomMot'},
                            {'Vis'}};

case 'Yan_Kong17'
    group_labels = ft_read_cifti(fullfile(getenv('CBIG_CODE_DIR'), 'stable_projects',...
    'brain_parcellation', 'Yan2023_homotopic', 'parcellations', 'HCP', 'fsLR32k', 'kong17',...
    [num2str(res) 'Parcels_Kong2022_17Networks.dlabel.nii']), 'mapname', 'array');
    
    parcelname = group_labels.dlabellabel';
    network_order_struct = {{'DefaultC', 'DefaultB', 'DefaultA'},
                            {'ContC', 'ContB', 'ContA'},
                            {'Language'},
                            {'SalVenAttnB', 'SalVenAttnA'},
                            {'DorsAttnB', 'DorsAttnA'},
                            {'Aud'},
                            {'SomMotB', 'SomMotA'}
                            {'VisualC', 'VisualB', 'VisualA'}};

case 'Yan_Yeo17'
    group_labels = ft_read_cifti(fullfile(getenv('CBIG_CODE_DIR'), 'stable_projects',...
    'brain_parcellation', 'Yan2023_homotopic', 'parcellations', 'HCP', 'fsLR32k', 'yeo17',...
    [num2str(res) 'Parcels_Yeo2011_17Networks.dlabel.nii']), 'mapname', 'array');

    parcelname = group_labels.dlabellabel';
    network_order_struct = {{'TempPar'},
                            {'DefaultC', 'DefaultB', 'DefaultA'},
                            {'ContC', 'ContB', 'ContA'},
                            {'LimbicA', 'LimbicB'},
                            {'SalVentAttnB', 'SalVentAttnA'},
                            {'DorsAttnB', 'DorsAttnA'},
                            {'SomMotB', 'SomMotA'},
                            {'VisPeri', 'VisCent'}};
case 'Yan_Yeo7'
    group_labels = ft_read_cifti(fullfile(getenv('CBIG_CODE_DIR'), 'stable_projects',...
    'brain_parcellation', 'Yan2023_homotopic', 'parcellations', 'HCP', 'fsLR32k', 'yeo7',...
    [num2str(res) 'Parcels_Yeo2011_7Networks.dlabel.nii']), 'mapname', 'array');

    parcelname = group_labels.dlabellabel';
    network_order_struct = {{'Default'},
                            {'Cont'},
                            {'Limbic'},
                            {'SalVentAttn'},
                            {'DorsAttn'},
                            {'SomMot'},
                            {'Vis'}};                        
end

CBIG_PlotCorrMat_general(res, corr_mat, parcelname, network_order_struct, scalelim);
    
end

