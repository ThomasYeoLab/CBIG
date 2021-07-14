function new_corr_mat = CBIG_TRBPC_rearrange_matrix(corr_mat, scale)

% new_corr_mat = CBIG_TRBPC_rearrange_matrix(corr_mat, scale)
% 
% This function will take a 419x419 matrix and rearrange it by the major
% networks and subcortical structures to generate a new matrix with the new
% ordering.
%
% This function assumes that the order of the subcortical ROIs in the
% input correlation matrices is in ascending order based on the labels in
% $FREESURFER_HOME/ASegStatsLUT.txt file.
%
% Ordering of major networks and subcortical regions from left to right:
%
%       Major network/Subcortical regions        Sub-network
%     
%     Cortical networks:
%
%     1)  Default:                                TempPar
%                                                 DefaultC 
%                                                 DefaultB
%                                                 DefaultA
%
%     2)  Control:                                ContC
%                                                 ContB 
%                                                 ContA
%
%     3)  Limbic
%
%     4)  SalVentAttn:                            SalVentAttnB
%                                                 SalVentAttnA
%
%     5)  DorsAttn:                               DorsAttnB
%                                                 DorsAttnA
%
%     6)  SomMot:                                 SomMotB
%                                                 SomMotA
%
%     7)  Visual:                                 VisPeri
%                                                 VisCent
%     
%     Subcortical structures:
%     The following 4 striatum structures are arranged together:
%     1)  Accumbens
%     2)  Caudate
%     3)  Pallidum
%     4)  Putamen
%
%     5)  Thalamus
%
%     6)  Amygdala
%     
%     7)  Hippocampus
%
%     8)  Brain Stem
%
%     9)  Diencephalon (Ventral)
%
%    10)  Cerebellum 
%
% Within each sub-network or subcortical region, correlation matrix entries
% start from left hemisphere, then right hemisphere entries (from left to
% right).
%
% Input:
%     - corr_mat:
%           419x419 matrix. The order of the subcortical ROIs in this
%           matrix is assumed to be in ascending order based on the labels
%           in $FREESURFER_HOME/ASegStatsLUT.txt file.
%     - scale:
%           integer, 7 or 17 for the resolution of cortical networks
%
% Written by Siyi Tang, Angela Tam & CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md


%% function starts here
% load original cortical networks labels
networks_path = fullfile(getenv('CBIG_CODE_DIR'),...
    'stable_projects/brain_parcellation/Schaefer2018_LocalGlobal/Parcellations/FreeSurfer5.3/fsaverage5/label/');
if scale == 7
    lh_annot_file = [networks_path 'lh.Schaefer2018_400Parcels_7Networks_order.annot'];
    rh_annot_file = [networks_path 'rh.Schaefer2018_400Parcels_7Networks_order.annot'];
else
    lh_annot_file = [networks_path 'lh.Schaefer2018_400Parcels_17Networks_order.annot'];
    rh_annot_file = [networks_path 'rh.Schaefer2018_400Parcels_17Networks_order.annot'];
end
[lh_vertex_labels, lh_colortable] = CBIG_read_annotation(lh_annot_file);
[rh_vertex_labels, rh_colortable] = CBIG_read_annotation(rh_annot_file);

lh_label = lh_colortable.struct_names(2:end);
rh_label = rh_colortable.struct_names(2:end);

% hard-coded, assuming the original input subcortical labels are in ascending order
subcor_labelname = {'Cerebellum_Left';'Thalamus_Left';'Caudate_Left';'Putamen_Left';'Pallidum_Left';...
    'BrainStem';'Hippocampus_Left';'Amygdala_Left';'Accumbens_Left'; ...
    'DiencephalonVentral_Left';'Cerebellum_Right';'Thalamus_Right';'Caudate_Right';'Putamen_Right';...
    'Pallidum_Right';'Hippocampus_Right';'Amygdala_Right';'Accumbens_Right';'DiencephalonVentral_Right'};

major_grid = [];
minor_grid = [];
subcor_grid = [];
%major_acc_index = [1, 4, 7, 8, 10, 12, 14, 16];
major_acc_index = [1, 4, 7, 9, 11, 13, 15, 17];
subcor_acc_index = 17:25;

lhrh_label = {lh_label{:} rh_label{:}}';
all_label = [lhrh_label; subcor_labelname];
all_netw = cell(numel(all_label),1);
all_subnet = cell(numel(all_label),1);

for i = 1:numel(all_label)
    % cortical
    if i <= numel(lhrh_label)
        netw = textscan(char(all_label(i,1)), '%s %s %s %s', 'delimiter', '_');
        tmp = netw{1,3};
        all_netw(i,1) = tmp;
        subnet = netw{1,4};
        if ~isempty(subnet)
            all_subnet(i,1) = subnet;
        end
        
    % subcortical    
    else
        netw = textscan(char(all_label(i,1)),'%s %s', 'delimiter', '_');
        all_netw(i,1) = netw{1,1};       
    end
end

% arrange new labels based on template order
if scale == 7
    tmplate = {'Default'; 'Cont'; 'Limbic'; 'SalVentAttn'; 'DorsAttn'; 'SomMot'; 'Vis';...
        'Accumbens'; 'Caudate'; 'Pallidum'; 'Putamen'; 'Thalamus'; 'Amygdala';...
        'Hippocampus'; 'BrainStem'; 'DiencephalonVentral'; 'Cerebellum'};
else
    tmplate = {'TempPar'; 'DefaultC'; 'DefaultB';'DefaultA'; 'ContC'; 'ContB'; 'ContA'; 'LimbicA'; 'LimbicB';...
        'SalVentAttnB'; 'SalVentAttnA'; 'DorsAttnB'; 'DorsAttnA'; 'SomMotB'; 'SomMotA'; 'VisPeri'; 'VisCent'; ...
        'Accumbens'; 'Caudate'; 'Pallidum'; 'Putamen'; 'Thalamus'; 'Amygdala'; 'Hippocampus'; 'BrainStem';...
        'DiencephalonVentral'; 'Cerebellum'};
    %tmplate2 = {'TempPole'; 'OFC'};
end

newlabel = [];

if scale == 7
    for j = 1:numel(tmplate)
        ind = find(strcmp(all_netw,tmplate(j)));
        newlabel = cat(1, newlabel,all_label(ind));
    end
else
    curInd = 0;
    for j = 1:numel(tmplate)
        ind = find(strcmp(all_netw,tmplate(j)));
        
        %    % for Limbic networks, further separate the networks into TempPole and OFC
        %    if sum(strcmp('Limbic', tmplate(j))) ~= 0
        %        ind2 = [];
        %        for s = 1:numel(tmplate2)
        %            if sum(strcmp(tmplate2(s), all_subnet(ind))) ~= 0
        %                ind2 = [ind2; ind(find(strcmp(tmplate2(s), all_subnet(ind))))];
        %                minor_grid = [minor_grid curInd+size(ind2,1)];
        %            end
        %
        %        end
        %        if numel(ind) ~= numel(ind2)
        %            disp('Wrong Index')
        %        end
        %        ind = ind2;
        %    end
        
        curInd = curInd+size(ind,1);
        if (j~=numel(tmplate))
            if (any(j==subcor_acc_index))
                subcor_grid = [subcor_grid curInd];
            else
                minor_grid = [minor_grid curInd];
            end
        end
        if (any(j==major_acc_index))
            major_grid = [major_grid curInd];
        end
        
        newlabel = cat(1, newlabel,all_label(ind));
    end
end

% create indexing from old labeling to new labeling
Index = zeros(numel(all_label),1);

for k = 1:numel(all_label)
    Index(k) = find(strcmp(newlabel(k),all_label));
end

fin_label = all_label(Index);

% rearrange the original matrix with new indices
new_corr_mat = corr_mat(Index,Index);
end

