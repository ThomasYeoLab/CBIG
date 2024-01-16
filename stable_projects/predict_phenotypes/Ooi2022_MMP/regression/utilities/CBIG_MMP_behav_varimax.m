function [varimax_1, varimax_2, varimax_3, loadings_varimax] = CBIG_MMP_behav_varimax(init_y) 

% [varimax_1, varimax_2, varimax_3, loadings_varimax] = CBIG_MMP_behav_varimax(init_y) 
% 
% This function performs PCA, extracts the first 3 components, and performs a varimax rotation.
% 
% Input:
% - init_y 
%   A matrix of #subjects x #behaviours unnormalized behavioural scores
%
% Output:
% - varimax_1 
%   Behaviour scores of first component after pca and varimax rotation.
%
% - varimax_2  
%   Behaviour scores of second component after pca and varimax rotation. 
%
% - varimax_3 
%   Behaviour scores of third component after pca and varimax rotation. 
%
% Written by Leon Ooi and CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md
    
    % z normalize
    init_y_norm = normalize(init_y,1);
    
    % PCA
    [loadings, components, ~, ~, var_exp, ~] = pca(init_y_norm);
    % print var explained and loadings
    pca_1 = components(:,1); 
    pca_2 = components(:,2); 
    pca_3 = components(:,3); 

    % EFA varimax
    loadings_varimax = rotatefactors(loadings(:,1:3));
    varimax_1 = init_y_norm*loadings_varimax(:,1);
    varimax_2 = init_y_norm*loadings_varimax(:,2);
    varimax_3 = init_y_norm*loadings_varimax(:,3);
    
end