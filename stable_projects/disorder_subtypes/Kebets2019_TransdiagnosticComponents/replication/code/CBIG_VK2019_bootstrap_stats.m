function [std_behav_boot,zscore_behav_boot,pvals_behav_boot,std_RSFC_boot,zscore_RSFC_boot,pvals_RSFC_boot] = ...
    CBIG_VK2019_bootstrap_stats(LC_behav_loadings,LC_behav_loadings_boot,LC_RSFC_loadings,LC_RSFC_loadings_boot,...
    nRois,signif_LC,scripts_dir,out_dir)
%
% This function computes the mean and std deviation across the RSFC &
% behavior loadings obtained with bootstrap resampling
% (CBIG_VK2019_bootstrap_loadings).
% Z-scores & p-values are then calculated by dividing the original loadings
% by their bootstrap-estimated 
% standard deviation. 
% For RSFC, both original & bootstrap-estimated loadings are first averaged
% within and between 18 networks.
%
% Inputs:
% - LC_behav_loadings        : B x L matrix, B is #behaviors, L is 
%                              #latent components (LCs), behavior loadings
% - LC_behav_loadings_boot   : B x S x P matrix, S is #significant LCs, 
%                              P is # bootstrap samples, bootstrapped 
%                              behavior loadings 
% - LC_RSFC_loadings         : M x L matrix, M is #FC, RSFC loadings 
%                              within/between 18 networks for significant LCs
% - LC_RSFC_loadings_boot    : M x S x P matrix, bootstrapped RSFC 
%                              loadings for significant LCs
% - nRois                    : number of ROIs in RSFC matrix (e.g. 419)
% - signif_LC                : significant latent components (e.g. [1,2])
% - scripts_dir              : directory where are saved .csv
%                              & .txt files used for re-indexing parcels 
%                              when averaging loadings within/between networks
%
% Outputs:
% std_behav_boot             : B x S matrix, standard deviation of behavior
%                              loadings across bootstrap samples for 
%                              significant LCs
% zscore_behav_boot          : B x S matrix, zscore of behavior
%                              loadings across bootstrap samples for significant LCs
% pvals_behav_boot           : B x S matrix, p-values of behavior loadings
% std_RSFC_boot              : 18 x 18 x S matrix, standard deviation of
%                              RSFC loadings across bootstrap samples for 
%                              significant LCs
% zscore_RSFC_boot           : 18 x 18 x S matrix, zscore of RSFC
%                              loadings across bootstrap samples for significant LCs
% pvals_RSFC_boot            : 18 x 18 x S matrix, p-values of RSFC 
%                              loadings for significant LCs
%
% Written by Valeria Kebets and CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md

nBehav = size(LC_behav_loadings,1);
nBootstraps = size(LC_behav_loadings_boot,3);

%%% Behavior loadings

% Compute mean, std deviation, z-scores, p-values of behavior loadings
for iter_lc = 1:size(signif_LC,1)   
    for iter_behav = 1:nBehav

        this_lc = signif_LC(iter_lc);
        
        % Std across samples
        std_behav_boot(iter_behav,iter_lc) = std(LC_behav_loadings_boot(iter_behav,iter_lc,:));
               
        % Z-score (original loading / std loading across samples)
        zscore_behav_boot(iter_behav,iter_lc) = ...
            LC_behav_loadings(iter_behav,this_lc) / std_behav_boot(iter_behav,iter_lc);
        
        % P-values for z-scores
        if zscore_behav_boot(iter_behav,iter_lc) >= 0
            pvals_behav_boot(iter_behav,iter_lc) = 1-cdf('norm',zscore_behav_boot(iter_behav,iter_lc),0,1);
        elseif zscore_behav_boot(iter_behav,iter_lc) < 0
            pvals_behav_boot(iter_behav,iter_lc) = cdf('norm',zscore_behav_boot(iter_behav,iter_lc),0,1);
        end
        
    end    
end

%%% RSFC loadings

% Convert vectorized RSFC loadings into 419x419 matrices 
for iter_lc = 1:size(signif_LC,1)
    this_lc = signif_LC(iter_lc);
    LC_RSFC_loadings_CM(:,:,iter_lc) = jVecToSymmetricMat(LC_RSFC_loadings(:,this_lc),nRois,1);
end

% Average RSFC loadings within/between 18 networks 
file1_in = fullfile(out_dir,'LC_RSFC_loadings_419x419.mat');
save(file1_in,'LC_RSFC_loadings_CM');
file1_out = fullfile(out_dir,'LC_RSFC_loadings_18x18.mat');
cd(scripts_dir);

a = system(['python CBIG_VK2019_avgCMbyNetworks.py ',file1_in,' ','LC_RSFC_loadings_CM ',scripts_dir,' ',file1_out]);

cd(scripts_dir);        
% Convert vectorized bootstrapped RSFC loadings into 419x419 matrices
% Files are saved separately for each LC because of size
for iter_lc = 1:size(signif_LC,1)
    this_lc = signif_LC(iter_lc);
    
    for iter_boot = 1:nBootstraps
        LC_RSFC_loadings_boot_CM(:,:,iter_boot) = ...
            jVecToSymmetricMat(LC_RSFC_loadings_boot(:,this_lc,iter_boot),nRois,1);
    end
    
        % Average bootstrapped RSFC loadings within/between 18 networks
        file2_in = fullfile(out_dir,['LC' num2str(this_lc) '_RSFC_loadings_boot_419x419.mat']);
        save(file2_in,'LC_RSFC_loadings_boot_CM');
        file2_out = fullfile(out_dir,['LC' num2str(this_lc) '_RSFC_loadings_boot_18x18.mat']);
        
        a = system(['python CBIG_VK2019_avgCMbyNetworks.py ',...
            file2_in,' ','LC_RSFC_loadings_boot_CM ',scripts_dir,' ',file2_out]);
        
        clear file2_out file2_in
end

% Load files with averaged RSFC loadings
load(file1_out);
LC_RSFC_loadings_avg = thisLC_RSFC_loadings_avg;

files = dir(fullfile(out_dir,'*_RSFC_loadings_boot_18x18.mat'));

for iter_lc = 1:size(files,1)
    load(fullfile(out_dir,files(iter_lc).name));
    LC_RSFC_loadings_boot_avg(:,:,iter_lc,:) = thisLC_RSFC_loadings_avg(:,:,:);
    clear thisLC_RSFC_loadings_avg
end

% Compute mean, std deviation, z-scores, p-values of RSFC loadings
nNetworks = size(LC_RSFC_loadings_boot_avg,1);

for iter_lc = 1:size(signif_LC,1)
    for iter_net1 = 1:nNetworks
        for iter_net2 = 1:nNetworks
            
            % Std across samples
            std_RSFC_boot(iter_net1,iter_net2,iter_lc) = ...
                std(LC_RSFC_loadings_boot_avg(iter_net1,iter_net2,iter_lc,:));
            
            % Z-score (original loading / std loading across samples)
            zscore_RSFC_boot(iter_net1,iter_net2,iter_lc) = ...
                LC_RSFC_loadings_avg(iter_net1,iter_net2,iter_lc) / std_RSFC_boot(iter_net1,iter_net2,iter_lc);
            
            % P-values for z-scores
            if zscore_RSFC_boot(iter_net1,iter_net2,iter_lc) >= 0
                pvals_RSFC_boot(iter_net1,iter_net2,iter_lc) = ...
                    1-cdf('norm',zscore_RSFC_boot(iter_net1,iter_net2,iter_lc),0,1);
            elseif zscore_RSFC_boot(iter_net1,iter_net2,iter_lc) < 0
                pvals_RSFC_boot(iter_net1,iter_net2,iter_lc) = ...
                    cdf('norm',zscore_RSFC_boot(iter_net1,iter_net2,iter_lc),0,1);
            end
        end
    end
end
