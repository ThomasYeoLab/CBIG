function [data_reg] = CBIG_VK2019_regrOutConfounds...
    (data_nonreg,names_regr,regr_dir,motion_dir,RSfile_stem,idx_subj,names_subj)
% This function regresses out the effects of confounds on a data matrix. 
% 
% Inputs:
%
% - data_nonreg       : N x X matrix, N is #subjects, X is #data variables
% - names_regr        : confound variables as saved in path_regressors
%                       For motion: 
%                       motion1 : only FDRMS
%                       motion2 : FDRMS & DVARS 
% - regr_dir          : directory where confounds are saved (.mat files)
% - motion_dir        : directory where motion files (e.g. FDRMS) are located
% - RSfile_stem       : stem of RS file (e.g. 'bld001_rest_skip4_stc')
% - idx_subj          : N x 1 vector, index of subjects included
% in the analysis (to select subjects from regressors files) (e.g. 1)
% - names_subj        : N x 1 string, names of subjects (e.g. 'sub-10159')
%
% Outputs:
% - data_reg          : regressed out data
%
% Written by Valeria Kebets and CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md

regressors = [];

for iter_reg = 1:numel(names_regr)
    
    a = strfind(names_regr{iter_reg},'motion');    
    if a == 1    
        
        b = strcmp(names_regr{iter_reg},'motion1');       
        if b == 1
            % FDRMS
            for iter_subj = 1:numel(names_subj)
                clear mot
                mot = load(fullfile(motion_dir,names_subj{iter_subj},'bold','mc',...
                    [names_subj{iter_subj} '_' RSfile_stem '_motion_outliers_FDRMS']));
                fdrms(iter_subj,:) = mean(mot);
            end
            
            regressors = [regressors fdrms];
            
        end
        
        c = strcmp(names_regr{iter_reg},'motion2'); % motion1=only FDRMS, motion2=FDRMS & DVARS       
        if c == 1
            % DVARS
            for iter_subj = 1:numel(names_subj)
                clear mot
                mot = load(fullfile(motion_dir,names_subj{iter_subj},'bold','mc',...
                    [names_subj{iter_subj} '_' RSfile_stem '_motion_outliers_DVARS']));
                dvars(iter_subj,:) = mean(mot);
            end
            
            regressors = [regressors dvars];
            
        end
        
        clear b c
        
    else
        
        % Other confounds
        load(fullfile(regr_dir,[names_regr{iter_reg} '.mat']));
        this_score = score(idx_subj);
        regressors = [regressors this_score];
        clear score scoreNames this_score
    end
    
end
        
% Regress out confounds
[data_reg, ~, ~, ~] = CBIG_glm_regress_matrix(data_nonreg, regressors, 0, []);

end


