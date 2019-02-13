function CBIG_LiGSR_explained_variance_HCP(restricted_csv, unrestricted_csv, trait_list, covariate_list, ...
    FD_file, DVARS_file, sub_list, rm_sub_list, FSM_file, outdir, ystem, m2_subdir, verbose)

% CBIG_LiGSR_explained_variance_HCP(restricted_csv, unrestricted_csv, trait_list, covariate_list, ...
%     FD_file, DVARS_file, sub_list, rm_sub_list, FSM_file, outdir, ystem, m2_subdir, verbose)
% 
% This function calls Morphometricity_1RandEffect.m (Sabuncu et al., 2016)
% to estimate the trait variance explained by RSFC for the Human Connectome
% Project (HCP) dataset.
% 
% Inputs:
%   - restricted_csv
%     A string. The full path of the restricted CSV file containing traits
%     and covariates of each subject, downloaded from the HCP website.
% 
%   - unrestricted_csv
%     A string. The full path of the unrestricted CSV file containing
%     traits and covariates of each subject, downloaded from the HCP website.
% 
%   - trait_list
%     A string. The full path to a text file of trait names. Each line
%     in this file corresponds to one trait measure. The trait names
%     should exist as headers in either "restricted_csv" or "unrestricted_csv".
% 
%   - covariate_list
%     A string. The full path to a text file containing the name of all
%     covariates. Each line corresponds to one covariate. Except for the
%     covariates FD and DVARS, all the other covariate names should exist
%     as headers in either "restricted_csv" or "unrestricted_csv".
% 
%   - FD_file (optional)
%     A string. The full path to a text file of the mean framewise
%     displacement (FD) of all subjects. The number of lines in "FD_file"
%     should be the same as the number of lines in "sub_list". If the user
%     wants to regress out FD, then covariate_list should contain 'FD'.
%     If the covariates do not include FD, this input variable is not
%     needed. The user can pass in 'NONE'.
%     
%   - DVARS_file (optional)
%     A string. The full path to a text file of the mean DVARS of all
%     subjects. The number of lines in "DVARS_file" should be the same as
%     the number of lines in "subject_list". If the user wants to regress
%     out 'DAVRS', then covariate_list should contain 'DVARS' (or 'DV').
%     If the covariates do not include DVARS, this input variable is not
%     needed. The user can pass in 'NONE'.
%
%   - sub_list
%     A string. The full path of a text file containing all subject IDs.
%     Each line is one subject ID.
% 
%   - rm_sub_list (optional)
%     A string. The full path of a text file containing the subjects that
%     need to be removed from the full list (i.e. "sub_list"). This
%     argument is only necessary for jackknife samples. If the algorithm
%     runs on the full subject list, pass in 'none'.
%
%   - FSM_file
%     A string. The full path of the functional similarity matrix (FSM)
%     file. It is assumed that a #subject x #subject matrix called "FSM" is
%     saved in this file. #subjects equals to the number of lines in "sub_list".
%
%   - outdir
%     A string. The full path of the output directory.
% 
%   - ystem
%     A string. The trait values will be read from "data_csv" and saved in
%     a .mat file: [outdir '/y_' ystem '.mat']. ystem is user-defined. For
%     example, if "trait_list" contains 13 cognitive behavioral names,
%     you can set ystem = '13cognitive'.
%   
%   - m2_subdir
%     A string. The subfolder name created to save the estimated "variance
%     explained" files. The "variance explained" files will be saved as 
%     [outdir '/' m2_subdir '/m2_QuantileNorm_' trait_name '.mat'] 
%     (for the traits with quantile normalization, i.e. non-binary traits), or
%     [outdir '/' m2_subdir '/m2_' trait_name '.mat'] 
%     (when the trait is NOT quantile normalized, i.e. binary trait (in our
%     case only "sex" is binary)).
% 
%   - verbose (optional)
%     0 (mute printout messages of this function) or 1 (print out every
%     message). Default is 1.
% 
% Written by Jingwei Li and CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md

addpath(fullfile(getenv('CBIG_CODE_DIR'), 'external_packages', 'LME'))

if(~exist('verbose', 'var') || isempty(verbose))
    verbose = 1;
end

if(~exist(outdir, 'dir'))
    mkdir(outdir)
end

%% Read the traits from CSV files
[y_names, num_y] = CBIG_text2cell(trait_list);
quantile_norm_flag = ones(num_y, 1);
for i = 1:num_y
    if(strcmp(y_names{i}, 'Gender'))
        y_types{i} = 'categorical';
        quantile_norm_flag(i) = 0;
    else
        y_types{i} = 'continuous';
    end
end

if(~isempty(ystem))
    ystem = ['_' ystem];
end
if(~exist([outdir '/y' ystem '.mat'], 'file'))
    CBIG_read_y_from_csv( {restricted_csv, unrestricted_csv}, 'Subject', y_names, y_types, ...
        sub_list, fullfile(outdir, ['y' ystem '.mat']), ',' );
end
load(fullfile(outdir, ['y' ystem '.mat']))

%% Read FSM
load(FSM_file);

%% Construct fixed-effects regressors
[cov_names, num_cov] = CBIG_text2cell(covariate_list);
for i = 1:num_cov
    if(strcmp(cov_names{i}, 'Gender') || strcmp(cov_names{i}, 'Race'))
        cov_types{i} = 'categorical';
    else
        cov_types{i} = 'continuous';
    end
end
cov_stem = ['_' strjoin(sort(cov_names), '_')];
if(~exist(fullfile(outdir, ['covariates' cov_stem '.mat']), 'file'))
    CBIG_generate_covariates_from_csv( {restricted_csv, unrestricted_csv}, ...
        'Subject', cov_names, cov_types, sub_list, FD_file, DVARS_file, ...
        fullfile(outdir, ['covariates' cov_stem '.mat']), ',' );
end
load(fullfile(outdir, ['covariates' cov_stem '.mat']))

%% Extract the subjects we care about (needed for jackknife)
if(~isempty(rm_sub_list) && ~strcmpi(rm_sub_list, 'none'))
    subject_name = CBIG_text2cell(sub_list);
    rm_sub_names = CBIG_text2cell(rm_sub_list);
    [sub_names_diff, I_not_rm] = setdiff(subject_name, rm_sub_names, 'stable');
    
    y = y(I_not_rm, :);
    covariates = covariates(I_not_rm, :);
    FSM = FSM(I_not_rm, I_not_rm);
end

% remove subjects with any empty trait data 
nan_mat = ~isnan(y);
prod_nan = prod(nan_mat, 2);
ind_pres = find(prod_nan ~= 0);

covariates = covariates(ind_pres, :);
y = y(ind_pres, :);
FSM = FSM(ind_pres, ind_pres);


%% Do real work: compute variance explained
% parameters: use default
alg = 0;
tol = 1e-4;
MaxIter = 100;
if(verbose == 1)
    fprintf('Parameters: alg = %d;   tol = %e;   MaxIter = %d\n\n', alg, tol, MaxIter);
end

mkdir(fullfile(outdir, m2_subdir))

if(sum(quantile_norm_flag) > 0)   
    % which means at least one trait needs to be quantile normalized
    eq_sam = linspace(0, 1, size(y, 1)+2);
    inv_sam = norminv(eq_sam, 0, 1);
    inv_sam = inv_sam(2:end-1);
end

% loop through each trait
for i = 1:length(y_names)
    if(verbose == 1)
        fprintf('Behavior: %s\n', y_names{i});
    end
    
    phenotype = y(:,i);
    
    if(quantile_norm_flag(i) == 1)
        [pheno_sort, pheno_sort_ind] = sort(phenotype, 'ascend');
        pheno_new(pheno_sort_ind) = inv_sam;
        if(size(pheno_new, 1) == 1)
            pheno_new = pheno_new';
        end
        phenotype = pheno_new;
    end
    
    [morpho.flag, morpho.m2, morpho.SE, morpho.Va, morpho.Ve, morpho.Lnew] = ...
        Morphometricity_1RandEffect(phenotype, covariates, FSM, alg, tol, MaxIter);
    
    if(verbose == 1)
        fprintf('\tVariance explained: %f;  SE: %f.\n\n', morpho.m2, morpho.SE);
    end
    
    if(quantile_norm_flag(i) == 1)
        save(fullfile(outdir, m2_subdir, ['m2_QuantileNorm_' y_names{i}]), 'morpho', '-v7.3');
    else
        save(fullfile(outdir, m2_subdir, ['m2_' y_names{i}]), 'morpho', '-v7.3');
    end
end

rmpath(fullfile(getenv('CBIG_CODE_DIR'), 'external_packages', 'LME'))

end

