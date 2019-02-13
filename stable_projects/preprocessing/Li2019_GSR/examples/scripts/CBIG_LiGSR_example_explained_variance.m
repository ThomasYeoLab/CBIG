function CBIG_LiGSR_example_explained_variance( data_csv1, data_csv2, ...
    sub_list, rm_sub_list, FSM_file, outdir, m2_subdir )

% CBIG_LiGSR_example_explained_variance( data_csv1, data_csv2, behavior_list, ...
%     covariate_list, sub_list, rm_sub_list, FSM_file, outdir, m2_subdir )
% 
% This is a wrapper function to call Morphometricity_1RandEffect.m to
% estimate the explained variance of the behaviors in the faked toy
% example. 
% Note that the wrapper function is supposed to be dataset specific. Hence
% this function is only applicable to our example data. For the wrapper
% functions of the Brain Genomics Superstruct Project (GSP) and the Human
% Connectome Project (HCP), which were utilized in our paper, please refer
% to
% ../../VarianceComponentModel/scripts/CBIG_LiGSR_explained_variance_GSP.m
% and
% ../../VarianceComponentModel/scripts/CBIG_LiGSR_explained_variance_HCP.m.
% 
% Input:
%   - data_csv1
%     In this toy example, we faked two CSV files. "data_csv1" is the
%     full-path name of the first CSV file. No compulsory ordering of these
%     two CSV files required.
%     In this example, data_csv1 contains the behavioral scores of each
%     fake subject.
% 
%   - data_csv2
%     In this toy example, we faked two CSV files. "data_csv2" is the
%     full-path name of the second CSV file.
%     In this example, data_csv2 contains the demographic measures of each
%     fake subject.
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
%   - m2_subdir
%     A string. The subfolder name created to save the estimated "variance
%     explained" files. The "variance explained" files will be saved as 
%     [outdir '/' m2_subdir '/m2_QuantileNorm_' behavior_name '.mat'] 
%     (for quantile normalization), or
%     [outdir '/' m2_subdir '/m2_' behavior_name '.mat'] 
%     (when the trait is NOT quantile normalized).
% 
% Written by Jingwei Li and CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md

addpath(fullfile(getenv('CBIG_CODE_DIR'), 'external_packages', 'LME'))

quantile_norm_flag = 1;
if(~exist(outdir, 'dir'))
    mkdir(outdir)
end

%% Read the traits from data_csv
y_names = {'Behavior_1', 'Behavior_2'};
y_types = {'continuous', 'continuous'};

if(~exist(fullfile(outdir, 'y.mat'), 'file'))
    CBIG_read_y_from_csv( {data_csv1, data_csv2}, 'Subject', y_names, y_types, ...
        sub_list, fullfile(outdir, 'y.mat'), ',' );
end
load(fullfile(outdir, 'y.mat'))

%% Read FSM
load(FSM_file);

%% Construct fixed-effects regressors
cov_names = {'Age', 'Sex'};
cov_types = {'continuous', 'categorical'};

if(~exist(fullfile(outdir, 'covariates.mat'), 'file'))
    CBIG_generate_covariates_from_csv( {data_csv1, data_csv2}, ...
        'Subject', cov_names, cov_types, sub_list, 'none', 'none', ...
        fullfile(outdir, 'covariates.mat'), ',' );
end
load(fullfile(outdir, 'covariates.mat'))

%% Extract the subjects we care about (needed for jackknife)
if(~isempty(rm_sub_list) && ~strcmpi(rm_sub_list, 'none'))
    subject_name = CBIG_text2cell(sub_list);
    rm_sub_names = CBIG_text2cell(rm_sub_list);
    [sub_names_diff, I_not_rm] = setdiff(subject_name, rm_sub_names, 'stable');
    
    y = y(I_not_rm, :);
    covariates = covariates(I_not_rm, :);
    FSM = FSM(I_not_rm, I_not_rm);
end

% remove subjects with any empty behavioral data 
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
fprintf('Parameters: alg = %d;   tol = %e;   MaxIter = %d\n\n', alg, tol, MaxIter);

mkdir([outdir '/' m2_subdir])

if(quantile_norm_flag == 1)
    eq_sam = linspace(0, 1, size(y, 1)+2);
    inv_sam = norminv(eq_sam, 0, 1);
    inv_sam = inv_sam(2:end-1);
end

% loop through each behavior
for i = 1:length(y_names)
    fprintf('Behavior: %s\n', y_names{i});
    
    phenotype = y(:,i);
    
    if(quantile_norm_flag == 1)
        [pheno_sort, pheno_sort_ind] = sort(phenotype, 'ascend');
        pheno_new(pheno_sort_ind) = inv_sam;
        if(size(pheno_new, 1) == 1)
            pheno_new = pheno_new';
        end
        phenotype = pheno_new;
    end
    
    [morpho.flag, morpho.m2, morpho.SE, morpho.Va, morpho.Ve, morpho.Lnew] = ...
        Morphometricity_1RandEffect(phenotype, covariates, FSM, alg, tol, MaxIter);
    
    fprintf('\tVariance explained: %f;  SE: %f.\n\n', morpho.m2, morpho.SE);
    
    if(quantile_norm_flag == 1)
        save(fullfile(outdir, m2_subdir, ['m2_QuantileNorm_' y_names{i}]), 'morpho', '-v7.3');
    else
        save(fullfile(outdir, m2_subdir, ['m2_' y_names{i}]), 'morpho', '-v7.3');
    end
end

rmpath(fullfile(getenv('CBIG_CODE_DIR'), 'external_packages', 'LME'))

end

