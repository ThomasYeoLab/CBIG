function CBIG_glm_regress_vol(fMRI_list, output_list, regressor_list, polynomial_fit, censor_list, per_run, save_beta)

% CBIG_glm_regress_vol(fMRI_list, output_list, regressor_list, polynomial_fit, censor_list, per_run)
%
% This function does GLM regression for fMRI volumetric data. It will:
%
% 1) read in fMRI volumetric data from fMRI data list (fMRI_list) and
%    reshape each input 4D volume to a TxN matrix, where T is number of time
%    points of each timecourse, N is the number of voxels.
%
% 2) do GLM regression for each TxN matrix via CBIG_glm_regress_mtx.m
%
% 3) reshape each GLM residual matrix back to 4D volume and write it
%    out as a new volumetric file. For example, a .nii.gz file.
%    each element of the output_cell is the volumetric file name.
%
% NOTE: For polynomial_fit and censor_list options, please <help
% CBIG_glm_regress_mtx> to see details. For the per_run option, if it's off,
% this function will first concatenate all fMRI data along the time dimension,
% then do the regression jointly. If it's on, this function will do
% regression for each fMRI data seperately.
%
% Input:
%     - fMRI_list:
%       a text file including fMRI data list of all runs, each line is the
%       full path of the fMRI file. For example:
%       fMRI_list = '/path_to_data_list/fMRI_data_list.txt', each line is
%       '/path_to_fMRI_data/rest_bld00?.nii.gz'.
%
%     - output_list:
%       a text file including output file list of all runs, each line is the
%       full path of the output file name. For example:
%       output_list = '/path_to_output_list/output_list.txt',
%       each line is '/path_to_output/rest_bld00?_resid.nii.gz'.
%
%     - regressor_list:
%       a text file incuding the regressor list of all runs, each line is
%       the full path of the regressor file. For the regressor file, each
%       column represents a regressor, each row is a time point.
%       For example:
%       regressor_list = '/path_to_regressor_list/regressor_list.txt',
%       each line is '/path_to_regressor/bld00?_regressors.dat'.
%
%     - polynomial_fit:
%       '-1'/'0'/'1'/ (default is '0'). If polynomial_fit is set to -1, we add
%       nothing in regressor matrix. If polynomial_fit is set to 0, we prepend a
%       Mx1 vector [1,1,1...1]' to the regressor matrix and end up with a
%       Mx(K+1) matrix. If polynomial_fit is set to 1, a Mx1 vector
%       [1,1,1...1]' will be added into the first column and a Mx1 matrix
%       linspace(-1, 1, M)' will be added into the second column of the
%       regressor matrix.
%
%     - censor_list:
%       a text file including the censor list of all runs, each line is
%       the full path of the censor file. For the censor file, each line is
%       a time point, 1 means this time point is not censored, 0 means the
%       time point is censored.
%
%     - per_run:
%       '0'/'1' (default is '1'). If per_run is set to '0', this function
%       will first concatenate all fMRI data, regressor, censor vector along
%       the time dimension, then do the regression jointly. If per_run is
%       set to '1', this function will do regression for each fMRI data
%       seperately.
%
%     - save_beta:
%       '0'/'1' (default is '0'). If save_beta is set to '1', this function
%       will save out the beta coefficients to a nifti file.
%
% Example:
% CBIG_glm_regress_vol(fMRI_list, output_list, regressor_list, '1', censor_list, '1')
%
% Written by CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md

%% default using per_run option
if (nargin < 6)
    per_run = 1;
end

% default not save beta
if (nargin < 7)
    save_beta = 0;
end

%% fMRI_list is a list including all fMRI file name
[fMRI_name, num_of_fMRI] = text2cell(fMRI_list);

%% output_list is a list including all output file name
[output_name, num_of_output] = text2cell(output_list);

if (num_of_output ~= num_of_fMRI)
    error('ERROR: number of output files is not even with the number of fMRI file.')
end

%% regressor_list is a list including all regressor's name
[regressor_name, num_of_regressor] = text2cell(regressor_list);

if (num_of_regressor ~= num_of_fMRI)
    error('ERROR: number of regressor files is not even with the number of fMRI file.')
end

%% polynomial_fit option is a string or not
if (ischar(polynomial_fit))
    polynomial_fit = str2num(polynomial_fit);
end

%% censor_file is a name of censor file or a censor vector
if (~isempty(censor_list))
    [censor_name, num_of_censor] = text2cell(censor_list);
    
    if (num_of_censor ~= num_of_fMRI)
        error('ERROR: number of censor files is not even with the number of fMRI file.')
    end
end



%% format input flags
if (ischar(per_run))
    per_run = str2num(per_run);
end

if (ischar(save_beta))
    save_beta = str2num(save_beta);
end

%% start doing regression
if (per_run == 0)
    % If per_run is set to '0', regress all runs jointly
    vol_2d_all = [];
    regressor_all = [];
    censor_all = [];
    tp_length_all = [];
    reg_name_str = regressor_name{1};
    for i = 1:num_of_fMRI
        if (isempty(strfind(fMRI_name{i}, '.dtseries.nii')))
            % if input volume is nifti file
            mri = MRIread(fMRI_name{i});
            vol_2d = single(mri.vol);
            mri.vol = [];
            mri_size = size(vol_2d);
            vol_2d = reshape(vol_2d, [mri_size(1)*mri_size(2)*mri_size(3) mri_size(4)]);
        else
            % if input volume is cifti file
            mri = ft_read_cifti(fMRI_name{i});
            mri_size = size(mri.dtseries);
            vol_2d = single(mri.dtseries);
            mri.dtseries = [];
        end
        % time points length
        tp_length = size(vol_2d, 2);
        tp_length_all = [tp_length_all tp_length];
        
        vol_2d_all = [vol_2d_all vol_2d];
        
        % concatenate all regressors
        regressor_mtx = load(regressor_name{i});
        if(i>1)
            reg_name_str = [reg_name_str ', ' regressor_name{i}];
        end
        regressor_all = [regressor_all; regressor_mtx];
        
        
        % concatenate all censor vectors
        if (~isempty(censor_list))
            censor_vec = load(censor_name{i});
            censor_all = [censor_all; censor_vec];
        end
    end
    vol_2d_all = transpose(vol_2d_all);
    if ~isempty(censor_all)
        zero_ind = find(sum(abs(regressor_all(censor_all == 1, :)), 1)==0);
        flag_zero = sum(abs(regressor_all(censor_all == 1)), 1)==0;
    else
        zero_ind = find(sum(abs(regressor_all), 1)==0);
        flag_zero = sum(abs(regressor_all), 1)==0;
    end
    if(~isempty(zero_ind))
        regressor_all(:, zero_ind) = [];
        zero_ind_str = num2str(zero_ind(1));
        for c = 2:length(zero_ind)
            zero_ind_str = [zero_ind_str ', ' num2str(zero_ind(c))];
        end
        warning('The columns %s in files %s are all zeros. They will be excluded from the regressors matrix.', ...
            zero_ind_str, reg_name_str);
    end
    [resid_mtx, beta, ~, ~] = CBIG_glm_regress_matrix(vol_2d_all, regressor_all, polynomial_fit, censor_all);
    
    
    for i=1:num_of_fMRI
        if (i == 1)
            res = resid_mtx(1:tp_length_all(1), :);
        elseif (i > 1)
            res = resid_mtx(sum(tp_length_all(1:i-1))+1:sum(tp_length_all(1:i)), :);
        end
        
        res = transpose(res);
        
        if (isempty(strfind(fMRI_name{i}, '.dtseries.nii')))
            % if output volume is nifti file
            mri.vol = reshape(res, mri_size);
            MRIwrite(mri, output_name{i});
        else
            % if output volume is cifti file
            output_name{i} = regexprep(output_name{i}, '.dtseries.nii', '');
            mri.dtseries = reshape(res, mri_size);
            ft_write_cifti(output_name{i}, mri, 'parameter', 'dtseries');
        end
    end
    
    if save_beta == 1
        beta = beta';
        beta_all_regressor = zeros(size(beta,1),length(flag_zero));
        % add beta of non zero regressors
        beta_all_regressor(:,~flag_zero) = beta(:,polynomial_fit+2:end);
        % add beta of polynomial trends
        beta_all_regressor = [beta(:,1:polynomial_fit+1) beta_all_regressor];
        
        % save beta to image file
        beta_4d = reshape(beta_all_regressor, [mri_size(1),mri_size(2),mri_size(3),size(beta_all_regressor,2)]);
        beta_file = strrep(output_name{i},'.nii','_beta.nii');
        if (isempty(strfind(fMRI_name{i}, '.dtseries.nii')))
            % if output volume is nifti file
            mri.vol = beta_4d;
            MRIwrite(mri, beta_file);
        else
            % if output volume is cifti file
            mri.dtseries = beta_4d;
            ft_write_cifti(beta_file, mri, 'parameter', 'dtseries');
        end
    end
    
    
    
elseif (per_run == 1)
    % If per_run is set to '1', regress each run seperately
    for i = 1:num_of_fMRI
        if isempty(strfind(fMRI_name{i}, '.dtseries.nii'))
            % if input volume is nifti file
            mri = MRIread(fMRI_name{i});
            vol_2d = single(mri.vol);
            mri.vol = [];
            mri_size = size(vol_2d);
            vol_2d = reshape(vol_2d, [mri_size(1)*mri_size(2)*mri_size(3) mri_size(4)]);
        else
            % if input volume is cifti file
            mri = ft_read_cifti(fMRI_name{i});
            mri_size = size(mri.dtseries);
            vol_2d = single(mri.dtseries);
            mri.dtseries = [];
        end
        
        % load regressor file
        regressor_mtx = load(regressor_name{i});
        
        censor_vec = [];
        
        % load censor file
        if (~isempty(censor_list))
            censor_vec = load(censor_name{i});
        end
        vol_2d = transpose(vol_2d);
        if ~isempty(censor_vec)
            zero_ind = find(sum(abs(regressor_mtx(censor_vec == 1, :)), 1)==0);
            flag_zero = sum(abs(regressor_mtx(censor_vec == 1, :)), 1)==0;
        else
            zero_ind = find(sum(abs(regressor_mtx), 1)==0);
            flag_zero = sum(abs(regressor_mtx), 1)==0;
        end
        if(~isempty(zero_ind))
            regressor_mtx(:, zero_ind) = [];
            zero_ind_str = num2str(zero_ind(1));
            for c = 2:length(zero_ind)
                zero_ind_str = [zero_ind_str ', ' num2str(zero_ind(c))];
            end
            warning('The columns %s in file %s are all zeros. They will be excluded from the regressors matrix.', ...
                zero_ind_str, regressor_name{i});
        end
        [resid_mtx, beta, ~, ~] = CBIG_glm_regress_matrix(vol_2d, regressor_mtx, polynomial_fit, censor_vec);
        resid_mtx = transpose(resid_mtx);
        
        if isempty(strfind(fMRI_name{i}, '.dtseries.nii'))
            % if input volume is nifti file
            mri.vol = reshape(resid_mtx, mri_size);
            MRIwrite(mri, output_name{i});
        else
            % if input volume is cifti file
            output_name{i} = regexprep(output_name{i}, '.dtseries.nii', '');
            mri.dtseries = reshape(resid_mtx, mri_size);
            ft_write_cifti(output_name{i}, mri, 'parameter', 'dtseries');
        end
        
        if save_beta == 1
            beta = beta';
            beta_all_regressor = zeros(size(beta,1),length(flag_zero));
            % add beta of non zero regressors
            beta_all_regressor(:,~flag_zero) = beta(:,polynomial_fit+2:end);
            % add beta of polynomial trends
            beta_all_regressor = [beta(:,1:polynomial_fit+1) beta_all_regressor];
            
            % save beta to image file
            beta_4d = reshape(beta_all_regressor, [mri_size(1),mri_size(2),mri_size(3),size(beta_all_regressor,2)]);
            beta_file = strrep(output_name{i},'.nii','_beta.nii');
            if (isempty(strfind(fMRI_name{i}, '.dtseries.nii')))
                % if output volume is nifti file
                mri.vol = beta_4d;
                MRIwrite(mri, beta_file);
            else
                % if output volume is cifti file
                mri.dtseries = beta_4d;
                ft_write_cifti(beta_file, mri, 'parameter', 'dtseries');
            end
        end
        
    end
end




function [cell_array, length] = text2cell(text)
% read text file line by line and write it into cell_array
length = 0;
fid = fopen(text);
while (~feof(fid))
    length = length + 1;
    cell_array{length} = fgetl(fid);
end
fclose(fid);

