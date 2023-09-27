function CBIG_bandpass_vol(input_vol, output_vol, low_f, high_f, detrend, retrend, censor_file, sample_period)

% CBIG_bandpass_vol(input_vol, output_vol, low_f, high_f, detrend, retrend, censor_file, sample_period)
%
% Given an input 4D volume file name (input_vol): nifti (.nii) or cifti 
% (.dtseries.nii), and sample period (sample_period), this function will:
%
% 1) reshape input 4D volume to a TxN matrix, where T is number of time 
%    points of each timecourse, N is the number of voxels.
%
% 2) demean and remove the linear trend of each timecourse by setting 
%    detrend to '1'. This step will be skipped by setting detrend 
%    to '0'. 
%
% 3) use CBIG_bpss_by_regression.m to do bandpass filtering. It will use
%    GLM regression to do bandpass filtering. For more details,
%    <help CBIG_bpss_by_regression>
%
% 4) add back the linear trend of each time course by setting retrend to
%    '1'. NOTE: retrend is not allowed to perform if detrend is off.
%
% Input:
%     - input_vol:
%       a 4D volume file name, it should be full path of the file. For example,
%       '/path_input_data/input_data.nii'.
%
%     - output_vol:
%       a 4D volume file name, it should be full path of the file. For example,
%       '/path_output_bp_data/output_bp_data.nii'.
%
%     - low_f:
%       a string, low cutoff frequency, e.g. '0.001'.
%
%     - high_f:
%       a string, high cutoff frequency, e.g. '0.08'.
%
%     - sample_period:
%       a string, the sample period (second), e.g. '2'. If sample_period 
%       is not given, the function will get sample period from TR of the 
%       input file.
%
%     - detrend:
%       '0'/'1'. If detrend is set to '1', it will demean and remove the 
%       linear trend of each time course before bandpass filtering. If
%       detrend is set to '0', this step will be skipped.
%
%     - retrend:
%       '0'/'1'. If retrend is set to '1', it will add back the linear 
%       trend of each time course after bandpass filtering. Retrend is not 
%       allowed to set to '1' if detrend is set to '0'.
%
%     - censor_file:
%       censor file name (default ''). If censor_file is '', this function will
%       not remove censored frames when it uses GLM to do bandpass filtering.
%       If censor_file is not empty, a vector with 1 or 0 will be loaded from 
%       censor_file. In the censor vector, 1 means kept frame, 0 means censored
%       frame. Then, we will only use kept frames to do GLM regression, and
%       then apply the coefficent to all frames.
%  
% Example:
% CBIG_bandpass_vol('/path_input_data/input_data.nii', '/path_output_bp_data/output_bp_data.nii', '0.001', '0.08', '1', '1', '', '2')
%
% Written by CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md

%% check whether the variable exists or not
if (~exist('input_vol'))
    error('ERROR: input_vol does not exist.')
end

if (~exist('output_vol'))
    error('ERROR: output_vol does not exist.')
end

if (exist('low_f'))
    low_f = str2num(low_f);    
else
    error('ERROR: low_f does not exist.')
end

if (exist('high_f'))
    high_f = str2num(high_f);
else
    error('ERROR: high_f does not exist.')
end

if (exist('detrend'))
    detrend = str2num(detrend);
else
    error('ERROR: detrend does not exist.')
end

if (exist('retrend'))
    retrend = str2num(retrend);
else
    error('ERROR: retrend does not exist.')
end


%% load fMRI image and reshape it into 2D

if (isempty(strfind(input_vol, '.dtseries.nii')))
    % if input_vol is a nifti file
    mri = MRIread(input_vol);
    mri_vol = single(mri.vol);
    mri.vol = [];
    mri_size = size(mri_vol);
    data = transpose(reshape(mri_vol, [mri_size(1)*mri_size(2)*mri_size(3) mri_size(4)]));
else
    % if input_vol is a cifti file
    mri = ft_read_cifti(input_vol);
    mri_size = size(mri.dtseries);
    data = single(transpose(mri.dtseries));
    mri.dtseries = [];
end

%% if censor_file is not set, censor is empty
if (~exist('censor_file')) 
    error('ERROR: censor_file does not exist.')
else 
    if (isempty(censor_file))
        censor = [];
    else
        censor = load(censor_file);
    end
end

%% if sample period is not set, assume TR is sample period
if (~exist('sample_period'))
    
    if (isempty(strfind(input_vol, '.dtseries.nii')))
        % if input_vol is a nifti file
        sample_period = mri.tr/1000;
    else
        % if input_vol is a cifti file
        sample_period = mri.time(2) - mri.time(1);
    end
else 
    sample_period = str2num(sample_period);
end

%% detrend each time course before bandpass
if (detrend)
    fprintf('Detrend each time course.\n')
    [data, ~, ~, retrend_mtx] = CBIG_glm_regress_matrix(data, [], 1, censor);
end

%% using GLM regression to do bandpass
data = CBIG_bpss_by_regression(data, low_f, high_f, sample_period, censor);

%% after bandpass, retrend each time course
if (~detrend && retrend)
    error('ERROR: retrend can not be done because detrend option is off');
elseif (detrend && retrend)
    fprintf('Retrend each time course.\n')
    data = data + retrend_mtx;
end

%% write the mri into a volume
data = transpose(data);
data = reshape(data, mri_size);

if (isempty(strfind(input_vol, '.dtseries.nii')))
    % if input_vol is a nifti file
    mri.vol = data;
    MRIwrite(mri, output_vol);
else
    % if input_vol is a cifti file
    output_vol = regexprep(output_vol, '.dtseries.nii', '');
    mri.dtseries = data;
    ft_write_cifti(output_vol, mri, 'parameter', 'dtseries');
end
