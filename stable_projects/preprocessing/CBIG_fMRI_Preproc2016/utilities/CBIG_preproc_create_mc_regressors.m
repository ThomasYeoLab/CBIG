function CBIG_preproc_create_mc_regressors(mc_par_list, output_list, detrend_method, mt_diff_flag)

% CBIG_preproc_create_mc_regressors(mc_par_list, output_list, detrend_method)
%
% This function creates motion regressors based on given mc.par produced by
% fsl mcflirt. Because all frames are registered to reference frame, if original 
% data have 120 frames, the mc.par will have 121 frames. This function will
%
% 1) read in mc.par for each run and remove the first row (reference frame)
%    to get mc_dat, 
%
% 2) demean each column of mc_dat to get mc_rdat, take the derivative of
%    mc_dat to get mc_ddat
% 
% 3) merge mc_rdat and mc_ddat column by column to get mc_rddat, then detrend
%    mc_rddat to get final motion regressor
%
% Input:
%     - mc_par_list:
%       a text file including mc.par list of all runs, each line is the 
%       full path of the par file. For example:
%       mc_par_list = '/path_to_data_list/par_data_list.txt', each line is
%       '/path_to_par_data/rest_bld00?_mc.par'.
%
%     - output_list:
%       a text file including output file list of all runs, each line is the
%       full path of the output file name. For example:
%       output_list = '/path_to_output_list/output_list.txt', 
%       each line is '/path_to_output/rest_bld00?_mc_regressor.txt'.   
%
%     - detrend_method:
%       'detrend' or 'trendout' (default 'detrend'). 'detrend' option will
%       remove constant and linear trend. 'trendout' will only remove the 
%       linear trend, this option is only used for debugging. Please use
%       'detrend' option. 
%
%
% 
% Example:
% CBIG_preproc_create_mc_regressors(mc_par_list, output_list)
%
% Written by CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md

%% default detrend_method is 'detrend'
if (nargin < 3)
    detrend_method = 'detrend';
end
if (nargin < 4)
    mt_diff_flag = 1;
end

if (ischar(mt_diff_flag))
    mt_diff_flag = str2num(mt_diff_flag);
end

%% get mc par names and output file names
[par_name, num_of_par] = text2cell(mc_par_list);
[output_name, num_of_output] = text2cell(output_list);

if (num_of_output ~= num_of_par)
    error('ERROR: number of output file is not even with number of par file')
end

%% create mc regressor by following Thomas's code
for i = 1:num_of_par
    par = load(par_name{i});
    
    % step1
    % reorder the column of the mc parameter and remove the first row,
    % corresponding to _mc.dat, purpose of reordering is only for debugging
    par_order = par(:, [4 5 6 1 2 3]);
    mc_dat = par_order(2:end, :);
    
    % step2
    % demean the mc parameters and remove the first row, corresponding to _mc.rdat
    par_order_demean = bsxfun(@minus, par_order, mean(par_order));
    mc_rdat = par_order_demean(2:end, :);
    
    % take the derivative for mc_dat and set first row as zero,
    % corresponding to _mc.ddat
    mc_ddat = [zeros(1,size(mc_dat,2)); diff(mc_dat)];
    
    % step3
    % concat two files together, corresponding to _mc.rddat
    if(mt_diff_flag == 1)
        mc_rddat = [mc_rdat mc_ddat];
    else
        mc_rddat = mc_rdat;
    end
        
    % detrend the regressors
    if (strcmp(detrend_method, 'trendout'))
        mc_out = CBIG_preproc_trendout(mc_rddat);   
    elseif (strcmp(detrend_method, 'detrend'))    
        [mc_out, ~, ~, ~] = CBIG_glm_regress_matrix(mc_rddat, [], 1, []);
    end
    
    % write into output file 
    dlmwrite(output_name{i}, mc_out, ' ');
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

