function CBIG_preproc_aCompCor_multipleruns(runs_list, ROI_file, out_files_list, nPCs, per_run, diff_flag)

% CBIG_preproc_aCompCor_multipleruns(runs_list, ROI_file, out_files_list, nPCs, per_run)
%
% Compute aCompCor regressors for all runs of one subjects. Detrend
% regressors. Write regressors into text files.
%
% Inputs:
%   - runs_list:
%     The filename (full path) of a list containing all input BOLD filenames.
%     E.g. '<sub_dir>/Sub0001_Ses1/bold/regression/fMRI_list.txt'
%
%     In this text file, each line is the filename of the BOLD file of one run. 
%     An example of this text file:
%     <sub_dir>/Sub0001_Ses1/bold/002/Sub0001_Ses1_bld002_rest_skip4_stc_mc.nii.gz
%     <sub_dir>/Sub0001_Ses1/bold/003/Sub0001_Ses1_bld003_rest_skip4_stc_mc.nii.gz
%
%   - ROI_file:
%     The filename of wm+ventricles mask (full path).
%     E.g. '<sub_dir>/Sub0001_Ses1/bold/mask/Sub0001_Ses1.func.wm.vent.nii.gz'
%
%   - out_files_list:
%     The filename (full path) of a list containing all output aCompCor
%     regressors filenames.
%     E.g. '<sub_dir>/Sub0001_Ses1/bold/regression/aCompCor_regressors_list.txt'
% 
%     If per_run == 1, this list has multiple lines. Each line is the 
%     regressor filename of one run. If per_run == 0, this list has only
%     one line, which is the filename of regressors of merged runs.
%     An example of this list:
%     <sub_dir>/Sub0001_Ses1/bold/regression/Sub0001_Ses1_bld002_aCompCor_regressor.txt
%     <sub_dir>/Sub0001_Ses1/bold/regression/Sub0001_Ses1_bld003_aCompCor_regressor.txt
%
%   - nPCs:
%     Number of principle components extracted from wm+ventricles signals,
%     e.g. '5' or 5. Default is 5.
%
%   - per_run:
%     0 (or '0') OR 1 (or '1'). If per_run == 1 (or '1'), generate aCompCor
%     regressors for each run separately; if per_run == 0 (or '0'),
%     concatenate all runs together, and generate aCompCor regressors of
%     merged runs. 
%
% Example:
%   CBIG_preproc_aCompCor_multipleruns('<sub_dir>/Sub0001_Ses1/bold/regression/fMRI_list.txt', ...
%          '<sub_dir>/Sub0001_Ses1/bold/mask/Sub0001_Ses1.func.wm.vent.nii.gz', ...
%          '<sub_dir>/Sub0001_Ses1/bold/regression/aCompCor_regressors_list.txt', ...
%          '5', '1')
%
% Written by Jingwei Li.
% Written by CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md


%%%%%%%%%%%%%%%%%%%%%%%%%
% Setup parameters
%%%%%%%%%%%%%%%%%%%%%%%%%
if (nargin < 6)
    diff_flag = 1;
end
if (nargin < 5)
    per_run = 1;
end
if (nargin < 4)
    nPCs = 5;               % default is to pick 5 components
end

if(ischar(diff_flag))
    diff_flag = str2num(diff_flag);
end
if(ischar(per_run))
    per_run = str2num(per_run);
end
per_run = logical(per_run);
if(ischar(nPCs))
    nPCs = str2num(nPCs);
end

[runs_file, nruns] = text2cell(runs_list);
[out_files, nout] = text2cell(out_files_list);


%%%%%%%%%%%%%%%%%%%%%%%%%
% Compute PCs and save out
%%%%%%%%%%%%%%%%%%%%%%%%%
if(per_run == 0)            % combine runs together and compute PCs
    % check the number of output files
    if(nout ~= 1)
        error('There should only be one output file if all runs are combined.');
    end
    
    % merge all runs
    merged_file = [dirname(runs_list) '/merged_func_before_regression.nii.gz'];
    merge_cmd = ['fslmerge -t ' merged_file ' ' sprintf('%s.nii.gz ', runs_file{:})];
    disp(merge_cmd);
    [merge_status merge_reslt] = system(merge_cmd);
    if (merge_status == 0)
        disp('Successfully concatenated.')
    else
        error('Unable to concatenate.')
    end
    
    % aCompCor for combined runs
    [ROI_PCs] = CBIG_preproc_aCompCor(merged_file, ROI_file, nPCs, diff_flag);
    
    % save out the aCompCor regressors for combined runs
    fid = fopen(out_files{1}, 'w+');
    if(isempty(find(isnan(ROI_PCs))))
        fprintf(fid, [repmat('%f\t', [1 size(ROI_PCs, 2)]) '\n'], ROI_PCs');
    else
        fprintf('WARNING: The wm + ventricles mask is empty. The aCompCor regressor is empty.\n');
    end
    fclose(fid);
    
else                        % multiple runs but compute PCs for each run OR only single run
    % check the number of input runs and output files
    if(nruns ~= nout)
        error('The number of runs should be equal to the number of output files if aCompCor is applied on each run separately.');
    end
    
    % loop through each run
    for i = 1:nruns
        % aCompCor for each run
        [ROI_PCs] = CBIG_preproc_aCompCor(runs_file{i}, ROI_file, nPCs, diff_flag);
        
        % save out the aCompCor regressors for each run
        fid = fopen(out_files{i}, 'w+');
        if(isempty(find(isnan(ROI_PCs))))
            fprintf(fid, [repmat('%f\t', [1 size(ROI_PCs, 2)]) '\n'], ROI_PCs');
        else
            fprintf('WARNING: The wm + ventricles mask is empty. The aCompCor regressor is empty.\n');
        end
        fclose(fid);
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
end