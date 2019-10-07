function CBIG_ASDf_checkJobStatus(progress_file, num_jobs, lag_time)
% CBIG_ASDf_checkJobStatus(progress_file, num_jobs, lag_time)
% 
% This is a helper function that holds off the next jobs until submitted jobs
% are finished. Useful when the next job is dependent on the submitted
% jobs.
%
% Input:
%     - progress_file:
%           Text file. When a job is finished, a line should be written to
%           this text file. When all submitted jobs are finished, the
%           number of text lines in progress_file should be the same as
%           num_jobs.
%     - num_jobs:
%           Integer or string. Number of jobs.
%     - lag_time:
%           Integer or string. Time lag between each check of finished job, in seconds. 
% 
% Example:
%	CBIG_ASDf_checkJobStatus('~/example_output/estimate/progressFile.txt',100,600)
%	This function will check whether the 100 jobs have finished by checking the output in progressFile.txt every 10min.
%
% Written by CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md

if ischar(num_jobs)
    num_jobs = str2double(num_jobs);
end

if ischar(lag_time)
    lag_time = str2double(lag_time);
end

fid = fopen(progress_file);
txts = textscan(fid,'%s','delimiter','\n');
fclose(fid);
num_jobs_done = length(txts{1});
while num_jobs_done < num_jobs
    fprintf('%d/%d jobs completed\n', num_jobs_done, num_jobs);
    pause(lag_time);
    fid = fopen(progress_file);
    txts = textscan(fid,'%s','delimiter','\n');
    fclose(fid);
    num_jobs_done = length(txts{1});
end
