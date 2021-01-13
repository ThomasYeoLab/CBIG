function status = CBIG_check_job_status(job, varargin)

% status = CBIG_check_job_status(job, varargin)
%
% This function checks whether your job is successful. It will first check
% whether your jobs are running or not. After they are all finished, it
% will check your log file or output file. You can check either, both or 
% neither of the two files (which we do not recommend).
% For log file check, if the log file includes 'FAILED' or does NOT include
% 'SUCCESS', then the job will be considered unsuccessful. To notify this 
% function about the status of your job, you need to manually print a 
% 'SUCCESS' message at the end of your code. If this is not applicable, you
% can choose to check the existance of the output file. If all jobs are 
% successfully completed, this function will return 0. Otherwise it will 
% return 1 and print an error message.
%
% This function uses varargin to handle optional inputs. Please use 'log'
% and 'out' to indicate your options. 
% If neither check is selected, this function will simply wait till all the
% jobs are completed. We suggest you choose at least one check to avoid
% unexpected job failure. 
%
% Input:
%     - job:
%       Name or ID of your job (char). If you choose job name, all jobs
%       with the same name will be checked. If you choose job id, you can
%       check one single job or multiple jobs. You can connect the job id 
%       with '|', for example, '100001|100002'. Or you can use a cell
%       array. 
%       Note that for our current scheduler, the length of the job name is
%       limited to 16 characters. Longer name will be cut down. Please be
%       careful.
%
%     - log (optional):
%       Path of the log file or a cell array of all log files. 
%
%     - out (optional):
%       Path of the output file or a cell array of all output files. 
%
% Output:
%     - status:
%       When job is successful, status = 0. Otherwise status = 1.
%
% Example:
%     status = CBIG_check_job_status('Unit_test_job', 'log', log_list);
%     status = CBIG_check_job_status('100001|100002', 'log', log_list, 'out', out_list);
%
% Written by XUE Aihuiping and CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md

% Handle optional inputs. Default values are empty. 
pnames = {'log' 'out'};
dflts =  {[] []};
[log_file, out_file] = internal.stats.parseArgs(pnames, dflts, varargin{:});

if(isempty(log_file) && isempty(out_file))
    warning('Both log file and output file are NOT checked.');
end

% Convert cell job id to string
if(iscell(job))
    temp = [];
    for i = 1:length(job)
        temp = [temp '|' job{i}];
    end
    job = temp(2:end);
end

% Wait till all jobs are completed
status = 0;
[~, whoami] = system('whoami');
whoami = whoami(1:end-1);
% Count the running job number
cmd = ['qstat -a | grep -E ''' job ''' | grep ' whoami ' | wc -l'];
[~, job_num] = system(cmd);
job_num = str2num(job_num(1: end-1));
while(job_num)
    pause(60);
    cmd = ['qstat -a | grep -E ''' job ''' | grep ' whoami ' | wc -l'];
    [~, job_num] = system(cmd);
    job_num = str2num(job_num(1: end-1));
end

% Check logs
if(~isempty(log_file))
    if(~iscell(log_file))
        temp{1} = log_file;
        log_file = temp;
    end
    for i=1:length(log_file)
        [~, error_messages] = system(['cat ' log_file{i} ' | grep FAILED']);
        [~, success_messages] = system(['cat ' log_file{i} ' | grep SUCCESS']);
        if(~isempty(error_messages) || isempty(success_messages))
            fprintf(['[FAILED] Job is unsuccesful. Log file: ' log_file{i} ' \n']);
            status = 1;
        end
    end
end

% Check output files
if(~isempty(out_file))
    if(~iscell(out_file))
        temp{1} = out_file;
        out_file = temp;
    end
    for i=1:length(out_file)
        if(~exist(out_file{i}, 'file'))
            fprintf(['[FAILED] Output file ' out_file{i} ' does not exist. \n']);
            status = 1;
        end
    end
end

end
