function CBIG_waitUntilFinished(progressFile, noJobs)

% CBIG_waitUntilFinished(progressFile, noJobs)
%
% This is a helper function that holds off the next jobs until submitted
% jobs are finished. This is useful when the next batch of jobs depend on
% this submitted batch.
%
% Input:
%     - progressFile:
%       Text file
%       When a job is completed, a line will be written to this file
%     - noJobs:
%       Integer
%   
% Example:
% See ../replicatePNAS/CBIG_LDA_wrapper.m
%
% Written by Xiuming Zhang and CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md

noJobsDone_prev = -1;
noJobsDone = 0;
while noJobsDone < noJobs
    if noJobsDone ~= noJobsDone_prev
        fprintf('%d/%d jobs completed\n', noJobsDone, noJobs);
    end
    pause(600); % pause for ten minutes
    noJobsDone_prev = noJobsDone;
    [~, noJobsDone] = system(['wc ' progressFile ' -l']);
end
