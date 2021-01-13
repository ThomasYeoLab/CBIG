function [pr_process, act_contribution] = CBIG_ComputePrProcessGivenAct(act, beta, threshold)

% [pr_process, act_contribution] = CBIG_ComputePrProcessGivenAct(act, beta, threshold)
%
% act  = 1 x 18715
% beta = # processes (K) x 18715 
%
% Written by B.T.Thomas Yeo and CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md

if size(act,1) ~= 1
    error('Input argument ''act'' should be a row vector');
end

if(nargin < 3)
   threshold = 1e-4; 
end

num_processes = size(beta, 1);
curr_pr_process = zeros(num_processes, 1)+1/num_processes; % K x 1

for i = 1:1000
    
    q = bsxfun(@times, beta, curr_pr_process);
    q = bsxfun(@times, q, 1./sum(q, 1));
    
    pr_process = sum(bsxfun(@times, q, act), 2);
    pr_process = pr_process./sum(pr_process);
    
    change = max(abs(curr_pr_process - pr_process));
    curr_pr_process = pr_process;
    
    if(change < threshold)
        disp(['iter ' num2str(i) ': ' num2str(change)]);
        break;
    end
end

% Compute activation contribution to pr_process
if(nargout > 1)
    curr_pr_process = zeros(num_processes, 1)+1/num_processes; % K x 1
    q = bsxfun(@times, beta, curr_pr_process);
    q = bsxfun(@times, q, 1./sum(q, 1));
    act_contribution = bsxfun(@times, q, act);
    act_contribution = act_contribution/sum(act_contribution(:));
end


