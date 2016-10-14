function pr_process_given_act = CBIG_EM_ComputePrProcessGivenAct(act, beta)

% pr_process_given_act = CBIG_EM_ComputePrProcessGivenAct(act, beta)
%
% act  = # studies x 18715
% beta = # processes (K) x 18715 
%
% Written by B.T.Thomas Yeo and CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md

num_processes = size(beta, 1);
curr_pr_process_given_act = zeros(num_processes, 1)+1/num_processes; % K x 1

act = sum(act, 1);

for i = 1:1000
    
    q = bsxfun(@times, beta, curr_pr_process_given_act);
    q = bsxfun(@times, q, 1./sum(q, 1));
    
    
    pr_process_given_act = sum(bsxfun(@times, q, act), 2);
    pr_process_given_act = pr_process_given_act./sum(pr_process_given_act);
    
    change = max(abs(curr_pr_process_given_act - pr_process_given_act));
    curr_pr_process_given_act = pr_process_given_act;
    
    if(change < 1e-4)
        disp(['iter ' num2str(i) ': ' num2str(change)]);
        break;
    end
end
