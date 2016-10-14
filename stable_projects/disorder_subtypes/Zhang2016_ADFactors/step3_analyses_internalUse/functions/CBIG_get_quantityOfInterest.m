function [q, q_str] = CBIG_get_quantityOfInterest(rid, Q, predementedOnly)

% [q, q_str] = CBIG_get_quantityOfInterest(rid, Q, predementedOnly)
% Quantity of interest
%
% Written by CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md

if Q == 1 % MMSE
    q = CBIG_get_mmse(rid, predementedOnly);
    q_str = 'mmse';
elseif Q == 2 % MEM
    q = CBIG_get_mem_ef(rid, predementedOnly);
    q = q(:, [1 2 3 4]);
    q_str = 'mem';
elseif Q == 3 % EF
    q = CBIG_get_mem_ef(rid, predementedOnly);
    q = q(:, [1 2 3 5]);
    q_str = 'ef';
else
    error('Wrong selection!');
end

disp(q_str);
