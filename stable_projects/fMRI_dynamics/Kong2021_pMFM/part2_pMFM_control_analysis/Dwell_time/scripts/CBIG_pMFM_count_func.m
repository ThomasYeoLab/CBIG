function num = CBIG_pMFM_count_func(array)

% This function is used to compute number of time points for high and low
% states
%
% Written by Kong Xiaolu and CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md

num = [];
count = 0;
for i = 1:length(array)
    if array(i) == 1
        count = count+1;
    elseif array(i) == 0 && count ~= 0
        num = [num;count];
        count = 0;
    else
        count = 0;
    end
        
end

if count ~= 0
    num = [num;count];
end

end