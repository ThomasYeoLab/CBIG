function traceAB = CBIG_MFMem_rfMRI_Trace_AXB(A,B)

%------------------------------------------------------------------------
% traceAB = CBIG_MFMem_rfMRI_Trace_AXB(A,B)
%
% This is a function to caculate y = trace(AB)
%
% Input:
%      - A: matrix A
%      - B: matrix B
%
% Output:
%      - traceAB: trace of AB
%
% Example:
%     y = CBIG_mfm_rfMRI_Trace_AXB(A,B)
%     suppose:  A = [1 2; 3 4]; B = [4 5; 6 7];
%     y = 69
%
% Written by Peng Wang and CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md
%------------------------------------------------------------------------

Num = size(A,1);
tr = zeros(1,Num);

for i = 1:Num
    temp = bsxfun(@times,A(i,:),B(:,i)');
    tr(i) = sum(temp);
end

traceAB = sum(tr);

end
