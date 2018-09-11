function j = CBIG_MFMem_rfMRI_diff_PC1(func,para,num)

%-----------------------------------------------------------------------------
% Function for Complex-step first order derivative approximation
%
%                       Im[f(x0+ih)]
%     diff(f(x0) =    --------------
%                             h
%
% Input:
%     - func:     model equation for derivation, f   
%     - para:     model parameter      
%     - num:      chosen parameter index, ith-parameter from para vector, 
%
% Output:
%     - j:        derivative f at chosen parameter 
%
% Example:
%    j = CBIG_MFMem_rfMRI_diff_PC1(func,para,houtput,num)
%    suppose: 
%          - model function Y with parameter vector X: Y = model(X)
%          - parameter vector X = [x1,x2,x3]
%          - funcY:  funcY = @(X) model(X); 
%    
%    then the derivative model Y at x1 (dY/dx1) is computed by:       
%          j_x1 = CBIG_mfm_rfMRI_diff_PC1(funcY,X,1)
%
% Reference: 
%    (Martins 2003), The complex-step derivative approximation.
%
% Written by Peng Wang and CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md
%----------------------------------------------------------------------------

    h1i = 1e-20*1i;
    h1 = 1e-20;

    para_1 = para;
    para_1(num) = para_1(num)+h1i;
    
    houtput_new = imag(func(para_1));
    
 
    j = houtput_new/(h1);
   
end











