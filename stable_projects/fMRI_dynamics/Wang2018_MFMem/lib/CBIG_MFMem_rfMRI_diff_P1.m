 function j = CBIG_MFMem_rfMRI_diff_P1(func,para,houtput,num)

%------------------------------------------------------------------------
% j = CBIG_mfm_rfMRI_diff_P1(func,para,houtput,num)
%
% Function for Newton-forwards first order derivative approximation
%
%                       f(x0+h)-f(x0)
%     diff(f(x0) =      --------------
%                             h
%
% Input:
%     - func:     model equation for derivation, f   
%     - para:     model parameter      
%     - houtput:  model output at "para" 
%     - num:      chosen parameter index, ith-parameter from para vector, 
%
% Output:
%     - j:        derivative f at chosen parameter, df(x0)/x0
%
% Example:
%    j = CBIG_MFMem_rfMRI_diff_P1(func,para,houtput,num)
%    suppose: 
%          - model function Y with parameter vector X: Y = model(X)
%          - parameter vector X = [x1,x2,x3]
%          - funcY:  funcY = @(X) model(X); 
%    
%    then the derivative model Y at x1 (dY/dx1) is computed by:       
%          j_x1 = CBIG_mfm_rfMRI_diff_P1(funcY,X,Y,1)
%
% Written by Peng Wang and CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md
%----------------------------------------------------------------------------

    Mesp = 2.2*1e-16;

    para_1 = para;
    if para(num) == 0
        h = sqrt(Mesp);
    else    
        h = sqrt(Mesp)*para(num);
    end
    para_1(num) = para_1(num)+h;
    
    houtput_new = func(para_1);
    
    j = (houtput_new - houtput)/(para_1(num)-para(num));
    
end











