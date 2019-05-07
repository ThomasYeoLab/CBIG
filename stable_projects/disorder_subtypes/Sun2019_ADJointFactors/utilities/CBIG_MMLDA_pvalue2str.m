function s = CBIG_MMLDA_pvalue2str(p)
% s = CBIG_MMLDA_pvalue2str(p)
%
% Convert p value to string 
%
% Input:
%   - p     : double number
%
% Output:
%   - s     : string
%
% Example:
%   s = CBIG_MMLDA_pvalue2str(0.0001);
%
% Written by Nanbo Sun and CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md

if p < 0.01
    s = num2str(p, '%.0e');
    s = regexprep(s, 'e-0(\d)', 'e-$1');
else 
    s = num2str(p, '%.2f');
end
