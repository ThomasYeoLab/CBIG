function [trans check full2sparse] = get_matlab_bgl_options(varargin)
% 
% Internal private function.
%
% Example:
%    Don't use this function!
%

%% History
%  2008-09-26: Changed to use merge_options instead
%%

doptions = set_matlab_bgl_default();
if nargin>0
    options = merge_options(doptions,varargin{:});
else
    options = doptions;
end

trans = ~options.istrans;
check = ~options.nocheck;
full2sparse = options.full2sparse;

