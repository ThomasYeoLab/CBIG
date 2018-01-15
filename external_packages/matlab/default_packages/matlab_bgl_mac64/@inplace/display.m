function display(ipa)
% INPLACE/DISPLAY Display an inplace array from the command line.
%
% Example:
%    ipa = inplace(ones(5));
%    ipa

disp([inputname(1), ' = ']);
disp(ipa.get_a());