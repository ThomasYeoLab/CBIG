function CBIG_SPGrad_save_dump(filename)

% This function save a dump file <filename>.dump
% Written by Ru(by) Kong and CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md

fclose(fopen([filename '.dump'],'w'));

end
