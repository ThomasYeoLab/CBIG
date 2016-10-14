function CBIG_WriteParcellationToAnnotation(labels, filename, color_mat, struct_names)

% CBIG_WriteParcellationToAnnotation(labels, filename, color_mat, struct_names)
%
% Write colors of parcellation into annotation file.
%
% Input arguments:
%     - labels    : labels of parcellation
%     - filename  : output file name
%     - color_mat : a N*3 color matrix, e.g. FreeSurfer's color matrix
%
% Written by CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md



if(nargin < 4)
    ct = CBIG_CreateCTfromColorMat(color_mat);
else
    ct = CBIG_CreateCTfromColorMat(color_mat, struct_names);
end
labels = ct.table(labels + 1, 5);
Write_Brain_Annotation(filename, 0:(length(labels)-1), labels, ct);










