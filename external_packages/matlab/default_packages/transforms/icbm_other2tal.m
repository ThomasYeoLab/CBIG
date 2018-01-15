function outpoints = icbm_other2tal(inpoints)
%
% This function converts coordinates from MNI space (normalized 
% using templates other than those contained in SPM and FSL) to 
% Talairach space using the icbm2tal transform developed and 
% validated by Jack Lancaster at the Research Imaging Center in 
% San Antonio, Texas.
%
% http://www3.interscience.wiley.com/cgi-bin/abstract/114104479/ABSTRACT
% 
% FORMAT outpoints = icbm_other2tal(inpoints)
% Where inpoints is N by 3 or 3 by N matrix of coordinates
% (N being the number of points)
%
% ric.uthscsa.edu 3/14/07

% find which dimensions are of size 3
dimdim = find(size(inpoints) == 3);
if isempty(dimdim)
  error('input must be a N by 3 or 3 by N matrix')
end

% 3x3 matrices are ambiguous
% default to coordinates within a row
if dimdim == [1 2]
  disp('input is an ambiguous 3 by 3 matrix')
  disp('assuming coordinates are row vectors')
  dimdim = 2;
end

% transpose if necessary
if dimdim == 2
  inpoints = inpoints';
end

% Transformation matrices, different for each software package
icbm_other = [0.9357 0.0029 -0.0072 -1.0423
			 -0.0065 0.9396 -0.0726 -1.3940
			  0.0103 0.0752  0.8967  3.6475
			  0.0000 0.0000  0.0000  1.0000];

% apply the transformation matrix
inpoints = [inpoints; ones(1, size(inpoints, 2))];
inpoints = icbm_other * inpoints;

% format the outpoints, transpose if necessary
outpoints = inpoints(1:3, :);
if dimdim == 2
  outpoints = outpoints';
end
