% Caret-toolbox for Matlab 
%   Functionality is provided for a file-interface with caret
%   Calculation of surface based statistics from metric files 
%   and Gaussian Field theory
% --------------------------------------------------------------
% Caret is a freely available surface-based mapping program developed 
% By David van Essen's lab at Washington University
% http://brainmap.wustl.edu/vanessen.html
% --------------------------------------------------------------
% v.1.0 author: Joern Diedrichsen
% v.1.1 5/01/06: added metriccalc and metric permute, added spmj_ECdensity
%       email: jdiedric@jhu.edu
% 
% 1. File-IO 
%  caret_load       - can load a metric,coord,topo,paint,surface_shape into a structure
%  caret_save       - saves a caret-structure into a file
%  caret_savemetric - saves a Metric-structure as file 
%  caret_savepaint  - saves a Paint-structure as file 
%  caret_struct     - generates a new structure to hold a file
% 2. Surface operations 
%  caret_getsurface - loads a .coord and .topo file into a surface structure 
%  caret_calcarea   - calulates the Euler Characteristic, circumference and Area of the search region
%  caret_calcneighbor - makes a list of neighborhood relationships (needed for caret_cluster)
%  caret_surfacesubset - makes a subset of surface for small-volume correction
%  caret_avrgPALSsurface - Generates average surces from PALS atlas
% 3. Metric file operations 
%  caret_regionstat - Extract the value from metric files based on ROIs specified in a paint file 
%  caret_mask       - Masks out areas from a metric file, by applying logical mask, specified in a different metric file
%  caret_threhold   - thresholds a metric file on a height- and clustersize-threshold 
%  caret_metriccalc - Calculator for metric files, works like spm_imcalc
%  caret_metricpermute - Makes N files with P columns to P files with N cols
% 4. Statistical Maps 
%  caret_getcSPM    - makes a statistical map from metric (using a GLM) 
%  caret_savecSPM   - saves a statistical map as a metric file
%  caret_estsmoothness - estimates the smoothness of a statistical map
%  caret_clusters   - finds clusters on statistical map
%  caret_P          - calculates corrected, and cluster-wise p-values
%  caret_list       - lists maxima, clusters, p-values
