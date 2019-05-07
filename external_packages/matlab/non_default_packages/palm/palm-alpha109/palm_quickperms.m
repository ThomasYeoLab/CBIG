function varargout = palm_quickperms(varargin)
% Create a set of permutations that is left in the Octave/Matlab workspace
% for later use, or to be exported to other programs.
%
% [Pset,VG] = palm_quickperms(M,EB,P,EE,ISE,CMCx,CMCp)
%
% Inputs (to skip an argument, use an empty array, []):
% M       : Design matrix. It can be the full, unpartitioned design, or
%           if there are nuisance, simply the part that contains the EVs
%           of interest. This distinction is only relevant if there are
%           discrete EVs of interest and nuisance variables that are
%           continuous. You may consider a partitioning as in the
%           function palm_partition.m.
%           If you have no idea what to use, it is in general, it is in
%           general safe to use as M simply a vector (1:N)'.
%           You can also simply leave it empty ([]) if EB is supplied, and
%           by default it will be (1:N)'. If an EB isn't supplied, you can
%           simply use N and by default it will be (1:N)'.
% EB      : Exchangeability blocks (can be multi-level). For freely
%           exchangeable data, use ones(N,1). You can also leave it
%           empty ([]) if a valid, non-empty M was supplied.
% P       : Desired number of permutations. The actual number may be
%           smaller if N is too small. Use 0 for exhaustive.
%           Default is 10000.
% EE      : True/False indicating whether to assume exchangeable errors,
%           which allow permutations.
% ISE     : True/False indicating whether to assume independent and
%           symmetric errors, which allow sign-flippings.
% CMCx    : True/False indicating whether repeated rows in the design
%           should be be ignored. Default is false.
% CMCp    : True/False indicating whether repeated permutations should
%           be ignored. Default is false.
%
% Outputs:
% Pset    : Permutation set. It contains one permutation per column.
% VG      : Variance groups (VGs), to be used with a statistic that is
%           robust to heteroscedasticity.
%
% _____________________________________
% Anderson M. Winkler
% FMRIB / University of Oxford
% Sep/2015
% http://brainder.org

% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
% PALM -- Permutation Analysis of Linear Models
% Copyright (C) 2015 Anderson M. Winkler
%
% This program is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% any later version.
%
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with this program.  If not, see <http://www.gnu.org/licenses/>.
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

% Take the inputs:
v.M       = [];
v.EB      = [];
v.P       = 10000;
v.EE      = true;
v.ISE     = false;
v.CMCx    = false;
v.CMCp    = false;
fields = fieldnames(v);
for a = 1:nargin,
    if ~ isempty(varargin{a}),
        v.(fields{a}) = varargin{a};
    end
end
if ~ v.EE && ~ v.ISE,
    error('At least one of EE or ISE must be given as "true".');
end
if v.P < 0,
    error('P must not be negative');
end
if isempty(v.M ), v.M  = 0; end
if isempty(v.EB), v.EB = 0; end 

% Number of subjects
N = max(size(v.M,1),size(v.EB,1));
if N == 1,
    N = max(v.M,v.EB);
elseif N == 0,
    error('Design matrix and exchangeability blocks cannot be both empty.');
end
if v.M  == 0, v.M  = N; end
if v.EB == 0, v.EB = N; end 
if ...
        (size(v.M,1)  > 1 && N ~= size(v.M,1))  || ...
        (isscalar(v.M)    && N ~= v.M)          || ...
        (size(v.EB,1) > 1 && N ~= size(v.EB,1)) || ...
        (isscalar(v.EB)   && N ~= v.EB),
    error('Design matrix and exchangeability blocks of incompatible sizes');
end

% Design matrix:
if isempty(v.M) || isscalar(v.M) || v.CMCx,
    v.M = (1:N)';
end

% Create the shufflings:
if ...
        isempty(v.EB)   || ...
        isscalar(v.EB)  || ...
        (isvector(v.EB) && numel(unique(v.EB)) == 1),
    
    % If no EBs were given, use simple shuffling:
    simpleshuf = true;
    Pset  = palm_shuffree(v.M,v.P,v.CMCp,v.EE,v.ISE,true);

else
    
    % Or use the multi-level blocks. Begin by reindexing the leaves:
    simpleshuf = false;
    v.EB  = palm_reindex(v.EB,'fixleaves');
    
    % Then create the permutation tree:
    Ptree = palm_tree(v.EB,v.M);
    
    % Then the set of permutations:
    Pset  = palm_shuftree(Ptree,v.P,v.CMCp,v.EE,v.ISE,true);
end
varargout{1} = Pset;

% Define the variance groups (for heteroscedasticity, if needed):
if nargout == 2,
    
    % Create the VGs:
    if simpleshuf,
        VG = ones(N,1);
    else
        VG = palm_ptree2vg(Ptree);
    end
    varargout{2} = VG;
end
