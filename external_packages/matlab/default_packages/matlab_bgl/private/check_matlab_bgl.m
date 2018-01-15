function check_matlab_bgl(A,options)
% CHECK_MATLAB_BGL Checks the input A for various properties
%
% check_matlab_bgl(A,options) throws an input error if...
%   A is not square
%   if options.values and A is not double valued
%   if options.sym = 1 and A is not symmetric
%   if options.flow_graph = 1 and A is not a flow graph data structure
%   if options.nosparse = 1 do not check if A is sparse
%   if options.nodefault = 1 then do not check default cases
%   if options.nodiag = 1 throw an error if A has any non-zero diagonal values

% David Gleich
% Copyright, Stanford University, 2006-2008

%% History
%  2007-04-20: Added nodefault option
%  2007-07-22: Fixed empty array error for noneg check
%  2008-09-23: Added no diagonal check, misc formatting fixes
%% 

if ~isfield(options, 'nodefault') || options.nodefault == 0
    if size(A,1) ~= size(A,2)
        error('matlab_bgl:invalidParameter', 'the matrix A must be square.');
    end
end

if isfield(options, 'values') && options.values == 1
    if ~isa(A,'double')
        error('matlab_bgl:invalidParameter', 'the matrix A must have double values.');
    end
end

if isfield(options, 'noneg') && options.noneg == 1
    v=min(min(A));
    if ~isempty(v) && v < 0
        error('matlab_bgl:invalidParameter', 'the matrix A must have non-negative values.');
    end
end

if isfield(options, 'sym') && options.sym == 1
    if ~isequal(A,A')
        error('matlab_bgl:invalidParameter', 'the matrix A must be symmetric.');
    end
end

if isfield(options, 'nosparse') && options.nosparse == 1
else
    if ~issparse(A)
        error('matlab_bgl:invalidParameter', 'the matrix A must be sparse.  (See set_matlab_bgl_default.)');
    end
end

if isfield(options,'nodiag') && options.nodiag == 1
    if any(diag(A))
        error('matlab_bgl:invalidParameter',...
            'the matrix A must not have any diagonal values')
    end
end

