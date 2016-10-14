function ccfs = clustering_coefficients(A,varargin)
% CLUSTERING_COEFFICIENTS Compute the clustering coefficients for vertices.
%
% ccfs = clustering_coefficients(A) returns the clustering coefficients for
% all vertices in A.  The clustering coefficient is the ratio of the number
% of edges between a vertex's neighbors to the total possible number of 
% edges between the vertex's neighbors. 
%
% This method works on directed or undirected graphs.
% The runtime is O(nd^2) where d is the maximum vertex degree.
%
% ... = clustering_coefficients(A,...) takes a set of
% key-value pairs or an options structure.  See set_matlab_bgl_options
% for the standard options. 
%   options.undirected: enable optimizations for undirected graphs [{0} | 1]
%   options.unweighted: an optional switch to perform the weighted 
%       computation [{0} | 1], see Note
%   options.edge_weight: a double array over the edges with an edge
%       weight for each node, see EDGE_INDEX and EXAMPLES/REWEIGHTED_GRAPHS
%       for information on how to use this option correctly
%       [{'matrix'} | length(nnz(A)) double vector]
%
% Note: Prior to version 3.0, this function did not depend on the values of
% the matrix A.  In version 3.0, the default computation changed to
% support weighted components.  To maintain backwards compatibility, the
% function will automatically set options.unweighted if the input matrix
% is not a double type.  That is, unless 'unweighted' is explictly set to
% 0, in which case we throw an error 
%
% Example:
%    load('graphs/clique-10.mat');
%    clustering_coefficients(A)

% David Gleich
% Copyright, Stanford University, 2006-2008

%%
%  2006-04-19: Initial version
%  2006-05-31: Added full2sparse check
%  2007-07-11: Added directed and weighted options
%  2007-07-12: Added non-negative edge check
%  2008-10-07: Changed options parsing
%%

[trans check full2sparse] = get_matlab_bgl_options(varargin{:});
if full2sparse && ~issparse(A), A = sparse(A); end

options = struct('edge_weight', 'matrix', 'undirected', 0, 'unweighted', 0);
options = merge_options(options,varargin{:});

% edge_weights is an indicator that is 1 if we are using edge_weights
% passed on the command line or 0 if we are using the matrix.
edge_weights = 0;
edge_weight_opt = 'matrix';

if strcmp(options.edge_weight, 'matrix')
    % do nothing if we are using the matrix weights
else
    edge_weights = 1;
    edge_weight_opt = options.edge_weight;
end

if check
    % possibly check the symmetry
    check_matlab_bgl(A,struct('sym',options.undirected));
    
    if options.unweighted ~= 1 && edge_weights~=1
        try
            check_matlab_bgl(A,struct('nodefault',1,'values',1,'noneg',1));
        catch
            % we got a value error, check if they specified unweighted
            if ~isempty(varargin) && ...
                isfield(varargin{1},'unweighted') && varargin{1}.('unweighted')~=1
                rethrow(lasterror);
            else
                % use the backward compatibility option
                options.unweighted = 1;
            end
        end
    elseif options.unweighted ~= 1 && edge_weights && any(edge_weight_opt < 0)
        error('matlab_bgl:invalidParameter', ...
                'the edge_weight array must be non-negative');
    end
end

% for the directed case, we don't transpose because the result is symmetric
% about the tranpose and we do something similar internally
if options.undirected && trans, A = A'; end

weight_arg = options.unweighted;
if ~weight_arg
    weight_arg = edge_weight_opt;
else
    weight_arg = 0;
end

ccfs=clustering_coefficients_mex(A,options.undirected,weight_arg);


