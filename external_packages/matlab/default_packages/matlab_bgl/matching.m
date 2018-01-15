function [m,max_card_matching]=matching(A,varargin)
% MATCHING Compute a maximum cardinality matching
%
% m=matching(A) returns a matching between the vertices of an undirected
% graph.  A matching is a subset of edges where each vertex in incident on
% only one edge.  A maximum cardinality matching is the largest possible
% set of edges with this property.  The output from matching is a vector
% where m(v) = u if vertex v is matched to vertex u or m(v) = 0 if vertex v
% is not matched.  Consequently, sum(m>0)/2 is the cardinality of the
% matching.  Symmetry implies if m(v) = u, then m(u) = v.  
%
% [m,v]=matching(A) returns an optional parameter v, where v=1 if the
% matching was verified as maximum cardinality and v=0 otherwise.  
%
% [M,v]=matching(A,struct('matrix_output',1)) encodes the matching as a 
% sparse matrix so that nnz(M)/2 is the cardinality and A.*M is the set of 
% edges picked in the matching.
%
% The matching function will throw a warning if the edmonds algorithm
% failed (that is, v=0, but v should be 1) and if verification is enabled.
% (Please don't disable verification, it really isn't that expensive.)
%
% In the Boost graph library, a matching algorithm has three components:
% 1.  an initial matching
% 2.  augmenting path extension
% 3.  verification
% which allow us to build our own matching algorithms from these three
% components.  All possibilities for these options are implemented in the
% optional parameters to the matching function in MatlabBGL. 
%
% For sparse graphs, the call matching(A,struct('initial_match','none'))
% may work better.
%
% The greedy and extra_greedy initial match options create maximal
% matchings, which cannot be increased in size by adding edges and are at
% least half the cardinality of the maximum matching.  The call 
% [m,v] = matching(A,struct('augmenting_path','none')) computes a maximal
% matching and sets v = 1 if the matching is maximum too.  
%
% This method works on undirected graphs (symmetric matrices) and ignores
% edge weights.  
% The runtime is O(mn*alpha(m,n)).  The alpha function is the inverse
% Ackermann function and is <= 4 for all valid Matlab inputs.
%
% ...=matching(A,...) takes a set of
% key-value pairs or an options structure.  See set_matlab_bgl_options
% for the standard options. 
%   options.initial_match: choice of initial matching 
%     ['none' | 'greedy' | {'extra_greedy'}]
%   options.augmenting_path: [{'edmonds'} | 'none']
%   options.verify: verify output is maximum cardinality matching [0 | {1}]
%   options.matrix_output: return the matching as a sparse matrix [{0} | 1]
%
% Example:
%   load graphs/matching_example.mat
%   m = matching(A) 
%   sum(m>0)/2                % matching cardinality, should be 8
%   [m,v] = matching(A,struct('augmenting_path','none'))% not maximum matching


% David Gleich
% Copyright, Stanford University, 2007-2008

%% History
%  2007-07-08: Initial version
%  2008-10-07: Changed options parsing
%%

[trans check full2sparse] = get_matlab_bgl_options(varargin{:});
if full2sparse && ~issparse(A), A = sparse(A); end
if trans, end % input must be symmetric, so no need to transpose

options = struct('initial_match', 'extra_greedy', ...
    'augmenting_path', 'edmonds', ...
    'verify', 1, ...
    'matrix_output', 0);
options = merge_options(options, varargin{:});    

if check, check_matlab_bgl(A,struct('sym',1)); end % make sure the matrix is symmetric

[m,max_card_matching] = matching_mex(A, ...
    options.verify, ...
    lower(options.initial_match), ...
    lower(options.augmenting_path));

if options.verify && ~max_card_matching && ...
        strcmpi(options.augmenting_path,'edmonds') && ...
        nargout < 2
    % display a warning if they requested verification, it wasn't a maximal
    % matching, and they tried a non-trivial algorithm
    warning('matlab_bgl:notAMaximumMatching', ...
        'edmonds algorithm did not return a maximum cardinality matching');
end

if options.matrix_output
    mask = m>0;
    i = (1:size(A,1)); i = i(mask);
    m = m(mask);
    m = sparse(i,m,true(1),size(A,1),size(A,1));
end
