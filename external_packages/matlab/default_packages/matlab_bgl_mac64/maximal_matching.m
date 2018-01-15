function m=maximal_matching(A,varargin)
% MAXIMAL_MATCHING Compute a maximal matching
%
% A maximal matching is a subset of edges where each vertex is incident on
% only one edge and the size of the matching cannot be increased by adding
% an edge left in the graph.
%
% m=maximal_matching(A) returns a maximal matching where m(v) = u and m(u)
% = v if edge (u,v) is in the matching and m(v) = 0 if v is not matched.
%
% There are two algorithms for computing a maximal matching: a greedy
% algorithm that simply adds an edge (u,v) to the matching unless m(u) or
% m(v) != 0; and an extra greedy algorithm that first sorts the edges in
% the graph by their degree.  See the options.algname option.
%
% This method works on undirected graphs (symmetric matrices) and ignores
% edge weights.  
% The runtime is O(M log N). 
%
% To change the algorithm, use the optional algname argument.
%   options.algname: the greedy algorithm to use ['greedy' | {'extra_greedy'}]
%
% See the MATCHING function for additional calling information.  This function 
% just calls matching(...,struct('initial_match',options.algname,...
% 'augmenting_path','none','verify',0)) and does _not_ verify the output,
% so there is no output v.
%
% See also MATCHING, TEST_MATCHING 
%
% Example:
%   load('graphs/matching_example.mat');
%   m=maximal_matching(A)
%   sum(m>0)/2                % maximal matching cardinality, should be < 8
%   mmax=matching(A);
%   sum(mmax>0)/2             % maximum matching cardinality, should be 8 

% David Gleich
% Copyright, Stanford University, 2007-2008

%% History
%  2007-07-09: Initial coding
%  2008-10-07: Changed options parsing
%%

options = merge_options(struct(),varargin{:});
options.augmenting_path = 'none';
options.verify = 0;
if isfield(options,'algname'), options.initial_match = options.algname; 
else options.initial_match='extra_greedy'; end

m = matching(A,options);
