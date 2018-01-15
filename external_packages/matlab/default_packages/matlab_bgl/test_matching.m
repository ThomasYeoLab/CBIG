function max_card_match=test_matching(A,m,varargin)
% TEST_MATCHING Tests if a matching is maximum cardinality
%
% max_card_match=test_matching(A,m) returns 1 if m is a maximum cardinality
% matching on A.  
%
% This method works on undirected graphs (symmetric matrices) and ignores
% edge weights.  
% The runtime is O(mn*alpha(m,n)).  The alpha function is the inverse
% Ackermann function and is <= 4 for all valid Matlab inputs.
%
% ... = test_matching(A,m,...) takes a set of
% key-value pairs or an options structure.  See set_matlab_bgl_options
% for the standard options. 
%   There are no additional options for this function.
%
% Example:
%   load('graphs/matching_example.mat');
%   [m_not_max,v] = matching(A,struct('augmenting_path','none'));
%   test_matching(A,m_not_max)

%
% David Gleich
% Copyright, Stanford University, 2007
%

%% History
%  2007-07-12: Initial version
%%


[trans check full2sparse] = get_matlab_bgl_options(varargin{:});
if full2sparse && ~issparse(A), A = sparse(A); end
if check, check_matlab_bgl(A,struct('sym',1)); end % ensure input is symmetric
% no trans check because the input is symmetric

max_card_match = test_matching_mex(A,m);
