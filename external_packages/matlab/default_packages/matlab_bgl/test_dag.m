function dag = test_dag(A,varargin)
% TEST_DAG Tests if a graph is directed and acyclic
%
% dag = test_dag(A) returns 1 if A is the adjacency matrix for a dag and
% returns 0 otherwise.  This function just calls topological order and
% checks the output to determine if the graph is a dag.  It is provided
% merely for expository purposes.
%
% Example:
%   n = 10; A = sparse(1:n-1, 2:n, 1, n, n); % construct a simple dag
%   test_dag(A)
%   A(10,1) = 1; % complete the cycle
%   test_dag(A)

% David Gleich
% Copyright, Stanford University, 2007-2008

%% History
%%

dag = 1;
if isempty(topological_order(A,varargin{:}))
    dag = 0;
end
