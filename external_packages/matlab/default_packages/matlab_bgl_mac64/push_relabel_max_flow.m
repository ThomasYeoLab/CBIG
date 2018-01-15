function [varargout]=push_relabel_max_flow(A,u,v,varargin)
% PUSH_RELABEL_MAX_FLOW Goldberg's push-relabel max flow algorithm
%
% This function calls Boost's implementation of Goldberg's push relabel
% algorithm to compute a maximum flow in a graph.
%
% See the max_flow function for calling information and return parameters.
% This function just calls max_flow(...,struct('algname','push_relabel'));
%
% Example:
%   load('graphs/max_flow_example.mat');
%   push_relabel_max_flow(A,1,8)

% David Gleich
% Copyright, Stanford University, 2007-2008

%% History
%  2007-07-07: Initial version
%  2008-10-07: Changed options parsing
%%

algname = 'push_relabel';
if ~isempty(varargin), 
    options = merge_options(struct(),varargin{:}); 
    options.algname= algname;
else options = struct('algname',algname); 
end

varargout = cell(1,max(nargout,1));
[varargout{:}] = max_flow(A,u,v,options);
