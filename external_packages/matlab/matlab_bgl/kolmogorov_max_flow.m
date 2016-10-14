function [varargout]=kolmogorov_max_flow(A,u,v,varargin)
% KOLMOGOROV_MAX_FLOW Kolmogorov's max flow algorithm
%
% Kolmogorov's algorithm implements a variation on the augmenting path idea
% to compute a maximum flow or minimum cut in a network.  
%
% See the max_flow function for calling information and return parameters.
% This function just calls max_flow(...,struct('algname','kolmogorov'));
%
% Example:
%   load('graphs/max_flow_example.mat');
%   kolmogorov_max_flow(A,1,8)

% David Gleich
% Copyright, Stanford University, 2007-2008

%% History
%  2007-07-07: Initial version
%  2008-10-07: Changed options parsing
%%

algname = 'kolmogorov';
if ~isempty(varargin), 
    options = merge_options(struct(),varargin{:}); 
    options.algname= algname;
else options = struct('algname',algname); 
end

varargout = cell(1,max(nargout,1));

[varargout{:}] = max_flow(A,u,v,options);
