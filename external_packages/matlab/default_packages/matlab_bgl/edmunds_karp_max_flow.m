function [varargout]=edmunds_karp_max_flow(A,u,v,varargin)
% EDMUNDS_KARP_MAX_FLOW Edmunds-Karp max flow algorithm
%
% The Edmunds-Karp algorithm implements the augmenting path idea
% to compute a maximum flow or minimum cut in a network.  
%
% See the max_flow function for calling information and return parameters.
% This function just calls max_flow(...,struct('algname','edmunds_karp'));
% The complexity is O(VE^2) or O(VEU) where U is the maximum capacity.
%
% Example:
%   load('graphs/max_flow_example.mat');
%   edmunds_karp_max_flow(A,1,8)

% David Gleich
% Copyright, Stanford University, 2007-2008

%% History
%  2007-07-07: Initial version
%  2008-10-07: Changed options parsing
%%

algname = 'edmunds_karp';

if ~isempty(varargin), 
    options = merge_options(struct(),varargin{:}); 
    options.algname= algname;
else options = struct('algname',algname); 
end

varargout = cell(1,max(nargout,1));

[varargout{:}] = max_flow(A,u,v,options);
