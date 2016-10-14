function X = gursoy_atun_layout(A,varargin)
% GURSOY_ATUN_LAYOUT Layout a graph by uniformly distributing vertices
% 
% Instead of trying to optimize an objective function of forces or springs
% in the layout, the Gursoy-Atun layout distributes vertices uniformly
% within a topology subject to keeping vertices nearby their neighbors.
%
% X = gursoy_atun_layout(A) compute a Gursoy-Atun layout.
%
% See
% http://www.boost.org/doc/libs/1_36_0/libs/graph/doc/gursoy_atun_layout.html
% for a description of all the parameters.
%
% ... = gursoy_atun_layout(A,...) takes a set of
% key-value pairs or an options structure.  See set_matlab_bgl_options
% for the standard options.  
%   options.topology: the topology spaced used for the layout
%     [{'square'} | 'heart' | 'sphere' | 'circle' | 'ballN'* | 'cubeN'*]
%     * for the ball and cube topolgy, N can be replaced by any number, so 
%       ball3 is the same as the sphere, cube2 is the same as the square,
%       but cube3 is the true cube and cube4 is a hypercube.  At the moment, 
%       only dimensions between 1 and 10 are implemented. 
%   options.iterations: the number of iterations [num_vertices(G) | integer]
%   options.diameter_range: The inital and final diameters when updating
%     [ {[sqrt(num_vertices(G)),1.0]} | [diameter_initial,diameter_final] ]
%     where all values are double.
%   options.learning_constant_range: The inital and final learning constants 
%     [ {[0.8,0.2]} | [learning_constant_initial,learning_constant_final] ]
%     where all values are double.
%
% Example:
%   G1 = cycle_graph(5000,struct('directed',0));
%   X1 = gursoy_atun_layout(G1,'topology','heart');
%   G2 = grid_graph(50,50);
%   X2 = gursoy_atun_layout(G2,'topology','square');
%   G3 = grid_graph(50,50);
%   X3 = gursoy_atun_layout(G3,'topology','circle');
%   subplot(1,3,1); gplot(G1,X1,'k'); subplot(1,3,2); gplot(G2,X2,'k');
%   subplot(1,3,3); gplot(G3,X3,'k');
%   
% See also KAMADA_KAWAI_SPRING_LAYOUT, 
% FRUCHTERMAN_REINGOLD_FORCE_DIRECTED_LAYOUT, LAYOUT

% David F. Gleich
% Copyright, Stanford University, 2008

%% History
%  2008-09-28: Initial coding
%%

[trans check full2sparse] = get_matlab_bgl_options(varargin{:});
if full2sparse && ~issparse(A), A = sparse(A); end

n = num_vertices(A);
options = struct('topology','square','iterations',n,...
    'diameter_range',[sqrt(n) 1.0],'learning_constant_range',[0.8 0.2],...
    'progressive',0,'edge_weight','matrix');
options = merge_options(options,varargin{:});
% edge_weights is an indicator that is 1 if we are using edge_weights
% passed on the command line or 0 if we are using the matrix, or -1 to use
% nothing
edge_weights = 0;
edge_weight_opt = 'matrix';

if strcmp(options.edge_weight, 'matrix'), % do nothing to use the matrix weights
elseif strcmp(options.edge_weight,'none'), edge_weights = -1; % don't use weights
else edge_weights = 1; edge_weight_opt = options.edge_weight;
end

if check
    % make sure the matrix is symmetric
    if ~edge_weights==1
        check_matlab_bgl(A,struct('values',edge_weights==0));
    else
        if nnz(A) ~= length(edge_weight_opt)
            error('matlab_bgl:invalidParameter', ...
             'the vector of edge weights must have length nnz(A)'); 
        end
    end
end

progressive_opt = [];
if ~isscalar(options.progressive), progressive_opt = options.progressive; end

X= gursoy_atun_mex(...
    A, options.topology, options.iterations, ...
    options.diameter_range(1), options.diameter_range(2),...
    options.learning_constant_range(1), options.learning_constant_range(2),...
    progressive_opt, edge_weights, edge_weight_opt);
