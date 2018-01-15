function X = random_graph_layout(G,box)
% RANDOM_GRAPH_LAYOUT Layout the vertices of a graph randomly
% 
% X = random_graph_layout(G) generate a random layout of the graph in the
% unit box [0,0],[1,1] with vertex points generated at random.
%
% X = random_graph_layout(G,box) gives the coordinates of the layout box
% [minX,minY,maxX,maxY].  If this parameter is empty, then the default 
% [0,0] to [1,1] box is used.  Alternatively, if this parameter is an
% integral type, then the vertices of the graph are chosen uniformly on the
% integer grid from [minX,minY] to [maxX,maxY].
%
% Example:
%   G = cycle_graph(1500);
%   X = random_graph_layout(G);
%   gplot(G,X); hold on; plot(X(:,1),X(:,2),'r.'); hold off;
%   % Layout on the grid
%   X = random_graph_layout(G,int32([0 0 5 5])); % random grid layout
%   gplot(G,X); grid on; hold on; plot(X(:,1),X(:,2),'r.'); hold off;

% David F. Gleich
% Copyright, Stanford University, 2008

%% History
%  2008-09-25: Initial coding
%%

if ~exist('box','var') || isempty(box), box = [0 0 1 1]; end;
n = num_vertices(G);
X = rand(n,2);
xrange = double(box(3)-box(1));
yrange = double(box(4)-box(2));
X(:,1) = xrange*X(:,1) + double(box(1));
X(:,2) = yrange*X(:,2) + double(box(2)); 
if isinteger(box), X = floor(X); end
    
