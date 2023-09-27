function [adj, sz_links] = CBIG_components_subgraphs(adj)

% Script to get components sub-graph. 
% 
%     function [adj, sz_links] = CBIG_components_subgraphs(adj)
%     Input:
%        adj:adjacency matrix with values going from 1 to number of components with
%            more than 1 node. The script assues that the adj matrix is a symmetric 
%            matrix.
%    Output:
%        adj:adjacency matrix 
%        sz_links: the size of components edge-links, that is, the number of edges 
%            it comprises.
%
% Written by CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md


adj = double(adj);
diff = adj - adj';
if(max(abs(diff(:))) ~= 0)
   error('adj not symmetric'); 
end

[a, sz] = components(adj);

%Convert size from number of nodes to number of edges
%Only consider components comprising more than one nodes (equivalent to at
%least one edge)
ind_sz = find(sz > 1);
sz_links = zeros(1, length(ind_sz));
for i=1:length(ind_sz)
    nodes = find(a == ind_sz(i)); 
    sz_links(i)=sum(sum(adj(nodes,nodes)))/2;
    adj(nodes, nodes) = adj(nodes, nodes) * (i+1);
end

adj(adj ~= 0) = adj(adj ~= 0) - 1;
adj = double(adj);
