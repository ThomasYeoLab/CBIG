function adj = palm_adjacency(fac,isvtx)
% Compute an adjacency matrix for vertexwise or facewise data.
% 
% Usage:
% adj = palm_adjacency(fac,isvtx)
% 
% - fac   : Face indices (see palm_srfread.m for details).
% - isvtx : Boolean indicating if the adjacency is for vertexwise
%           or facewise data.
% - adj   : Sparse adjacency matrix.
% 
% _____________________________________
% Anderson M. Winkler
% FMRIB / Univ. of Oxford
% Jul/2015
% http://brainder.org

% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
% PALM -- Permutation Analysis of Linear Models
% Copyright (C) 2015 Anderson M. Winkler
% 
% This program is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% any later version.
% 
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
% 
% You should have received a copy of the GNU General Public License
% along with this program.  If not, see <http://www.gnu.org/licenses/>.
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

% Make sure the face indices aren't stored as int,
% otherwise sparse won't work.
fac = double(fac);

% Number of vertices and faces
nV = max(fac(:));
nF = size(fac,1);

if isvtx,
    
    % For vertexwise data:
    % Two vertices are connected if they share an edge.
    % Note difference compared to the previous palm_vtxlabel.m
    % in which two vertices are connected if they are all in
    % the same face.
    adj  = sparse( ...
        [ ...
        fac(:,1); fac(:,1);  ...
        fac(:,2); fac(:,2);  ...
        fac(:,3); fac(:,3)], ...
        [ ...
        fac(:,2); fac(:,3);  ...
        fac(:,1); fac(:,3);  ...
        fac(:,1); fac(:,2)],1,nV,nV);
    
else
    
    % For facewise data:
    % Reparameterize fac in terms of edges:
    edg  = [fac(:,[1 2]); fac(:,[2 3]); fac(:,[1 3])];
    edg  = sort(edg,2);
    [edg,~,fec] = unique(edg,'rows');
    nE   = size(edg,1);
    fec  = reshape(fec,size(fac));
    
    % Two faces are connected if they share an edge:
    fidx = repmat(1:nF,[3 1])';
    ef   = unique([fec(:),fidx(:)],'rows');
    adj  = reshape(ef(:,2),[2 nE])';
    adj  = sparse( ...
        [adj(:,1); adj(:,2)],...
        [adj(:,2); adj(:,1)],1,nF,nF);
    
%     % It's also possible to have a weaker adjacency, in which two
%     % faces are connected if they share a vertex. However, this is
%     % too slow to run in practice.
%     fidx = repmat(1:size(fac,1),[3 1])';
%     vf   = unique([fac(:),fidx(:)],'rows');
%     v    = unique(vf(:,1));
%     adj = sparse([],[],[],nF,nF,nF*13);
%     for vv = 1:numel(v),
%         vidx = vf(:,1) == vv;
%         adj(vf(vidx,2),vf(vidx,2)) = 1;
%     end
end

% Add the otherwise missing diagonal.
adj = (adj > 0) + speye(size(adj));
