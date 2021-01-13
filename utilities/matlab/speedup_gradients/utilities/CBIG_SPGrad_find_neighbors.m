function K_neighbors = CBIG_SPGrad_find_neighbors(neighbors,K)

% K_neighbors = CBIG_SPGrad_find_neighbors(neighbors,K)
%
% This script find the K nearest neighbors for a given neighborhood structure
%
% Input:
%     - neighbors:
%       the #vertice x 6 neighborhood of a mesh structure 
%
%     - K:
%       the number of nearest neighbors.
%
% Output:
%     - K_neighbors:
%       the K nearest neighbors for a mesh.
%
% Written by Ru(by) Kong and CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md

neighbors(isnan(neighbors)) = size(neighbors,1)+1;
neighbors = [neighbors; ones(1,size(neighbors,2))*(size(neighbors,1)+1)];
K_neighbors = neighbors(:,2:end);
if(K > 1)
    curr_K = 1; 
    while(curr_K ~= K)
        if(curr_K == 1)
            curr_K_neighbors = neighbors(:,2:end);
            curr_K_neighbors = sort(curr_K_neighbors,2);
        end
        num_neigh = size(curr_K_neighbors,2);
        for nn = 1:num_neigh
        
            curr_neigh = neighbors(curr_K_neighbors(:,nn),2:end);
            
            for i = 1:size(K_neighbors,2)
                tmp = bsxfun(@minus, curr_neigh, K_neighbors(:,i));
                curr_neigh(tmp == 0) = size(neighbors,1);
            end
            %rm repeated neighbors
            curr_neigh = sort(curr_neigh,2);
            repeat_neigh = min(curr_neigh) == size(neighbors,1);
            curr_neigh(:,repeat_neigh) = [];
            if(nn == 1)
                tmp_curr_K_neighbors = curr_neigh;
            else
                
                for i = 1:size(tmp_curr_K_neighbors,2)
                    tmp = bsxfun(@minus, curr_neigh, tmp_curr_K_neighbors(:,i));
                    curr_neigh(tmp == 0) = size(neighbors,1);
                end
                %rm repeated neighbors
                curr_neigh = sort(curr_neigh,2);
                repeat_neigh = min(curr_neigh) == size(neighbors,1);
                curr_neigh(:,repeat_neigh) = [];
                tmp_curr_K_neighbors = [tmp_curr_K_neighbors curr_neigh];
            end
            K_neighbors = sort([K_neighbors curr_neigh],2);
            %find new neighbors
            repeat_neigh = min(K_neighbors) == size(neighbors,1);
            K_neighbors(:,repeat_neigh) = [];
           
            
            
        end
        curr_K_neighbors = tmp_curr_K_neighbors; 
        
       
        curr_K = curr_K + 1;
    end

        
end
%exclude the node itself and append it in the first column
tmp = bsxfun(@minus, K_neighbors, neighbors(:,1));
K_neighbors(tmp == 0) = size(neighbors,1);
K_neighbors = sort(K_neighbors,2);
repeat_neigh = min(K_neighbors) == size(neighbors,1);
K_neighbors(:,repeat_neigh) = [];
K_neighbors = [neighbors(:,1) K_neighbors];

K_neighbors(K_neighbors == size(neighbors,1)) = NaN;
K_neighbors(end,:) = []; 
   
end
