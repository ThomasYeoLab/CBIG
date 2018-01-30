function [lh_grad_matrix,rh_grad_matrix]=CBIG_gwMRF_gradient_vertices_to_matrix(lh_grad,rh_grad,fsaverage)

%% this function creates two  6xN GradientMatrices from two 1xN GradientVectors
   % - Input
    %   - lh_grad = left hemisphere vector
    %   - rh_grad = right hemisphere gradient vector
    %   - fsaverage = the resolution of the underlying fsaverage space
    %
    % Ouput
    %   - lh_grad_matrix = connected components for the left hemisphere
    %   - rh_grad_matrix = connected components for the right hemisphere
    %
    % Example
    %   - [lh_grad_matrix,rh_grad_matrix]=CBIG_gwMRF_gradient_vertices_to_matrix(lh_grad,rh_grad,'fsaverage6')
    %Written by Alexander Schaefer and CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md

    lh_avg_mesh6 = CBIG_ReadNCAvgMesh('lh', fsaverage, 'white', 'cortex'); %which fsaverage is given from outside
    rh_avg_mesh6 = CBIG_ReadNCAvgMesh('rh', fsaverage, 'white', 'cortex');


    for i=1:max(size(lh_avg_mesh6.vertexNbors))
        for j=1:min(size(lh_avg_mesh6.vertexNbors))
            if (lh_avg_mesh6.vertexNbors(j,i)>0)
                lh_grad_matrix(i,j)=(lh_grad(i)+lh_grad(lh_avg_mesh6.vertexNbors(j,i)))/2;
                rh_grad_matrix(i,j)=(rh_grad(i)+rh_grad(rh_avg_mesh6.vertexNbors(j,i)))/2;
            end
        end
    end

end
