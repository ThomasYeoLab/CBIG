function [lh_ci,rh_ci,lh_sizes,rh_sizes]=CBIG_gwMRF_generate_components(lh_avg_mesh,rh_avg_mesh,lh_labels,rh_labels)

    % This function computes connected components for a labeled vector.
    %
    % - Input
    %   - lh_avg_mesh = left hemisphere mesh
    %   - rh_avg_mesh = right hemisphere mesh
    %   - lh_labels = label vector for left hemisphere for corresponding mesh
    %   - rh_labels = label vector for right hemisphere for corresponding mesh
    %
    % Ouput
    %   - lh_ci = connected components for the left hemisphere
    %   - rh_ci = connected components for the right hemisphere
    %   - lh_sizes = size of the corressponding components for the left hemisphere
    %   - rh_sizes = size of the corressponding components for the right hemisphere
    %
    % Example
    %   - [lh_ci,rh_ci,lh_sizes,rh_sizes]=CBIG_gwMRF_generate_components(lh_avg_mesh,rh_avg_mesh,lh_labels,rh_labels)
    %Written by Alexander Schaefer and CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md

    tic
    lh_ini_labels = lh_labels;
    n = size(lh_avg_mesh.vertices, 2);
    neighbors=min(size(lh_avg_mesh.vertexNbors));%usually 6
    b = zeros(n*neighbors, 1);
    c = zeros(n*neighbors, 1);
    d = zeros(n*neighbors, 1);
    for i = 1:12 % vertices with fewer neighbors
        b((i-1)*neighbors+1:i*neighbors-1)=i;
        c((i-1)*neighbors+1:i*neighbors-1)=lh_avg_mesh.vertexNbors(1:neighbors-1, i);
        d((i-1)*neighbors+1:i*neighbors-1)=(lh_ini_labels(lh_avg_mesh.vertexNbors(1:neighbors-1, i))==lh_ini_labels(i));
    end
    for i = 13:n % vertices with regual neighbors
        b((i-1)*neighbors+1:i*neighbors)=i;
        c((i-1)*neighbors+1:i*neighbors)=lh_avg_mesh.vertexNbors(:, i);
        d((i-1)*neighbors+1:i*neighbors)=(lh_ini_labels(lh_avg_mesh.vertexNbors(:, i))==lh_ini_labels(i));
    end
    for i = 12:-1:1% account for zero vertices
        b(i*neighbors)=[];
        c(i*neighbors)=[];
        d(i*neighbors)=[];
    end
    g = sparse(b,c,d);
    toc
    [lh_ci, lh_sizes] = components(g);
    
       
    clear g;
    tic
    rh_ini_labels = rh_labels;
    neighbors=min(size(rh_avg_mesh.vertexNbors));%usually 6
    b = zeros(n*neighbors, 1);
    c = zeros(n*neighbors, 1);
    d = zeros(n*neighbors, 1);
    for i = 1:12 % vertices with fewer neighbors
        b((i-1)*neighbors+1:i*neighbors-1)=i;
        c((i-1)*neighbors+1:i*neighbors-1)=rh_avg_mesh.vertexNbors(1:neighbors-1, i);
        d((i-1)*neighbors+1:i*neighbors-1)=(rh_ini_labels(rh_avg_mesh.vertexNbors(1:neighbors-1, i))==rh_ini_labels(i));
    end
    for i = 13:n % vertices with regual neighbors
        b((i-1)*neighbors+1:i*neighbors)=i;
        c((i-1)*neighbors+1:i*neighbors)=rh_avg_mesh.vertexNbors(:, i);
        d((i-1)*neighbors+1:i*neighbors)=(rh_ini_labels(rh_avg_mesh.vertexNbors(:, i))==rh_ini_labels(i));
    end
    for i = 12:-1:1% account for zero vertices
        b(i*neighbors)=[];
        c(i*neighbors)=[];
        d(i*neighbors)=[];
    end
    
    g = sparse(b,c,d);
    toc
    [rh_ci, rh_sizes] = components(g);
end
