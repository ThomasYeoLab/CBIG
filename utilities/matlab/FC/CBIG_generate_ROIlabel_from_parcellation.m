function CBIG_generate_ROIlabel_from_parcellation(lh_label, rh_label, labelname, mesh, outdir)

% CBIG_generate_ROIlabel_from_parcellation(lh_label, rh_label, labelname, mesh, outdir)
% this function will generate a .label file from two vectors lh_label
% and rh_label which are assumed to be of size 1xN
% for e.g. lh_label and rh_label is of size 1x163842
% labelname should be like 'craddock' or 'shen'
%
% Written by CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md

    
    if(size(lh_label, 1) ~= 1)
        error('assume lh_label is 1 x N'); 
    end

    if(size(rh_label, 1) ~= 1)
        error('assume rh_label is 1 x N'); 
    end

    if(size(lh_label, 2) ~= size(rh_label, 2))
        error('size of lh_label not equal rh_label');
    end

    lh_highest_value = max(lh_label);
    rh_highest_value = max(rh_label);
    
    % read in fsaverage mesh struct
    if (~isempty(strfind(mesh, 'fs_LR')))
        lh_avg_mesh = CBIG_read_fslr_surface('lh',mesh,'inflated');
        rh_avg_mesh = CBIG_read_fslr_surface('rh',mesh,'inflated');
    else
        lh_avg_mesh = CBIG_ReadNCAvgMesh('lh',mesh,'inflated','cortex');
        rh_avg_mesh = CBIG_ReadNCAvgMesh('rh',mesh,'inflated','cortex');
    end
    % write out .label file for each region in left hemishpere
    file_lh = fopen(strcat('lh.', labelname, '_components.txt'),'w');
    for i = 0:lh_highest_value
       lindex =  find(lh_label==i); % vertex number
       lxyz = lh_avg_mesh.vertices(:,lindex)'; % coordinates of the vertex
       lvals = zeros(1,length(lindex));
       labelfile = [outdir,'/',strcat('lh.',labelname,'_parcel_', num2str(i), '.label')];
       
       if (~isempty(lindex))
            write_label(lindex-1,lxyz,lvals,labelfile);
            % write to txt file
            fprintf(file_lh, '%s\n', labelfile);
       end
       
       
    end
    
    % write out .label file for each region in right hemisphere
    file_rh = fopen(strcat('rh.', labelname, '_components.txt'),'w');
    for i = 0:rh_highest_value
       lindex =  find(rh_label==i);
       lxyz = rh_avg_mesh.vertices(:,lindex)';
       lvals = zeros(1,length(lindex));
       labelfile = [outdir,'/',strcat('rh.',labelname,'_parcel_', num2str(i), '.label')];
       
       if (~isempty(lindex))
            write_label(lindex-1,lxyz,lvals,labelfile);
            % write to txt file
            fprintf(file_rh, '%s\n', labelfile);
       end 
    end
    
end
