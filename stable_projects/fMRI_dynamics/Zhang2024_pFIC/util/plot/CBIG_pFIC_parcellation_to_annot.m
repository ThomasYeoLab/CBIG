function CBIG_pFIC_parcellation_to_annot(lh_label, rh_label, lh_filename, rh_filename, colortable)
% CBIG_pFIC_parcellation_to_annot(lh_label, rh_label, lh_filename, rh_filename, colortable)
% This function writes the parcellation file into a annot file which can be
% used for plotting figures
% Input:
%   - lh_label: 163842-by-1 parcellation label for left hemisphere
%   - rh_label: 163842-by-1 parcellation label for right hemisphere
%   - lh_filename: filename used to save left hemisphere annot file
%   - rh_filename: filename used to save right hemisphere annot file
%   - colortable: a n-by-3 matrix, n is the number colors, and 3 columns
%   correpsond to the RGB value 
% Written by CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md
    
    if(size(lh_label, 1) ~= 1)
        lh_label = lh_label';
    end
    
    if(size(rh_label, 1) ~= 1)
        rh_label = rh_label';
    end

    if(size(lh_label, 2) ~= size(rh_label, 2))
        error('size of lh_label not equal rh_label');
    end
    
    if(~exist('colortable', 'var'))
        colortable = 'jet'; 
    end
    
    label = [lh_label rh_label];
    vertices=0:(length(lh_label)-1);
    
    ct.numEntries = double(max(label)+1); %assume there is a 0 label, which is considered background/medial wall
    
    if(ischar(colortable))
        xcolormap = int32([30 10 30; feval(colortable, ct.numEntries-1)*255]);
        ct.orig_tab = colortable;
    else
        if(size(colortable, 2) ~= 3)
            error('If input colortable is a matrix, it should be should have 3 columns containing RGB values.');
        end
        xcolormap = int32(colortable); 
        ct.orig_tab = 'manual';
    end
    
    %! in annot files, medial wall is annotated as 1
    ct.table = ...
        [xcolormap, zeros([ct.numEntries 1],'int32'), xcolormap(:,1) + xcolormap(:, 2)*2^8 + xcolormap(:, 3)*2^16]; 
    
    new_lh_label = lh_label;
    new_rh_label = rh_label;
    
    for i=0:(ct.numEntries-1)
        if(i == 0)
            ct.struct_names{i+1} = 'unknown';
        else
            ct.struct_names{i+1}= ['parcel' num2str(i)];
        end
        new_lh_label(lh_label == i) = ct.table(i+1, 5);  
        new_rh_label(rh_label == i) = ct.table(i+1, 5); 
    end

    write_annotation(lh_filename, vertices, new_lh_label, ct)
    write_annotation(rh_filename, vertices, new_rh_label, ct)
end
