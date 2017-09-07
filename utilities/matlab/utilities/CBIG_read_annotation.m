function vertex_label = CBIG_read_annotation(annot_file)

% vertex_label = CBIG_read_annotation(annot_file)
%
% This function takes in a .annot file and converts it to a Nx1 column vector vertex_label
% where N is the number of vertices.
% Each entry in vertex_label represents the label of the corresponding vertex.
% If you only wish to convert .annot file to colortable, please use the
% original read_annotation.m function.
%
% Input:
%     - annot_file:
%       A .annot file.
%
% Output:
%     - vertex_label
%       A Nx1 column vector where N is the number of vertices. Each entry is the label of the corresponding vertex.
%
% Example:
% vertex_label = CBIG_read_annotation('CBIG_private/stable_projects/brain_parcellation/Schaefer2018_LocalGlobal/Parcellations/FreeSurfer5.3/fsaverage5/label/lh.Schaefer2018_400Parcels_17Networks_order.annot')
%
% Original Author: Bruce Fischl
% CVS Revision Info:
%    $Author: nicks $
%    $Date: 2007/01/10 22:55:09 $
%    $Revision: 1.4 $
%
% Copyright (C) 2002-2007,
% The General Hospital Corporation (Boston, MA).
% All rights reserved.
%
% Distribution, usage and copying of this software is covered under the
% terms found in the License Agreement file named 'COPYING' found in the
% FreeSurfer source code root directory, and duplicated here:
% https://surfer.nmr.mgh.harvard.edu/fswiki/FreeSurferOpenSourceLicense
%
% General inquiries: freesurfer@nmr.mgh.harvard.edu
% Bug reports: analysis-bugs@nmr.mgh.harvard.edu
%
% Written by Siyi Tang and CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md

%% get colortable and annotation labels from .annot file
fp = fopen(annot_file, 'r', 'b');

if(fp < 0)
    disp('Annotation file cannot be opened');
    return;
end

A = fread(fp, 1, 'int');

tmp = fread(fp, 2*A, 'int');
label = tmp(2:2:end);

bool = fread(fp, 1, 'int');

% no colortable
if(isempty(bool))
    disp('No Colortable found.');
    colortable = struct([]);
    fclose(fp);
    return;
end

% read colortable
if(bool)
    
    numEntries = fread(fp, 1, 'int');
    
    if(numEntries > 0)
        
        disp('Reading from Original Version');
        colortable.numEntries = numEntries;
        len = fread(fp, 1, 'int');
        colortable.orig_tab = fread(fp, len, '*char')';
        colortable.orig_tab = colortable.orig_tab(1:end-1);
        
        colortable.struct_names = cell(numEntries,1);
        colortable.table = zeros(numEntries,5);
        for i = 1:numEntries
            len = fread(fp, 1, 'int');
            colortable.struct_names{i} = fread(fp, len, '*char')';
            colortable.struct_names{i} = colortable.struct_names{i}(1:end-1);
            colortable.table(i,1) = fread(fp, 1, 'int');
            colortable.table(i,2) = fread(fp, 1, 'int');
            colortable.table(i,3) = fread(fp, 1, 'int');
            colortable.table(i,4) = fread(fp, 1, 'int');
            colortable.table(i,5) = colortable.table(i,1) + colortable.table(i,2)*2^8 + colortable.table(i,3)*2^16 + colortable.table(i,4)*2^24;
        end
        disp(['colortable with ' num2str(colortable.numEntries) ' entries read (originally ' colortable.orig_tab ')']);
        
    else
        version = -numEntries;
        if(version~=2)
            disp(['Error! Does not handle version ' num2str(version)]);
        else
            disp(['Reading from version ' num2str(version)]);
        end
        numEntries = fread(fp, 1, 'int');
        colortable.numEntries = numEntries;
        len = fread(fp, 1, 'int');
        colortable.orig_tab = fread(fp, len, '*char')';
        colortable.orig_tab = colortable.orig_tab(1:end-1);
        
        colortable.struct_names = cell(numEntries,1);
        colortable.table = zeros(numEntries,5);
        
        numEntriesToRead = fread(fp, 1, 'int');
        for i = 1:numEntriesToRead
            structure = fread(fp, 1, 'int')+1;
            if (structure < 0)
                disp(['Error! Read entry, index ' num2str(structure)]);
            end
            if(~isempty(colortable.struct_names{structure}))
                disp(['Error! Duplicate Structure ' num2str(structure)]);
            end
            len = fread(fp, 1, 'int');
            colortable.struct_names{structure} = fread(fp, len, '*char')';
            colortable.struct_names{structure} = colortable.struct_names{structure}(1:end-1);
            colortable.table(structure,1) = fread(fp, 1, 'int');
            colortable.table(structure,2) = fread(fp, 1, 'int');
            colortable.table(structure,3) = fread(fp, 1, 'int');
            colortable.table(structure,4) = fread(fp, 1, 'int');
            colortable.table(structure,5) = colortable.table(structure,1) + colortable.table(structure,2)*2^8 + colortable.table(structure,3)*2^16 + colortable.table(structure,4)*2^24;
        end
        disp(['colortable with ' num2str(colortable.numEntries) ' entries read (originally ' colortable.orig_tab ')']);
        
        % check if colortable is valid
        if numel(unique(colortable.table(:,5))) ~= size(colortable.table,1)
            disp('Error! Colortable has repeatd values!');
        end
        
    end
else
    disp('Error! Should not be expecting bool = 0');
end

fclose(fp);

%% match annotation labels with labels of vertices
vertex_label = zeros(size(label));
for i = 1:numel(label)
    vertex_label(i) = find(colortable.table(:,5)==label(i));
end

end


