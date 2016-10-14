function s = struct(this)
% Struct method for GIfTI objects
% FORMAT s = struct(this)
% this   -  GIfTI object
% s      -  a structure containing public fields of the object
%__________________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% Guillaume Flandin
% $Id: struct.m 4716 2012-04-19 11:01:32Z guillaume $

names = fieldnames(this);
names = unique(names);
values = cell(length(names), length(this(:)));

for i=1:length(names)
    [values{i,:}] = subsref(this(:), substruct('.',names{i}));
end
s  = reshape(cell2struct(values,names,1),size(this));