function [bins,binNames] = discretize(inputData,edges)
%DISCRETIZE Create a categorical array by discretizing numeric data.
%   BINS = DISCRETIZE(A,EDGES) returns the bin numbers obtained by discretizing
%   the numeric array A into categories.  EDGES is a numeric vector that
%   contains bin edges defining the categories.  EDGES has one more element than
%   B has categories, and the I-th category in B includes values in A where
%   EDGES(I) <= X < EDGES(I+1).  The uppermost category includes values in A
%   equal to the rightmost edge in EDGES.
%
%   [BINS,BINNAMES] = DISCRETIZE(A,EDGES) returns bins names of the form
%    "[A,B)", where A and B are consecutive values from EDGES.
%
%   See also CATEGORICAL, HISTC.

if ~isnumeric(inputData) || ~isreal(inputData)
    error(message('MATLAB:categorical:discretize:NonnumericData'));
elseif ~isnumeric(edges) || ~isreal(edges) || ~isvector(edges) || length(edges) < 2
    error(message('MATLAB:categorical:discretize:InvalidEdges'));
elseif length(edges)-1 > categorical.maxNumCategories
    error(message('MATLAB:categorical:MaxNumCategoriesExceeded',categorical.maxNumCategories));
end
nbins = length(edges) - 1;

if (nargout > 1)
    % Create names from the edges
    binNames = cell(1,nbins);
    for i = 1:nbins
        binNames{i} = sprintf('[%0.5g, %0.5g)',edges(i),edges(i+1));
    end
    binNames{end}(end) = ']';
    if length(unique(binNames)) < length(binNames)
        error(message('MATLAB:categorical:discretize:CantCreateCategoryNames'));
    end
end

% HISTC includes a rightmost "x==edges(end)" bin, remove it
[~,bins] = histc(inputData(:),edges);
bins(bins==nbins+1) = nbins;
bins = reshape(bins,size(inputData));
