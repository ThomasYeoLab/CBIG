%
% KDTREEIDX Find closest points using a k-D tree.
% 
%  IDX = KDTREEIDX( REFERENCE, MODEL ) finds the closest points in
%  REFERENCE for each point in MODEL. The search is performed in an
%  efficient manner by building a k-D tree from the datapoints in
%  REFERENCE, and querying the tree for each datapoint in
%  MODEL. 
%
%  Input :
%    REFERENCE is an NxD matrix, where each row is a D-dimensional
%    point. MODEL is an MxD matrix, where each row is a D-dimensional
%    query point. 
%
%  Output:
%    IDX is a vector of length M. The i-th value of IDX is the row
%    index into the matrix REFERENCE, which is the closest point to
%    the i-th row (point) of MODEL. The "closest"  metric is
%    defined as the D-dimensional Euclidean (2-norm) distance.
%    The closest point values can be found by: CP = REFERENCE(IDX,:)
%
%  
%  [IDX, DIST] = KDTREEIDX( ... ) returns the distances between
%  each row of MODEL and its closest point match from the k-D tree
%  in the vector DIST. DIST(i) corresponds to the i-th row (point)
%  of MODEL.
%
%  The default behavior of the function is that the k-D tree is
%  destroyed when the function returns. If you would like to save
%  the k-D tree in memory for use at a later time for additional
%  queries on the same REFERENCE data, then call the function with
%  an additional output:
%
%      [IDX, DIST, ROOT] = KDTREEIDX(REFERENCE, MODEL) where ROOT
%      receives a pointer to the root of the k-D tree.
%
%  Subsequently, use the following call to pass the k-D tree back
%  into the mex function:
%
%      [IDX, DIST, ROOT] = KDTREEIDX([], MODEL, ROOT)
% 
%  Note that ROOT is again an output, preventing the tree from
%  being removed from memory. 
%
%  Ultimately, to clear the k-D tree from memory, pass ROOT as
%  input, but do not receive it as output:
%
%      KDTREEIDX([], [], ROOT)
%
%
%  See also KDTREE and KDRANGEQUERY.
% 
%  Written by / send comments or suggestions to :
%     Guy Shechter
%     guy at jhu dot edu
%     June 2004
%