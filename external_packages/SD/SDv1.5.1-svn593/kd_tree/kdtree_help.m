%
% KDTREE Find closest points using a k-D tree.
% 
%  CP = KDTREE( REFERENCE, MODEL ) finds the closest points in
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
%    CP is the same dimension as MODEL. There is a one-to-one
%    relationship between the rows of MODEL and the rows of CP. The
%    i-th row (point) of CP is a row (point) from REFERENCE which
%    is closest to the i-th row (point) of MODEL. The "closest"
%    metric is defined as the D-dimensional Euclidean (2-norm)
%    distance.
%
%  
%  [CP, DIST] = KDTREE( ... ) returns the distances between
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
%      [CP, DIST, ROOT] = KDTREE(REFERENCE, MODEL) where ROOT
%      receives a pointer to the root of the k-D tree.
%
%  Subsequently, use the following call to pass the k-D tree back
%  into the mex function:
%
%      [CP, DIST, ROOT] = KDTREE([], MODEL, ROOT)
% 
%  Note that ROOT is again an output, preventing the tree from
%  being removed from memory. 
%
%  Ultimately, to clear the k-D tree from memory, pass ROOT as
%  input, but do not receive it as output:
%
%      KDTREE([], [], ROOT)
%
%  New since June 2004: This k-D tree library now handles points
%  with dimension greater than 3.
%
%  See also KDTREEIDX and KDRANGEQUERY.
% 
%  Written by / send comments or suggestions to :
%     Guy Shechter
%     guy at jhu dot edu
%     June 2004
%
