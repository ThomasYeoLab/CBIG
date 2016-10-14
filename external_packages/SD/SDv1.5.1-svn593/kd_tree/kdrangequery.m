%
% KDRANGEQUERY Find all points within a local neighborhood.
% 
%  PTS = KDRANGEQUERY( ROOT, QUERYPT, DISTLIM ) finds all the
%  points stored in the k-D tree ROOT that are within DISTLIM units
%  from the QUERYPT. Proximity is quantified usng a D-dimensional
%  Euclidean (2-norm) distance.
%
%  Input :
%    ROOT is a pointer to a k-D tree which must be constructed with
%    the function KDTREE or KDTREEIDX. QUERYPT is a 1xD vector
%    representing a point in D-dimensional space. DISTLIM is a
%    scalar which specifies the radius of the neighborhood around
%    QUERYPT.
%
%  Output:
%    PTS is an NxD matrix, where each row is a datapoint from the
%    k-D tree ROOT. Each of these datapoints is found within a
%    distance DISTLIM from QUERYPT.
%
%  [PTS, DIST] = KDRANGEQUERY( ... ) returns the distances between
%  each row of PTS and QUERYPT in the Nx1 vector DIST.
%
%  [PTS, DIST, IDX ] = KDRANGEQUERY( ... ) returns the index value
%  for each row (point) of PTS. The index value maps back to a row
%  from the matrix REFERENCE used to build the k-D tree (see the 
%  KDTREE or KDTREEIDX functions).
%
%  Limitations: 
%    QUERYPT must be a 1xD dimensional array meaning that the range
%    query can be performed for one point at a time.
%
%  See also KDTREE and KDTREEIDX.
%
%  Written by / send comments or suggestions to :
%     Guy Shechter
%     guy at jhu dot edu
%     June 2004
%