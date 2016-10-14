%% Demo of the KDTREE and KDTREEIDX function.
%%
%%   Guy Shechter
%%   guy at jhu dot edu
%%   June 2004

%% Generate a set of 10000 reference datapoints in R^3
ReferencePts = rand(100000,3); 

%%% Build the k-d Tree once from the reference datapoints.
[tmp, tmp, TreeRoot] = kdtree( ReferencePts, []);

for Iteration = 1:10

  disp(sprintf('Iteration %3d: ',Iteration));
  
  %% Generate a random test set of 100 points in R^3
  TestPoints = rand(10000, 3);

  %% Find the closest point in ReferencePts for each random TestPoint
  [ ClosestPts, DistA, TreeRoot ] = kdtree([], TestPoints, ...
					       TreeRoot);

  %% Find the row index of the closest point in ReferencePts for
  %each random TestPoint
  [ ClosestPtIndex, DistB, TreeRoot ] = kdtreeidx([], TestPoints, ...
						  TreeRoot);
    
  %% Are the two solutions equivalent?
  IndexedPts = ReferencePts(ClosestPtIndex,:);
  if sum(sum(abs([ IndexedPts - ClosestPts ]))) == 0
    disp('KDTREE and KDTREEIDX found the same set of closest points');
  end
  
end

%%% Free the k-D Tree from memory.
kdtree([],[],TreeRoot);



