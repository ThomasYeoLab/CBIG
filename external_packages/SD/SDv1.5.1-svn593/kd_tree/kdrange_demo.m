%% Demo of the KDRANGEQUERY function.
%%
%%   Guy Shechter
%%   guy at jhu dot edu
%%   June 2004



%% Generate a set of 1000 reference datapoints in R^2
ReferencePts = rand(1000,2); 

%%% Build the k-d Tree once from the reference datapoints.
[tmp, tmp, TreeRoot] = kdtree( ReferencePts, []);

%%% and find all the points in the k-d tree that are within 0.4
%%% units (D-dimensional Euclidean, 2-norm, distance) from the origin
[ PtsInNeighborhood, Dist ] = kdrangequery( TreeRoot, [0 0], 0.4 );

%%% Free the k-D Tree from memory.
kdtree([],[],TreeRoot);


figure; clf; hold on; axis equal

%% Plot all the points in the k-D tree
plot(ReferencePts(:,1), ReferencePts(:,2),'k.');

%% Draw a red circle around every point found in the neighborhood
plot(PtsInNeighborhood(:,1), PtsInNeighborhood(:,2), 'ro');

%% Show the arc with radius=0.4 centered at the origin.
t=0:.01:pi/2;
plot(0.4*cos(t), 0.4*sin(t),'g-','LineWidth',2);


