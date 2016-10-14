%% Example 1 to visualize topic model results 
% This example shows how to visualize topics in a 2D map. 
%
% load a document-topic count matrix saved for the nips dataset
load 'ldasingle_nips';
load 'words_nips';

%%
% extract the topics in a cell array of strings
[S]=WriteTopics( WP,BETA,WO,5,0.6 );

fprintf( 'Please wait while calculating visualization...\n' );
drawnow;

%%
% visualize these topics in a 2D map. Have each topic by displayed
% vertically
VisualizeTopics( DP,ALPHA,S,'vertical');