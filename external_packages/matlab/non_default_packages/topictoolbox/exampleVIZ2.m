%% Example 2 to visualize topic model results
%
% This example shows how to visualize documents in a 2D map.
% The documents are located such that neighboring documents have similar
% topic distributions. Documents are shown by their corresponding titles

%%
% Load a document-topic count matrix saved for the nips dataset
load 'ldasingle_nips';

%%
% Get the title descriptions of the papers
load 'titles_nips';

%%
% We cannot visualize every paper, take a random subset
SEED = 5;
rand( 'state' , SEED ); subs = randperm( size( DP,1));

NS = 200; 

%%
% These are the indices of the random subset of documents
subs = subs( 1:200 );

%%
% Visualize the documents of this subset
NCHARS   = 30; % maximum number of characters on one line
MAXLINES = 3; % maximum number of lines

fprintf( 'Please wait while calculating visualization...\n' );
drawnow;

VisualizeDocs( DP(subs,:) , ALPHA , titles(subs) , NCHARS , MAXLINES )