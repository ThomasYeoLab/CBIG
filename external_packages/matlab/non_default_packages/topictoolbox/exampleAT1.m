%% Example 1 of applying author topic model (AT)
%
% This example shows how to run the AT Gibbs sampler on a small dataset to
% extract a set of topics. This code will produce a single |WP| and |AT| count matrix
% and will show the most likely words and authors per topic. It also writes the results
% to a file

load 'bagofwords_nips'; % load the nips word document dataset
load 'authordoc_nips'; % load the author-document pairings for nips
load 'words_nips'; % Load the vocabulary
load 'authors_nips'; % Load the author names

%%
% The text file to show the topic-word and topic-author distributions
filename = 'topics_nips_at.txt';
 
%%
% Set the number of topics
T = 50; 

%%
% Set the hyperparameters
BETA  = 0.01;
ALPHA = 50/T;

%%
% The number of iterations
N = 50; 

%%
% The random seed
SEED = 3;

%%
% What output to show (0=no output; 1=iterations; 2=all output)
OUTPUT = 1;

%%
% This function might require 30-45 minutes of compute time
fprintf( 'The following computation might take 30-45 minutes...\n' );
tic
[ WP, AT , Z , X ] = GibbsSamplerAT( WS , DS , AD , T , N , ALPHA , BETA , SEED , OUTPUT );
toc

%%
% Just in case, save this WP sample and associated information 
save 'temp' WP AT Z X ALPHA BETA SEED N;

WPM{1} = WP; WPM{2} = AT; 
BETAM(1)=BETA; BETAM(2) = ALPHA;
WOM{1}=WO; WOM{2}=AN;

%%
% Write the word topic and author topic distributions to a text file
[ SM ] = WriteTopicsMult( WPM , BETAM , WOM , 7 , 0.7 , 4 , filename );

%%
% Show the most likely words in the first ten topics
SM{1}(1:10)

%%
% Show the most likely authors in the first ten topics
SM{2}(1:10) 
