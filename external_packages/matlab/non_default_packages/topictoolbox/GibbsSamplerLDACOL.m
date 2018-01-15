%% function GibbsSamplerLDACOL 

%%
% Runs the Gibbs sampler for the LDA-COL (Collocation) model. 
%
% |[ WP,DP,WC,C,Z ] = 
%      GibbsSamplerLDACOL( WS , DS , SI , WW , T , N , ALPHA ,
%                          BETA , GAMMA0, GAMMA1 , DELTA , SEED , OUTPUT
%                          );|
% 
%%
% will run the Gibbs sampler for the LDA-COL model on a word stream 
% specified by word |WS|, document |DS|, and status |SI| indices. All these
% vectors are of size 1 x |n| where |n| is the number of word tokens. In
% |SI|, a status value of 1 at position |k| indicates that the word at
% position |k| can form a collocation with the word at position |k|-1. This
% indicator is necessary to flag (with value SI=0) word positions that
% cross document boundaries or that cross word positions that originally
% had stop words in between. |WW| is a |W| x |W| sparse matrix where |W| is
% the number of words in the vocabulary. |WW(i,j)| contains the count of
% the number of times that word |i| follows word |j| in the word stream.
%
%% 
% The output matrices |WP| and |DP| contain the word-topic and
% document-topic counts for all words that were assigned to the topic
% route. |WC( i )| counts the number of times that word |i| led to
% a collocation with the next word in the word stream. |C| is of size 1 x |n| where |C( i )|=0 if word token |i|
% was assigned to the topic route, and |C( i )|=1 otherwise. |Z| is of size
% 1 x |n| and contains the topic assignments for all word tokens
%
% |[ WP,DP,WC,C,Z ] = 
%      GibbsSamplerLDACOL( WS , DS , SI , WW , T , N , ALPHA ,
%                          BETA , GAMMA0, GAMMA1 , DELTA , SEED , OUTPUT ,
%                          C_IN , Z_IN );|
%
% will run the model without random initialization. The initial conditions
% are specified by |C_IN| and |Z_IN|. This allows the model to continue
% from a previous state and provides the possibility of extracting multiple
% Gibbs samples from a single Markov chain.  
%
%
% Notes
%
%
% |N| determines the number of iterations to run the Gibbs sampler
%
% |ALPHA|, |BETA|, |DELTA|, |GAMMA0| and |GAMMA1| are the hyperparameters
% in the model.  
% 
% |SEED| sets the seed for the random number generator
%
% |OUTPUT| determines the screen output by the sampler
%   0 = no output provided 1 = show the iteration number only 2 = show all
%   output
%
% The sampler uses its own random number generator and setting the seed for
% this function will not influence the random number seed for Matlab
% functions
% 
%% 
% References
%%
% Tom Griffiths, Technical Report, July 18, 2005.

