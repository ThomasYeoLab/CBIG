%% Function GibbsSamplerHMMLDA
%
%%
% Runs the Gibbs sampler for the HMM-LDA model 
%
%%
% |[WP,DP,MP,Z,X]=GibbsSamplerHMMLDA(WS,DS,T,NS,N,ALPHA,BETA,GAMMA,SEED,OUT
% PUT)| will run the Gibbs sampler for the HMM-LDA model on a sequence of
% word indices |WS| and document indices |DS|, both of size |1| x |n| where
% |n| is the number of word tokens. The word stream |WS| contains word
% indices in order of occurence with |WS|=0 representing the
% end-of-sentence-end marker. The document indices |DS| contains all
% document indices and |max(DS)| = |D| = number of documents. |T| is the
% number of topics, |NS| is the number of syntactic states, and |N| is the
% number of iterations for the Gibbs sampler. |ALPHA|, |BETA|, |GAMMA| are
% hyperparameters of the generative model (see reference for an explanation
% of these parameters). The first output is a sparse matrix |WP|, of size
% |W| x |T| where |WP(i,j)| contains the number of times word |i| has been
% assigned to topic |j|. The second output is a sparse matrix |DP|, a |D| x
% |T| matrix, where |DP(i,j)| contains the number of times a word in
% document |d| has been assigned to topic |j|.  The third output is a
% sparse matrix |MP| of size |W x NS|, where |MP(i,j)| contains the number
% of times word |i| has been assigned to syntactic state |j|. The vectors
% |Z| and |X| are both of size 1 x |n| containing the topic and hmm-state
% assignments respectively.
%
% |[WP,DP,MP,Z,X]=GibbsSamplerHMMLDA(WS,DS,T,NS,N,ALPHA,BETA,GAMMA,SEED,OUT
% PUT,ZIN,XIN)| will run the sampler from a previous state as specified by
% |ZIN| and |XIN|
%%
% NOTES
%
% |WS| and |DS| should be double precision vectors of indices. |WS(k)=0|
% when the kth position in the text is the end-of-sentence marker.
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
%%
% REFERENCES
%%
% * Griffiths, T.L., & Steyvers, M.,  Blei, D.M., & Tenenbaum, J.B. (2004). 
% Integrating Topics and Syntax. In: Advances in Neural Information
% Processing Systems, 17. 
%

