%% Function GibbsSamplerLDA
% Runs the Gibbs sampler for the Latent Dirichlet Allocation (LDA) model
%
% |[ WP,DP,Z ] = GibbsSamplerLDA( WS,DS,T,N,ALPHA,BETA,SEED,OUTPUT )| will
% run the Gibbs sampler for the LDA model on a bag of words data provided
% by the vectors |WS| and |DS|. |WS(k)| and |DS(k)| contain the word and
% document indices for the kth token. The maximum of |WS| is |W|, the
% vocabulary size. The maximum of |DS| is |D|, the number of documents. |T|
% is the number of topics. The first output is the sparse matrix |WP|, of
% size |W| x |T|, where |WP(i,j)| contains the number of times word |i| has
% been assigned to topic |j|.   The second output is |DP|, a sparse |D| x
% |T| matrix, where |DP(i,j)| contains the number of times a word in
% document |d| has been assigned to topic |j|. The third output |Z|
% contains the topic assignments; |Z(k)| contains the topic assignment for
% token k.    
%
% |[ WP,DP,Z ] = GibbsSamplerLDA( WS,DS,T,N,ALPHA,BETA,SEED,OUTPUT,ZIN )|
% will run the sampler from starting state |ZIN|, where |ZIN(k)| contains
% the topic assignment for token k, saved from a previous sample.
%
% NOTES
%
% |WS| and |DS| should be in double precision
% |N| determines the number of iterations to run the Gibbs sampler.
% |ALPHA| and |BETA| are the hyperparameters on the Dirichlet priors for 
% the topic distributions (|theta|) and the topic-word distributions (|phi|)
% respectively
% 
% |SEED| sets the seed for the random number generator
%
% |OUTPUT| determines the screen output by the sampler
%   0 = no output provided
%   1 = show the iteration number only
%   2 = show all output
%
% The time to complete the procedure scales linearly with the number of
% topics and the number of iterations. The memory requirements scale
% linearly with the number of topics and documents.
%
% A good setting for the number of iterations will depend on the number of
% topics and the complexity of problem. For most problems, 500 to 2000
% iterations will suffice.
%
% Appropriate values for |ALPHA| and |BETA| depend on the number of topics and
% the number of words in vocabulary. For most applications, good results
% can be obtained by setting |ALPHA = 50 / T| and |BETA = 200 / W|
%
% The sampler uses its own random number generator and setting the seed for
% this function will not influence the random number seed for Matlab
% functions
%
%%
% REFERENCES
%%
% * Griffiths, T., & Steyvers, M. (2004). Finding Scientific Topics. 
% Proceedings of the National Academy of Sciences, 101 (suppl. 1), 5228-5235. 
% * Steyvers, M. & Griffiths, T. (in press). Probabilistic topic models. In T. 
% Landauer, D. McNamara, S. Dennis & W. Kintsch Eds. Latent Semantic Anlysis: 
% A road to meaning. Erlbaum
% * D. Blei, A. Ng, and M. Jordan (2003). Latent Dirichlet allocation.
% Journal of Machine Learning Research, 3:993-1022.


