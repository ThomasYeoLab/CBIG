%% function GibbsSamplerAT 

%%
% Runs the AT Gibbs sampler for the Author-Topic model. 
%
% |[ WP , AT, Z, X ] = GibbsSamplerAT( WD,AD,T,N,ALPHA,BETA,SEED,OUTPUT )|
% will run the Gibbs sampler for the AT model on a bag of words sparse count
% matrix |WD| of size |W| x |D| (|W|=vocabulary size, |D|=number of documents) and a
% author document sparse matrix |AD| of size |A| x |D| (|A|=number of authors).
% |WD(i,j)| contains the frequency of word |i| in document |j|. |AD(i,j)|
% contains binary values indicating the presence and absence of author |i|
% on document |j|. The first output is matrix |WP|, a |W| x |T| matrix,
% where |WP(i,j)| contains the number of times word |i| has been assigned
% to topic |j|. |T| is the number of topics. The second output is matrix |AT|, 
% a matrix of size |A| x |T|, where |AT(i,j)| contains the number of times a 
% word in a document with author |i| has been assigned to topic |j|. The outputs
% Z and X are vectors that contain the topic and author assignments
% respectively.
%
% |[ WP , AT, Z, X ] = GibbsSamplerAT( WD,AD,T,N,ALPHA,BETA,SEED,OUTPUT,ZIN,XIN )|
% will run the sampler from starting state |ZIN|, and |XIN| where |ZIN(k)| contains
% the topic assignment for token k, and |XIN(k)| contains the author assignment 
% for token k, saved from a previous sample
% 
%
%% 
% Notes
%
% |WD| and |AD| should be sparse matrices in double precision
%
% |N| determines the number of iterations to run the Gibbs sampler
%
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
% topics, the average number of authors on documents, and the number of
% iterations. The memory requirements scale linearly with the number of
% topics, documents, and authors. 
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
% Equivalence of LDA and AT model
%
% As shown in Rosen-Zvi et al. (2004), the AT model is identical to the LDA
% model if each document has a single unique author. Therefore, if the AD
% matrix is diagonal, the output of the AT model and the LDA model are
% identical

%% 
% References
%%
% * Steyvers, M., Smyth, P., Rosen-Zvi, M., & Griffiths, T. (2004).
% Probabilistic Author-Topic Models for Information Discovery. The 
% Tenth ACM SIGKDD International Conference on Knowledge Discovery and 
% Data Mining. Seattle, Washington.    
% * Rosen-Zvi, M., Griffiths T., Steyvers, M., & Smyth, P. (2004). The
% Author-Topic Model for Authors and Documents. In 20th Conference on
% Uncertainty in Artificial Intelligence. Banff, Canada    
