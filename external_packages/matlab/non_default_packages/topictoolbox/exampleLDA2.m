%% Example 2 of running basic topic model -- getting multiple samples
%
% This example shows how to run the LDA Gibbs sampler on a small dataset to
% produce *multiple* Gibbs samples from the same and different Markov
% chains

%%
% Choose the dataset
dataset = 1; % 1 = psych review abstracts 2 = NIPS papers

if (dataset == 1)
    % load the psych review data in bag of words format
    load 'bagofwords_psychreview'; 
    % Load the psych review vocabulary
    load 'words_psychreview' 
elseif (dataset == 2)
    % load the nips dataset
    load 'bagofwords_nips'; 
    % load the nips vocabulary
    load 'words_nips' 
end

%%
% Set the number of topics
T=50; 

%%
% Set the hyperparameters
BETA=0.01;
ALPHA=50/T;

%%
% What output to show (0=no output; 1=iterations; 2=all output)
OUTPUT = 1;

%%
% The number of iterations
BURNIN   = 100; % the number of iterations before taking samples
LAG      = 10; % the lag between samples
NSAMPLES = 2; % the number of samples for each chain
NCHAINS  = 2; % the number of chains to run 

%% 
% The starting seed number
SEED = 1;

for c=1:NCHAINS
    SEED = SEED + 1; 
    N = BURNIN;  
    fprintf( 'Running Gibbs sampler for burnin\n' );
    [ WP,DP,Z ] = GibbsSamplerLDA( WS , DS , T , N , ALPHA , BETA , SEED , OUTPUT );
    
    fprintf( 'Continue to run sampler to collect samples\n' );
    for s=1:NSAMPLES
        filename = sprintf( 'lda_chain%d_sample%d' , c , s );
        fprintf( 'Saving sample #%d from chain #%d: filename=%s\n' , s , c , filename );
        comm = sprintf( 'save ''%s'' WP DP Z ALPHA BETA SEED N Z T s c' , filename );
        eval( comm );
        
        WPM{ s , c } = WP; 
        
        if (s < NSAMPLES)
           N = LAG;
           SEED = SEED + 1; % important -- change the seed between samples !!
           [ WP,DP,Z ] = GibbsSamplerLDA( WS , DS , T , N , ALPHA , BETA , SEED , OUTPUT , Z );
        end
    end
end

%%
% Inspect the first few topics of a few samples
for c=1:NCHAINS
    for s=1:NSAMPLES
       
      [S] = WriteTopics( WPM{s,c} , BETA , WO ); 

      fprintf( 'Example topics of chain %d sample %d\n' , c , s );
      S(1:5)
    end
end



