%% Example 2 of running HMM-LDA topic model
%
% This example shows how to collect multiple samples from the HMM-LDA Gibbs sampler 
% from the same chain and different chains 

%% 
% Choose the dataset
dataset = 1; % 1 = psych review; 2 = nips papers

if (dataset == 1)
    % Load the psych review word stream
    load 'psychreviewstream';

    % Set the parameters for the model
    T      = 50;     % number of topics
    NS     = 12;     % number of syntactic states
    ALPHA  = 50 / T; % ALPHA hyperparameter
    BETA   = 0.01;   % BETA hyperparameter
    GAMMA  = 0.1;    % GAMMA hyperparameter    
end

if (dataset == 2)
    % Load the nips paper word stream
    load 'nips_stream';

    % Set the parameters for the model
    T      = 50;     % number of topics
    NS     = 16;     % number of syntactic states
    ALPHA  = 50 / T; % ALPHA hyperparameter
    BETA   = 0.01;   % BETA hyperparameter
    GAMMA  = 0.1;    % GAMMA hyperparameter
end

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
    [WP,DP,MP,Z,X]=GibbsSamplerHMMLDA( WS,DS,T,NS,N,ALPHA,BETA,GAMMA,SEED,OUTPUT);
    
    fprintf( 'Continue to run sampler to collect samples\n' );
    for s=1:NSAMPLES
        filename = sprintf( 'ldahmm_chain%d_sample%d' , c , s );
        fprintf( 'Saving sample #%d from chain #%d: filename=%s\n' , s , c , filename );
        comm = sprintf( 'save ''%s'' MP WP DP Z ALPHA BETA SEED N Z T s c' , filename );
        eval( comm );
        
        WPM{ s , c } = WP;
        MPM{ s , c } = MP;
        
        if (s < NSAMPLES)
           N = LAG;
           SEED = SEED + 1; % important -- change the seed between samples !!
           [WP,DP,MP,Z,X]=GibbsSamplerHMMLDA( WS,DS,T,NS,N,ALPHA,BETA,GAMMA,SEED,OUTPUT,Z,X);
        end
    end
end

%%
% Inspect the first few topics of a few samples
for c=1:NCHAINS
    for s=1:NSAMPLES
       
      [S] = WriteTopics( WPM{s,c} , BETA , WO , 7 , 0.8 ); 

      fprintf( '\n\nExample topic-word distributions of chain %d sample %d\n' , c , s );
      S(1:5)
      
      [S] = WriteTopics( MPM{s,c} , BETA , WO , 7 , 0.8 ); 

      fprintf( '\nExample hmm state-word distributions of chain %d sample %d\n' , c , s );
      S(1:5)
    end
end

