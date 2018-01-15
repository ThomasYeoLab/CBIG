%% Example 2 of applying author topic model (AT)
%
% This example shows how to run the AT Gibbs sampler and extract samples
% from the same and different chains. 

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
BURNIN   = 100; % the number of iterations before taking samples
LAG      = 10; % the lag between samples
NSAMPLES = 2; % the number of samples for each chain
NCHAINS  = 2; % the number of chains to run 

%%
% What output to show (0=no output; 1=iterations; 2=all output)
OUTPUT = 1;

%% 
% The starting seed number
SEED = 1;

%%
% This function might require 30-45 minutes of compute time
fprintf( 'The following computation might take 30-45 minutes...\n' );

for c=1:NCHAINS
    SEED = SEED + 1; 
    N = BURNIN;  
    fprintf( 'Running Gibbs sampler for burnin\n' );
    [ WP, AT , Z , X ] = GibbsSamplerAT( WS , DS , AD , T , N , ALPHA , BETA , SEED , OUTPUT );
        
    fprintf( 'Continue to run sampler to collect samples\n' );
    for s=1:NSAMPLES
        filename = sprintf( 'at_chain%d_sample%d' , c , s );
        fprintf( 'Saving sample #%d from chain #%d: filename=%s\n' , s , c , filename );
        comm = sprintf( 'save ''%s'' WP AT Z X ALPHA BETA SEED N T s c' , filename );
        eval( comm );
        
        if (s < NSAMPLES)
           N = LAG;
           SEED = SEED + 1; % important -- change the seed between samples !!
           [ WP, AT , Z , X ] = GibbsSamplerAT( WS , DS , AD , T , N , ALPHA , BETA , SEED , OUTPUT , Z , X );
        end
    end
end

%%
% Load in the samples and show some topic-word and author-topic
% distributions
for c=1:NCHAINS
    for s=1:NSAMPLES
        filename = sprintf( 'lda_chain%d_sample%d' , c , s );
        fprintf( 'Loading sample #%d from chain #%d: filename=%s\n' , s , c , filename );
        comm = sprintf( 'load ''%s''' , filename );
        eval( comm );

        WPM{1} = WP; WPM{2} = AT;
        BETAM(1)=BETA; BETAM(2) = ALPHA;
        WOM{1}=WO; WOM{2}=AN;
        %%
        % Write the word topic and author topic distributions to a text file
        [ SM ] = WriteTopicsMult( WPM , BETAM , WOM , 7 , 0.7 , 4 , filename );

        fprintf( 'The first ten topic-word distributions:\n' );
        SM{1}(1:10)

        fprintf( 'The first ten author-topic distributions:\n' );
        SM{2}(1:10)
    end
end
