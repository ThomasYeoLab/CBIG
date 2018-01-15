%% Running LDA model on bag of words data from stream for use with 
% SWR model

%%
% Choose the dataset
dataset = 2; % 1 = psych review abstracts 2 = NIPS papers 5 = NYTimes

if (dataset == 1)
    load 'psychreviewbagofwordsfromstream'; % WO DS WS ORIGWSPOS;
    T      = 50; % Set the number of topics
    BETA   = 0.01;
    ALPHA  = 0.1;

    BURNIN   = 500; % the number of iterations before taking samples
    LAG      = 1; % the lag between samples
    NSAMPLES = 1; % the number of samples for each chain
    NCHAINS  = 10; % the number of chains to run

elseif (dataset == 2)
    load 'nipsbagofwordsfromstream'; % WO DS WS ORIGWSPOS;
    T      = 200; % Set the number of topics
    BETA   = 0.01;
    ALPHA  = 0.1;
    
    BURNIN   = 500; % the number of iterations before taking samples
    LAG      = 1; % the lag between samples
    NSAMPLES = 1; % the number of samples for each chain
    NCHAINS  = 10; % the number of chains to run
elseif (dataset == 5)
    load '..\kdd\NYTimesCollocation'; % load in variables: WW WO DS WS SI
    T      = 200; % Set the number of topics
    BETA   = 0.01;
    ALPHA  = 0.1;
    
    BURNIN   = 500; % the number of iterations before taking samples
    LAG      = 1; % the lag between samples
    NSAMPLES = 1; % the number of samples for each chain
    NCHAINS  = 10; % the number of chains to run
end

OUTPUT = 1; % What output to show (0=no output; 1=iterations; 2=all output)
SEED   = 1; % The starting seed number


mkdir( sprintf( 'swrsamples_%d' , dataset )); 

for c=1:NCHAINS
    SEED = SEED + 1; 
    N = BURNIN;  
    fprintf( 'Running Gibbs sampler for burnin\n' );
    [ WP,DP,Z ] = GibbsSamplerLDA( WS , DS , T , N , ALPHA , BETA , SEED , OUTPUT );
    
    fprintf( 'Continue to run sampler to collect samples\n' );
    for s=1:NSAMPLES
        filename = sprintf( 'swrsamples_%d\\lda_chain%d_sample%d' , dataset , c , s );
        fprintf( 'Saving sample #%d from chain #%d: filename=%s\n' , s , c , filename );
        comm = sprintf( 'save ''%s'' WP DP Z ALPHA BETA SEED N Z T s c' , filename );
        eval( comm );
                
        if (s < NSAMPLES)
           N = LAG;
           SEED = SEED + 1; % important -- change the seed between samples !!
           [ WP,DP,Z ] = GibbsSamplerLDA( WS , DS , T , N , ALPHA , BETA , SEED , OUTPUT , Z );
        end
    end
end


