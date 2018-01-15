%% Example 2 of Running Collocation Topic Model (LDACOL)

%%
% How to save multiple samples from the same or different chains in the
% LDACOL model

%%
% Load a dataset
dataset = 1; % 1 = psych review; 2 = nips

if (dataset == 1)
    fprintf( 'Loading Psych Review Abstracts - Collocation Data\n' );
    load 'psychreviewcollocation'; % load in variables: WW WO DS WS SI
elseif (dataset == 2 )
    fprintf( 'Loading NIPS papers - Collocation Data\n' );
    load 'nipscollocation'; % load in variables: WW WO DS WS SI
end

%%
% The number of topics
T = 50;

%%
% What output to show (0=no output; 1=iterations; 2=all output)
OUTPUT = 1;

%%
% Set the hyperparameters of the model
BETA   = 0.01;
ALPHA  = 50/T;
GAMMA0 = 0.1;
GAMMA1 = 0.1;
DELTA  = 0.1;

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
    [ WP,DP,WC,C,Z ] = GibbsSamplerLDACOL( WS , DS , SI , WW , T , N , ALPHA , BETA , GAMMA0, GAMMA1 , DELTA , SEED , OUTPUT  );

    fprintf( 'Continue to run sampler to collect samples\n' );
    for s=1:NSAMPLES
        filename = sprintf( 'ldacol_chain%d_sample%d' , c , s );
        fprintf( 'Saving sample #%d from chain #%d: filename=%s\n' , s , c , filename );
        comm = sprintf( 'save ''%s'' WP DP WC C ALPHA BETA GAMMA0 GAMMA1 DELTA SEED N Z T s c' , filename );
        eval( comm );
        
        if (s < NSAMPLES)
           N = LAG;
           SEED = SEED + 1; % important -- change the seed between samples !!
           [ WP,DP,WC,C,Z ] = GibbsSamplerLDACOL( WS , DS , SI , WW , T , N , ALPHA , BETA , GAMMA0, GAMMA1 , DELTA , SEED , OUTPUT , C , Z );
        end
    end
end



    




