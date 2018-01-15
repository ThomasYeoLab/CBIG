%% Example 1 of Running Collocation Topic Model (LDACOL)

%%
% Load a dataset
dataset = 1; % 1 = psych review; 2 = nips

if (dataset == 1)
    fprintf( 'Loading Psych Review Abstracts - Collocation Data\n' );
    load 'psychreviewcollocation'; % load in variables: WW WO DS WS SI
    filenm = 'topics_psychreview_col.txt';
elseif (dataset == 2 )
    fprintf( 'Loading NIPS papers - Collocation Data\n' );
    load 'nipscollocation'; % load in variables: WW WO DS WS SI
    filenm = 'topics_nips_col.txt';
end

%%
% The number of topics
T = 100;

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

MAXC   = 4;   % maximum collocation length (in post-processing topics)

%%
% The number of iterations
N = 100;

%% 
% The seed number
SEED = 1;

%%
% This function might need a few minutes to finish
tic
[ WP,DP,WC,C,Z ] = GibbsSamplerLDACOL( WS , DS , SI , WW , T , N , ALPHA , BETA , GAMMA0, GAMMA1 , DELTA , SEED , OUTPUT  );
toc

% convert topics to include collocations
[ WPNEW , DPNEW , WONEW ] = CreateCollocationTopics( C , Z , WO , DS , WS , T , MAXC );

fprintf( 'Writing collocation topics to file: %s\n' , filenm );

%%
% Post-process the vocabulary to include collocations as separate entries. Recalculate word-topic distributions with expanded vocabulary 
if (dataset == 1)
    S = WriteTopics( WPNEW , BETA , WONEW , 20 , 0.7 , 4 , filenm );
    save 'ldacol_psychreview' WPNEW DPNEW WONEW WP DP WC C ALPHA BETA GAMMA0 GAMMA1 DELTA SEED N Z T;
elseif (dataset ==2)
    S = WriteTopics( WPNEW , BETA , WONEW , 40 , 0.7 , 4 , filenm );
    save 'ldacol_nips' WPNEW DPNEW WONEW WP DP WC C ALPHA BETA GAMMA0 GAMMA1 DELTA SEED N Z T;
end

%%
% Show some topics
S{1}( 1:90 )
S{2}( 1:90 )
S{3}( 1:90 )
S{4}( 1:90 )
S{5}( 1:90 )
S{6}( 1:90 )
S{7}( 1:90 )
S{8}( 1:90 )
S{9}( 1:90 )
S{10}( 1:90 )



