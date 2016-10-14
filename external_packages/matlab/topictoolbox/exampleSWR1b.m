%% Precompute distributions for special word retrieval model (SWR)
%

%%
% Choose the dataset
dataset = 1; % 1 = psych review abstracts 2 = NIPS papers

if (dataset == 1)
    % Set the hyperparameters 
    BETA1 = 0.000001;
    BETA2 = 0.000001;
        
    % load the psych review data in bag of words format   
    load 'psychreviewbagofwordsfromstream'; % WO DS WS ORIGWSPOS;

    NSAMPLES = 10; % number of topic samples
elseif (dataset == 2)
    % Set the hyperparameters 
    BETA1 = 0.000001;
    BETA2 = 0.000001;
    
    % load the psych review data in bag of words format   
    load 'nipsbagofwordsfromstream'; % WO DS WS ORIGWSPOS;

    NSAMPLES = 6; % number of topic samples
end

% load in word-topic probabilities from different samples
for ss=1:NSAMPLES
    filename = sprintf( 'swrsamples_%d\\lda_chain%d_sample%d' , dataset , ss , 1 );
    comm = sprintf( 'load %s' , filename );
    eval( comm );
    
    fprintf( 'Loading sample: %s\n' , filename );
    
    W = size( WP , 1 );
    T = size( WP , 2 );
    PWZ = full( WP + BETA );
    for j=1:T
        PWZ( : , j ) = PWZ( : , j ) / sum( PWZ( : , j ));
    end

    % create document-topic probabilities
    PTD = full( DP + ALPHA );
    D = size( DP , 1 );
    for d=1:D
        PTD( d , : ) = PTD( d , : ) / sum( PTD( d , : ));
    end
    
    % create the NQ x D word probabilities according to topic model
    if (ss==1)
        PQ0_ALLW = PWZ * PTD'; 
    else
        PQ0_ALLW = PQ0_ALLW + PWZ * PTD';  
    end
end

PQ0_ALLW = PQ0_ALLW / NSAMPLES;

% create document word sparse count matrix
WD = sparse( WS , DS , ones( size( DS )));

% create the NQ x D word probabilities according to document model

%PQ1_ALLW = log( full( WD )+1 ) + BETA1; % !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

PQ1_ALLW = full( WD + BETA1 );


for d=1:D
    PQ1_ALLW( : , d ) = PQ1_ALLW( : , d ) / sum( PQ1_ALLW( : , d ));
end

PQ2_ALLW = full( sum( WD , 2 ))' + BETA2;
PQ2_ALLW = PQ2_ALLW / sum( PQ2_ALLW );

comm = sprintf( 'save swrsamples_%d\\swrdata PQ0_ALLW PQ1_ALLW PQ2_ALLW' , dataset );
eval( comm );



