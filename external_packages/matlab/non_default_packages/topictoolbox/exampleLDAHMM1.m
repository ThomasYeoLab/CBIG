%% Example 1 of running HMM-LDA topic model
%
% This example shows how to run the HMM-LDA Gibbs sampler on a small
% dataset to extract a set of topics and a set of syntactic states. Unlike
% the LDA and AT topic models, there is no need to exclude the stop words
% from a corpus of text. Also, this model differs from the LDA and AT topic
% models by utilizing the word order information of the text. The output of
% this code is a count matrix |WP| with the number of times words are 
% assigned to topics and a count matrix |MP| that contains the number
% of times words are assigned to syntactic states.

%% 
% Choose the dataset
dataset = 1; % 1 = psych review; 2 = nips papers

if (dataset == 1)
    % Load the psych review word stream
    load 'psychreviewstream';

    % Set the parameters for the model
    T      = 50;     % number of topics
    NS     = 12;     % number of syntactic states
    N      = 200;    % number of iterations
    ALPHA  = 50 / T; % ALPHA hyperparameter
    BETA   = 0.01;   % BETA hyperparameter
    GAMMA  = 0.1;    % GAMMA hyperparameter
    SEED    = 2;     % random SEED
    
    filename1 = 'topics_psychreview_hmmlda.txt'; % text file showing topic-word distributions
    filename2 = 'states_psychreview_hmmlda.txt'; % text file showing hmm state-word distributions
end

if (dataset == 2)
    % Load the nips paper word stream
    load 'nips_stream';

    % Set the parameters for the model
    T      = 50;     % number of topics
    NS     = 16;     % number of syntactic states
    N      = 400;    % number of iterations
    ALPHA  = 50 / T; % ALPHA hyperparameter
    BETA   = 0.01;   % BETA hyperparameter
    GAMMA  = 0.1;    % GAMMA hyperparameter
    SEED    = 2;     % random SEED
    
    filename1 = 'topics_nips_hmmlda_2.txt'; % text file showing topic-word distributions
    filename2 = 'states_nips_hmmlda_2.txt'; % text file showing hmm state-word distributions
end

%%
% What output to show (0=no output; 1=iterations; 2=all output)
OUTPUT = 1;

%%
% Run the HMM-LDA Gibbs sampler
tic
[WP,DP,MP,Z,X]=GibbsSamplerHMMLDA( WS,DS,T,NS,N,ALPHA,BETA,GAMMA,SEED,OUTPUT);

fprintf( 'Elapsed time = %5.0f seconds\n' , toc );

%%
% Save the results to a file
if (dataset==1)
   save 'hmmldasingle_psychreview' WP DP MP Z X ALPHA BETA GAMMA N;
elseif (dataset==2)
   save 'hmmldasingle_nips' WP DP MP Z X ALPHA BETA GAMMA N; 
end

%%
% Calculate the most likely words in each topic and write to a cell array
% of strings
[S] = WriteTopics( WP , BETA , WO , 7 , 0.8 , 4 , filename1 );


%%
% Show the most likely words in the topics
fprintf( '\n\nMost likely words in the topics:\n' );
S( 1:T )  

%%
% Calculate the most likely words in each syntactic state and write to a
% cell array of strings
[S] = WriteTopics( MP , BETA , WO , 7 , 0.8 , 4 , filename2 );

%%
% Show the most likely words in the syntactic states
fprintf( '\n\nMost likely words in the syntactic states:\n' );
S( 1:NS ) 
