%% Example 3 of Running Collocation Topic Model (LDACOL)
%%

%%
% How to convert stream data into collocation stream data

%%
% Choose a dataset
dataset = 1; % 1 = psych review; 2 = nips

if (dataset==1)
    % example 1: load in PSYCH REVIEW word stream and convert to
    % collocation data
    load 'psychreviewstream';
    [ WW , WO , WS , DS , SI]=stream_to_collocation_data(DS,WS,WO,'stopwordlist.txt');
    save 'psychreviewcollocation' WW WO DS WS SI;
end

if (dataset==2)
    % example 2: load in NIPS word stream and convert to collocation data
    load 'nips_stream';
    [ WW , WO , WS , DS , SI]=stream_to_collocation_data(DS,WS,WO,'stopwordlist.txt');
    save 'nipscollocation' WW WO DS WS SI;
end