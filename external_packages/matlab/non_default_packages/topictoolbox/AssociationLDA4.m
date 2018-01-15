function [ probs ] = AssociationLDA4( WP , BETA , W1 , pz )
%% Function WriteTopics
%
% works on single samples as well as multiple samples
% TAKE THE PRIOR TOPIC PROBABILITY INTO ACCOUNT
if iscell( WP )
    % we are dealing with multiple Gibbs samples. Need to average
    % predictions over these samples. Note that the WP samples do not have
    % to contain the same number of topics
    
    NS = length( WP ); % number of samples
    NW = size( WP{ 1 } , 1 ); % number of words
    NW1 = length( W1 ); % number of W1 cues
    probs = zeros( NW , NW1 );
    
    for s=1:NS
        T  = size( WP{ s } , 2 ); % number of topics
        
        %pz = full( sum( WP{ s } , 1 ));
        %pz = pz / sum( pz );

        pwz = WP{ s } + BETA;
        for i=1:T
            pwz( : , i ) = pwz( : , i ) / sum( pwz( : , i ));
        end

        pzw = zeros( NW , T );
        for i=1:NW
            tprobs = pwz( i , : ) .* pz{ s }( : )';
            pzw( i , : ) = tprobs / sum( tprobs );
        end

        probs = probs + pwz * pzw( W1 , : )';
    end
    
    probs = probs / NS;
else
    % we are dealing with just one Gibbs sample
    %pz = full( sum( WP , 1 ));
    %pz = pz / sum( pz );
    
    NW = size( WP , 1 ); % number of words
    T  = size( WP , 2 ); % number of topics
    NW1 = length( W1 ); % number of W1 cues

    pwz = WP + BETA;
    for i=1:T
        pwz( : , i ) = pwz( : , i ) / sum( pwz( : , i ));
    end
    
    pzw = zeros( NW , T );
    for i=1:NW
        tprobs = pwz( i , : ) .* pz( : )';
        pzw( i , : ) = tprobs / sum( tprobs );
    end
    
    probs = pwz * pzw( W1 , : )'; 
end
