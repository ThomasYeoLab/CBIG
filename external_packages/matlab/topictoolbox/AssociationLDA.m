function [ probs_sorted , index , probs ] = AssociationLDA( WP , BETA , W1 )
%% Function WriteTopics
%
% works on single samples as well as multiple samples
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
        
        sumt  = full( sum( WP{ s } , 1 )) + NW * BETA;
        beta2 = BETA ./ sumt;

        WP2   = WP{ s }; % note WP2 stays sparse
        for z=1:T
            WP2( : , z ) = WP2( : , z ) / sumt( z );
        end

        % take transpose of WP2 for faster processing
        WP3 = WP2( W1 , : )';

        % this should also be sparse
        crossprods = WP2 * WP2( W1 , : )';

        sum4 = sum( beta2 .* beta2 );
        sumbeta = sum( beta2 );
        for i=1:NW1
            sum1 = crossprods( : , i ) ;
            sum2 = sum( WP3( : , i )' .* beta2 );
            sum3 = WP2 * beta2';

            sum5 = sum( WP3( : , i )) + sumbeta;

            sumt = ( sum1 + sum2 + sum3 + sum4 ) / sum5;
            
            sumt( W1( i )) = 0;
            sumt = sumt / sum( sumt );
            
            probs( : , i ) = probs( : , i ) + sumt;
        end
    end
    
    probs = probs / NS;
else
    % we are dealing with just one Gibbs sample
    
    NW = size( WP , 1 ); % number of words
    T  = size( WP , 2 ); % number of topics
    NW1 = length( W1 ); % number of W1 cues

    sumt  = full( sum( WP , 1 )) + NW * BETA;
    beta2 = BETA ./ sumt;

    WP2   = WP; % note WP2 stays sparse
    for z=1:T
        WP2( : , z ) = WP2( : , z ) / sumt( z );
    end

    % take transpose of WP2 for faster processing
    WP3 = WP2( W1 , : )';

    % this should also be sparse
    crossprods = WP2 * WP2( W1 , : )';

    probs = zeros( NW , NW1 );
    sum4 = sum( beta2 .* beta2 );
    sumbeta = sum( beta2 );
    for i=1:NW1
        sum1 = crossprods( : , i ) ;
        sum2 = sum( WP3( : , i )' .* beta2 );
        sum3 = WP2 * beta2';

        sum5 = sum( WP3( : , i )) + sumbeta;

        sumt = ( sum1 + sum2 + sum3 + sum4 ) / sum5;
        
        sumt( W1( i )) = 0;
        sumt = sumt / sum( sumt );
            
        probs( : , i ) = sumt;
    end
end

[ probs_sorted , index ] = sort( probs , 1 , 'descend' );

% TEST CODE
% pwz = WP + BETA;
% for i=1:T
%    pwz( : , i ) = pwz( : , i ) / sum( pwz( : , i )); 
% end
% probs2 = zeros( NW , NW1 );
% for i=1:NW1
%    w1 = W1( i );
%    p = zeros( NW , 1 );
%    for z=1:T
%       p = p + pwz( : , z ) * pwz( w1 , z ) / sum( pwz( w1 , : ));
%    end
%     
%    probs2( : , i ) =  p;
% end
% 
