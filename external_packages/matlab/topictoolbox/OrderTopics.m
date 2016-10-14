function [ Order ] = OrderTopics( DP , ALPHA );
%% Function OrderTopics
%%
% Orders topics according to similarity in topic distributions
%%
% In this procedure, for each topic, the probability distribution over
% documents is calculated. For each topic pair (i,j), the (symmetrized) KL
% distance between topic distributions i and j is calculated. The KL
% distances are then subjected multidimensional scaling (MDS). The first
% dimension of the MDS solution is then used to get an ordering for topics.
% 
%%
% Notes:
% The variable |ALPHA| is the hyperparameter on the topic distributions needed to
% calculate the topic-document distributions

% calculate prob distribution over documents for each topic
D = size( DP , 1 );
T = size( DP , 2 );

sumdp = full( sum( DP , 1 )) + ALPHA * D;

topicdist = zeros( T , T );
for i1=1:T
    dp1 = ( DP( : , i1 ) + ALPHA ) / sumdp( i1 );
    
    for i2=i1+1:T
        dp2 = ( DP( : , i2 ) + ALPHA ) / sumdp( i2 );
        
        % calculate KL distances both ways       
        KL1 = sum( dp1 .* log2( dp1 ./ dp2 ));    
        KL2 = sum( dp2 .* log2( dp2 ./ dp1 )); 
        
        topicdist( i1,i2 ) = (KL1+KL2)/2;
        topicdist( i2,i1 ) = (KL1+KL2)/2;
    end
end

Coords = mdscale(topicdist,1);
%[ Coords,e ] = cmdscale( topicdist );

[ temp,Order]=sort( Coords(:,1) );

