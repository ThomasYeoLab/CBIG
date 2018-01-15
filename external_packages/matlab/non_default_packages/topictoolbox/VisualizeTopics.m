function VisualizeTopics( DP , ALPHA , S , VIZMODE );
%% Function VisualizeTopics
%%
% Visualizes topics in a 2D map
%%
% |visualizetopics( DP,ALPHA,S,VIZMODE )| takes the topic strings in |S| and
% the document-topic count matrix |DP| where |DP(i,j)| contains the number of
% times a word in document |i| has been assigned to topic |j|. These counts are
% transformed to probability distributions over documents for each topic.
% For each topic pair, the symmetrized Kullback Leibler distance between
% document distributions is calculated. Classical multidimensional scaling
% is used to visualize all pairwise topic distances in a 2D map.
%
%%
% Notes:
%%
% The variable |VIZMODE| can be set to 'horizontal' or 'vertical' which
% determines the way the topics are displayed.
% 
%%
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

[ Coords,e ] = cmdscale( topicdist );

figure;
clf;

if (strcmp( VIZMODE , 'horizontal' )) 
    str2 = S;
else
    for i=1:T
        str = S{ i };
        str = [ ' ' str ' ' ];
        nspace = sum( str == ' ' );
        newstr = cell( 1,nspace-1 );
        whspace = find( str == ' ' );
        for k=1:nspace-1
            newstr{ k } = str( whspace( k )+1:whspace( k+1)-1 );
        end

        str2{ i } = newstr;
    end
end

% show the topics in order of probability. Most likely topics are drawn
% last
[ temp , index ] = sort( sumdp );

for i=1:T
    wh = index( i );
    backcol = 'b';
    text( Coords(wh,1) , Coords(wh,2) , str2{ wh } , 'backgroundcolor' , 'w' , 'EdgeColor' , backcol , ...
          'FontSize' , 9 , 'HorizontalAlignment' , 'center' , 'LineWidth' , 1.5 );
end
xlim( [ min( Coords(:,1)) max( Coords(:,1)) ] );
ylim( [ min( Coords(:,2)) max( Coords(:,2)) ] ); 
axis off;
