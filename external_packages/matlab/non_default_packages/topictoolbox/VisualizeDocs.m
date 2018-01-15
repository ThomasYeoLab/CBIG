function VisualizeDocs( DP , ALPHA , S , NCHARS , MAXLINES );
%% Function VisualizeDocs
%%
% Visualize documents in a 2D map
%%
% |visualizedocuments( DP,ALPHA,S,NCHARS,MAXLINES )| takes the document
% strings in |S| and the document-topic count matrix |DP| where |DP(i,j)|
% contains the number of times a word in document |i| has been assigned to
% topic |j|. These counts are transformed to probability distributions over
% topics for each document. For each document pair, the symmetrized
% Kullback Leibler distance between topic distributions is calculated.
% Classical multidimensional scaling is used to visualize all pairwise
% document distances in a 2D map. 
%
% Example
%%
% look at exampleVIZ2.m for an example script utilizing this function
%
% Notes
%%
% The variable |NCHARS| limits the number of characters on one line within
% a text box 
% 
% The variable |MAXLINES| limited the number of lines within a text box. 

% calculate prob distribution over documents for each topic
D = size( DP , 1 );
T = size( DP , 2 );

sumdp = full( sum( DP , 2 )) + ALPHA * T;

docdist = zeros( D , D );
for i1=1:D
    dp1 = ( DP( i1, : ) + ALPHA ) / sumdp( i1 );
    
    for i2=i1+1:D
        dp2 = ( DP( i2 , : ) + ALPHA ) / sumdp( i2 );
        
        % calculate KL distances both ways       
        KL1 = sum( dp1 .* log2( dp1 ./ dp2 ));    
        KL2 = sum( dp2 .* log2( dp2 ./ dp1 )); 
        
        docdist( i1,i2 ) = (KL1+KL2)/2;
        docdist( i2,i1 ) = (KL1+KL2)/2;
    end
end

[ Coords,e ] = cmdscale( docdist );

figure;
clf;

if NCHARS == 0
    str2 = S;
else
    for i=1:D
        str = strtrim( S{ i } );
        
        str = [ ' ' str ' ' ];
        whspace = find( str == ' ' );
        nspaces = length( whspace );
        
        nlines = 1;
        nchars = 0;
        newstr{ nlines } = '';
        for k=1:nspaces-1
           whstart = whspace( k   )+1;
           whend   = whspace( k+1 )-1;
           word    = str( whstart:whend );
           wordlen = length( word );
           
           
           if nchars > NCHARS
               nchars = 0;
               nlines = nlines + 1;
           end
           
           if (nchars==0)
               newstr{ nlines } = word;
           else
               newstr{ nlines } = [ newstr{ nlines } ' ' word ];
           end  
           
           nchars = nchars + wordlen + 1;
        end
       
        str2{ i } = newstr( 1:min( MAXLINES,nlines ));
    end
end

for i=1:D
    wh = i;
    backcol = 'b';
    text( Coords(wh,1) , Coords(wh,2) , str2{ wh } , 'backgroundcolor' , 'w' , 'EdgeColor' , backcol , ...
          'FontSize' , 9 , 'HorizontalAlignment' , 'center' , 'LineWidth' , 1.5 );
end
xlim( [ min( Coords(:,1)) max( Coords(:,1)) ] );
ylim( [ min( Coords(:,2)) max( Coords(:,2)) ] ); 
axis off;
