function [ S ] = WriteTopics( WP , BETA , WO , K , E , M , FILENAME )
%% Function WriteTopics
%%
% |[ S ] = WriteTopics( WP , BETA , WO )| writes the most likely entities
% (e.g. words, authors) per topic to a cell array of strings. |WP| is a
% sparse |W| x |T| count matrix such that |WP(i,j)| contains the number of times
% entity |i| has been assigned to topic |j| (|W|=number of entities, |T|=number of
% topics). |BETA| is the relevant hyperparameter used in the Gibbs sampling
% routine. |WO| is a cell structure of strings such that |WO{i}| contains the ith 
% entity string (e.g., word or author name)
%
% |[ S ] = WriteTopics( WP , BETA , WO , K , E , M , FILENAME )|
% writes the |K| most likely entities per topic to a text file |FILENAME|  
% with |M| columns. |E| is a threshold on the topic listings in |S|. Only
% entities that do not exceed this cumulative probability are listed.
% 
%%
% |WriteTopics( WP , BETA , WO , K , E , M , FILENAME )| writes the |K| 
% most likely entities per topic without producing a cell array of strings
% 
%
%%
% |[ S ] = WriteTopics( WP , BETA , WO  )| uses default values |K=5| and
% |E=1|
% 
%%
% Example
%
%%
%     |load 'ldasingle_psychreview'|
%     |[ S ] = WriteTopics( WP , BETA , WO )|
%
%%
% will list the most likely words per topic in the saved |WP| sample 
%
%%
%     |WriteTopics( WP , BETA , WO , 10 , 1.0 , 4 ,
%     'topics50_psychreview.txt' )|
%%
% will write 10 most likely words per topic to a four column text file

if (nargin>7)
    error( 'Too many input arguments' );
elseif (nargin<3)
    error( 'Too few input arguments' );
end

if (nargin==3)
    K = 5;
    E   = 1.0;
    dofile = 0;
end

if (nargin==4)
    E  = 1.0;
    dofile = 0;
end

if (nargin==5)
    dofile = 0;
end

if (nargin==7)
    dofile = 1;
end

if (nargout>1)
    error( 'Too many output arguments' );
end
    
W = size( WP , 1 );
T = size( WP , 2 );
K = min( [ W K ] );

sumWP = sum( WP , 1 ) + BETA*W;
probtopic = sumWP / sum( sumWP ); 

Sorted_P_w_z = zeros( K , T );
Index_P_w_z = zeros( K , T );

for t=1:T
   [ temp1 , temp2 ] = sort( -WP( : , t ) );
   Sorted_P_w_z( : , t )  = ( full( -temp1( 1:K )) + BETA ) ./ ( repmat( sumWP( t ) , K , 1 ));
   Index_P_w_z( : , t )   = temp2( 1:K );
end

if (nargout==1)
    S = cell( T,1 );
    for t=1:T
        cumprob = 0;
        for r=1:K
            index = Index_P_w_z( r , t );
            prob  = Sorted_P_w_z( r , t );
            cumprob = cumprob + prob;
            word = WO{ index };

            if (cumprob < E) || (r==1)
                if (r==1)
                    str = word;
                else
                    str = [ str ' ' word ];
                end
            end
        end
        S{ t } = str;
    end
end

% !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
largewidth = 1;

if (dofile==1)
    M = min( [ T M ] );
    
    startt = 1;
    endt   = M;
    fid = fopen( FILENAME , 'W' );
    while startt < T
        for c=startt:endt
            tempstr = sprintf( '%s_%d' , 'TOPIC' , c );
            
            if (largewidth==1) 
                fprintf( fid , '%25s\t%6.5f\t' , tempstr , probtopic( c ) ); else
                fprintf( fid , '%15s\t%6.5f\t' , tempstr , probtopic( c ) ); 
            end
        end
        fprintf( fid , '\r\n\r\n' );

        for r=1:K
            for c=startt:endt
                index = Index_P_w_z( r , c );
                prob = Sorted_P_w_z( r , c );

                if (largewidth==1)
                    fprintf( fid , '%25s\t%6.5f\t' , WO{ index } , prob ); else
                    fprintf( fid , '%15s\t%6.5f\t' , WO{ index } , prob );
                end
            end
            fprintf( fid , '\r\n' );
        end
        fprintf( fid , '\r\n\r\n' );

        startt = endt + 1;
        endt   = min( T , startt + M - 1 );
    end
    fclose( fid );
end
