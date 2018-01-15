function [ SM ] = WriteTopicsMult( WPM , BETAM , WOM , K , E , M , FILENAME )
%% Function WriteTopicsMult
%%
% |[ SM ] = WriteTopicsMult( WPM , BETAM , WOM )| writes the most likely entities
% (e.g. words, authors) per topic to a cell array of strings |SM|. The
% difference with function WriteTopics is that for each topics, the
% distributions over multiple entities are calculated 
%
% Note that |WPM| is a cell array of matrices and |BETAM| is an array of
% hyperparameter constants
%
% Each entry of |WPM| is a sparse |W| x |T| count matrix such that
% |WPM{k}(i,j)| contains the number of times entity |i| of type |k| has been assigned to
% topic |j|. |BETA(k)| is the relevant hyperparameter used in the Gibbs sampling
% routine. |WOM| is a cell structure of strings such that |WOM{k}{i}| contains the ith 
% entity string of type |k|.
%
% |[ SM ] = WriteTopics( WPM , BETAM , WOM , K , E , M , FILENAME )|
% writes the |K| most likely entities per topic to a text file |FILENAME|
% with |M| columns. |E| is a threshold on the topic listings in |SM|. Only
% entities that do not exceed this cumulative probability are listed.
%  
% |WriteTopicsMult( WPM , BETAM , WOM , K , E , M , FILENAME )| writes the |K| most
% likely entities per topic without producing a cell array of strings 
%
% |[ SM ] = WriteTopicsMult( WPM , BETAM , WOM  )| uses default values
% |K=5| and |E=1|
% 
%%
% Example
%%
%     |load 'atsingle_nips'|
%     |load 'words_nips'|
%     |load 'authors_nips'|
%     |WPM{1} = WP; WPM{2} = AT;| 
%     |BETAM(1)=BETA; BETAM(2) = ALPHA;|
%     |WOM{1}=WO; WOM{2}=AN;|
%     |[ SM ] = WriteTopicsMult( WPM , BETAM , WOM );|
%     |SM{1}(1:10)|
%     |SM{2}(1:10)| 
%
%%
% will list the most likely words and authors for the first 10 topics 
% 
%%
%     |WriteTopicsMult( WPM , BETAM , WOM , 10 , 1.0 , 4 ,
%     'topics_nips.txt' )|
%
% will write 10 most likely words and authors per topic to a four column text file

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
  
NK = length( WPM );

if NK ~= length( WOM )
   error( 'WPM and WOM must have equal number of entity types' ); 
end

if NK ~= length( BETAM )
   error( 'WPM and BETAM must have equal number of entity types' ); 
end

for k=1:NK
    WP = WPM{k};
    BETA = BETAM( k );
    
    W = size( WP , 1 );
    T = size( WP , 2 );
  
    WM( k ) = W;
    
    sumWP = sum( WP , 1 ) + BETA*W;
    probtopic = sumWP / sum( sumWP );

    Sorted_P_w_z{k} = zeros( K , T );
    Index_P_w_z{k} = zeros( K , T );

    for t=1:T
        [ temp1 , temp2 ] = sort( -WP( : , t ) );
        Sorted_P_w_z{k}( : , t )  = ( full( -temp1( 1:K )) + BETA ) ./ ( repmat( sumWP( t ) , K , 1 ));
        Index_P_w_z{k}( : , t )   = temp2( 1:K );
    end
end

if (nargout==1)
    for k=1:NK
        W = WM( k );
        WO = WOM{ k };
        
        S = cell( T,1 );
        for t=1:T
            cumprob = 0;
            for r=1:K
                index = Index_P_w_z{k}( r , t );
                prob  = Sorted_P_w_z{k}( r , t );
                cumprob = cumprob + prob;
                word = WO{ index };

                if (cumprob < E)
                    if (r==1)
                        str = word;
                    else
                        str = [ str ' ' word ];
                    end
                end
            end
            S{ t } = str;
        end
        
        SM{ k } = S;
    end
end

if (dofile==1)
    startt = 1;
    endt   = M;
    fid = fopen( FILENAME , 'W' );

    while startt < T
        for c=startt:endt
            tempstr = sprintf( '%s %d' , 'TOPIC' , c );
            %fprintf( fid , '%15s\t%6.5f\t' , tempstr , probtopic( c ) );
            fprintf( fid , '%35s\t%6.5f\t' , tempstr , probtopic( c ) );
        end
        fprintf( fid , '\r\n\r\n' );

        for k=1:NK
            for r=1:K
                for c=startt:endt
                    index = Index_P_w_z{k}( r , c );
                    prob = Sorted_P_w_z{k}( r , c );

                    %fprintf( fid , '%15s\t%6.5f\t' , WOM{k}{ index } , prob );
                    fprintf( fid , '%35s\t%6.5f\t' , WOM{k}{ index } , prob );
                end
                fprintf( fid , '\r\n' );
            end
            fprintf( fid , '\r\n' );
        end
        fprintf( fid , '\r\n' );

        startt = endt + 1;
        endt   = min( T , startt + M - 1 );
    end
    
    fclose( fid );
end

