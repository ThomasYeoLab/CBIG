function [ WW , WO , WS , DS , SI , ORIGWSPOS ] = stream_to_collocation_data( DS_IN , WS_IN , WO_IN , stopwordsfile , freqcutoff , doexcludenumbers )
%% Function stream_to_collocation_data

%%
% A utility to convert the stream data used by HMM-LDA model into stream
% data useful for the LDA-COL model.

%%
% |[ WW , WO , WS , DS , SI ] = stream_to_collocation_data( DS_IN , WS_IN ,
% WO_IN , stopwordsfile )|
%
% creates the new stream of words data. The input document and word indices
% are provided by |DS_IN| and |WS_IN| respectively. The vocabulary is given
% in the cell array |WO_IN|. The string |stopwordsfile| contains the name
% of the file with stop words. The output |WW| is a matrix where
% |WW(w2,w1)| gives frequency with which |w2| follows |w1|. The status
% vector |SI| has values where |SI(k)| = 1 for positions |k| where previous
% word |k|-1 can be considered as part of a collocation with current word.
% For |SI(|k|)| = 0, the previous word |k|-1 did not precede word k in the
% original text (because there were stop words in between) or because the
% word ended a previous document. 

if (nargin<=4)
    freqcutoff = 0;
end

if (nargin<=5)
    doexcludenumbers = 0;
end

n  = length( WS_IN );
W  = length( WO_IN ) + 1;
D  = max( DS_IN );

if (nargout > 5)
   ORIGWSPOS = 1:n; 
end

F = hist( WS_IN , 1:(W-1) ); % calculate the frequencies of all words

w1 = WS_IN( 1:n-1 ) + 1; % add one to get sentence marker on index 1
w2 = WS_IN( 2:n   ) + 1;

% find document boundaries
DB = diff( DS_IN );
notboundary = find( DB == 0 );
w1 = w1( notboundary );
w2 = w2( notboundary );

WW = sparse( w2 , w1 , ones( size( w1 )) , W , W );

% remove the sentence marker from these matrices
WW( 1,: ) = [];
WW( :,1 ) = [];

% load in stopwords file
[ stopwords ] = lower( textread( stopwordsfile , '%s' ));

WO_IN = lower( WO_IN );

% find word indices that are not stop words
[ temp , ind1 ] = setdiff( WO_IN , stopwords );

% find word indices for words above frequency cutoff
ind2 = find( F >= freqcutoff );

% now take the intersection
ind = intersect( ind1 , ind2 );

% do we remove any strings with numbers?
if (doexcludenumbers==1)
    nind = length( ind );
    oklist = ones( 1,nind );
    for i=1:nind
        word = WO_IN{ ind( i )};
        if any( isstrprop( word , 'digit' ))
            oklist( i ) = 0;
        end
    end
    ind = ind( find( oklist ));
end

WO = WO_IN( ind );

WW = WW( ind,ind );

[ SI , WS ] = ismember( WS_IN , ind );
whok = find( SI );

% mark the locations where collocation cannot be formed
SI( 2:end ) = SI( 2:end ) .* SI( 1:end-1); 

% also exclude document boundaries
DB = diff( DS_IN );
whboundary = find( DB == 1 );
SI( whboundary+1 ) = 0;


SI = double( SI( whok ));

WS = WS( whok );
DS = DS_IN( whok );

if (nargout > 5)
   ORIGWSPOS = ORIGWSPOS( whok );
end


