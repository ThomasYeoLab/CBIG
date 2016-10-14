function [ sims , index , csims ] = AssociationTFIDF( WD , Q )

NW = size( WD , 1 ); % number of words
ND = size( WD , 2 ); % number of docs
NQ = length( Q ); % number of queries

%---------------------------------
%  Calculate TF-IDF weight
% --------------------------------
df = full( sum( WD>0 , 1 ));
idf = log2( NW ./ df );

%----------------------------------
%  Convert counts to TF-IDF weights
% ---------------------------------
for i=1:ND
   WD( : , i ) = WD( : , i ) * idf( i );  
end

% calculate squared terms for cosine calculation
squaredterms = full( sum( WD.^2 , 2 ));

csims = zeros( NW , NQ );
for i=1:NQ
    wh = Q( i );
    queryvector = full( WD( wh , : )' );
    innerprods = WD * queryvector;
    normfactor = sqrt( squaredterms( wh ) * squaredterms( : ) );
    
    cosvalues = innerprods ./ normfactor;
    cosvalues( wh ) = min( cosvalues );
    cosvalues = cosvalues / sum( cosvalues ); % normalize to probabilities
    
    csims( : , i ) = cosvalues;
end

%---------------------------------
%  Sort the outcome
% --------------------------------
[ sims , index ] = sort( csims , 1 , 'descend' );
