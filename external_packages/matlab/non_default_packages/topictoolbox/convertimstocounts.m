function [ WS , DS , n ] = convertimstocounts( ims )

% count the total number of pixels activated
n = sum( sum( ims > 0 ));

w  = zeros( n , 1 );
d  = zeros( n , 1 );
c  = zeros( n , 1 );

nims = size( ims , 1 );

ii = 0;
for di=1:nims
   imnow = ims( di , : );
   
   % find activated pixels
   wh = find( imnow > 0 );
   nz = length( wh );
   
   % what are the intensities at those pixels?
   counts = imnow( wh );
   
   w( ii+1:ii+nz ) = wh;
   d( ii+1:ii+nz ) = di;
   c( ii+1:ii+nz ) = counts;
   
   ii = ii + nz;
end

WD = sparse( w , d , c );

[ WS , DS ] = SparseMatrixtoCounts( WD );

%maxcount = max( c );
%fprintf( 'Maximum count for a token is %d\n' , maxcount );