function [ imcollage , labels ] = createcollage( ims , nx , ny , nimsx , nimsy , gapx , gapy )

% ims is a K x W array, where ims( k,w ) has the intensity of pixel w of
% image k

nims = size( ims , 1 );

if (nims > nimsx*nimsy)
    warning( 'Not all images fit in this collage' );
end

nxtotal = nimsx * (nx + gapx) - gapx;
nytotal = nimsy * (ny + gapy) - gapy;
imcollage = zeros( nytotal , nxtotal ) * NaN;
ii = 0;
for x=1:nimsx
    for y=1:nimsy
        starty = 1+(y-1)*(ny+gapy);
        startx = 1+(x-1)*(nx+gapx);     
        
        if (ii < nims)
            ii = ii + 1;
            imnow = reshape( ims( ii , : ) , ny , nx );
            imcollage( starty:starty+ny-1 , startx:startx+nx-1 ) = imnow;
            labels( ii , 1 ) = x;
            labels( ii , 2 ) = y;
        end
    end
end
