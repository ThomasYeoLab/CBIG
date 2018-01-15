%% Example 2 of applying topic models to images
% This example illustrates how to apply topic models to a set of handwritten
% digits and characters. The dataset is based on Sam Roweis's "binaryalphadigs.mat" at  
% http://www.cs.toronto.edu/~roweis/data.html
% It contains binary 20x16 digits of "0" through "9" and capital "A"
% through "Z" where there are 39 examples of each class. 
%
% The images are represented by matrix |ims| where |ims( i,j )| contains the
% the value of binary pixel |j| of image |i|. Each image |i| is represented by
% a flattened vector representation |ims( i,: )|. The variables |nx| and
% |ny| indicate the size of the images (20x16). The variables |nimsx| and
% |nimsy| indicate how the examples are organized -- there are |nimsx|=39 examples of
% each digit, and there are |nimsy|=36 image types total
%
% The LDA model is applied to the images by representing each binary pixel
% as a word and each image as a document. The resulting |WD| sparse matrix is
% a |W| x |D| matrix (|W| = number of pixels = |nx| x |ny|, |D| = number of
% images) where |WD(i,j)| contains the binary pixel value of pixel |i| of image
% |j|   

load 'binaryalphabet'; % loading these variables: ims nx ny nimsx nimsy W D;
nims = size( ims , 1 );

%%
% Parameters
ALPHA = 0.5;   % ALPHA hyperparameter
BETA  = 1;   % BETA hyperparameter
NITER = [ 2 5 10 25 50 ]; % the number of iterations per run
seed  = 4; % seed for the random number generator
T     = 20; % the number of topics is determined by the number of rows and columns of the images

NEXAMPLES = 10; % the number of examples to illustrate decomposition into topics
NTOPICSEXAMPLES = 4; % the number of most likely topics to show for each example image

rand( 'state' , seed );
randn( 'state' , seed );

gapx  = 5;
gapy  = 5;

fprintf( 'This script will take approx 1 minute to complete ....\n' );

%%
% Show the images
[ imcollage , labels ] = createcollage( ims , nx , ny , nimsx , nimsy , gapx , gapy );
figure; clf;
colormap( gray );
imagesc( imcollage );
axis off;
axis image;
title( 'Input Images' );

%%
% Convert data to sparse counts
[ WS , DS , n ] = convertimstocounts( ims );
 
%%
% Run the topic model
figure;
colormap( gray );

wh = 0;
for N=NITER
    [WP,DP,Z] = GibbsSamplerLDA( WS , DS , T , N , ALPHA , BETA , seed , 0 );
    
    % create a collage of these topics
    topicims = full( WP )';
    
    [ topiccollage , labels ] = createcollage( topicims , nx , ny , T , 1 , gapx , 0 );
    
    wh = wh + 1;
    subplot( length(NITER) , 1 , wh );
    maxcount = max( topiccollage(:));
    imagesc( topiccollage , [ 0 maxcount ] ); 
    
    axis image;
    axis off;
    
    title( sprintf( 'Topics after %d iterations' , N )); 
    drawnow;
end

%%
% For some example images, show the most likely topics
imsdecompose = zeros( NEXAMPLES * ( NTOPICSEXAMPLES + 1 ), W );

maxintensity = max( ims(:) );
maxwp        = max( WP(:) );

wh = 0;
for i=1:NEXAMPLES
   whx = unidrnd(nimsx);
   why = unidrnd(nimsy);
  
   %whx = 1;
   %why = i;
   
   whim = (whx-1) * nimsy + why;
   
   whim = unidrnd( nims );
   
   wh = i; 
   imsdecompose( wh , : ) = ims( whim , : ) / maxintensity;
   
   dp = DP( whim , : ) + ALPHA;
   dp = dp / sum( dp );
   [ prob , sortindex ] = sort( -dp ); prob = -prob;
   
   for t=1:NTOPICSEXAMPLES
      topic = sortindex( t );
      
      wh = i + t * NEXAMPLES;
      imsdecompose( wh , : ) = full( WP( : , topic ))' / maxwp;
   end
end

[ imscollage2 , labels ] = createcollage( imsdecompose , nx , ny ,  NTOPICSEXAMPLES+1 , NEXAMPLES , gapx , gapy );
   
figure
colormap( gray );
imagesc( imscollage2 );
axis on;
axis image;
   
set(gca,'XTickLabel',{})
set(gca,'YTickLabel',{})

xlabel( sprintf( 'Best %d topics' , NTOPICSEXAMPLES ));
ylabel( sprintf( 'Example images' ));
