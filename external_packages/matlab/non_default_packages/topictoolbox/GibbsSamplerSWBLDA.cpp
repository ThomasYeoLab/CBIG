#include "mex.h"
#include "cokus.cpp"

// Syntax
//   [ WP , DP , Z , X ] = GibbsSamplerSWBLDA( WS , DS , T , N , ALPHA , BETA0 , BETA1, BETA2 , GAMMA0, GAMMA1, GAMMA2, SEED , OUTPUT , ZIN , XIN )

// The Special Words plus Background (SWB) LDA model


void GibbsSampler( double ALPHA, double BETA0, double BETA1, double BETA2, double GAMMA0, double GAMMA1, double GAMMA2, 
                   int W, int T, int D, int NN, int OUTPUT, int n, int *z, int *x, int *d, int *w, int *wp, int *dp, 
                   int *sumdp, int *ztot, int *xcounts0, int *xcounts1, int *xcounts2, int *order, double *probs, int startcond )
{
  int wi,di,xi,i,ii,j,topic, rp, temp, iter, wioffset, dioffset, T2, sumx1, sumx2;
  double totprob, WBETA0, WBETA1, WBETA2, TALPHA, r, max, normfact1;

  T2 = T + D + 1;
  
  sumx1 = 0;
  sumx2 = 0;
  
  if (startcond == 1) {
      /* start from previously saved state */
      for (i=0; i<n; i++)
      {
          wi    = w[ i ];
          di    = d[ i ];
          topic = z[ i ];
          xi    = x[ i ];
          
          if (xi==1) sumx1 += 1;
          if (xi==2) sumx2 += 1;
          
          sumdp[ di ]++;
          
          if (xi==0) 
          { // regular topic assignment
              xcounts0[ di ]++;
              wp[ wi*T2 + topic ]++; // increment wp count matrix
              dp[ di*(T+2) + topic ]++; // increment dp count matrix
              ztot[ topic ]++; // increment ztot matrix
          } else if (xi==1) 
          { // special topic
              xcounts1[ di ]++;
              wp[ wi*T2 + T + di ]++; // increment wp count matrix
              dp[ di*(T+2) + T ]++; // increment count for special topic
              ztot[ T + di ]++; // increment ztot matrix
          } else if (xi==2) 
          { // background distribution
              xcounts2[ di ]++;
              wp[ wi*T2 + T + D ]++; // increment wp count matrix    
              dp[ di*(T+2) + T + 1 ]++; // increment count for special topic
              ztot[ T + D ]++; // increment ztot matrix
          }
              
      }
  }
  
  if (startcond == 0) {
  /* random initialization */
      if (OUTPUT==2) mexPrintf( "Starting Random initialization\n" );
      for (i=0; i<n; i++)
      {
          xi = 0;
          wi = w[ i ];
          di = d[ i ];
          // pick a random topic 0..T-1
          topic = (int) ( (double) randomMT() * (double) T / (double) (4294967296.0 + 1.0) );
          z[ i ] = topic; // assign this word token to this topic
          x[ i ] = xi;
          
          sumdp[ di ]++;
          
          if (xi==0) { // regular topic assignment
              xcounts0[ di ]++;
              wp[ wi*T2 + topic ]++; // increment wp count matrix
              dp[ di*(T+2) + topic ]++; // increment dp count matrix
              ztot[ topic ]++; // increment ztot matrix
          } else if (xi==1)
          {
              sumx1++;
              xcounts1[ di ]++;
              wp[ wi*T2 + T + di ]++; // increment wp count matrix
              dp[ di*(T+2) + T ]++; // increment count for special topic
              ztot[ T + di ]++; // increment ztot matrix
          } else if (xi==2)
          {
              sumx2++;
              xcounts2[ di ]++;
              wp[ wi*T2 + T + D ]++; // increment wp count matrix
              dp[ di*(T+2) + T + 1]++; // increment count for special topic
              ztot[ T + D ]++; // increment ztot matrix
          }
      }
  }
  
  if (OUTPUT==2) mexPrintf( "Determining random order update sequence\n" );
  
  for (i=0; i<n; i++) order[i]=i; // fill with increasing series
  for (i=0; i<(n-1); i++) {
      // pick a random integer between i and nw
      rp = i + (int) ((double) (n-i) * (double) randomMT() / (double) (4294967296.0 + 1.0));
      
      // switch contents on position i and position rp
      temp = order[rp];
      order[rp]=order[i];
      order[i]=temp;
  }
  
  //for (i=0; i<n; i++) mexPrintf( "i=%3d order[i]=%3d\n" , i , order[ i ] );
  WBETA0 = (double) W * (double) BETA0;
  WBETA1 = (double) W * (double) BETA1;
  WBETA2 = (double) W * (double) BETA2;
  TALPHA = (double) T * (double) ALPHA;
  
  for (iter=0; iter<NN; iter++) {
      if (OUTPUT >=1) {
          if ((iter % 10)==0) mexPrintf( "\tIteration %d of %d  (words to special topics and background: %4.2f%% %4.2f%%)\n" , iter , NN , (double) sumx1 * (double) 100.0 / (double) n , (double) sumx2 * (double) 100.0 / (double) n);
          if ((iter % 10)==0) mexEvalString("drawnow;");
      }
      for (ii = 0; ii < n; ii++) {
          i = order[ ii ]; // current word token to assess
          
          wi    = w[i]; // current word index
          di    = d[i]; // current document index
          xi    = x[i]; // current status
          topic = z[i]; // current topic assignment to word token
          
          wioffset = wi*T2;
          dioffset = di*(T+2);
          
          sumdp[ di ]--;
          
          if (xi==0) {
              xcounts0[ di ]--;
              ztot[topic]--;  
              wp[wioffset+topic]--;
              dp[dioffset+topic]--;
          } else if (xi==1)
          {
              xcounts1[ di ]--;
              sumx1--;
              ztot[ T + di ]--;  
              wp[wioffset + T + di ]--;
              dp[dioffset + T ]--;
          } else if (xi==2)
          {
              xcounts2[ di ]--;
              sumx2--;
              ztot[ T + D ]--;  
              wp[wioffset + T + D ]--;
              dp[dioffset + T + 1 ]--;
          }
              
          //mexPrintf( "(1) Working on ii=%d i=%d wi=%d di=%d topic=%d wp=%d dp=%d\n" , ii , i , wi , di , topic , wp[wi+topic*W] , dp[wi+topic*D] );
          
          // LOOP OVER ALL REGULAR TOPICS
          totprob = (double) 0;
          
          normfact1 = ((double) GAMMA0  + (double) xcounts0[ di ] ) / ((double) GAMMA0  + (double) GAMMA1 + (double) GAMMA2 + (double) xcounts0[ di ]  + (double) xcounts1[ di ] + (double) xcounts2[ di ]);
          
          for (j = 0; j < T; j++) {
              probs[j] =  normfact1 * 
                         ((double) wp[ wioffset+j ] + (double) BETA0)/( (double) ztot[j]+ (double) WBETA0) * 
                         ((double) dp[ dioffset+ j ] + (double) ALPHA) / ((double) xcounts0[ di ] + (double) TALPHA );
                            
              totprob += probs[j];
          }
          
          // CALCULATE PROB. FOR SPECIAL WORD TOPIC
          probs[T] = ((double) GAMMA1  + (double) xcounts1[ di ]) / ( ((double) GAMMA0  + (double) GAMMA1 + (double) GAMMA2 + (double) xcounts0[ di ])  + (double) xcounts1[ di ] + (double) xcounts2[ di ]) *
                     ((double) wp[ wioffset + T + di ] + (double) BETA1) / ( (double) xcounts1[ di ] + (double) WBETA1);
          
          totprob += probs[T];
          
          // CALCULATE PROB. FOR BACKGROUND DISTRIBUTION
          probs[T+1] = ((double) GAMMA2  + (double) xcounts2[ di ]) / ( ((double) GAMMA0  + (double) GAMMA1 + (double) GAMMA2 + (double) xcounts0[ di ])  + (double) xcounts1[ di ] + (double) xcounts2[ di ]) *
                     ((double) wp[ wioffset + T + D ] + (double) BETA2) / ( (double) sumx2 + (double) WBETA2);
          
          totprob += probs[T+1];
          
          // sample a topic from the distribution
          r = (double) totprob * (double) randomMT() / (double) 4294967296.0;
          max = probs[0];
          topic = 0;
          while (r>max) {
              topic++;
              max += probs[topic];
          }
           
          xi = 0;
          if (topic == T)   xi = 1;
          if (topic == T+1) xi = 2;
          
          x[i] = xi;
          
          z[i] = topic; // assign current word token i to topic j
          
          sumdp[ di ]++;
          
          if (xi==0) {   
             xcounts0[ di ]++;
             wp[wioffset + topic ]++; // and update counts
             dp[dioffset + topic ]++;
             ztot[topic]++;
          } else if (xi==1)
          {
             sumx1++;
             xcounts1[ di ]++;
             ztot[ T + di ]++;  
             wp[wioffset + T + di ]++;
             dp[dioffset + T ]++; 
          } else if (xi==2)
          {
             sumx2++;
             xcounts2[ di ]++;
             ztot[ T + D ]++;  
             wp[wioffset + T + D ]++;
             dp[dioffset + T + 1 ]++; 
          }
              
          //mexPrintf( "(2) Working on ii=%d i=%d wi=%d di=%d topic=%d wp=%d dp=%d\n" , ii , i , wi , di , topic , wp[wi+topic*W] , dp[wi+topic*D] );
      }
  }
}

// Syntax
//   [ WP , DP , Z , X ] = GibbsSamplerSWBLDA( WS , DS , T , N , ALPHA, BETA0 , BETA1, BETA2 , GAMMA0, GAMMA1, GAMMA2, SEED , OUTPUT , ZIN , XIN )


void mexFunction(int nlhs, mxArray *plhs[], int nrhs,
                 const mxArray *prhs[])
{
  double *srwp, *srdp, *probs, *Z, *X, *WS, *DS, *ZIN, *XIN;
  double ALPHA, BETA0, BETA1, BETA2, GAMMA1, GAMMA0, GAMMA2;
  int *irwp, *jcwp, *irdp, *jcdp;
  int *z,*d,*w, *x, *xcounts0, *xcounts1, *xcounts2, *order, *wp, *dp, *sumdp, *ztot;
  int W,T,T2,D,NN,SEED,OUTPUT, nzmax, nzmaxwp, nzmaxdp, ntokens;
  int i,j,c,n,nt,wi,di, startcond;
  
  /* Check for proper number of arguments. */
  if (nrhs < 13) {
    mexErrMsgTxt("At least 13 input arguments required");
  } else if (nlhs < 4) {
    mexErrMsgTxt("4 output arguments required");
  }
  
  startcond = 0;
  if (nrhs >= 14) startcond = 1;
  
  /* process the input arguments */
  if (mxIsDouble( prhs[ 0 ] ) != 1) mexErrMsgTxt("WS input vector must be a double precision matrix");
  if (mxIsDouble( prhs[ 1 ] ) != 1) mexErrMsgTxt("DS input vector must be a double precision matrix");
  
  // pointer to word indices
  WS = mxGetPr( prhs[ 0 ] );
     
  // pointer to document indices
  DS = mxGetPr( prhs[ 1 ] );
  
  // get the number of tokens
  ntokens = mxGetM( prhs[ 0 ] ) * mxGetN( prhs[ 0 ] );
  
  
  if (ntokens == 0) mexErrMsgTxt("WS vector is empty"); 
  if (ntokens != ( mxGetM( prhs[ 1 ] ) * mxGetN( prhs[ 1 ] ))) mexErrMsgTxt("WS and DS vectors should have same number of entries");
  
  T    = (int) mxGetScalar(prhs[2]);
  if (T<=0) mexErrMsgTxt("Number of topics must be greater than zero");
  
  NN    = (int) mxGetScalar(prhs[3]);
  if (NN<0) mexErrMsgTxt("Number of iterations must be positive");
  
  ALPHA = (double) mxGetScalar(prhs[4]);
  if (ALPHA<=0) mexErrMsgTxt("ALPHA0 must be greater than zero");

  BETA0 = (double) mxGetScalar(prhs[5]);
  if (BETA0<=0) mexErrMsgTxt("BETA0 must be greater than zero");

  BETA1 = (double) mxGetScalar(prhs[6]);
  if (BETA1<=0) mexErrMsgTxt("BETA1 must be greater than zero");

  BETA2 = (double) mxGetScalar(prhs[7]);
  if (BETA2<=0) mexErrMsgTxt("BETA2 must be greater than zero");

  GAMMA0 = (double) mxGetScalar(prhs[8]);
  if (GAMMA0<=0) mexErrMsgTxt("GAMMA0 must be greater than zero");

  GAMMA1 = (double) mxGetScalar(prhs[9]);
  if (GAMMA1<=0) mexErrMsgTxt("GAMMA1 must be greater than zero");

  GAMMA2 = (double) mxGetScalar(prhs[10]);
  if (GAMMA2<=0) mexErrMsgTxt("GAMMA2 must be greater than zero");
  
  SEED = (int) mxGetScalar(prhs[11]);
  
  OUTPUT = (int) mxGetScalar(prhs[12]);
  
  if (startcond == 1) {
      ZIN = mxGetPr( prhs[ 13 ] );
      if (ntokens != ( mxGetM( prhs[ 13 ] ) * mxGetN( prhs[ 13 ] ))) mexErrMsgTxt("WS and ZIN vectors should have same number of entries");
      
      XIN = mxGetPr( prhs[ 14 ] );
      if (ntokens != ( mxGetM( prhs[ 14 ] ) * mxGetN( prhs[ 14 ] ))) mexErrMsgTxt("WS and XIN vectors should have same number of entries");
  }
  
  // seeding
  seedMT( 1 + SEED * 2 ); // seeding only works on uneven numbers
  
   
  
  /* allocate memory */
  z  = (int *) mxCalloc( ntokens , sizeof( int ));
  x  = (int *) mxCalloc( ntokens , sizeof( int ));
  
  if (startcond == 1) {
     for (i=0; i<ntokens; i++) z[ i ] = (int) ZIN[ i ] - 1;
     for (i=0; i<ntokens; i++) x[ i ] = (int) XIN[ i ]; // 0 = normal topic assignment 1 = special document topic 2 = background distribution
  }
  
  d  = (int *) mxCalloc( ntokens , sizeof( int ));
  w  = (int *) mxCalloc( ntokens , sizeof( int ));
  order  = (int *) mxCalloc( ntokens , sizeof( int ));  
  
  
  // copy over the word and document indices into internal format
  for (i=0; i<ntokens; i++) {
     w[ i ] = (int) WS[ i ] - 1;
     d[ i ] = (int) DS[ i ] - 1;
  }
  
  n = ntokens;
  
  W = 0;
  D = 0;
  for (i=0; i<n; i++) {
     if (w[ i ] > W) W = w[ i ];
     if (d[ i ] > D) D = d[ i ];
  }
  W = W + 1;
  D = D + 1;
  
  // NOTE: the wp matrix has T+D+1 topics for T regular topics, D special word topics and 1 background distribution
  T2 = T + D + 1;
  wp  = (int *) mxCalloc( T2*W , sizeof( int ));
  
  // NOTE: the last two topic probabilities are for the special word topic and background topic
  dp  = (int *) mxCalloc( (T+2)*D , sizeof( int ));
  
  sumdp  = (int *) mxCalloc( D , sizeof( int ));
  
  ztot  = (int *) mxCalloc( T2 , sizeof( int ));
  probs  = (double *) mxCalloc( T+2 , sizeof( double ));
  xcounts0 = (int *) mxCalloc( D , sizeof( int ));
  xcounts1 = (int *) mxCalloc( D , sizeof( int ));
  xcounts2 = (int *) mxCalloc( D , sizeof( int ));
  
  //mexPrintf( "N=%d  T=%d W=%d D=%d\n" , ntokens , T , W , D );
  
  if (OUTPUT==2) {
      mexPrintf( "Running SWB LDA Gibbs Sampler\n" );
      if (startcond==1) mexPrintf( "Starting from previous state ZIN\n" );
      mexPrintf( "Arguments:\n" );
      mexPrintf( "\tNumber of words      W = %d\n"    , W );
      mexPrintf( "\tNumber of docs       D = %d\n"    , D );
      mexPrintf( "\tNumber of topics     T = %d\n"    , T );
      mexPrintf( "\tNumber of iterations N = %d\n"    , NN );
      mexPrintf( "\tHyperparameter   ALPHA = %4.4f\n" , ALPHA );
      mexPrintf( "\tHyperparameter   BETA0 = %4.4f\n" , BETA0 );
      mexPrintf( "\tHyperparameter   BETA1 = %4.4f\n" , BETA1 );
      mexPrintf( "\tHyperparameter   BETA2 = %4.4f\n" , BETA2 );
      mexPrintf( "\tHyperparameter  GAMMA0 = %4.4f\n" , GAMMA0 );
      mexPrintf( "\tHyperparameter  GAMMA1 = %4.4f\n" , GAMMA1 );
      mexPrintf( "\tHyperparameter  GAMMA2 = %4.4f\n" , GAMMA2 );
      mexPrintf( "\tSeed number            = %d\n"    , SEED );
      mexPrintf( "\tNumber of tokens       = %d\n"    , ntokens );
  }
  
  /* run the model */
  GibbsSampler( ALPHA, BETA0 , BETA1, BETA2, GAMMA0, GAMMA1, GAMMA2, 
                W, T, D, NN, OUTPUT, n, z, x, d, w, wp, dp, sumdp, ztot, 
                xcounts0, xcounts1, xcounts2, order, probs, startcond );
  
  /* convert the full wp matrix into a sparse matrix */
  nzmaxwp = 0;
  for (i=0; i<W; i++) {
     for (j=0; j<T2; j++)
         nzmaxwp += (int) ( *( wp + j + i*T2 )) > 0;
  }  
  /*if (OUTPUT==2) {
      mexPrintf( "Constructing sparse output matrix wp\n" );
      mexPrintf( "Number of nonzero entries for WP = %d\n" , nzmaxwp );
  }*/
  
  // MAKE THE WP SPARSE MATRIX
  plhs[0] = mxCreateSparse( W,T2,nzmaxwp,mxREAL);
  srwp  = mxGetPr(plhs[0]);
  irwp = mxGetIr(plhs[0]);
  jcwp = mxGetJc(plhs[0]);  
  n = 0;
  for (j=0; j<T2; j++) {
      *( jcwp + j ) = n;
      for (i=0; i<W; i++) {
         c = (int) *( wp + i*T2 + j );
         if (c >0) {
             *( srwp + n ) = c;
             *( irwp + n ) = i;
             n++;
         }
      }    
  }  
  *( jcwp + T2 ) = n;    
   
  // MAKE THE DP SPARSE MATRIX
  nzmaxdp = 0;
  for (i=0; i<D; i++) {
      for (j=0; j<(T+2); j++)
          nzmaxdp += (int) ( *( dp + j + i*(T+2) )) > 0;
  }  
  /*if (OUTPUT==2) {
      mexPrintf( "Constructing sparse output matrix dp\n" );
      mexPrintf( "Number of nonzero entries for DP = %d\n" , nzmaxdp );
  } */ 
  plhs[1] = mxCreateSparse( D,T+2,nzmaxdp,mxREAL);
  srdp  = mxGetPr(plhs[1]);
  irdp = mxGetIr(plhs[1]);
  jcdp = mxGetJc(plhs[1]);
  n = 0;
  for (j=0; j<T+2; j++) {
      *( jcdp + j ) = n;
      for (i=0; i<D; i++) {
          c = (int) *( dp + i*(T+2) + j );
          if (c >0) {
              *( srdp + n ) = c;
              *( irdp + n ) = i;
              n++;
          }
      }
  }
  *( jcdp + (T+2) ) = n;
  
  plhs[ 2 ] = mxCreateDoubleMatrix( 1,ntokens , mxREAL );
  Z = mxGetPr( plhs[ 2 ] );
  for (i=0; i<ntokens; i++) Z[ i ] = (double) z[ i ] + 1;
  
  plhs[ 3 ] = mxCreateDoubleMatrix( 1,ntokens , mxREAL );
  X = mxGetPr( plhs[ 3 ] );
  for (i=0; i<ntokens; i++) X[ i ] = (double) x[ i ];
}
