#include "mex.h"
#include "cokus.cpp"

// Syntax
//   [ WP , DP , Z , X ] = GibbsSamplerSWLDA( WS , DS , T , N , ALPHA , BETA , GAMMA0, GAMMA1, SEED , OUTPUT , ZIN , XIN )

// The Special Words (SW) model


void GibbsSampler( double ALPHA, double BETA, double GAMMA0, double GAMMA1, int W, int T, int D, int NN, int OUTPUT, int n, int *z, int *x, int *d, int *w, int *wp, int *dp, int *sumdp, int *ztot, int *xcounts0, int *xcounts1, int *order, double *probs, int startcond )
{
  int wi,di,xi,i,ii,j,topic, rp, temp, iter, wioffset, dioffset, T2, sumx;
  double totprob, WBETA, TALPHA, r, max;

  T2 = T + D;
  
  sumx = 0;
  
  if (startcond == 1) {
      /* start from previously saved state */
      for (i=0; i<n; i++)
      {
          wi    = w[ i ];
          di    = d[ i ];
          topic = z[ i ];
          xi    = x[ i ];
          
          sumx += xi;
          
          sumdp[ di ]++;
          
          if (xi==0) { // regular topic assignment
              xcounts0[ di ]++;
              wp[ wi*T2 + topic ]++; // increment wp count matrix
              dp[ di*(T+1) + topic ]++; // increment dp count matrix
              ztot[ topic ]++; // increment ztot matrix
          } else {
              xcounts1[ di ]++;
              wp[ wi*T2 + T + di ]++; // increment wp count matrix
              dp[ di*(T+1) + T ]++; // increment count for special topic
              ztot[ T + di ]++; // increment ztot matrix
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
          sumx += xi;
          
          sumdp[ di ]++;
          
          if (xi==0) { // regular topic assignment
              xcounts0[ di ]++;
              wp[ wi*T2 + topic ]++; // increment wp count matrix
              dp[ di*(T+1) + topic ]++; // increment dp count matrix
              ztot[ topic ]++; // increment ztot matrix
          } else {
              xcounts1[ di ]++;
              wp[ wi*T2 + T + di ]++; // increment wp count matrix
              dp[ di*(T+1) + T ]++; // increment count for special topic
              ztot[ T + di ]++; // increment ztot matrix
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
  WBETA = (double) W * (double) BETA;
  TALPHA = (double) T * (double) ALPHA;
  
  for (iter=0; iter<NN; iter++) {
      if (OUTPUT >=1) {
          if ((iter % 10)==0) mexPrintf( "\tIteration %d of %d  (sumx=%4.2f%%)\n" , iter , NN , (double) sumx * (double) 100.0 / (double) n );
          if ((iter % 10)==0) mexEvalString("drawnow;");
      }
      for (ii = 0; ii < n; ii++) {
          i = order[ ii ]; // current word token to assess
          
          wi    = w[i]; // current word index
          di    = d[i]; // current document index
          xi    = x[i]; // current status
          topic = z[i]; // current topic assignment to word token
          
          wioffset = wi*T2;
          dioffset = di*(T+1);
          
          sumdp[ di ]--;
          
          if (xi==0) {
              xcounts0[ di ]--;
              ztot[topic]--;  
              wp[wioffset+topic]--;
              dp[dioffset+topic]--;
          } else
          {
              xcounts1[ di ]--;
              sumx--;
              ztot[ T + di ]--;  
              wp[wioffset + T + di ]--;
              dp[dioffset + T ]--;
          }
              
          //mexPrintf( "(1) Working on ii=%d i=%d wi=%d di=%d topic=%d wp=%d dp=%d\n" , ii , i , wi , di , topic , wp[wi+topic*W] , dp[wi+topic*D] );
          
          // LOOP OVER ALL REGULAR TOPICS
          totprob = (double) 0;
          for (j = 0; j < T; j++) {
              probs[j] = ((double) GAMMA0  + (double) xcounts0[ di ]) / ( ((double) GAMMA0  + (double) GAMMA1 + (double) xcounts0[ di ])  + (double) xcounts1[ di ]) * 
                         ((double) wp[ wioffset+j ] + (double) BETA)/( (double) ztot[j]+ (double) WBETA) * 
                         ((double) dp[ dioffset+ j ] + (double) ALPHA) / ((double) xcounts0[ di ] + (double) TALPHA );
                            
              totprob += probs[j];
          }
          
          // CALCULATE PROB. FOR SPECIAL TOPIC
          probs[T] = ((double) GAMMA1  + (double) xcounts1[ di ]) / ( ((double) GAMMA0  + (double) GAMMA1 + (double) xcounts0[ di ])  + (double) xcounts1[ di ]) *
                     ((double) wp[ wioffset + T + di ] + (double) BETA) / ( (double) xcounts1[ di ] + (double) WBETA);
                    
          totprob += probs[T];
          
          // sample a topic from the distribution
          r = (double) totprob * (double) randomMT() / (double) 4294967296.0;
          max = probs[0];
          topic = 0;
          while (r>max) {
              topic++;
              max += probs[topic];
          }
           
          xi = 0;
          if (topic == T) xi = 1;
          
          x[i] = xi;
          
          z[i] = topic; // assign current word token i to topic j
          
          sumdp[ di ]++;
          
          if (xi==0) {   
             xcounts0[ di ]++;
             wp[wioffset + topic ]++; // and update counts
             dp[dioffset + topic ]++;
             ztot[topic]++;
          } else
          {
             sumx++;
             xcounts1[ di ]++;
             ztot[ T + di ]++;  
             wp[wioffset + T + di ]++;
             dp[dioffset + T ]++; 
          }
              
          //mexPrintf( "(2) Working on ii=%d i=%d wi=%d di=%d topic=%d wp=%d dp=%d\n" , ii , i , wi , di , topic , wp[wi+topic*W] , dp[wi+topic*D] );
      }
  }
}

// Syntax
//   [ WP , DP , Z , X ] = GibbsSamplerSWLDA( WS , DS , T , N , ALPHA, BETA , GAMMA0, GAMMA1, SEED , OUTPUT , ZIN , XIN )


void mexFunction(int nlhs, mxArray *plhs[], int nrhs,
                 const mxArray *prhs[])
{
  double *srwp, *srdp, *probs, *Z, *X, *WS, *DS, *ZIN, *XIN;
  double ALPHA, BETA, GAMMA1, GAMMA0;
  int *irwp, *jcwp, *irdp, *jcdp;
  int *z,*d,*w, *x, *xcounts0, *xcounts1, *order, *wp, *dp, *sumdp, *ztot;
  int W,T,T2,D,NN,SEED,OUTPUT, nzmax, nzmaxwp, nzmaxdp, ntokens;
  int i,j,c,n,nt,wi,di, startcond;
  
  /* Check for proper number of arguments. */
  if (nrhs < 10) {
    mexErrMsgTxt("At least 10 input arguments required");
  } else if (nlhs < 4) {
    mexErrMsgTxt("4 output arguments required");
  }
  
  startcond = 0;
  if (nrhs >= 11) startcond = 1;
  
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

  BETA = (double) mxGetScalar(prhs[5]);
  if (BETA<=0) mexErrMsgTxt("BETA must be greater than zero");

  GAMMA0 = (double) mxGetScalar(prhs[6]);
  if (GAMMA0<=0) mexErrMsgTxt("GAMMA0 must be greater than zero");

  GAMMA1 = (double) mxGetScalar(prhs[7]);
  if (GAMMA1<=0) mexErrMsgTxt("GAMMA1 must be greater than zero");

  SEED = (int) mxGetScalar(prhs[8]);
  
  OUTPUT = (int) mxGetScalar(prhs[9]);
  
  if (startcond == 1) {
      ZIN = mxGetPr( prhs[ 10 ] );
      if (ntokens != ( mxGetM( prhs[ 10 ] ) * mxGetN( prhs[ 10 ] ))) mexErrMsgTxt("WS and ZIN vectors should have same number of entries");
      
      XIN = mxGetPr( prhs[ 11 ] );
      if (ntokens != ( mxGetM( prhs[ 11 ] ) * mxGetN( prhs[ 11 ] ))) mexErrMsgTxt("WS and XIN vectors should have same number of entries");
  }
  
  // seeding
  seedMT( 1 + SEED * 2 ); // seeding only works on uneven numbers
  
   
  
  /* allocate memory */
  z  = (int *) mxCalloc( ntokens , sizeof( int ));
  x  = (int *) mxCalloc( ntokens , sizeof( int ));
  
  if (startcond == 1) {
     for (i=0; i<ntokens; i++) z[ i ] = (int) ZIN[ i ] - 1;
     for (i=0; i<ntokens; i++) x[ i ] = (int) XIN[ i ]; // 0 = normal topic assignment 1 = idiosyncratic topic
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
  
  // NOTE: the wp matrix has T+D topics where the last D topics are idiosyncratic
  T2 = T + D;
  wp  = (int *) mxCalloc( T2*W , sizeof( int ));
  
  // NOTE: the last topic probability is the special topic probability
  dp  = (int *) mxCalloc( (T+1)*D , sizeof( int ));
  
  sumdp  = (int *) mxCalloc( D , sizeof( int ));
  
  ztot  = (int *) mxCalloc( T2 , sizeof( int ));
  probs  = (double *) mxCalloc( T+1 , sizeof( double ));
  xcounts0 = (int *) mxCalloc( D , sizeof( int ));
  xcounts1 = (int *) mxCalloc( D , sizeof( int ));
  
  //mexPrintf( "N=%d  T=%d W=%d D=%d\n" , ntokens , T , W , D );
  
  if (OUTPUT==2) {
      mexPrintf( "Running LDA Gibbs Sampler -- with special topics\n" );
      if (startcond==1) mexPrintf( "Starting from previous state ZIN\n" );
      mexPrintf( "Arguments:\n" );
      mexPrintf( "\tNumber of words      W = %d\n"    , W );
      mexPrintf( "\tNumber of docs       D = %d\n"    , D );
      mexPrintf( "\tNumber of topics     T = %d\n"    , T );
      mexPrintf( "\tNumber of iterations N = %d\n"    , NN );
      mexPrintf( "\tHyperparameter   ALPHA = %4.4f\n" , ALPHA );
      mexPrintf( "\tHyperparameter    BETA = %4.4f\n" , BETA );
      mexPrintf( "\tHyperparameter  GAMMA0 = %4.4f\n" , GAMMA0 );
      mexPrintf( "\tHyperparameter  GAMMA1 = %4.4f\n" , GAMMA1 );
      mexPrintf( "\tSeed number            = %d\n"    , SEED );
      mexPrintf( "\tNumber of tokens       = %d\n"    , ntokens );
      //mexPrintf( "Internal Memory Allocation\n" );
      //mexPrintf( "\tw,d,z,order indices combined = %d bytes\n" , 4 * sizeof( int) * ntokens );
      //mexPrintf( "\twp (full) matrix = %d bytes\n" , sizeof( int ) * W * T  );
      //mexPrintf( "\tdp (full) matrix = %d bytes\n" , sizeof( int ) * D * T  );
      //mexPrintf( "Checking: sizeof(int)=%d sizeof(long)=%d sizeof(double)=%d\n" , sizeof(int) , sizeof(long) , sizeof(double));
  }
  
  /* run the model */
  GibbsSampler( ALPHA, BETA, GAMMA0, GAMMA1, W, T, D, NN, OUTPUT, n, z, x, d, w, wp, dp, sumdp, ztot, xcounts0, xcounts1, order, probs, startcond );
  
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
      for (j=0; j<(T+1); j++)
          nzmaxdp += (int) ( *( dp + j + i*(T+1) )) > 0;
  }  
  /*if (OUTPUT==2) {
      mexPrintf( "Constructing sparse output matrix dp\n" );
      mexPrintf( "Number of nonzero entries for DP = %d\n" , nzmaxdp );
  } */ 
  plhs[1] = mxCreateSparse( D,T+1,nzmaxdp,mxREAL);
  srdp  = mxGetPr(plhs[1]);
  irdp = mxGetIr(plhs[1]);
  jcdp = mxGetJc(plhs[1]);
  n = 0;
  for (j=0; j<T+1; j++) {
      *( jcdp + j ) = n;
      for (i=0; i<D; i++) {
          c = (int) *( dp + i*(T+1) + j );
          if (c >0) {
              *( srdp + n ) = c;
              *( irdp + n ) = i;
              n++;
          }
      }
  }
  *( jcdp + (T+1) ) = n;
  
  plhs[ 2 ] = mxCreateDoubleMatrix( 1,ntokens , mxREAL );
  Z = mxGetPr( plhs[ 2 ] );
  for (i=0; i<ntokens; i++) Z[ i ] = (double) z[ i ] + 1;
  
  plhs[ 3 ] = mxCreateDoubleMatrix( 1,ntokens , mxREAL );
  X = mxGetPr( plhs[ 3 ] );
  for (i=0; i<ntokens; i++) X[ i ] = (double) x[ i ];
}
