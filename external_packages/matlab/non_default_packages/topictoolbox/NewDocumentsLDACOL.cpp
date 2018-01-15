#include "mex.h"
#include "cokus.cpp"
      
// Syntax                                       0       1          2       3       4   5    6         7      8       9      10        11      12      13      14       15       16       17     18 
//   [ DP_NEW,C,Z,PROB ] = NewDocumentsLDACOL( WS_NEW , DS_NEW , SI_NEW , WW_NEW , T , N , NSAMPLE , ALPHA , BETA , GAMMA0, GAMMA1 , DELTA , SEED , OUTPUT , WW_OLD , WP_OLD , WC_OLD , WS_OLD, UPDATECOLS );
     
// work on C, Z and PROB output
  
void GibbsSamplerLDACOL( double ALPHA, double BETA, double GAMMA0, double GAMMA1, double DELTA, int W, int T, int D, int NITER, 
                            int NSAMPLE, int OUTPUT, int n, int *z, int *d, int *w, int *s, int *c, int *sc, int *wp, int *dp, int *dp2, 
                            int *ztot, int *wtot, int *order, double *probs, int *wc, double *PROBA, int updatecols )
{
  int wi,wipre,di,si,i,ii,j,topic,count, route,rp, temp, iter, wioffset, dioffset;
  int startindex,endindex,currentrow,prevroute, ncol, k;
  double totprob, WBETA, WDELTA,r, max, C0, C1;

  if (OUTPUT==2) mexPrintf( "Starting Random initialization\n" );
  
  ncol = 0;
  for (i=0; i<n; i++)
  {
      wi = w[ i ]; // word index
      di = d[ i ]; // document index
      si = s[ i ]; // status of word
      
      topic = (int) ( (double) randomMT() * (double) T / (double) (4294967296.0 + 1.0) ); // pick a random topic 0..T-1
      
      z[ i ] = topic; // assign this word token to this topic
      dp[ di*T + topic ]++; // increment dp count matrix
            
      // initially, set the route to 0 for topic route
      route = 0;
      c[ i ] = route;
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
  WBETA = (double) (W*BETA);
  WDELTA = (double) (W*DELTA);
  
  //mexPrintf( "Nsample = %d\n" , NSAMPLE );
  
  for (iter=0; iter<(NITER+NSAMPLE); iter++) {
      if (OUTPUT >=1) {
          if (iter == NITER) mexPrintf( "Starting oversampling DP matrix\n" );
          if ((iter % 10)==0) mexPrintf( "\tIteration %d of %d;   Number of tokens in collocation = %d\n" , iter , NITER+NSAMPLE , ncol );
          if ((iter % 10)==0) mexEvalString("drawnow;");
      }
      
      if (iter==NITER) {
          // copy over current DP counts into DP2 matrix -- this matrix will be oversampled   
          for (i=0; i<(T*D); i++) dp2[ i ] = dp[ i ];
      }
      
      for (ii = 0; ii < n; ii++) {
          i = order[ ii ]; // current word token to assess
          
          wi  = w[i]; // current word index
          di  = d[i]; // current document index  
          si  = s[i]; // current status
          route = c[i]; // current route assignment
          topic = z[i]; // current topic assignment to word token
          
          /* ---------------------------------------
              sample a topic for this word token
            ---------------------------------------- */
          wioffset = wi*T;
          dioffset = di*T; 
                  
          if (route==0) { // ROUTE = 0; TOPICS
              dp[dioffset+topic]--;
              totprob = (double) 0;
              for (j = 0; j < T; j++) {
                  probs[j] = ((double) wp[ wioffset+j ] + BETA)/( (double) ztot[j]+ WBETA)*( (double) dp[ dioffset+ j ] + ALPHA);
                  totprob += probs[j];
              }
          } else { // ROUTE = 1; COLLOCATION
              totprob = (double) 0;
              for (j = 0; j < T; j++) {
                  probs[j] = ( (double) dp[ dioffset+ j ] + ALPHA);
                  totprob += probs[j];
              }
          }
          
          // sample a topic from the distribution
          r = (double) totprob * (double) randomMT() / (double) 4294967296.0;
          max = probs[0];
          topic = 0;
          while (r>max) {
              topic++;
              max += probs[topic];
          }
          
          z[i] = topic; // assign current word token i to topic j
          
           /* ---------------------------------------
              sample a route for this word token
            ---------------------------------------- */
          
          if ((si == 1) && (i>0)) { // ony consider a route change for words that can be part of a collocation with a *previous* word
              wipre = w[ i-1 ];
              
              prevroute = route;
              
              if (updatecols==1) wtot[ wipre ]--;
              if (prevroute==1) {
                  if (updatecols==1) wc[ wipre ]--; // if the current word token was assigned to a collocation, take this out of the count for the previous word
                  ncol--;
              }
         
              
              // calculate the (unnormalized) probability for route = 0 (topic)
              C0 = ((double) wp[ wioffset+topic ] + (double) BETA)/( (double) ztot[topic]+ (double) WBETA) * ((double) (wtot[ wipre ] - wc[ wipre ]) + GAMMA0) / ((double) wtot[ wipre ] + GAMMA0 + GAMMA1);
              
              // calculate the (unnormalized) probability for route = 1 (collocation)
              C1 = ((double) (sc[ i ]-1)  + (double) DELTA)/( (double) wtot[ wipre ] + (double) WDELTA) * ((double) wc[ wipre ] + GAMMA1) / ((double) wtot[ wipre ] + GAMMA0 + GAMMA1);
              
              C0 = C0 / (C0 + C1);
              
              r = (double) randomMT() / (double) 4294967296.0;
              route = 0;
              if (r > C0) route = 1;
              
              if (updatecols==1) wtot[ wipre ]++;
              
              // is the word part of a collocation?
              if (route==1) {
                  ncol++;
                  if (updatecols==1) wc[ wipre ]++; // increment collocation count for this word
              } else
              {
                  dp[dioffset + topic ]++;
                  if (iter >= NITER) dp2[dioffset + topic ]++;
              }
            
              c[ i ] = route;
          } else {
            // add counts to topic route
            dp[dioffset + topic ]++;
            if (iter >= NITER) dp2[dioffset + topic ]++;
          }
      }
  }
  
  if (OUTPUT >=1) mexPrintf( "Calcalating most probable assignments\n" );
  
  // Calculate most probable assignments
  for (i = 0; i < n; i++) {
      
      wi  = w[i]; // current word index
      di  = d[i]; // current document index
      si  = s[i]; // current status
      wioffset = wi*T;
      dioffset = di*T;
      
      // sample most likely topic assignment, if this were a topic route      
      totprob = (double) 0;
      for (j = 0; j < T; j++) {
          // PICK DP2 !!!!!!!!!!!1111
          probs[j] = ((double) wp[ wioffset+j ] + BETA)/( (double) ztot[j]+ WBETA)*( (double) dp2[ dioffset+ j ] + ALPHA);
          totprob += probs[j];
      }    
      max   = 0;
      topic = 0;
      for (k=0; k<T; k++) {
          if (probs[ k ] > max) {
              max = probs[ k ];
              topic = k;
          }
      }
      max = max / totprob;
      
      z[i] = topic; // assign current word token i to topic j
      PROBA[ i ] = max;
      
      //   calculate most likely route
      route = 0;
      if ((si == 1) && (i>0)) { // ony consider a route change for words that can be part of a collocation with a *previous* word
          wipre = w[ i-1 ];
          
          // calculate the (unnormalized) probability for route = 0 (topic)
          C0 = ((double) wp[ wioffset+topic ] + (double) BETA)/( (double) ztot[topic]+ (double) WBETA) * ((double) (wtot[ wipre ] - wc[ wipre ]) + GAMMA0) / ((double) wtot[ wipre ] + GAMMA0 + GAMMA1);
          
          // calculate the (unnormalized) probability for route = 1 (collocation)
          C1 = ((double) (sc[ i ]-1)  + (double) DELTA)/( (double) wtot[ wipre ] + (double) WDELTA) * ((double) wc[ wipre ] + GAMMA1) / ((double) wtot[ wipre ] + GAMMA0 + GAMMA1);
          
          C0 = C0 / (C0 + C1);
          C1 = 1 - C0;
          
          if (C1 > C0) {
              route = 1;
              PROBA[ i ] = C1;
          }
      }
      c[ i ] = route;
   }
}

void mexFunction(int nlhs, mxArray *plhs[], int nrhs,
                 const mxArray *prhs[])
{
  double *srww_old, *srww_new,*srwp_old,*srdp_new, *probs, *DS_NEW, *WS_NEW, *WS_OLD, *SI_NEW, *WC_OLD, *PROBA, *C, *Z;
  double ALPHA,BETA, GAMMA0, GAMMA1, DELTA;
  int *irww_old, *irww_new, *jcww_old, *jcww_new, *irwp_old, *jcwp_old, *irdp_new, *jcdp_new;
  int *z,*d,*w, *s, *c, *sc, *order, *wp_old, *dp_new, *dp_new2, *ztot_old, *wtot_new, *wc_new;
  int W_NEW, W_OLD,T_NEW, T_OLD,D_NEW,NITER,SEED,OUTPUT, nzmax_new, nzmax_old, nzmaxdp_new, ntokens_new, ntokens_old;
  int i,j,cc,nt,n,wi,wipre,di, startindex, endindex, currentrow, count;
  int NSAMPLE, updatecols;
  
  /* Check for proper number of arguments. */
  if (nrhs != 19) {
    mexErrMsgTxt("19 input arguments required");
  } else if (nlhs != 4) {
    mexErrMsgTxt("4 output arguments required");
  }

// Syntax                                       0       1          2       3       4   5    6         7      8       9      10        11      12      13      14       15       16       17     18 
//   [ DP_NEW,C,Z,PROB ] = NewDocumentsLDACOL( WS_NEW , DS_NEW , SI_NEW , WW_NEW , T , N , NSAMPLE , ALPHA , BETA , GAMMA0, GAMMA1 , DELTA , SEED , OUTPUT , WW_OLD , WP_OLD , WC_OLD , WS_OLD, UPDATECOLS );
  
  
  /* process the input arguments */
  if (mxIsDouble( prhs[ 0 ] ) != 1)
      mexErrMsgTxt("WS_NEW must be double precision"); 
  WS_NEW = mxGetPr(prhs[0]);
  ntokens_new = mxGetM( prhs[ 0 ] ) * mxGetN( prhs[ 0 ] );
  
  if (mxIsDouble( prhs[ 1 ] ) != 1)
      mexErrMsgTxt("DS_NEW must be double precision"); 
  DS_NEW = mxGetPr(prhs[1]);
  
  if (mxIsDouble( prhs[ 2 ] ) != 1)
      mexErrMsgTxt("SI_NEW must be double precision");
  SI_NEW = mxGetPr(prhs[2]);
  
  if ((mxIsSparse( prhs[ 3 ] ) != 1) || (mxIsDouble( prhs[ 3 ] ) != 1))
      mexErrMsgTxt("WW_NEW collocation matrix must be a sparse double precision matrix");

  /* dealing with sparse array WW_NEW */
  srww_new  = mxGetPr(prhs[3]);
  irww_new  = mxGetIr(prhs[3]);
  jcww_new  = mxGetJc(prhs[3]);
  nzmax_new = mxGetNzmax(prhs[3]);  
  W_NEW     = mxGetM( prhs[3] );
  
  if (mxGetN( prhs[3] ) != W_NEW) mexErrMsgTxt("WW_NEW matrix should be square");
  
  D_NEW    = 0;
  for (i=0; i<ntokens_new; i++) {
      if (DS_NEW[ i ] > D_NEW) D_NEW = (int) DS_NEW[ i ];
      if (WS_NEW[ i ] > W_NEW) mexErrMsgTxt("Some word tokens in WS_NEW stream exceed number of word types in WW_NEW matrix");
  }
   
  T_NEW    = (int) mxGetScalar(prhs[4]);
  if (T_NEW<=0) mexErrMsgTxt("Number of topics must be greater than zero");
  
  NITER    = (int) mxGetScalar(prhs[5]);
  if (NITER<0) mexErrMsgTxt("Number of iterations must be greater than zero");
  
  NSAMPLE    = (int) mxGetScalar(prhs[6]);
  if (NSAMPLE<0) mexErrMsgTxt("Number of samples must be greater than zero");
  
  ALPHA = (double) mxGetScalar(prhs[7]);
  if (ALPHA<=0) mexErrMsgTxt("ALPHA must be greater than zero");
  
  BETA = (double) mxGetScalar(prhs[8]);
  if (BETA<=0) mexErrMsgTxt("BETA must be greater than zero");
  
  GAMMA0 = (double) mxGetScalar(prhs[9]);
  if (GAMMA0<=0) mexErrMsgTxt("GAMMA0 must be greater than zero");
  
  GAMMA1 = (double) mxGetScalar(prhs[10]);
  if (GAMMA1<=0) mexErrMsgTxt("GAMMA1 must be greater than zero");
  
  DELTA = (double) mxGetScalar(prhs[11]);
  if (DELTA<=0) mexErrMsgTxt("DELTA must be greater than zero");
  
  SEED = (int) mxGetScalar(prhs[12]);
  // set the seed of the random number generator
  
  OUTPUT = (int) mxGetScalar(prhs[13]);
  
  // dealing with sparse array WW_OLD 
  if ((mxIsSparse( prhs[ 14 ] ) != 1) || (mxIsDouble( prhs[ 14 ] ) != 1))
      mexErrMsgTxt("WW_OLD collocation frequency matrix must be a sparse double precision matrix");  
  srww_old  = mxGetPr(prhs[14]);
  irww_old  = mxGetIr(prhs[14]);
  jcww_old  = mxGetJc(prhs[14]);
  nzmax_old = mxGetNzmax(prhs[14]);
  W_OLD     = mxGetM( prhs[14] );
  if (W_OLD != W_NEW) mexErrMsgTxt("WW_OLD and WW_NEW matrices have different dimensions");
  if (mxGetN( prhs[14] ) != W_NEW) mexErrMsgTxt("WW_OLD matrix should be square");
  
  // dealing with sparse array WP_OLD 
  if ((mxIsSparse( prhs[ 15 ] ) != 1) || (mxIsDouble( prhs[ 15 ] ) != 1))
      mexErrMsgTxt("WP_OLD topic-word matrix must be a sparse double precision matrix");  
  srwp_old  = mxGetPr(prhs[15]);
  irwp_old  = mxGetIr(prhs[15]);
  jcwp_old  = mxGetJc(prhs[15]);
  if (mxGetM( prhs[15] ) != W_NEW) mexErrMsgTxt("Number of words in WP_OLD matrix does not match WW_NEW number of words");
  if (mxGetN( prhs[15] ) != T_NEW) mexErrMsgTxt("Number of topics in WP_OLD matrix does not match given number of topics");
 
  WC_OLD = mxGetPr( prhs[ 16 ]);
  if ((mxGetM( prhs[16] ) * mxGetN( prhs[16] )) != W_NEW ) mexErrMsgTxt("Number of words in WC_OLD matrix does not match WW_NEW number of words");   
  
  if (mxIsDouble( prhs[ 17 ] ) != 1) mexErrMsgTxt("WS_OLD must be double precision"); 
  WS_OLD = mxGetPr(prhs[ 17 ]);
  ntokens_old = mxGetM( prhs[ 17 ] ) * mxGetN( prhs[ 17 ] );
  
  updatecols = (int) mxGetScalar(prhs[18]);
  
  // seeding
  seedMT( 1 + SEED * 2 ); // seeding only works on uneven numbers
  
  /* allocate memory */
  z  = (int *) mxCalloc( ntokens_new , sizeof( int ));
  d  = (int *) mxCalloc( ntokens_new , sizeof( int ));
  w  = (int *) mxCalloc( ntokens_new , sizeof( int ));
  s  = (int *) mxCalloc( ntokens_new , sizeof( int ));
  c  = (int *) mxCalloc( ntokens_new , sizeof( int ));
  sc  = (int *) mxCalloc( ntokens_new , sizeof( int ));
 
  for (i=0; i<ntokens_new; i++) w[ i ] = (int) (WS_NEW[ i ] - 1); // Matlab indexing not zero based
  for (i=0; i<ntokens_new; i++) d[ i ] = (int) (DS_NEW[ i ] - 1); // Matlab indexing not zero based
  for (i=0; i<ntokens_new; i++) s[ i ] = (int) SI_NEW[ i ];
  
  order   = (int *) mxCalloc( ntokens_new , sizeof( int ));
  wp_old  = (int *) mxCalloc( T_NEW*W_NEW , sizeof( int ));
  dp_new  = (int *) mxCalloc( T_NEW*D_NEW , sizeof( int ));
  dp_new2  = (int *) mxCalloc( T_NEW*D_NEW , sizeof( int ));
  wc_new  = (int *) mxCalloc( W_NEW , sizeof( int ));
  ztot_old  = (int *) mxCalloc( T_NEW , sizeof( int ));
  wtot_new  = (int *) mxCalloc( W_NEW , sizeof( int ));
  probs  = (double *) mxCalloc( T_NEW , sizeof( double ));
  
  // FILL IN WTOT_NEW with OLD COUNTS
  for (i=0; i<ntokens_old; i++) {
     j = (int) WS_OLD[ i ] - 1;
     wtot_new[ j ]++; // calculate number of words of each type
  }
  
  // FILL IN WTOT_NEW with NEW COUNTS
  if (updatecols==1)
  for (i=0; i<ntokens_new; i++) {
     j = w[ i ];
     wtot_new[ j ]++; // calculate number of words of each type
  }
  
  // FILL IN WC WITH COUNTS FROM OLD SET
  for (i=0; i<W_NEW; i++) {
     wc_new[ i ] = (int) WC_OLD[ i ];   
  }
  
  // FILL IN wp_old and ztot_old
  for (j=0; j<T_NEW; j++) {
      startindex = *( jcwp_old + j );
      endindex   = *( jcwp_old + j + 1 ) - 1;
      
      for (i=startindex; i<=endindex; i++) {
          currentrow = *( irwp_old + i );
          count      = (int) *( srwp_old + i );
          wp_old[ j + currentrow*T_NEW ] = count;
          ztot_old[ j ] += count;
      }    
  }
  
  // FILL IN SC COUNTS 
  for (i=1; i<ntokens_new; i++) // start with second item
  {
      wi = w[ i ]; // word index
      
      count = 0;
      
      // what is the previous word?
      wipre = w[ i-1 ];
      
      // calculate how many times the current word follows the previous word IN THE OLD DOCUMENT SET    
      startindex = *( jcww_old + wipre ); // look up start and end index for column wipre in WW_OLD matrix
      endindex   = *( jcww_old + wipre + 1 ) - 1;   
      for (j=startindex; j<=endindex; j++) {
          currentrow = *( irww_old + j );
          if (currentrow == wi) count = (int) *( srww_old + j );
      }
      
      // calculate how many times the current word follows the previous word IN THE NEW DOCUMENT SET
      if (updatecols==1) {
          startindex = *( jcww_new + wipre ); // look up start and end index for column wipre in WW_NEW matrix
          endindex   = *( jcww_new + wipre + 1 ) - 1;
          for (j=startindex; j<=endindex; j++) {
              currentrow = *( irww_new + j );
              if (currentrow == wi) count += (int) *( srww_new + j ); // add the count to what was there previously !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
          }
      }
      
      sc[ i ] = count;
  }
  
  
  
  /*for (i=0; i<10; i++) {
     mexPrintf( "i=%4d w[i]=%3d d[i]=%d z[i]=%d s[i]=%d\n" , i , w[i] , d[i] , z[i] , s[i] );    
  }*/
  
  if (OUTPUT==2) {
      mexPrintf( "Running NEW DOCUMENT LDA COL Gibbs Sampler Version 1.0\n" );
      mexPrintf( "Arguments:\n" );
      mexPrintf( "\tNumber of words      W = %d\n" , W_NEW );
      mexPrintf( "\tNumber of docs       D = %d\n" , D_NEW );
      mexPrintf( "\tNumber of topics     T = %d\n" , T_NEW );
      mexPrintf( "\tNumber of iterations N = %d\n" , NITER );
      mexPrintf( "\tNumber of samples    S = %d\n" , NSAMPLE );
      mexPrintf( "\tHyperparameter   ALPHA = %4.4f\n" , ALPHA );
      mexPrintf( "\tHyperparameter    BETA = %4.4f\n" , BETA );
      mexPrintf( "\tHyperparameter  GAMMA0 = %4.4f\n" , GAMMA0 );
      mexPrintf( "\tHyperparameter  GAMMA1 = %4.4f\n" , GAMMA1 );
      mexPrintf( "\tHyperparameter  DELTA  = %4.4f\n" , DELTA );
      mexPrintf( "\tSeed number            = %d\n" , SEED );
      mexPrintf( "\tNumber of tokens       = %d\n" , ntokens_new );
      mexPrintf( "\tUpdating collocation counts with new documents? = %d\n" , updatecols );
  }
  
  /* ---------------------------------------------
     create the PROB vector
   -----------------------------------------------*/
  plhs[ 3 ] = mxCreateDoubleMatrix( ntokens_new , 1 , mxREAL );
  PROBA = mxGetPr( plhs[ 3 ] );
  
  /* run the model */  
  GibbsSamplerLDACOL( ALPHA, BETA, GAMMA0,GAMMA1,DELTA, W_NEW, T_NEW, D_NEW, NITER, NSAMPLE , OUTPUT, ntokens_new, z, d, w, s, c, sc, wp_old, dp_new, dp_new2, ztot_old, wtot_new, order, probs, wc_new , PROBA , updatecols );
  

  /* ---------------------------------------------
   convert the full DP matrix into a sparse matrix 
   -----------------------------------------------*/
  nzmaxdp_new = 0;
  for (i=0; i<D_NEW; i++) {
      for (j=0; j<T_NEW; j++)
          nzmaxdp_new += (int) ( *( dp_new2 + j + i*T_NEW )) > 0; // !!!!!!!!!!!!!!!!!! copy from dp_new2
  }
  
  if (OUTPUT==2) {
      mexPrintf( "Constructing sparse output matrix dp_new\n" );
      mexPrintf( "Number of nonzero entries for DP = %d\n" , nzmaxdp_new );
  }
  
  plhs[0] = mxCreateSparse( D_NEW,T_NEW,nzmaxdp_new,mxREAL);
  srdp_new  = mxGetPr(plhs[0]);
  irdp_new = mxGetIr(plhs[0]);
  jcdp_new = mxGetJc(plhs[0]);
  n = 0;
  for (j=0; j<T_NEW; j++) {
      *( jcdp_new + j ) = n;
      for (i=0; i<D_NEW; i++) {
          cc = (int) *( dp_new2 + i*T_NEW + j ); // !!!!!!!!!!!!!!!!!! copy from dp_new2
          if (cc >0) {
              *( srdp_new + n ) = cc;
              *( irdp_new + n ) = i;
              n++;
          }
      }
  } 
  *( jcdp_new + T_NEW ) = n;
  
  /* ---------------------------------------------
     create the C route vector
   -----------------------------------------------*/
  plhs[ 1 ] = mxCreateDoubleMatrix( ntokens_new , 1 , mxREAL );
  C = mxGetPr( plhs[ 1 ] );
  for (i=0; i<ntokens_new; i++) C[ i ] = (double) c[ i ];
  
  /* ---------------------------------------------
     create the topic assignment vector
   -----------------------------------------------*/
  plhs[ 2 ] = mxCreateDoubleMatrix( ntokens_new , 1 , mxREAL );
  Z = mxGetPr( plhs[ 2 ] );
  for (i=0; i<ntokens_new; i++) Z[ i ] = (double) z[ i ] + 1;

}
