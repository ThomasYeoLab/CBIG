#include "mex.h"
#include "cokus.cpp"

// Syntax
//   [ WP,DP,WC,C,Z ] = GibbsSamplerLDACOL( WS , DS , SI , WW , T , N , ALPHA , BETA , GAMMA0, GAMMA1 , DELTA , SEED , OUTPUT );
    
// updated for 64 bit

void GibbsSamplerLDACOL( double ALPHA, double BETA, double GAMMA0, double GAMMA1, double DELTA, int W, int T, int D, int NN, 
                         int OUTPUT, int n, int *z, int *d, int *w, int *s, int *c, int *sc, int *wp, int *dp, int *ztot, 
                         int *wtot, int *order, double *probs, double *srww, 
                         mwIndex *irww, mwIndex *jcww, 
                         int *wc, int startcond )
{
  int wi,wipre,di,si,i,ii,j,topic,count, route,rp, temp, iter, wioffset, dioffset;
  int startindex,endindex,currentrow,prevroute, ncol;
  double totprob, WBETA, WDELTA,r, max, C0, C1;

  /* random initialization */
  if (startcond==1) {
      ncol = 0;
      for (i=0; i<n; i++)
      {
          wi     = w[ i ]; // word index
          di     = d[ i ]; // document index
          si     = s[ i ]; // status of word
          route  = c[ i ];
          topic  = z[ i ]; 
          
          wtot[ wi ]++; // increment word count for this type
          
          // only update topic counts for words assigned to LDA model
          if (route == 0) {              
              wp[ wi*T + topic ]++; // increment wp count matrix
              dp[ di*T + topic ]++; // increment dp count matrix
              ztot[ topic ]++; // increment ztot matrix
          }
          
          // calculate how many times the current word follows the previous word (over the whole corpus)
          count = 0;
          if (i>0) {           
              // what is the previous word?
              wipre = w[ i-1 ];
              
              // look up start and end index for column wipre in WW matrix
              startindex = *( jcww + wipre );
              endindex   = *( jcww + wipre + 1 ) - 1;
              
              for (j=startindex; j<=endindex; j++) {
                  currentrow = *( irww + j );
                  if (currentrow == wi) count = (int) *( srww + j );
              }
          }
          sc[ i ] = count;
           
          ncol += route;
          if ((route==1) && (i>0)) {
              wipre = w[ i-1 ];
              wc[ wipre ]++; // increment the route count for this word type  
          }
      }
  } 
  
  if (startcond==0) {
      if (OUTPUT==2) mexPrintf( "Starting Random initialization\n" );
      
      ncol = 0;
      for (i=0; i<n; i++)
      {
          wi = w[ i ]; // word index
          di = d[ i ]; // document index
          si = s[ i ]; // status of word
          
          wtot[ wi ]++; // increment word count for this type
          
          topic = (int) ( (double) randomMT() * (double) T / (double) (4294967296.0 + 1.0) ); // pick a random topic 0..T-1
          
          z[ i ] = topic; // assign this word token to this topic
          wp[ wi*T + topic ]++; // increment wp count matrix
          dp[ di*T + topic ]++; // increment dp count matrix
          ztot[ topic ]++; // increment ztot matrix
          
          // calculate how many times the current word follows the previous word (over the whole corpus)
          count = 0;
          if (i>0) {
              
              // what is the previous word?
              wipre = w[ i-1 ];
              
              // look up start and end index for column wipre in WW matrix
              startindex = *( jcww + wipre );
              endindex   = *( jcww + wipre + 1 ) - 1;
              
              for (j=startindex; j<=endindex; j++) {
                  currentrow = *( irww + j );
                  if (currentrow == wi) count = (int) *( srww + j );
              }
          }
          sc[ i ] = count;
          
          // initially, set the route to 0 for topic route
          route = 0;          
          c[ i ] = route;
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
  WBETA = (double) (W*BETA);
  WDELTA = (double) (W*DELTA);
  
  for (iter=0; iter<NN; iter++) {
      if (OUTPUT >=1) {
          if ((iter % 10)==0) mexPrintf( "\tIteration %d of %d;   Number of tokens in collocation = %d\n" , iter , NN , ncol );
          if ((iter % 10)==0) mexEvalString("drawnow;");
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
              ztot[topic]--;  // substract this from counts, only if going with topic route
              wp[wioffset+topic]--;
              dp[dioffset+topic]--;

              //mexPrintf( "(1) Working on ii=%d i=%d wi=%d di=%d topic=%d wp=%d dp=%d\n" , ii , i , wi , di , topic , wp[wi+topic*W] , dp[wi+topic*D] );
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
              
              wtot[ wipre ]--;
              if (prevroute==1) {
                  wc[ wipre ]--; // if the current word token was assigned to a collocation, take this out of the count for the previous word
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
              
              wtot[ wipre ]++;
              
              // is the word part of a collocation?
              if (route==1) {
                  ncol++;
                  wc[ wipre ]++; // increment collocation count for this word
              } else
              {
                  wp[wioffset + topic ]++; // and subtract this token from topic counts
                  dp[dioffset + topic ]++;
                  ztot[topic]++;
              }
            
              c[ i ] = route;
          } else {
            // add counts to topic route
            wp[wioffset + topic ]++; // and update counts
            dp[dioffset + topic ]++;
            ztot[topic]++;   
          }
          
          //mexPrintf( "(2) Working on ii=%d i=%d wi=%d di=%d topic=%d wp=%d dp=%d\n" , ii , i , wi , di , topic , wp[wi+topic*W] , dp[wi+topic*D] );
      }
  }
  
  /*sumc = 0;
  for (i=0; i<n; i++) sumc += c[ i ];
  mexPrintf( "(3) SumC=%d  ncol=%d\n" , sumc , ncol );
   */
}

void mexFunction(int nlhs, mxArray *plhs[], int nrhs,
                 const mxArray *prhs[])
{
  double *srww,*srwp,*srdp, *probs, *DS, *WS, *SI, *WC, *WWC, *ZZ, *CIN, *ZIN;
  double ALPHA,BETA, GAMMA0, GAMMA1, DELTA;
  mwIndex *irww, *jcww, *irwp, *jcwp, *irdp, *jcdp;
  int *z,*d,*w, *s, *c, *sc, *order, *wp, *dp, *ztot, *wtot, *wc;
  int W,T,D,NN,SEED,OUTPUT, nzmax, nzmaxwp, nzmaxdp, ntokens;
  int i,j,cc,n,nt,wi,di;
  int startcond;
  
  /* Check for proper number of arguments. */
  if (nrhs < 13) {
    mexErrMsgTxt("At least 13 input arguments required");
  } else if (nlhs != 5) {
    mexErrMsgTxt("5 output arguments required");
  }
  
  //[ WP,DP,WC,C,Z ] = GibbsSamplerLDACOL( WS , DS , SI , WW , T , N , ALPHA , BETA , GAMMA0, GAMMA1 , DELTA , SEED , OUTPUT , CIN , ZIN );
  
  /* process the input arguments */
  if (mxIsDouble( prhs[ 0 ] ) != 1)
      mexErrMsgTxt("WS must be double precision"); 
  WS = mxGetPr(prhs[0]);
  
  if (mxIsDouble( prhs[ 1 ] ) != 1)
      mexErrMsgTxt("DS must be double precision"); 
  DS = mxGetPr(prhs[1]);
  
  if (mxIsDouble( prhs[ 2 ] ) != 1)
      mexErrMsgTxt("SI must be double precision");
  SI = mxGetPr(prhs[2]);
  
  if ((mxIsSparse( prhs[ 3 ] ) != 1) || (mxIsDouble( prhs[ 3 ] ) != 1))
      mexErrMsgTxt("WW collocation matrix must be a sparse double precision matrix");

  /* dealing with sparse array WW */
  srww = mxGetPr(prhs[3]);
  irww  = mxGetIr(prhs[3]);
  jcww  = mxGetJc(prhs[3]);
  nzmax= (int) mxGetNzmax(prhs[3]);
  
  W    = (int) mxGetM( prhs[3] );
  ntokens = (int) mxGetM( prhs[ 0 ] ) * (int) mxGetN( prhs[ 0 ] );
  
  D    = 0;
  for (i=0; i<ntokens; i++) {
      if (DS[ i ] > D) D = (int) DS[ i ];
  }
   
  T    = (int) mxGetScalar(prhs[4]);
  if (T<=0) mexErrMsgTxt("Number of topics must be greater than zero");
  
  NN    = (int) mxGetScalar(prhs[5]);
  if (NN<0) mexErrMsgTxt("Number of iterations must be greater than zero");
  
  ALPHA = (double) mxGetScalar(prhs[6]);
  if (ALPHA<=0) mexErrMsgTxt("ALPHA must be greater than zero");
  
  BETA = (double) mxGetScalar(prhs[7]);
  if (BETA<=0) mexErrMsgTxt("BETA must be greater than zero");
  
  GAMMA0 = (double) mxGetScalar(prhs[8]);
  if (GAMMA0<=0) mexErrMsgTxt("GAMMA0 must be greater than zero");
  
  GAMMA1 = (double) mxGetScalar(prhs[9]);
  if (GAMMA1<=0) mexErrMsgTxt("GAMMA1 must be greater than zero");
  
  DELTA = (double) mxGetScalar(prhs[10]);
  if (DELTA<=0) mexErrMsgTxt("DELTA must be greater than zero");
  
  SEED = (int) mxGetScalar(prhs[11]);
  // set the seed of the random number generator
  
  OUTPUT = (int) mxGetScalar(prhs[12]);
  
  // assume that we start a new chain
  startcond = 0;
  
  if (nrhs > 13) {
      startcond = 1;
      CIN = mxGetPr( prhs[ 13 ]);
      ZIN = mxGetPr( prhs[ 14 ]);
  }
  
  // seeding
  seedMT( 1 + SEED * 2 ); // seeding only works on uneven numbers
  
  /* allocate memory */
  z  = (int *) mxCalloc( ntokens , sizeof( int ));
  d  = (int *) mxCalloc( ntokens , sizeof( int ));
  w  = (int *) mxCalloc( ntokens , sizeof( int ));
  s  = (int *) mxCalloc( ntokens , sizeof( int ));
  c  = (int *) mxCalloc( ntokens , sizeof( int ));
  sc  = (int *) mxCalloc( ntokens , sizeof( int ));
 
  if (startcond==1) {
      for (i=0; i<ntokens; i++) c[ i ] = (int) CIN[ i ];
      for (i=0; i<ntokens; i++) z[ i ] = (int) ZIN[ i ] - (int) 1;
  }
  
  for (i=0; i<ntokens; i++) w[ i ] = (int) (WS[ i ] - 1); // Matlab indexing not zero based
  for (i=0; i<ntokens; i++) d[ i ] = (int) (DS[ i ] - 1); // Matlab indexing not zero based
  for (i=0; i<ntokens; i++) s[ i ] = (int) SI[ i ];
  
  order  = (int *) mxCalloc( ntokens , sizeof( int ));
  wp  = (int *) mxCalloc( T*W , sizeof( int ));
  dp  = (int *) mxCalloc( T*D , sizeof( int ));
  wc  = (int *) mxCalloc( W , sizeof( int ));
  ztot  = (int *) mxCalloc( T , sizeof( int ));
  wtot  = (int *) mxCalloc( W , sizeof( int ));
  probs  = (double *) mxCalloc( T , sizeof( double ));
  
  n = ntokens;
  
  /*for (i=0; i<10; i++) {
     mexPrintf( "i=%4d w[i]=%3d d[i]=%d z[i]=%d s[i]=%d\n" , i , w[i] , d[i] , z[i] , s[i] );    
  }*/
  
  if (OUTPUT==2) {
      mexPrintf( "Running LDA COL Gibbs Sampler Version 1.0\n" );
      if (startcond==1) mexPrintf( "Starting with C and Z vector given as input\n" );
      mexPrintf( "Arguments:\n" );
      mexPrintf( "\tNumber of words      W = %d\n" , W );
      mexPrintf( "\tNumber of docs       D = %d\n" , D );
      mexPrintf( "\tNumber of topics     T = %d\n" , T );
      mexPrintf( "\tNumber of iterations N = %d\n" , NN );
      mexPrintf( "\tHyperparameter   ALPHA = %4.4f\n" , ALPHA );
      mexPrintf( "\tHyperparameter    BETA = %4.4f\n" , BETA );
      mexPrintf( "\tHyperparameter  GAMMA0 = %4.4f\n" , GAMMA0 );
      mexPrintf( "\tHyperparameter  GAMMA1 = %4.4f\n" , GAMMA1 );
      mexPrintf( "\tHyperparameter  DELTA  = %4.4f\n" , DELTA );
      mexPrintf( "\tSeed number            = %d\n" , SEED );
      mexPrintf( "\tNumber of tokens       = %d\n" , ntokens );
      mexPrintf( "Internal Memory Allocation\n" );
      mexPrintf( "\tw,d,z,s,c,order indices combined = %d bytes\n" , 7 * sizeof( int) * ntokens );
      mexPrintf( "\twp (full) matrix = %d bytes\n" , sizeof( int ) * W * T  );
      mexPrintf( "\tdp (full) matrix = %d bytes\n" , sizeof( int ) * D * T  );
      //mexPrintf( "Checking: sizeof(int)=%d sizeof(long)=%d sizeof(double)=%d\n" , sizeof(int) , sizeof(long) , sizeof(double));
  }
  
  /* run the model */  
  GibbsSamplerLDACOL( ALPHA, BETA, GAMMA0,GAMMA1,DELTA, W, T, D, NN, OUTPUT, n, z, d, w, s, c, sc, wp, dp, ztot, wtot, order, probs , srww, irww, jcww, wc, startcond );
  
  /* ---------------------------------------------
   convert the full wp matrix into a sparse matrix 
   -----------------------------------------------*/
  nzmaxwp = 0;
  for (i=0; i<W; i++) {
     for (j=0; j<T; j++)
         nzmaxwp += (int) ( *( wp + j + i*T )) > 0;
  }
  
  if (OUTPUT==2) {
      mexPrintf( "Constructing sparse output matrix wp\n" );
      mexPrintf( "Number of nonzero entries for WP = %d\n" , nzmaxwp );
  }
  
  plhs[0] = mxCreateSparse( W,T,nzmaxwp,mxREAL);
  srwp  = mxGetPr(plhs[0]);
  irwp = mxGetIr(plhs[0]);
  jcwp = mxGetJc(plhs[0]);
  
  n = 0;
  for (j=0; j<T; j++) {
      *( jcwp + j ) = n;
      for (i=0; i<W; i++) {
         cc = (int) *( wp + i*T + j );
         if (cc >0) {
             *( srwp + n ) = cc;
             *( irwp + n ) = i;
             n++;
         }
      }    
  }
  
  *( jcwp + T ) = n;
 
  /* ---------------------------------------------
   convert the full DP matrix into a sparse matrix 
   -----------------------------------------------*/
  nzmaxdp = 0;
  for (i=0; i<D; i++) {
      for (j=0; j<T; j++)
          nzmaxdp += (int) ( *( dp + j + i*T )) > 0;
  }
  
  if (OUTPUT==2) {
      mexPrintf( "Constructing sparse output matrix dp\n" );
      mexPrintf( "Number of nonzero entries for DP = %d\n" , nzmaxdp );
  }
  
  plhs[1] = mxCreateSparse( D,T,nzmaxdp,mxREAL);
  srdp  = mxGetPr(plhs[1]);
  irdp = mxGetIr(plhs[1]);
  jcdp = mxGetJc(plhs[1]);
  
  n = 0;
  for (j=0; j<T; j++) {
      *( jcdp + j ) = n;
      for (i=0; i<D; i++) {
          cc = (int) *( dp + i*T + j );
          if (cc >0) {
              *( srdp + n ) = cc;
              *( irdp + n ) = i;
              n++;
          }
      }
  }
  
  *( jcdp + T ) = n;

  /* ---------------------------------------------
     create the WC count matrix
   -----------------------------------------------*/
  plhs[ 2 ] = mxCreateDoubleMatrix( W , 1 , mxREAL );
  WC = mxGetPr( plhs[ 2 ] );
  for (i=0; i<W; i++) WC[ i ] = (double) wc[ i ];
  
  /* ---------------------------------------------
     create the C route vector
   -----------------------------------------------*/
  plhs[ 3 ] = mxCreateDoubleMatrix( ntokens , 1 , mxREAL );
  WWC = mxGetPr( plhs[ 3 ] );
  for (i=0; i<ntokens; i++) WWC[ i ] = (double) c[ i ];
  
  /* ---------------------------------------------
     create the topic assignment vector
   -----------------------------------------------*/
  plhs[ 4 ] = mxCreateDoubleMatrix( ntokens , 1 , mxREAL );
  ZZ = mxGetPr( plhs[ 4 ] );
  for (i=0; i<ntokens; i++) ZZ[ i ] = (double) z[ i ] + 1;
}
