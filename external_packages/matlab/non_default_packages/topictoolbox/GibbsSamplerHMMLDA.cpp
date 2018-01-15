#include "mex.h"
#include "cokus.cpp"

// Syntax
//   [ WP , DP , MP , Z , X ] = GibbsSamplerHMMLDA( WS , DS , T , S , N , ALPHA , BETA , GAMMA , SEED , OUTPUT )

// Syntax
//   [ WP , DP , MP , Z , X ] = GibbsSamplerHMMLDA( WS , DS , T , S , N , ALPHA , BETA , GAMMA , SEED , OUTPUT , ZIN , XIN )

// updated for 64 bit

void GibbsSamplerHMMLDA( double ALPHA, double BETA, double GAMMA, int W, int T, int S, int D, int NN, int OUTPUT, int n, int *z, int *x, double *d, double *w, int *wp, int *mp, int *dp, int *ztot, int *mtot, int *stot, int *sp, int *first, int *second, int *third, double *probs, int startcond )
{
  int     i,j,current,prev,preprev, dioffset,wioffset; 
  int     S2,S3,di,wi,xi,topic, iter;
  double  TALPHA,WBETA,SGAMMA,max,best,totprob,r; 

  // some constants
  S2 = S * S;
  S3 = S * S * S;
  WBETA = (double) W * BETA;
  TALPHA = (double) T * ALPHA;
  SGAMMA = (double) S * GAMMA;
  
  //mexPrintf( "W * BETA = %5.5f\n" , (double) (W * BETA) );
  //mexPrintf( "W * BETA = %5.5f\n" , (double) W * BETA );
  //mexErrMsgTxt( "exit..\n" );
  
  if (startcond == 1) {
      // Starting from a previous state specified by ZIN and XIN\n" );
      
      current = 0;
      prev = 0;
      preprev = 0;
      first[0] = 0;
      second[0] = 0;
      third[0] = 0;
      
      // why do we start at i=1 ?????????????????????????????????????
      for (i = 1; i < n; i++) {
          wi = (int) w[ i ] - 1;
          
          if (wi == -1) {
              // sentence marker
              sp[ preprev*S3 + prev*S2 + current*S + 0 ]++;
              first[i]  = current;
              second[i] = prev;
              third[i]  = preprev;
              preprev   = 0;
              prev      = 0;
              current   = 0;
          } else {
              di = (int) d[ i ] - 1;
              dioffset = di*T;
              wioffset = wi*T;
              
              if (x[i] == 1) {
                  topic = z[ i ];
                  wp[ wioffset + topic]++;
                  dp[ dioffset + topic]++;
                  ztot[ topic ]++;
              } else {
                  mp[ wi*S + x[i] ]++;
                  
                  //mexPrintf( "wi=%d xi=%d mp[wi*S+xi]=%d\n" , wi , xi , mp[wi*S+x[i]] );
              }
              xi = x[ i ];
              mtot[ xi ]++;
              stot[ prev*S2 + current*S + xi ]++;
              sp[ preprev*S3 + prev*S2 + current*S + xi ]++;
              first[i] = current;
              second[i] = prev;
              third[i] = preprev;
              preprev = prev;
              prev = current;
              current = xi;
          }
      }
  }
  
  if (startcond == 0) {
      if (OUTPUT==2) mexPrintf( "Starting Random initialization\n" );
      
      // online (Anderson) initialization
      current = 0;
      prev = 0;
      preprev = 0;
      first[0] = 0;
      second[0] = 0;
      third[0] = 0;
      
      for (i = 1; i < n; i++) {
          wi = (int) w[ i ] - 1;
          
          if (wi == -1) {
              // sentence marker
              sp[ preprev*S3 + prev*S2 + current*S + 0 ]++;
              first[i]  = current;
              second[i] = prev;
              third[i]  = preprev;
              preprev   = 0;
              prev      = 0;
              current   = 0;
          } else {
              di = (int) d[ i ] - 1;
              dioffset = di*T;
              wioffset = wi*T;
              
              max = 0; best = 0; totprob = 0;
              x[i] = (int) (((double) randomMT() / (double) (4294967296.0)) > 0.5 );
              
              for (j = 0; j < T; j++) {
                  probs[j] = (double) dp[ dioffset + j ] + ALPHA;
                  if (x[i]==1) {
                      probs[j] *= ( (double) wp[wioffset + j] + BETA)/( (double) ztot[j]+WBETA);
                  }
                  totprob += probs[j];
              }
              r = ((double) randomMT() / (double) (4294967296.0))*totprob;
              max = probs[0];
              j = 0;
              while (r>max) {
                  j++;
                  max += probs[j];
              }
              z[i] = j;
              max = 0; best = 0;
              if (x[i]==0) {
                  probs[1] = ((double) wp[ wioffset + j] + BETA) / ( (double) ztot[j]+WBETA)
                  *( (double) sp[ preprev*S3 + prev*S2 + current*S + 1]+GAMMA);
                  totprob = probs[1];
                  for (j = 2; j < S; j++) {
                      probs[j] = ( (double) mp[ wi*S + j ] + BETA ) / ( (double) mtot[j] + WBETA)
                      *( (double) sp[ preprev*S3 + prev*S2 + current*S + j ]+GAMMA);
                      totprob += probs[j];
                  }
                  r = ((double) randomMT() / (double) (4294967296.0))*totprob;
                  max = probs[1];
                  j = 1;
                  while (r>max) {
                      j++;
                      max += probs[j];
                  }
                  x[i] = j;
              }
              if (x[i] == 1) {
                  topic = z[ i ];
                  wp[ wioffset + topic]++;
                  dp[ dioffset + topic]++;
                  ztot[ topic ]++;
              } else {
                  mp[ wi*S + x[i] ]++;
                  
                  //mexPrintf( "wi=%d xi=%d mp[wi*S+xi]=%d\n" , wi , xi , mp[wi*S+x[i]] );
              }
              xi = x[ i ];
              mtot[ xi ]++;
              stot[ prev*S2 + current*S + xi ]++;
              sp[ preprev*S3 + prev*S2 + current*S + xi ]++;
              first[i] = current;
              second[i] = prev;
              third[i] = preprev;
              preprev = prev;
              prev = current;
              current = xi;
          }
      }
  }
    
  /*mexPrintf("Initial assignments:\n\n");
  tot = 0;
  for (i = 0; i < T; i++) {
    tot += ztot[i];
    mexPrintf(" T%d=%d", i , ztot[i]);
  }
  mexPrintf("\n");
  mexPrintf("LDA TOTAL: %d\n", tot);
  tot = 0;
  for (i = 0; i < S; i++) {
    tot += mtot[i];
    mexPrintf(" S%d=%d", i,mtot[i]);
  }
  mexPrintf("\n");
  mexPrintf("HMM TOTAL: %d\n\n", tot);
  */
  
  probs[0] = 0;
  
  //for (i=0; i<n; i++) mexPrintf( "i=%4d w[i]=%3d d[i]=%d z[i]=%d x[i]=%d\n" , i , (int) w[i] , (int) d[i] , z[i] , x[i] );  
  //mexPrintf( "\n\n" );
  
  for (iter=0; iter<NN; iter++) {
      if (OUTPUT >=1) {
          if ((iter % 10)==0) mexPrintf( "\tIteration %d of %d\n" , iter , NN );
          if ((iter % 10)==0) mexEvalString("drawnow;");
      }
      
      current = 0;
      prev = 0;
      preprev = 0;
      for (i = 1; i < n-3; i++)
      {
          //mexPrintf( "i=%4d w[i]=%3d d[i]=%d z[i]=%d x[i]=%d\n" , i , (int) w[i] , (int) d[i] , z[i] , x[i] );  
          
          wi = (int) w[ i ] - 1; 
          
          if (wi == -1) // sentence marker
          {
              sp[third[i]*S3+second[i]*S2+first[i]*S+0]--;
              sp[preprev*S3+prev*S2+current*S+0]++;
              first[i] = current;
              second[i] = prev;
              third[i] = preprev;
              current = 0;
              prev = 0;
              preprev = 0;
          } else
          {
              di = (int) d[ i ] - 1; 
              dioffset = di*T;
              wioffset = wi*T;
              
              sp[third[i]*S3+second[i]*S2+first[i]*S+x[i]]--;
              if (x[i] == 1)
              {
                  topic = z[ i ];
                  wp[ wioffset + topic ]--;
                  dp[ dioffset + topic ]--;
                  ztot[ topic ]--;
              } else
              {
                  mp[ wi*S + x[i] ]--;
              }
              mtot[x[i]]--;
              stot[second[i]*S2+first[i]*S+x[i]]--;
              max = 0; best = 0; totprob = 0;
              for (j = 0; j < T; j++)
              {
                  probs[j] = (double) dp[ dioffset + j ]+ALPHA;
                  if (x[i] == 1)
                  {
                      probs[j] *= ((double) wp[ wioffset + j ]+BETA)/((double) ztot[j]+WBETA);
                  }
                  totprob += probs[j];
              }
              r = ((double) randomMT() / (double) (4294967296.0))*totprob;
              max = probs[0];
              j = 0;
              while (r>max)
              {
                  j++;
                  max += probs[j];
              }
              z[i] = j;
              probs[1]=((double) wp[wioffset+j]+BETA)/((double) ztot[j]+WBETA)
              *((double) sp[preprev*S3+prev*S2+current*S+1]+GAMMA)
              *((double) sp[prev*S3+current*S2+1*S+x[i+1]]+GAMMA)/((double) stot[prev*S2+current*S+1]+SGAMMA)
              *((double) sp[current*S3+1*S2+x[i+1]*S+x[i+2]]+GAMMA)/((double) stot[current*S2+1*S+x[i+1]]+SGAMMA)
              *((double) sp[1*S3+x[i+1]*S2+x[i+2]*S+x[i+3]]+GAMMA)/((double) stot[1*S2+x[i+1]*S+x[i+2]]+SGAMMA);
              totprob = probs[1];
              for (j = 2; j < S; j++)
              {
                  probs[j]=((double) mp[wi*S+j]+BETA)/((double) mtot[j]+WBETA)
                  *((double) sp[preprev*S3+prev*S2+current*S+j]+GAMMA)
                  *((double) sp[prev*S3+current*S2+j*S+x[i+1]]+GAMMA)/((double) stot[prev*S2+current*S+j]+SGAMMA)
                  *((double) sp[current*S3+j*S2+x[i+1]*S+x[i+2]]+GAMMA)/((double) stot[current*S2+j*S+x[i+1]]+SGAMMA)
                  *((double) sp[j*S3+x[i+1]*S2+x[i+2]*S+x[i+3]]+GAMMA)/((double) stot[j*S2+x[i+1]*S+x[i+2]]+SGAMMA);
                  totprob += probs[j];
              }
              r = ((double) randomMT() / (double) (4294967296.0))*totprob;
              max = probs[1];
              j = 1;
              while (r>max)
              {
                  j++;
                  max += probs[j];
              }
              x[i] = j;
              sp[preprev*S3+prev*S2+current*S+x[i]]++;
              if (x[i] == 1)
              {
                  topic = z[ i ];
                  wp[wioffset+topic]++;
                  dp[dioffset+topic]++;
                  ztot[topic]++;
              } else
              {
                  mp[wi*S+x[i]]++;
              }
              mtot[x[i]]++;
              stot[prev*S2+current*S+x[i]]++;
              first[i] = current;
              second[i] = prev;
              third[i] = preprev;
              preprev = prev;
              prev = current;
              current = x[i];
          }
      }
      

      //for (i=0; i<n; i++) mexPrintf( "i=%4d w[i]=%3d d[i]=%d z[i]=%d x[i]=%d\n" , i , (int) w[i] , (int) d[i] , z[i] , x[i] );  
  }
}

void mexFunction(int nlhs, mxArray *plhs[], int nrhs,
                 const mxArray *prhs[])
{
  double *srwp, *srdp, *srmp, *probs, *d, *w, *ZIN, *XIN, *Z, *X;
  double ALPHA,BETA, GAMMA;
  mwIndex *irwp, *jcwp, *irdp, *jcdp, *irmp, *jcmp;
  int *z,*x, *wp, *mp, *dp, *ztot, *mtot, *stot, *sp, *first, *second, *third;
  int W,T,S,SINPUT,D,NN,SEED,OUTPUT, nzmax, nzmaxwp,nzmaxmp, nzmaxdp, ntokens;
  int NWS,NDS,i,j,c,n,wi,di,i1,i2,i3,i4,S2,S3, startcond;
  
  // Syntax
  //   [ WP , DP , MP , Z , X ] = GibbsSamplerHMMLDA( WS , DS , T , S , N , ALPHA , BETA , GAMMA , SEED , OUTPUT , ZIN , XIN )

  /* Check for proper number of arguments. */
  if (nrhs < 10) {
    mexErrMsgTxt("At least 10 input arguments required");
  } else if (nlhs != 5) {
    mexErrMsgTxt("5 output arguments required");
  }
  
  startcond = 0;
  if (nrhs > 10) startcond = 1;
  
  if (sizeof( int ) != sizeof( mxINT32_CLASS ))
    mexErrMsgTxt("Problem with internal integer representation -- contact programmer" );
  
  w   = mxGetPr( prhs[0] );
  NWS = mxGetN( prhs[0] );
  
  d   = mxGetPr( prhs[1] );
  NDS = mxGetN( prhs[1] );
  
  if ((mxIsSparse( prhs[ 0 ] ) == 1) || (mxIsDouble( prhs[ 0 ] ) != 1) || (mxGetM( prhs[0]) != 1))
      mexErrMsgTxt("WS input stream must be a one-dimensional, non-sparse, double precision vector of word indices");
 
  if ((mxIsSparse( prhs[ 1 ] ) == 1) || (mxIsDouble( prhs[ 1 ] ) != 1) || (mxGetM( prhs[1]) != 1))
      mexErrMsgTxt("DS input stream must be a one-dimensional, non-sparse, double precision vector of document indices");
  
  if (NWS != NDS)
      mexErrMsgTxt("WS and DS input streams must have equal dimensions" );
    
  ntokens = NWS;
  
  T    = (int) mxGetScalar(prhs[2]);
  if (T<=0) mexErrMsgTxt("Number of topics must be greater than zero");
  
  SINPUT    = (int) mxGetScalar(prhs[3]);
  if (SINPUT<=0) mexErrMsgTxt("Number of syntactic states must be greater than zero");
  
  // need one syntactic state for sentence start 
  // need one syntactic state for topic state
  S = SINPUT + 2;
  
  NN    = (int) mxGetScalar(prhs[4]);
  if (NN<0) mexErrMsgTxt("Number of iterations must be greater than zero");
  
  ALPHA = (double) mxGetScalar(prhs[5]);
  if (ALPHA<=0) mexErrMsgTxt("ALPHA must be greater than zero");
  
  BETA = (double) mxGetScalar(prhs[6]);
  if (BETA<=0) mexErrMsgTxt("BETA must be greater than zero");
  
  GAMMA = (double) mxGetScalar(prhs[7]);
  if (GAMMA<=0) mexErrMsgTxt("GAMMA must be greater than zero");
  
  SEED = (int) mxGetScalar(prhs[8]);
  // set the seed of the random number generator
  
  OUTPUT = (int) mxGetScalar(prhs[9]);
  
  if (startcond==1) {
      if ((mxGetN( prhs[ 10 ] )*mxGetM( prhs[ 10 ] )) != ntokens) mexErrMsgTxt( "The ZIN vector should have have the same number of tokens as WS" );
      if ((mxGetN( prhs[ 11 ] )*mxGetM( prhs[ 11 ] )) != ntokens) mexErrMsgTxt( "The ZIN vector should have have the same number of tokens as WS" );
      
      ZIN = mxGetPr( prhs[ 10 ]);
      XIN = mxGetPr( prhs[ 11 ]);
  }
  
  // seeding
  seedMT( 1 + SEED * 2 ); // seeding only works on uneven numbers
  
  W       = 0;
  D       = 0;
  for (i=0; i<ntokens; i++) {
      // in Matlab WI=1...W word occurence, WI=0 sentence marker  
      wi = (int) w[ i ] - 1; // word indices are not zero based in Matlab
      di = (int) d[ i ] - 1; // document indices are not zero based in Matlab
      
      //mexPrintf( "i=%d wi=%d di=%d\n" , i , wi , di );
      
      if (wi > W) W = wi;
      if (di > D) D = di;
      
      if (wi < -1) mexErrMsgTxt("Unrecognized code in word stream (<-2)");
      if (di <  0) mexErrMsgTxt("Unrecognized code in document stream (<0)");
  }
  W = W + 1;
  D = D + 1;
 
  if ( wi != -1 )
       mexErrMsgTxt("Word stream should end with sentence marker");
  
  n = ntokens;
  //mexPrintf( "W=%d D=%d n=%d w[n-1]=%d\n" , W , D , n , (int) w[ n-1 ] );
  
  
  // allocate memory 
  z  = (int *) mxCalloc( ntokens , sizeof( int ));
  x  = (int *) mxCalloc( ntokens , sizeof( int ));
  
  if (startcond==1) {
     for (i=0; i<ntokens; i++) {
        z[ i ] = (int) ZIN[ i ] - 1;
        x[ i ] = (int) XIN[ i ] - 1;
     }
  }
  
  first   = (int *) mxCalloc( ntokens , sizeof( int ));
  second  = (int *) mxCalloc( ntokens , sizeof( int ));
  third   = (int *) mxCalloc( ntokens , sizeof( int ));
  mp  = (int *) mxCalloc( W*S , sizeof( int ));
  wp  = (int *) mxCalloc( T*W , sizeof( int ));
  dp  = (int *) mxCalloc( T*D , sizeof( int ));
  ztot  = (int *) mxCalloc( T , sizeof( int ));
  mtot  = (int *) mxCalloc( S , sizeof( int ));
  stot  = (int *) mxCalloc( S * S * S , sizeof( int ));
  sp    = (int *) mxCalloc( S * S * S * S , sizeof( int ));
  probs  = (double *) mxCalloc( T+S , sizeof( double ));
  
  
  //for (i=0; i<n; i++) mexPrintf( "i=%4d w[i]=%3d d[i]=%d\n" , i , w[i] , d[i] , z[i] );    
  
  if (OUTPUT==2) {
      mexPrintf( "Running HMM-LDA Gibbs Sampler Version 1.0\n" );
      if (startcond==1) mexPrintf( "Starting sampler from previous state\n" );
      mexPrintf( "Arguments:\n" );
      mexPrintf( "\tNumber of words            W = %d\n" , W );
      mexPrintf( "\tNumber of docs             D = %d\n" , D );
      mexPrintf( "\tNumber of topics           T = %d\n" , T );
      mexPrintf( "\tNumber of syntactic states S = %d\n" , S );
      mexPrintf( "\tNumber of iterations       N = %d\n" , NN );
      mexPrintf( "\tHyperparameter         ALPHA = %4.4f\n" , ALPHA );
      mexPrintf( "\tHyperparameter          BETA = %4.4f\n" , BETA );
      mexPrintf( "\tHyperparameter         GAMMA = %4.4f\n" , GAMMA );
      mexPrintf( "\tSeed number             SEED = %d\n" , SEED );
      mexPrintf( "Properties of WS stream\n" );
      mexPrintf( "\tNumber of tokens       NZ = %d\n" , ntokens );
      mexPrintf( "Internal Memory Allocation\n" );
      mexPrintf( "\tz,x,first,second,third indices combined = %d bytes\n" , 5 * sizeof( int) * ntokens );
      mexPrintf( "\twp   matrix = %d bytes\n" , sizeof( int ) * W * T  );
      mexPrintf( "\tmp   matrix = %d bytes\n" , sizeof( int ) * W * S  );
      mexPrintf( "\tdp   matrix = %d bytes\n" , sizeof( int ) * D * T  );
      mexPrintf( "\tstot matrix = %d bytes\n" , sizeof( int ) * S * S * S );
      mexPrintf( "\tsp   matrix = %d bytes\n" , sizeof( int ) * S * S * S * S );
      //mexPrintf( "Checking: sizeof(int)=%d sizeof(long)=%d sizeof(double)=%d\n" , sizeof(int) , sizeof(long) , sizeof(double));
  }
  
  // run the model 
  GibbsSamplerHMMLDA( ALPHA, BETA, GAMMA, W, T, S , D, NN, OUTPUT, n, z, x, d, w, wp, mp, dp, ztot, mtot, stot , sp , first,second,third, probs, startcond );
  
  // convert the full wp matrix into a sparse matrix 
  nzmaxwp = 0;
  for (i=0; i<W; i++) {
     for (j=0; j<T; j++)
         nzmaxwp += (int) ( *( wp + j + i*T )) > 0;
  }
  if (OUTPUT==2) mexPrintf( "Constructing sparse output matrix WP  nnz=%d\n" , nzmaxwp );

  plhs[0] = mxCreateSparse( W,T,nzmaxwp,mxREAL);
  srwp  = mxGetPr(plhs[0]);
  irwp = mxGetIr(plhs[0]);
  jcwp = mxGetJc(plhs[0]); 
  n = 0;
  for (j=0; j<T; j++) {
      *( jcwp + j ) = n;
      for (i=0; i<W; i++) {
         c = (int) *( wp + i*T + j );
         if (c >0) {
             *( srwp + n ) = c;
             *( irwp + n ) = i;
             n++;
         }
      }    
  }
  
  *( jcwp + T ) = n;    
   
  // processing DP as sparse output matrix 
  nzmaxdp = 0;
  for (i=0; i<D; i++) {
      for (j=0; j<T; j++)
          nzmaxdp += (int) ( *( dp + j + i*T )) > 0;
  }  
  if (OUTPUT==2) mexPrintf( "Constructing sparse output matrix DP  nnz=%d\n" , nzmaxdp );
  
  plhs[1] = mxCreateSparse( D,T,nzmaxdp,mxREAL);
  srdp  = mxGetPr(plhs[1]);
  irdp = mxGetIr(plhs[1]);
  jcdp = mxGetJc(plhs[1]);
  
  n = 0;
  for (j=0; j<T; j++) {
      *( jcdp + j ) = n;
      for (i=0; i<D; i++) {
          c = (int) *( dp + i*T + j );
          if (c >0) {
              *( srdp + n ) = c;
              *( irdp + n ) = i;
              n++;
          }
      }
  }
  *( jcdp + T ) = n;
 
  // processing MP as sparse output matrix 
  nzmaxmp = 0;
  for (i=0; i<W; i++) {
      for (j=0; j<S; j++) {
          nzmaxmp += (int) ( mp[ j + i*S ] > 0 );
      }    
  }  
  if (OUTPUT==2) mexPrintf( "Constructing sparse output matrix MP nnz=%d\n" , nzmaxmp );
  
  plhs[2] = mxCreateSparse( W,SINPUT,nzmaxmp,mxREAL);
  srmp  = mxGetPr(plhs[2]);
  irmp = mxGetIr(plhs[2]);
  jcmp = mxGetJc(plhs[2]);
  
  n = 0;
  for (j=0; j<SINPUT; j++) {
      *( jcmp + j ) = n;
      for (i=0; i<W; i++) {
          c = mp[ i*S + (j+2) ];
          if (c >0) {
              *( srmp + n ) = c;
              *( irmp + n ) = i;
              n++;
          }
      }
  }
  *( jcmp + SINPUT ) = n;
  
  
  
  /* ---------------------------------------------
     create the topic assignment vector
   -----------------------------------------------*/
  plhs[ 3 ] = mxCreateDoubleMatrix( ntokens , 1 , mxREAL );
  Z = mxGetPr( plhs[ 3 ] );
  for (i=0; i<ntokens; i++) Z[ i ] = (double) z[ i ] + 1;

   /* ---------------------------------------------
     create the HMM state assignment vector
   -----------------------------------------------*/
  plhs[ 4 ] = mxCreateDoubleMatrix( ntokens , 1 , mxREAL );
  X = mxGetPr( plhs[ 4 ] );
  for (i=0; i<ntokens; i++) X[ i ] = (double) x[ i ] + 1;
}
