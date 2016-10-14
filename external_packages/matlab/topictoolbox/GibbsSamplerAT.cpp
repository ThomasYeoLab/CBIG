#include "mex.h"
#include "cokus.cpp"

// Syntax
//   [ WP , AT , Z , X ] = GibbsSamplerAT( WS , DS , AD , T , N , ALPHA , BETA , SEED , OUTPUT )

// Syntax
//   [ WP , AT , Z , X ] = GibbsSamplerAT( WS , DS , AD , T , N , ALPHA , BETA , SEED , OUTPUT , ZIN , XIN )

// Updated for 64 bit compiler

int NAMAX = 100; // maximum number of authors on a document

void GibbsSamplerAT( double ALPHA, double BETA, int W, int T, int D, 
                    int NN, int OUTPUT, int n, 
                    int *z, int *d, int *w, int *x, int *wp, int *at, 
                     int *ztot, int *atot, int *order, double *probs, mwIndex *irad, mwIndex *jcad, int startcond )
{
  int wi,di,i,ii,j,k, topic, rp, temp, iter, wioffset, aioffset, i_start, i_end, i_pick, author;
  int nauthors, kk;
  double totprob, WBETA, KALPHA, r, max;

  if (startcond==1) {
      /* start from previous state */
      for (i=0; i<n; i++)
      {
          wi     = w[ i ];
          di     = d[ i ];
          topic  = z[ i ];

          wp[ wi*T + topic ]++; // increment wp count matrix
          ztot[ topic ]++; // increment ztot matrix

          // i_start and i_end are the start and end indices of the sparse AD matrix for document di
          i_start = jcad[ di   ];
          i_end   = jcad[ di+1 ];
          nauthors = i_end - i_start;

          if (nauthors==1)
              author = (int) irad[ i_start ];
          else  {
              // pick a random number between i_start and i_end-1
              i_pick = i_start + (int) ( (double) randomMT() * (double) nauthors / (double) (4294967296.0 + 1.0) );

              // which author does this correspond to?
              author = (int) irad[ i_pick ];
          }

          author = x[ i ];

          // update counts for this author
          at[ author*T + topic ]++; // increment at count matrix
          atot[ author ]++;

          //mexPrintf( "For doc=%d nauthors=%d i_start=%d i_end=%d i_pick=%d author=%d\n" , di , nauthors , i_start , i_end , i_pick , author );
      }
  }

  if (startcond==0) {
      /* random initialization */
      if (OUTPUT==2) mexPrintf( "Starting Random initialization\n" );
      for (i=0; i<n; i++)
      {
          wi = w[ i ];
          di = d[ i ];
          // pick a random topic 0..T-1
          topic = (int) ( (double) randomMT() * (double) T / (double) (4294967296.0 + 1.0) );
          z[ i ] = topic; // assign this word token to this topic
          wp[ wi*T + topic ]++; // increment wp count matrix
          ztot[ topic ]++; // increment ztot matrix

          // i_start and i_end are the start and end indices of the sparse AD matrix for document di
          i_start = jcad[ di   ];
          i_end   = jcad[ di+1 ];
          nauthors = i_end - i_start;

          if (nauthors==1)
              author = (int) irad[ i_start ];
          else  {
              // pick a random number between i_start and i_end-1
              i_pick = i_start + (int) ( (double) randomMT() * (double) nauthors / (double) (4294967296.0 + 1.0) );

              // which author does this correspond to?
              author = (int) irad[ i_pick ];
          }

          x[ i ] = author;

          // update counts for this author
          at[ author*T + topic ]++; // increment at count matrix
          atot[ author ]++;

          //mexPrintf( "For doc=%d nauthors=%d i_start=%d i_end=%d i_pick=%d author=%d\n" , di , nauthors , i_start , i_end , i_pick , author );
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
  KALPHA = (double) T * (double) ALPHA;

  for (iter=0; iter<NN; iter++) {
      if (OUTPUT >=1) {
          if ((iter % 10)==0) mexPrintf( "\tIteration %d of %d\n" , iter , NN );
          if ((iter % 10)==0) mexEvalString("drawnow;");
      }
      for (ii = 0; ii < n; ii++) {
          i = order[ ii ]; // current word token to assess

          wi     = w[i]; // current word index
          di     = d[i]; // current document index
          topic  = z[i]; // current topic assignment for word token
          author = x[i]; // current author assignment for word token

          i_start  = jcad[ di   ];
          i_end    = jcad[ di+1 ];
          nauthors = i_end - i_start;

          ztot[topic]--;  // substract this from counts
          atot[author]--;

          wioffset = wi*T;
          aioffset = author*T;

          wp[wioffset+topic]--;
          at[aioffset+topic]--;

          //mexPrintf( "(1) Working on ii=%d i=%d wi=%d di=%d topic=%d author=%d wp=%d at=%d\n" , ii , i , wi , di , topic , author , wp[wioffset+topic] , at[aioffset+topic] );

          totprob = (double) 0;
          kk = 0;
          for (k=0; k<nauthors; k++)
          {
             author = (int) irad[ i_start + k ]; // current author index under consideration
			 aioffset = author * T;

			 for (j = 0; j < T; j++)
			 {
				  // probs[j][k] contains the (unnormalized) probability of assigning this word token to topic j and author k
			      probs[kk] =

					  ((double) wp[wioffset+j] + (double) BETA)  /( (double) ztot[j]  + (double) WBETA) *
				      ((double) at[aioffset+j] + (double) ALPHA) /( (double) atot[author] + (double) KALPHA );

			      totprob += probs[kk];
                  kk++;
			  }
		  }

          // sample a topic from the distribution
          r = (double) totprob * (double) randomMT() / (double) 4294967296.0;
          max = probs[0];
          if (nauthors==1) {
              topic = 0;
              while (r>max) {
                  topic++;
                  max += probs[topic];
              }
          } else
          {
              kk = 0;
              k  = 0;
              topic = 0;
              while (r>max) {
                  kk++;
                  topic++;
                  if (topic == T) {
                      k++;
                      topic = 0;
                  }
                  max += probs[kk];
              }

              author = (int) irad[ i_start + k ]; // sampled author
          }

          z[i] = topic; // assign current word token i to topic j
          x[i] = author;
          wp[wioffset + topic ]++; // and update counts
          ztot[topic]++;

          aioffset = author*T;
          atot[ author ]++;
          at[aioffset+topic]++;

          //mexPrintf( "(2) Working on ii=%d i=%d wi=%d di=%d topic=%d author=%d wp=%d at=%d\n" , ii , i , wi , di , topic , author , wp[wioffset+topic] , at[aioffset+topic] );

      }
  }

}

void mexFunction(int nlhs, mxArray *plhs[], int nrhs,
                 const mxArray *prhs[])
{
  double *srwp, *srat, *srad, *probs, *WS, *DS, *ZIN, *XIN, *Z, *X;
  double ALPHA,BETA;
  //int *irwp, *jcwp, *irat, *jcat, *irad, *jcad;
  mwIndex *irwp, *jcwp, *irat, *jcat, *irad, *jcad;
  int *z,*d,*w, *x, *order, *wp, *at, *ztot, *atot;
  int W,T,D,A,NN,SEED,OUTPUT, nzmax, nzmaxad, nzmaxwp, nzmaxat, ntokens;
  int i,j,c,n,nt,wi,ntcount,di;
  int i_start, i_end, a, nauthors, startcond;

// Syntax
//   [ WP , AT , Z , X ] = GibbsSamplerAT( WS , DS , AD , T , N , ALPHA , BETA , SEED , OUTPUT )

// Syntax
//   [ WP , AT , Z , X ] = GibbsSamplerAT( WS , DS , AD , T , N , ALPHA , BETA , SEED , OUTPUT , ZIN , XIN )


  /* Check for proper number of arguments. */
  if (nrhs < 9) {
    mexErrMsgTxt("At least 9 input arguments required");
  } else if (nlhs != 4) {
    mexErrMsgTxt("Four output arguments required");
  }

  startcond = 0;
  if (nrhs > 9) startcond = 1;

  /* process the input arguments */
  if (mxIsDouble( prhs[ 0 ] ) != 1) mexErrMsgTxt("WS input vector must be double precision");
  if (mxIsDouble( prhs[ 1 ] ) != 1) mexErrMsgTxt("DS input vector must be double precision");

  WS = mxGetPr( prhs[ 0 ] );
  DS = mxGetPr( prhs[ 1 ] );

  // get the number of tokens
  ntokens = (int) mxGetM( prhs[ 0 ] ) * (int) mxGetN( prhs[ 0 ] );
  if (ntokens == 0) mexErrMsgTxt("WS vector is empty");
  if (ntokens != ( mxGetM( prhs[ 1 ] ) * mxGetN( prhs[ 1 ] ))) mexErrMsgTxt("WS and DS vectors should have same number of entries");

  d  = (int *) mxCalloc( ntokens , sizeof( int ));
  w  = (int *) mxCalloc( ntokens , sizeof( int ));

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

  if ((mxIsSparse( prhs[ 2 ] ) != 1) || (mxIsDouble( prhs[ 2 ] ) != 1))
      mexErrMsgTxt("AT input matrix must be a sparse double precision matrix");

  /* dealing with sparse array AD */
  srad = mxGetPr(prhs[2]);
  irad  = mxGetIr(prhs[2]);
  jcad  = mxGetJc(prhs[2]);
  nzmaxad= (int) mxGetNzmax(prhs[2]);
  A    = (int) mxGetM( prhs[2] );
  if (mxGetN( prhs[2] ) != D) mexErrMsgTxt("The number of documents in DS vector and AD matrix must be equal" );

  T    = (int) mxGetScalar(prhs[3]);
  if (T<=0) mexErrMsgTxt("Number of topics must be greater than zero");

  NN    = (int) mxGetScalar(prhs[4]);
  if (NN<0) mexErrMsgTxt("Number of iterations must be greater than zero");

  ALPHA = (double) mxGetScalar(prhs[5]);
  if (ALPHA<=0) mexErrMsgTxt("ALPHA must be greater than zero");

  BETA = (double) mxGetScalar(prhs[6]);
  if (BETA<=0) mexErrMsgTxt("BETA must be greater than zero");

  SEED = (int) mxGetScalar(prhs[7]);
  // set the seed of the random number generator

  OUTPUT = (int) mxGetScalar(prhs[8]);

  if (startcond == 1) {
      ZIN = mxGetPr( prhs[ 9 ] );
      if (ntokens != ( mxGetM( prhs[ 9 ] ) * mxGetN( prhs[ 9 ] ))) mexErrMsgTxt("WS and ZIN vectors should have same number of entries");

      XIN = mxGetPr( prhs[ 10 ] );
      if (ntokens != ( mxGetM( prhs[ 10 ] ) * mxGetN( prhs[ 10 ] ))) mexErrMsgTxt("WS and XIN vectors should have same number of entries");
  }

  // seeding
  seedMT( 1 + SEED * 2 ); // seeding only works on uneven numbers

  // check entries of AD matrix
  for (i=0; i<nzmaxad; i++) {
      nt = (int) srad[ i ];
      if ((nt<0) || (nt>1)) mexErrMsgTxt("Entries in AD matrix can only be 0 or 1");
  }

  /* allocate memory */
  z  = (int *) mxCalloc( ntokens , sizeof( int ));
  x  = (int *) mxCalloc( ntokens , sizeof( int ));

  if (startcond == 1) {
     for (i=0; i<ntokens; i++) {
         z[ i ] = (int) ZIN[ i ] - 1;
         x[ i ] = (int) XIN[ i ] - 1;
     }
  }

  order  = (int *) mxCalloc( ntokens , sizeof( int ));
  wp     = (int *) mxCalloc( T*W , sizeof( int ));
  at     = (int *) mxCalloc( T*A , sizeof( int ));
  ztot   = (int *) mxCalloc( T , sizeof( int ));
  atot   = (int *) mxCalloc( A , sizeof( int ));
  probs  = (double *) mxCalloc( T*NAMAX , sizeof( double ));

  /* check that every document has some authors */
  for (j=0; j<D; j++) {
      i_start = (int) jcad[j];
      i_end   = (int) jcad[j+1];
      nauthors = i_end - i_start;
      if (nauthors==0) mexErrMsgTxt("there are some documents without authors in AD matrix ");
      if (nauthors>NAMAX) mexErrMsgTxt("Too many authors in some documents ... reached the NAMAX limit");
  }

  /*
  for (j=0; j<D; j++) {
      i_start = jcad[j];
      i_end   = jcad[j+1]-1;
      for (i=i_start; i<=i_end; i++) {
          a = irad[ i ];
          mexPrintf( "Document j=%d i_start=%d i_end=%d author=%d\n" , j , i_start , i_end , a );
      }
  }
  */

  /*for (i=0; i<n; i++) {
     mexPrintf( "i=%4d w[i]=%3d d[i]=%d z[i]=%d\n" , i , w[i] , d[i] , z[i] );
  }
  */

  if (OUTPUT==2) {
      mexPrintf( "Running AT Gibbs Sampler Version 1.0\n" );
      if (startcond==1) mexPrintf( "Starting from a previous state\n" );
      mexPrintf( "Arguments:\n" );
      mexPrintf( "\tNumber of words       W = %d\n" , W );
      mexPrintf( "\tNumber of docs        D = %d\n" , D );
      mexPrintf( "\tNumber of docs        A = %d\n" , A );
      mexPrintf( "\tNumber of topics      T = %d\n" , T );
      mexPrintf( "\tNumber of iterations  N = %d\n" , NN );
      mexPrintf( "\tHyperparameter    ALPHA = %4.4f\n" , ALPHA );
      mexPrintf( "\tHyperparameter     BETA = %4.4f\n" , BETA );
      mexPrintf( "\tSeed number             = %d\n" , SEED );
      mexPrintf( "\tNumber of tokens        = %d\n" , ntokens );
      mexPrintf( "Internal Memory Allocation\n" );
      mexPrintf( "\tw,d,z,x,order indices combined  = %d bytes\n" , 5 * sizeof( int) * ntokens );
      mexPrintf( "\twp (full) matrix                = %d bytes\n" , sizeof( int ) * W * T  );
      mexPrintf( "\tdp (full) matrix                = %d bytes\n" , sizeof( int ) * D * T  );
      //mexPrintf( "Checking: sizeof(int)=%d sizeof(long)=%d sizeof(double)=%d\n" , sizeof(int) , sizeof(long) , sizeof(double));
  }

  /* run the model */
  GibbsSamplerAT( ALPHA, BETA, W, T, D, NN, OUTPUT, n, z, d, w, x, wp, at, ztot, atot, order, probs, irad, jcad, startcond );

  /* convert the full wp matrix into a sparse matrix */
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
         c = (int) *( wp + i*T + j );
         if (c >0) {
             *( srwp + n ) = c;
             *( irwp + n ) = i;
             n++;
         }
      }
  }

  *( jcwp + T ) = n;

  // create sparse matrix AT
  nzmaxat = 0;
  for (i=0; i<A; i++) {
      for (j=0; j<T; j++)
          nzmaxat += (int) ( *( at + j + i*T )) > 0;
  }

  if (OUTPUT==2) {
      mexPrintf( "Constructing sparse output matrix dp\n" );
      mexPrintf( "Number of nonzero entries for AT = %d\n" , nzmaxat );
  }

  plhs[1] = mxCreateSparse( A,T,nzmaxat,mxREAL);
  srat  = mxGetPr(plhs[1]);
  irat = mxGetIr(plhs[1]);
  jcat = mxGetJc(plhs[1]);

  n = 0;
  for (j=0; j<T; j++) {
      *( jcat + j ) = n;
      for (i=0; i<A; i++) {
          c = (int) *( at + i*T + j );
          if (c >0) {
              *( srat + n ) = c;
              *( irat + n ) = i;
              n++;
          }
      }
  }

      *( jcat + T ) = n;


  plhs[ 2 ] = mxCreateDoubleMatrix( 1,ntokens , mxREAL );
  Z = mxGetPr( plhs[ 2 ] );
  for (i=0; i<ntokens; i++) Z[ i ] = (double) z[ i ] + 1;

  plhs[ 3 ] = mxCreateDoubleMatrix( 1,ntokens , mxREAL );
  X = mxGetPr( plhs[ 3 ] );
  for (i=0; i<ntokens; i++) X[ i ] = (double) x[ i ] + 1;
}
