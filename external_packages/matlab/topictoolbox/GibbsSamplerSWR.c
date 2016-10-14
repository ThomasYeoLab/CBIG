#include "mex.h"
#include "cokus.cpp"

// Syntax
//  [ X , XD , PQ ] = GibbsSamplerSWR( PQ0 , PQ1 , PQ2 , GAMMA , BURNIN , NS , LAG , SEED , OUTPUT );

// INPUTS
// [0]  PQ0:    NQ x D matrix of query word probabilities for each document according to topic model
// [1]  PQ1:    NQ x D matrix of query word probabilities for each document according to document model
// [2]  PQ2:    1  x NQ matrix of query word probabilities for each document according to corpus model
// [3]  GAMMA:  1 x 3 array for hyperparameters for three routes: [0] = topics [1] document words [2] corpus level words
// [4]  BURNIN: number of burnin iterations for Gibbs sampler
// [5]  NS:     number of samples to take
// [6]  LAG:    lag between samples
// [7]  SEED:   seed for random number generator
// [8]  OUTPUT: verbosity flag

// OUTPUT
// [0]  X:     D x 3 array of route probabilities (for each docs)
// [1]  XD:    D x 3 x NQ array of route probabilities (for each doc)
// [2]  PQ:    NQ x D array of retrieval probabilities

void mexFunction(int nlhs, mxArray *plhs[], int nrhs,
                 const mxArray *prhs[])
{
  double *PQ0, *PQ1, *PQ2, *GAMMA, *X, *XD, *PQ, *probs;
  double totweight, sumlogprob, pwz, pwd, pwc, proute0, proute1, proute2, prob, rn, max;
  int *dims, *x, *xtot;
  int NQ, D, BURNIN, S, NS, NN, LAG, SEED, OUTPUT, i, iter, d, route, r, oldroute, newroute;
    
    
  // Check for proper number of arguments.
  if (nrhs < 9) {
    mexErrMsgTxt("At least 9 input arguments required");
  } else if (nlhs < 3) {
    mexErrMsgTxt("3 output arguments required");
  }
  
  // process the input arguments
  if (mxIsDouble( prhs[ 0 ] ) != 1) mexErrMsgTxt("PQ0 input matrix must be a double precision matrix");
  if (mxIsDouble( prhs[ 1 ] ) != 1) mexErrMsgTxt("PQ1 input matrix must be a double precision matrix");
  if (mxIsDouble( prhs[ 2 ] ) != 1) mexErrMsgTxt("PQ2 input matrix must be a double precision matrix");
  if (mxIsDouble( prhs[ 3 ] ) != 1) mexErrMsgTxt("GAMMA input vector must be a double precision matrix");
 
  // pointer to PQ0 matrix
  PQ0 = mxGetPr( prhs[ 0 ] );
  NQ = mxGetM( prhs[ 0 ] );
  D  = mxGetN( prhs[ 0 ] );
  
  // pointer to PQ1 matrix
  PQ1 = mxGetPr( prhs[ 1 ] );
  if ( mxGetM( prhs[ 1 ] ) != NQ ) mexErrMsgTxt("PQ1 matrix should have same dimensions as PQ0 matrix");
  if ( mxGetN( prhs[ 1 ] ) != D )  mexErrMsgTxt("PQ1 matrix should have same dimensions as PQ0 matrix");
  
  // pointer to PQ2 matrix
  PQ2 = mxGetPr( prhs[ 2 ] );
  if ( mxGetM( prhs[ 2 ] ) != 1 )  mexErrMsgTxt("PQ2 matrix should be a row vector");
  if ( mxGetN( prhs[ 2 ] ) != NQ )  mexErrMsgTxt("PQ2 matrix should same number of query words as PQ0 matrix");
    
  // pointer to gamma hyperparameter vector
  GAMMA = mxGetPr( prhs[ 3 ] );
  if ( mxGetM( prhs[ 3 ] ) != 1 ) mexErrMsgTxt("GAMMA input vector must be a row vector");
  if ( mxGetN( prhs[ 3 ] ) != 3 ) mexErrMsgTxt("GAMMA input vector must have three entries");
  
  BURNIN    = (int) mxGetScalar(prhs[4]);
  if (BURNIN<1) mexErrMsgTxt("Number of burnin iterations must be positive");

  NS    = (int) mxGetScalar(prhs[5]);
  if (NS<1) mexErrMsgTxt("Number of samples must be greater than zero");
  
  LAG    = (int) mxGetScalar(prhs[6]);
  if (LAG<1) mexErrMsgTxt("Lag must be greater than zero");
  
  // seeding
  SEED = (int) mxGetScalar(prhs[7]);
  seedMT( 1 + SEED * 2 ); // seeding only works on uneven numbers
  
  OUTPUT = (int) mxGetScalar(prhs[8]);
  
  // create output matrices
  plhs[ 0 ] = mxCreateDoubleMatrix( D , 3, mxREAL );
  X = mxGetPr( plhs[ 0 ] );

  dims = (int *) mxCalloc( 3 , sizeof( int ));
  dims[ 0 ] = D;
  dims[ 1 ] = 3;
  dims[ 2 ] = NQ;
  
  plhs[ 1 ] = mxCreateNumericArray( 3,dims,mxDOUBLE_CLASS,mxREAL);
  XD = mxGetPr( plhs[ 1 ] );
  
  plhs[ 2 ] = mxCreateDoubleMatrix( NQ , D, mxREAL );
  PQ = mxGetPr( plhs[ 2 ] );
  
  if (OUTPUT==2) 
  {
      mexPrintf( "Running Special Word Retrieval Gibbs Sampler Version 1.0\n" );
      mexPrintf( "Arguments:\n" );
      mexPrintf( "\tNumber of docs          D = %d\n" , D );
      mexPrintf( "\tNumber of query words  NQ = %d\n" , NQ );        
      mexPrintf( "\tHyperparameter     GAMMA0 = %4.4f\n" , GAMMA[0] );
      mexPrintf( "\tHyperparameter     GAMMA1 = %4.4f\n" , GAMMA[1] );
      mexPrintf( "\tHyperparameter     GAMMA2 = %4.4f\n" , GAMMA[2] );
      mexPrintf( "\tSeed number          SEED = %d\n" , SEED );
      mexPrintf( "\tBurnin             BURNIN = %d\n" , BURNIN );
      mexPrintf( "\tNumber of samples      NS = %d\n" , NS );
      mexPrintf( "\tLag between samples   LAG = %d\n" , LAG );
  }
  
  // allocate memory for sampling arrays
  x  = (int *) mxCalloc( NQ * D , sizeof( int ));
  xtot = (int *) mxCalloc( 3 * D , sizeof( int ));
  probs = (double *) mxCalloc( 3 , sizeof( double ));
  
  // initialize the sampler
  for (d=0; d<D; d++)
  {
      for (i=0; i<NQ; i++)
      {
          // pick a random route between 0 and 2
          route = (int) ( (double) randomMT() * (double) 3 / (double) (4294967296.0 + 1.0) );
          
          // assign this query word to this route
          x[ i + d * NQ ] = route;
          
          // update total
          xtot[ route + d * 3 ] += 1;
          
          //if (d<5) mexPrintf( "d=%d i=%d route=%d\n" , d , i , route );
      }
  }
  
  // RUN THE GIBBS SAMPLER
  NN = BURNIN + NS * LAG - LAG + 1;
  S = 0;
  for (iter=0; iter<NN; iter++)
  {
      if (OUTPUT >=1)
      {
          if ((iter % 10)==0) mexPrintf( "\tIteration %d of %d\n" , iter , NN );
          if ((iter % 10)==0) mexEvalString("drawnow;");
      }
      
      // do an iteration
      for (d=0; d<D; d++) // loop over all docs
      {
          for (i=0; i<NQ; i++) // loop over all words in query
          {
              oldroute = x[ i + d * NQ ];
              xtot[ oldroute + d * 3 ]--; // subtract this from the counts
                      
              pwz = PQ0[ i + d * NQ ];
              pwd = PQ1[ i + d * NQ ];
              pwc = PQ2[ i ];
        
              probs[ 0 ] = ((double) xtot[ 0 + d * 3 ] + (double) GAMMA[ 0 ] ) * pwz;  
              probs[ 1 ] = ((double) xtot[ 1 + d * 3 ] + (double) GAMMA[ 1 ] ) * pwd;
              probs[ 2 ] = ((double) xtot[ 2 + d * 3 ] + (double) GAMMA[ 2 ] ) * pwc;
              
              totweight = probs[ 0 ] + probs[ 1 ] + probs[ 2 ];
              
              // sample a route from this distribution
              rn = (double) totweight * (double) randomMT() / (double) 4294967296.0;
              max = probs[0];
              newroute = 0;
              while (rn>max) {
                  newroute++;
                  max += probs[newroute];
              }

              x[ i + d * NQ ] = newroute;
              xtot[ newroute + d * 3 ]++; // add this to the counts
          }
      }
      
      
      if ((iter >= BURNIN) && (((iter - BURNIN) % LAG) == 0))
      {
          S++;
          
          if (OUTPUT >=1)
          {
              //mexPrintf( "\tDrawing sample %d of %d\n" , S , NS );
              //mexEvalString("drawnow;");
          }
          
          // update the return variables with the new counts
          for (d=0; d<D; d++)
          {
              
              X[ d + 0 * D ] += (double) xtot[ 0 + d * 3 ] + GAMMA[ 0 ];
              X[ d + 1 * D ] += (double) xtot[ 1 + d * 3 ] + GAMMA[ 1 ];
              X[ d + 2 * D ] += (double) xtot[ 2 + d * 3 ] + GAMMA[ 2 ];
              
              for (i=0; i<NQ; i++)
              {
                  route = x[ i + d * NQ ]; // current route assignment
                  
                  for (r=0; r<3; r++)
                  {
                      if (r==route)
                      {   
                          XD[ d + r*D + i*D*3 ] += (double) 1 + GAMMA[ r ];  // D x 3 x NQ
                      } else
                      {
                          XD[ d + r*D + i*D*3 ] += (double) GAMMA[ r ];  // D x 3 x NQ
                      }
                  } 
                           
              }
          }
      }
  }
  
  // NOW CALCULATE PROBABILITY DISTRIBUTIONS
  for (d=0; d<D; d++)
  {      
     totweight = 0;
     for (route=0; route<3; route++) totweight += X[ d + route * D ];
     for (route=0; route<3; route++) X[ d + route * D ] /= totweight; 
     
     for (i=0; i<NQ; i++)
     {
        totweight = 0;
        for (route=0; route<3; route++) totweight += XD[ d + route*D + i*D*3 ];
        for (route=0; route<3; route++) XD[ d + route*D + i*D*3 ] /= totweight;
     }
  }
  
  // NOW CALCULATE RETRIEVAL PROBABILITY
  for (d=0; d<D; d++)
  {  
     // the following is correct with a model that has lambda outside the plate -- one route probability for each query 
     proute0 = X[ d + 0 * D ];
     proute1 = X[ d + 1 * D ];
     proute2 = X[ d + 2 * D ];
      
     for (i=0; i<NQ; i++)
     {
        // the following is correct with a model that has lambda inside the plate -- a separate probability for each word       
        //proute0 = XD[ d + 0 * D + i * D * 3 ];
        //proute1 = XD[ d + 1 * D + i * D * 3 ];
        //proute2 = XD[ d + 2 * D + i * D * 3 ];
          
        pwz = PQ0[ i + d * NQ ];
        pwd = PQ1[ i + d * NQ ];
        pwc = PQ2[ i ];
        
        prob = proute0 * pwz + proute1 * pwd + proute2 * pwc;
        
        PQ[ i + d * NQ ] = prob;
        //if (d<5) mexPrintf( "d=%d i=%d prob=%4.7f logprob=%4.7f p0=%4.5f p1=%4.5f p2=%4.5f  pwz=%4.5f  pwd=%4.5f  pwc=%4.5f\n" , d , i , prob , sumlogprob , proute0 , proute1 , proute2 , pwz , pwd , pwc );
     }
     
     //if (d<5) mexPrintf( "d=%d  logprob=%4.4f\n" , d , sumlogprob );
          
  }
   
}
