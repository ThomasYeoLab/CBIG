// Guy Shechter
// June 2004


// ********
// Edit kdtree_common.h to set architecture specific pathnames.
// ********

//
// You can turn on some of the following definitions to debug this file.
//
#undef DEBUG             //General input-output debugging information
#undef DEBUG_RUN_RANGE_SEARCH  //Follow the tree querying process

#undef TIME              //Find out how long it takes to build the k-D tree
                         //and to perform queries.


//
// Standard includes 
//

#include <math.h>
#include <stdio.h>


#include "kdtree_common.h"
#include "kdtree_common.cc"


void mexFunction( int nlhs, mxArray **plhs, int nrhs, const mxArray **prhs){
  Tree         *tree;
  double       *distlim, *querypt, *tmppt;
  int          *index, i,j;
  unsigned int  N, M;
  double       *found_pts, *distances, *pointer_to_tree, *pIndices;
  unsigned int *indices;
  unsigned int L;
  int          SkipQueries=0;
  char buf[100];

  if (nrhs <3 ) {
    mexErrMsgTxt("Must have at least 3 inputs (ROOT, QUERYPT, DISTLIM).\n"); //, [MinNumberOfNeighbors].");
  }
  
#ifdef DEBUG
  mexPrintf("Mex function called with %d inputs and %d explicit outputs\n",nrhs,nlhs);
#endif

  if (   (pointer_to_tree = mxGetPr(prhs[0])) == NULL )
    mexErrMsgTxt("First argument must be a valid pointer to a k-D tree\n");

  if ( (tree = (Tree *) ((long) pointer_to_tree[0]))== NULL )
    mexErrMsgTxt("First argument must be a valid pointer to a k-D tree\n");
  

  querypt = mxGetPr(prhs[1]);
  if (  (M = mxGetM(prhs[1])) == 0)  // There are no points to query
    SkipQueries=1;    

  if ( mxGetN(prhs[1]) != tree->dims){
    sprintf(buf,"The dimension of the k-D tree (%d) and queried point (%d) do not match",tree->dims, mxGetN(prhs[1]));
    mexErrMsgTxt(buf);
  }

  if (M > 1)         // There are too manypoints to query
    mexErrMsgTxt("Only one search may be performed. QUERYPT must be 1xD array");


  distlim = mxGetPr(prhs[2]);
  if ( (mxGetM(prhs[2]) != 1)  || (mxGetN(prhs[2])!=1))
    mexErrMsgTxt("Distance limit (DISTLIM) must be a scalar value >= 0.");
  if ( distlim[0] < 0 )
    mexErrMsgTxt("Distance limit (DISTLIM) must be a scalar value >= 0.");

  if (!SkipQueries) {
    
#ifdef TIME
    gettimeofday(&tv1,&tz);
#endif  
    

#ifdef DEBUG
    mexPrintf("-----------------------\n");
    mexPrintf("Running Range Search...\n");
#endif
    
    run_range_search(tree->rootptr,querypt,M,tree->dims,distlim[0],
		     &found_pts, &L, &indices);
    
#ifdef DEBUG
    mexPrintf("Done Range Search \n");
    mexPrintf("-----------------------\n");
#endif
    
#ifdef TIME
    gettimeofday(&tv2,&tz);
    if (tv2.tv_usec - tv1.tv_usec < 0) {
      tv2.tv_sec--;
      tv2.tv_usec += 1000000;
    }    
    mexPrintf("Time for Range Search : %f\n", 
	      (tv2.tv_sec - tv1.tv_sec + 
	       (tv2.tv_usec-tv1.tv_usec) /1000000.0 )/(double)M);
#endif
    
    //
    // Copy the points to the output matrix
    //
    if (nlhs >=0){
      plhs[0] = mxCreateDoubleMatrix(L,tree->dims,mxREAL);
      tmppt = mxGetPr(plhs[0]);
      for (i=0; i < L*(tree->dims); i++) {
	tmppt[i]=found_pts[i];
      }
      if (nlhs >=2) {
	plhs[1] = mxCreateDoubleMatrix(L,1,mxREAL);
	distances = mxGetPr(plhs[1]);	
	for (i=0; i < L; i++) {
	  distances[i]=0;
	  for (j=0; j < (tree->dims); j++) 
	    distances[i]+= (tmppt[j*L+i]-querypt[j])*(tmppt[j*L+i]-querypt[j]);
	}
	for (i=0; i < L; i++) 
	  distances[i]=sqrt(distances[i]);
      }
      if (nlhs >=3) {
	plhs[2] = mxCreateDoubleMatrix(L,1,mxREAL);
	pIndices = mxGetPr(plhs[2]);
	for (i=0; i < L; i++) {
	  pIndices[i]=indices[i]+1; //Add 1 because of Matlab indexing 1..N
	}
      }
    }
    if(L) {
      free(indices);
      free(found_pts);
    }
  }
  
#ifdef DEBUG
  mexPrintf("Mex function has exited normally.\n");
#endif
}


void kdrangequery_main() {
}







