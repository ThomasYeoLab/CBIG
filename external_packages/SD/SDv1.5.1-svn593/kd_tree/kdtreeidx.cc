// Guy Shechter
// June 2004


// ********
// Edit kdtree_common.h to set architecture specific pathnames.
// ********

//
// You can turn on some of the following definitions to debug this file.
//
#undef DEBUG            //General input-output debugging information
#undef DEBUG_BUILD_TREE  //Follow how the k-D tree is being built.
#undef DISPLAY_TREE      //Output the tree in a depth first traversal
#undef DEBUG_RUN_QUERIES //Follow the tree querying process

#undef TIME             //Find out how long it takes to build the k-D tree
                         //and to perform queries.


//
// Standard includes 
//
#include <math.h>
#include <stdio.h>

//
// Core functions are located in the following files
//
#include "kdtree_common.h"
#include "kdtree_common.cc"


void mexFunction( int nlhs, mxArray **plhs, int nrhs, const mxArray **prhs){
  Tree         *tree;
  double       *reference, *model;
  int          *index, i;
  unsigned int  N, D, M;
  double       *closest_pts, *distances, *pointer_to_tree;
  int          SkipQueries=0;
  
  if (nrhs <2 ){
    mexErrMsgTxt("Must have at least two input arrays.");
  }
  
#ifdef DEBUG
  mexPrintf("Mex function called with %d inputs and %d explicit outputs\n",nrhs,nlhs);
#endif
  
  reference = mxGetPr(prhs[0]);
  N = mxGetM(prhs[0]);
  D = mxGetN(prhs[0]);
  
  if ((!N || !D ) && ( nrhs < 3) )
    mexErrMsgTxt("You have to supply some reference points to build a k-D tree.");
  
#ifdef TIME
  gettimeofday(&tv1,&tz);
#endif
  
  //
  //
  // If the tree is not passed in as a third input, we must build it
  //
  //
  if (nrhs < 3 ){   
    
#ifdef DEBUG
    mexPrintf("----------------------\n");
    mexPrintf("Building k-D Tree ...\n");
#endif
    
    index = (int*) malloc( sizeof(int) * N);
    for (i=0; i < N; i++) index[i]=i;  
    if ( (tree = build_kdtree(reference,N,D,index,N,0))==NULL ){
      free(index);
      mexErrMsgTxt("Not enough free memory to build k-D tree\n");
    } else {
      tree->dims = D;
      free(index);
    }

#ifdef DEBUG
    mexPrintf("Done Building k-D Tree\n");
    mexPrintf("----------------------\n");
#endif
    
  } else {
    
    //
    // The tree was built previously, and is now being passed in to the function.
    //
    if (   (pointer_to_tree = mxGetPr(prhs[2])) == NULL )
      mexErrMsgTxt("Third argument is not a valid pointer to a k-D tree\n");

    if ( (tree = (Tree *) ((long) pointer_to_tree[0]))== NULL )
      mexErrMsgTxt("Third argument is not a valid pointer to a k-D tree\n");

    
  }
  
#ifdef TIME
  gettimeofday(&tv2,&tz);
  if (tv2.tv_usec - tv1.tv_usec < 0) {
    tv2.tv_sec--;
    tv2.tv_usec += 1000000;
  }    
  mexPrintf("Time to Build Tree : %f\n", tv2.tv_sec -tv1.tv_sec+(tv2.tv_usec-tv1.tv_usec)/1000000.0);
#endif
  

#ifdef DISPLAY_TREE
  mexPrintf("\nDepth first traversal of the k-D tree\n");
  mexPrintf("-------------------------------------\n");
  display_tree(tree->rootptr,D);
  mexPrintf("-------------------------------------\n");
#endif

  
  //
  //  Query section
  //
  //
  
  model = mxGetPr(prhs[1]);
  M = mxGetM(prhs[1]);


  if (!model && !M) { 

    // There are no points to query
    SkipQueries=1;

  } else {

    // Check that the model points are of the same dimension as the 
    // reference points in the k-d tree.
    if (mxGetN(prhs[1]) != tree->dims)
      mexErrMsgTxt("Reference and Model Vectors must be of the same dimension");
  }

  if (nlhs >=0){
    plhs[0] = mxCreateDoubleMatrix(M,1,mxREAL);
    closest_pts = mxGetPr(plhs[0]);
  }
  else{ 
    closest_pts = (double *) malloc (sizeof(double) *M);
  }

  if (nlhs >=2) {
    plhs[1] = mxCreateDoubleMatrix(M,1,mxREAL);
    distances = mxGetPr(plhs[1]);
  }
  else {
    distances = (double *) malloc (sizeof(double)*M);
  }

  if (nlhs >=3) {
    plhs[2] = mxCreateDoubleMatrix(1,1,mxREAL);
    pointer_to_tree = mxGetPr(plhs[2]);
    pointer_to_tree[0] = (long) tree;
  }
  
  if (!SkipQueries) {

#ifdef TIME
    gettimeofday(&tv1,&tz);
#endif  
  
#ifdef DEBUG
    mexPrintf("--------------------\n");
    mexPrintf("Running Queries...\n");
#endif
    
    run_queries(tree->rootptr, model, M, tree->dims, closest_pts,
		distances, RETURN_INDEX);
    
    //Since MATLAB is a 1..N indexing language, we add 1 to each 
    //index value
    for (i=0; i < M ; i++)
      closest_pts[i]++;


    
#ifdef DEBUG
    mexPrintf("Done Running Queries\n");
    mexPrintf("--------------------\n");
#endif

#ifdef TIME
    gettimeofday(&tv2,&tz);
    if (tv2.tv_usec - tv1.tv_usec < 0) {
      tv2.tv_sec--;
      tv2.tv_usec += 1000000;
    }    
    mexPrintf("Time per Search : %f\n", 
	      (tv2.tv_sec - tv1.tv_sec + 
	       (tv2.tv_usec-tv1.tv_usec) /1000000.0 )/(double)M);
#endif
  }

  
  if (nlhs<3) {
#ifdef DEBUG
    mexPrintf("-------------------------------------\n");
    mexPrintf("Removing k-D Tree from system memory.\n");
#endif
    free_tree(tree->rootptr);
    free(tree);
#ifdef DEBUG
    mexPrintf("Done.\n");
    mexPrintf("-------------------------------------\n");
#endif
    if (nlhs < 2){
      free(distances);
    }
  } else{
#ifdef DEBUG
    mexPrintf("--------------------------------\n");
    mexPrintf("k-D Tree saved in system memory.\n");
    mexPrintf("Don't forget to remove it later.\n");
    mexPrintf("--------------------------------\n");
#endif
  }
    
#ifdef DEBUG
  mexPrintf("Mex function has exited normally.\n");
#endif
}



//
//  k-D Tree Index main function 
//
void kdtreeidx_main() {
}
