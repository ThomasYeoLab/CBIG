#include "mex.h"

void revord(char *input_buf, int buflen, char *output_buf)
{
  int   i;

  /* reverse the order of the input string */
  for(i=0;i<buflen-1;i++) 
    *(output_buf+i) = *(input_buf+buflen-i-2);
}

void mexFunction( int nlhs, mxArray *plhs[],
                  int nrhs, const mxArray *prhs[])
{
    char *entry, *string;
    int   buflen,status;
    int   N, index_now, outcome, low_index, high_index;
    const mxArray *cell_element;
    const mxArray *cell_array, *compare_string;
    double *outcomed;
    
    /* check for proper number of arguments */
    if(nrhs!=2) 
      mexErrMsgTxt("TWO inputs required.");
    else if(nlhs != 1) 
      mexErrMsgTxt("ONE output reuired.");

    /* input 1 must be a string */
    if ( mxIsChar(prhs[0]) != 1)
      mexErrMsgTxt("Input must be a string");
    
    compare_string = prhs[0];
    buflen = mxGetN( compare_string ) + 1;
    string =mxCalloc(buflen, sizeof(char));
    status = mxGetString( compare_string , string, buflen);
    if(status != 0) mexWarnMsgTxt("Not enough space. String is truncated.");
    
    /* input 2 must be a cell array */
    if ( mxIsCell(prhs[1]) != 1)
      mexErrMsgTxt("Input must be a cell array");

    if ( mxGetNumberOfDimensions(prhs[1]) != 2)
      mexErrMsgTxt("Cell array must be unidimensional.");
    
    cell_array = prhs[ 1 ];
    N = mxGetNumberOfElements( cell_array );
    
    low_index  = 0;
    high_index = N-1;
    
    while (low_index <= high_index) {
        index_now = (low_index + high_index) / (int) 2;
        
        cell_element = mxGetCell( cell_array , index_now );
        
        buflen = (mxGetN( cell_element )) + 1;
        entry  = mxCalloc(buflen, sizeof(char));
        status = mxGetString( cell_element , entry, buflen);
        if(status != 0) mexWarnMsgTxt("Not enough space. String is truncated.");
                
        outcome = strcmp(string , entry);
                
        if (outcome < 0)
            high_index = index_now - 1;
        else if (outcome > 0)
            low_index = index_now + 1;
        else break;
    }
       
    index_now++;
    if (outcome != 0) index_now = 0;
       
    plhs[0] = mxCreateDoubleMatrix(1,1, mxREAL);
    outcomed = mxGetPr(plhs[0]);
    
    outcomed[0] = (double) index_now;
    
    return;
}

