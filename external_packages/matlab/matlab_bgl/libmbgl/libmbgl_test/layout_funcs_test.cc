

#include <include/matlab_bgl.h>
#include <stdio.h>
#include "libmbgl_funcs_test.h"

int test_layout_1() {
  const mbglIndex n = 9;
  mbglIndex rp[] = {0,3,7,10,14,19,23,26,30,33};
  mbglIndex ci[] = {0,1,3,0,1,2,4,1,2,5,0,3,4,6,1,3,4,5,7,2,4,5,8,3,6,7,4,6,7,8,5,7,8};
  double w[33], X[2*n], dists[n*n], springs[n*n];
  int rval;
  for (int i=0; i<33; i++) { w[i]= 1.0; }
  rval= kamada_kawai_spring_layout(n, ci, rp, w, 0.0001, 100, 1.0, 0, 1.0,
                 X, dists, springs);
  if (rval!=0) {
    errstr = "function error";
    return -1;
  }
  return 0;
}

int layout_funcs_test() 
{
  int nfail= 0, ntotal= 0, rval;
  const char* name;

  name= "kamada_kawai 3-grid";
  rval= test_layout_1(); ntotal++;
  if (rval!= 0) { nfail++; printf("%20s  %50s\n", name, errstr); }
  else { printf("%20s  success\n", name); }
   
  printf("\n");
  printf("Total tests  : %3i\n", ntotal);
  printf("Total failed : %3i\n", nfail);

  return nfail!=0;
}

