#ifndef MATLAB_BGL_EXPAND_MACROS_H
#define MATLAB_BGL_EXPAND_MACROS_H

/** @file planar.cc
 * @author David F. Gleich
 * @date 2007-02-22
 * @copyright Stanford University, 2007-2008
 * Write a few simple ``macro''-ish functions to expand compressed
 * quantities to their original type.  This permits us to save on
 * duplicating space and allocating additional RAM.
 */

/** History
 *  2007-02-22: Initial version
 */

void expand_int_to_double(int* src, double* dst, mwSize len, double offset)
{
  int i;
  for (i = len - 1; i >= 0; i--) {
    dst[i] = (double) src[i] + offset;
  }
}

void expand_degree_to_double(mbglDegreeType* src, double* dst, mwSize len,
     double offset) 
{
  while (len>0) {
    --len;
    dst[len] = (double) src[len] + offset;
  }
}

void expand_index_to_double(mwIndex* src, double* dst, mwSize len,
    double offset)
{
  int i;
  for (i = len - 1; i >= 0; i--) {
    dst[i] = (double) src[i] + offset;
  }
}

/**
 * This function is used to update something like a predecessor array
 * where the output in the double array should be 0 if the input array
 * is equal to its index.
 */
void expand_index_to_double_zero_equality(mwIndex* src, double* dst,
    mwSize len, double offset)
{
  int i;
  for (i = len - 1; i >= 0; i--) {
    if (src[i] == i) {
      dst[i] = 0.0;
    } else {
      dst[i] = (double) src[i] + offset;
    }
  }
}

/**
 * This function is used to update something like an initialized
 * array where there is one special value used to initialize the
 * structure and we want to do something different with that
 * special value on output (like set it to 0 :-))
 */
void expand_index_to_double_zero_special(mwIndex* src, double* dst, mwSize len,
    double offset, mwIndex special)
{
  int i;
  for (i = len - 1; i >= 0; i--) {
    if (src[i] == special) {
      dst[i] = 0.0;
    } else {
      dst[i] = (double) src[i] + offset;
    }
  }
}

#endif /* MATLAB_BGL_EXPAND_MACROS_H */

