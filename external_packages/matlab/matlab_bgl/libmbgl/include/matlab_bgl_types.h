#ifndef MATLAB_BGL_TYPES_H
#define MATLAB_BGL_TYPES_H

/*
 * David Gleich
 * 19 February 2007
 * Copyright, Stanford University, 2007
 */

/**
 * @file matlab_bgl_types.h
 * Implement a series of types for the matlab_bgl project to support
 * large sparse matrices on 64-bit platforms.
 */
 
/** History
 *
 *  2007-07-30: Added mbglDegreeType to be a type that is a 64-bit 
 *    integer when for large graphs and a 32-bit unsigned int for
 *    small graphs.
 * 2007-08-27: Added proper comments for C files.
 * 2008-09-19: Fixed comments
 */

#ifdef MATLAB_BGL_LARGE_ARRAYS
#include <stdlib.h>
typedef size_t mbglIndex;
typedef mbglIndex mbglDegreeType;
#else
typedef int mbglIndex;
typedef unsigned int mbglDegreeType;
#endif /* MATLAB_BGL_LARGE_ARRAYS */

#endif /* MATLAB_BGL_TYPES_H */
