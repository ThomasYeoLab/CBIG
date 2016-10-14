/*
 * David Gleich
 * Copyright, Stanford University, 2006-2007
 */
 
 /**
 * @file libmbgl_test.cc
 * Implement a unit test for libmbgl
 */
 
/**History

 * 31 July 2007
 * Initial version
 */

#include "matlab_bgl.h"

int main(int argc, char **argv)
{
    mbglIndex ai[] = {0};
    mbglIndex aj[] = {0};
    mbglIndex ei[6], ej[6];
    mbglIndex nedges = 0;

    triangulate_graph(0, aj, ai, 1, 1, 0, ei, ej, &nedges);

    return 0;
}