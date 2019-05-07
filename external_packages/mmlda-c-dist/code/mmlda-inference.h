// (C) Copyright 2018, Nanbo Sun (nanbosun@u.nus.edu)

// This file is part of MMLDA-C.

// MMLDA-C is free software; you can redistribute it and/or modify it under
// the terms of the GNU General Public License as published by the Free
// Software Foundation; either version 2 of the License, or (at your
// option) any later version.

// MMLDA-C is distributed in the hope that it will be useful, but WITHOUT
// ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
// FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
// for more details.

// You should have received a copy of the GNU General Public License
// along with this program; if not, write to the Free Software
// Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
// USA

#ifndef MMLDA_INFERENCE_H
#define MMLDA_INFERENCE_H

#include <math.h>
#include <float.h>
#include <assert.h>
#include "mmlda.h"
#include "utils.h"

float VAR_CONVERGED;
int VAR_MAX_ITER;

double mmlda_inference(document*, document*, mmlda_model*, double*, double**, double**);
double compute_likelihood(document*, document*, mmlda_model*, double**, double**, double*);

double mmlda_inference1(document* doc1, mmlda_model* model, double* var_gamma, double** phi1);
double compute_likelihood1(document* doc1, mmlda_model* model, double** phi1, double* var_gamma);

double mmlda_inference2(document* doc2, mmlda_model* model, double* var_gamma, double** phi2);
double compute_likelihood2(document* doc2, mmlda_model* model, double** phi2, double* var_gamma);

#endif
