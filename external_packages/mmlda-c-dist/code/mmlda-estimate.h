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

#ifndef MMLDA_ESTIMATE_H
#define MMLDA_ESTIMATE_H

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <float.h>
#include <string.h>
#include <time.h>

#include "mmlda.h"
#include "mmlda-data.h"
#include "mmlda-inference.h"
#include "mmlda-model.h"
#include "mmlda-alpha.h"
#include "utils.h"
#include "cokus.h"

int LAG = 5;

float EM_CONVERGED;
int EM_MAX_ITER;
int ESTIMATE_ALPHA;
double INITIAL_ALPHA;
int NTOPICS;

double doc_e_step(document* doc1,
                  document* doc2,
                  double* gamma,
                  double** phi1,
                  double** phi2,
                  mmlda_model* model,
                  mmlda_suffstats* ss);

void save_gamma(char* filename,
                double** gamma,
                int num_docs,
                int num_topics);

void run_em(char* start,
            char* directory,
            corpus* corpus1,
            corpus* corpus2);

void read_settings(char* filename);

void infer(char* model_root,
           char* save,
           corpus* corpus1,
           corpus* corpus2);

#endif


