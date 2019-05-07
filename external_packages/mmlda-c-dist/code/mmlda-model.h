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

#ifndef MMLDA_MODEL_H
#define MMLDA_MODEL_H

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "mmlda.h"
#include "mmlda-alpha.h"
#include "cokus.h"

#define myrand() (double) (((unsigned long) randomMT()) / 4294967296.)
#define NUM_INIT 1

void free_mmlda_model(mmlda_model*);
void save_mmlda_model(mmlda_model*, char*);
mmlda_model* new_mmlda_model(int, int, int);
mmlda_suffstats* new_mmlda_suffstats(mmlda_model* model);
void corpus_initialize_ss(mmlda_suffstats* ss, mmlda_model* model, corpus* c1, corpus* c2);
void random_initialize_ss(mmlda_suffstats* ss, mmlda_model* model);
void zero_initialize_ss(mmlda_suffstats* ss, mmlda_model* model);
void mmlda_mle(mmlda_model* model, mmlda_suffstats* ss, int estimate_alpha);
mmlda_model* load_mmlda_model(char* model_root);

#endif
