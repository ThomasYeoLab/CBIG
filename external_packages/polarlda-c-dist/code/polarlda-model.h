// (C) Copyright 2018, Nanbo Sun (nanbosun@u.nus.edu), Xiuming Zhang

// This file is part of POLARLDA-C.

// POLARLDA-C is free software; you can redistribute it and/or modify it under
// the terms of the GNU General Public License as published by the Free
// Software Foundation; either version 2 of the License, or (at your
// option) any later version.

// POLARLDA-C is distributed in the hope that it will be useful, but WITHOUT
// ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
// FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
// for more details.

// You should have received a copy of the GNU General Public License
// along with this program; if not, write to the Free Software
// Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
// USA

#ifndef POLARLDA_MODEL_H
#define POLARLDA_MODEL_H

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <assert.h>
#include "polarlda.h"
#include "polarlda-alpha.h"
#include "cokus.h"

#define myrand() (double) (((unsigned long) randomMT()) / 4294967296.)
#define NUM_INIT 1

void free_polarlda_model(polarlda_model*);
void save_polarlda_model(polarlda_model*, char*);
polarlda_model* new_polarlda_model(int, int);
polarlda_suffstats* new_polarlda_suffstats(polarlda_model* model);
void corpus_initialize_ss(polarlda_suffstats* ss, polarlda_model* model, corpus* c);
void random_initialize_ss(polarlda_suffstats* ss, polarlda_model* model);
void zero_initialize_ss(polarlda_suffstats* ss, polarlda_model* model);
void polarlda_mle(polarlda_model* model, polarlda_suffstats* ss, int estimate_alpha);
polarlda_model* load_polarlda_model(char* model_root);

#endif
