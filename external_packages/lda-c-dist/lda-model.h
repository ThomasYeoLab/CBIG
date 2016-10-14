#ifndef LDA_MODEL_H
#define LDA_MODEL

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "lda.h"
#include "lda-alpha.h"
#include "cokus.h"

#define myrand() (double) (((unsigned long) randomMT()) / 4294967296.)
#define NUM_INIT 1

void free_lda_model(lda_model*);
void save_lda_model(lda_model*, char*);
lda_model* new_lda_model(int, int);
lda_suffstats* new_lda_suffstats(lda_model* model);
void corpus_initialize_ss(lda_suffstats* ss, lda_model* model, corpus* c);
void random_initialize_ss(lda_suffstats* ss, lda_model* model);
void zero_initialize_ss(lda_suffstats* ss, lda_model* model);
void lda_mle(lda_model* model, lda_suffstats* ss, int estimate_alpha);
lda_model* load_lda_model(char* model_root);

#endif
