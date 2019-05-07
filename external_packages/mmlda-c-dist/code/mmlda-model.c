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

#include "mmlda-model.h"

/*
 * compute MLE mmlda model from sufficient statistics
 *
 */

void mmlda_mle(mmlda_model* model, mmlda_suffstats* ss, int estimate_alpha)
{
    int k; int w;

    for (k = 0; k < model->num_topics; k++)
    {
        for (w = 0; w < model->num_terms1; w++)
        {
            if (ss->class_word1[k][w] > 0)
            {
                model->log_prob_w1[k][w] =
                    log(ss->class_word1[k][w]) -
                    log(ss->class_total1[k]);
            }
            else
                model->log_prob_w1[k][w] = -100;
        }
    }

    for (k = 0; k < model->num_topics; k++)
    {
        for (w = 0; w < model->num_terms2; w++)
        {
            if (ss->class_word2[k][w] > 0)
            {
                model->log_prob_w2[k][w] =
                    log(ss->class_word2[k][w]) -
                    log(ss->class_total2[k]);
            }
            else
                model->log_prob_w2[k][w] = -100;
        }
    }
    if (estimate_alpha == 1)
    {
        model->alpha = opt_alpha(ss->alpha_suffstats,
                                 ss->num_docs,
                                 model->num_topics);

        printf("new alpha = %5.5f\n", model->alpha);
    }
}

/*
 * allocate sufficient statistics
 *
 */

mmlda_suffstats* new_mmlda_suffstats(mmlda_model* model)
{
    int num_topics = model->num_topics;
    int num_terms1 = model->num_terms1;
    int num_terms2 = model->num_terms2;
    int i,j;

    mmlda_suffstats* ss = malloc(sizeof(mmlda_suffstats));
    ss->class_total1 = malloc(sizeof(double)*num_topics);
    ss->class_word1 = malloc(sizeof(double*)*num_topics);
    ss->class_total2 = malloc(sizeof(double)*num_topics);
    ss->class_word2 = malloc(sizeof(double*)*num_topics);
    for (i = 0; i < num_topics; i++)
    {
	ss->class_total1[i] = 0;
	ss->class_word1[i] = malloc(sizeof(double)*num_terms1);
	for (j = 0; j < num_terms1; j++)
	{
	    ss->class_word1[i][j] = 0;
	}
    }

    for (i = 0; i < num_topics; i++)
    {
    ss->class_total2[i] = 0;
    ss->class_word2[i] = malloc(sizeof(double)*num_terms2);
    for (j = 0; j < num_terms2; j++)
    {
        ss->class_word2[i][j] = 0;
    }
    }
    return(ss);
}


/*
 * various intializations for the sufficient statistics
 *
 */

void zero_initialize_ss(mmlda_suffstats* ss, mmlda_model* model)
{
    int k, w;
    for (k = 0; k < model->num_topics; k++)
    {
        ss->class_total1[k] = 0;
        for (w = 0; w < model->num_terms1; w++)
        {
            ss->class_word1[k][w] = 0;
        }
    }

    for (k = 0; k < model->num_topics; k++)
    {
        ss->class_total2[k] = 0;
        for (w = 0; w < model->num_terms2; w++)
        {
            ss->class_word2[k][w] = 0;
        }
    }
    ss->num_docs = 0;
    ss->alpha_suffstats = 0;
}


void random_initialize_ss(mmlda_suffstats* ss, mmlda_model* model)
{
    int num_topics = model->num_topics;
    int num_terms1 = model->num_terms1;
    int num_terms2 = model->num_terms2;
    int k, n;
    for (k = 0; k < num_topics; k++)
    {
        for (n = 0; n < num_terms1; n++)
        {
            ss->class_word1[k][n] += 1.0/num_terms1 + myrand();
            ss->class_total1[k] += ss->class_word1[k][n];
        }
    }

    for (k = 0; k < num_topics; k++)
    {
        for (n = 0; n < num_terms2; n++)
        {
            ss->class_word2[k][n] += 1.0/num_terms2 + myrand();
            ss->class_total2[k] += ss->class_word2[k][n];
        }
    }
}


void corpus_initialize_ss(mmlda_suffstats* ss, mmlda_model* model, corpus* c1, corpus* c2)
{
    int num_topics = model->num_topics;
    int i, k, d, n;
    document* doc1;
    document* doc2;

    for (k = 0; k < num_topics; k++)
    {
        for (i = 0; i < NUM_INIT; i++)
        {
            d = floor(myrand() * c1->num_docs);
            printf("initialized with document %d\n", d);
            doc1 = &(c1->docs[d]);
            for (n = 0; n < doc1->length; n++)
            {
                ss->class_word1[k][doc1->words[n]] += doc1->counts[n];
            }
        }
        for (n = 0; n < model->num_terms1; n++)
        {
            ss->class_word1[k][n] += 1.0;
            ss->class_total1[k] = ss->class_total1[k] + ss->class_word1[k][n];
        }
    }

    for (k = 0; k < num_topics; k++)
    {
        for (i = 0; i < NUM_INIT; i++)
        {
            d = floor(myrand() * c2->num_docs);
            printf("initialized with document %d\n", d);
            doc2 = &(c2->docs[d]);
            for (n = 0; n < doc2->length; n++)
            {
                ss->class_word2[k][doc2->words[n]] += doc2->counts[n];
            }
        }
        for (n = 0; n < model->num_terms2; n++)
        {
            ss->class_word2[k][n] += 1.0;
            ss->class_total2[k] = ss->class_total2[k] + ss->class_word2[k][n];
        }
    }
}

/*
 * allocate new mmlda model
 *
 */

mmlda_model* new_mmlda_model(int num_terms1, int num_terms2, int num_topics)
{
    int i,j;
    mmlda_model* model;

    model = malloc(sizeof(mmlda_model));
    model->num_topics = num_topics;
    model->num_terms1 = num_terms1;
    model->num_terms2 = num_terms2;
    model->alpha = 1.0;

    model->log_prob_w1 = malloc(sizeof(double*)*num_topics);
    for (i = 0; i < num_topics; i++)
    {
	model->log_prob_w1[i] = malloc(sizeof(double)*num_terms1);
	for (j = 0; j < num_terms1; j++)
	    model->log_prob_w1[i][j] = 0;
    }

    model->log_prob_w2 = malloc(sizeof(double*)*num_topics);
    for (i = 0; i < num_topics; i++)
    {
    model->log_prob_w2[i] = malloc(sizeof(double)*num_terms2);
    for (j = 0; j < num_terms2; j++)
        model->log_prob_w2[i][j] = 0;
    }
    return(model);
}


/*
 * deallocate new mmlda model
 *
 */

void free_mmlda_model(mmlda_model* model)
{
    int i;

    for (i = 0; i < model->num_topics; i++)
    {
	free(model->log_prob_w1[i]);
    }
    free(model->log_prob_w1);

    for (i = 0; i < model->num_topics; i++)
    {
    free(model->log_prob_w2[i]);
    }
    free(model->log_prob_w2);
}


/*
 * save an mmlda model
 *
 */

void save_mmlda_model(mmlda_model* model, char* model_root)
{
    char filename[200];
    FILE* fileptr;
    int i, j;

    sprintf(filename, "%s.beta1", model_root);
    fileptr = fopen(filename, "w");
    for (i = 0; i < model->num_topics; i++)
    {
	for (j = 0; j < model->num_terms1; j++)
	{
	    fprintf(fileptr, " %5.10f", model->log_prob_w1[i][j]);
	}
	fprintf(fileptr, "\n");
    }
    fclose(fileptr);

    sprintf(filename, "%s.beta2", model_root);
    fileptr = fopen(filename, "w");
    for (i = 0; i < model->num_topics; i++)
    {
    for (j = 0; j < model->num_terms2; j++)
    {
        fprintf(fileptr, " %5.10f", model->log_prob_w2[i][j]);
    }
    fprintf(fileptr, "\n");
    }
    fclose(fileptr);

    sprintf(filename, "%s.other", model_root);
    fileptr = fopen(filename, "w");
    fprintf(fileptr, "num_topics %d\n", model->num_topics);
    fprintf(fileptr, "num_terms1 %d\n", model->num_terms1);
    fprintf(fileptr, "num_terms2 %d\n", model->num_terms2);
    fprintf(fileptr, "alpha %5.10f\n", model->alpha);
    fclose(fileptr);
}


mmlda_model* load_mmlda_model(char* model_root)
{
    char filename[200];
    FILE* fileptr;
    int i, j, num_terms1, num_terms2, num_topics;
    float x, alpha;

    sprintf(filename, "%s.other", model_root);
    printf("loading %s\n", filename);
    fileptr = fopen(filename, "r");
    fscanf(fileptr, "num_topics %d\n", &num_topics);
    fscanf(fileptr, "num_terms1 %d\n", &num_terms1);
    fscanf(fileptr, "num_terms2 %d\n", &num_terms2);
    fscanf(fileptr, "alpha %f\n", &alpha);
    fclose(fileptr);

    mmlda_model* model = new_mmlda_model(num_terms1, num_terms2, num_topics);
    model->alpha = alpha;

    sprintf(filename, "%s.beta1", model_root);
    printf("loading %s\n", filename);
    fileptr = fopen(filename, "r");
    for (i = 0; i < num_topics; i++)
    {
        for (j = 0; j < num_terms1; j++)
        {
            fscanf(fileptr, "%f", &x);
            model->log_prob_w1[i][j] = x;
        }
    }
    fclose(fileptr);

    sprintf(filename, "%s.beta2", model_root);
    printf("loading %s\n", filename);
    fileptr = fopen(filename, "r");
    for (i = 0; i < num_topics; i++)
    {
        for (j = 0; j < num_terms2; j++)
        {
            fscanf(fileptr, "%f", &x);
            model->log_prob_w2[i][j] = x;
        }
    }
    fclose(fileptr);
    return(model);
}
