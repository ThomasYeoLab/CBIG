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

#include "polarlda-model.h"

/*
 * compute MLE polarlda model from sufficient statistics
 *
 */

void polarlda_mle(polarlda_model* model, polarlda_suffstats* ss, int estimate_alpha)
{
    int k, w;

    for (k = 0; k < model->num_topics; k++)
    {
        for (w = 0; w < model->num_terms; w++)
        {
            if (ss->class_word[k][w] > 0) // word w appears for topic k
            {
                model->log_prob_w[k][w] = log(ss->class_word[k][w]) - log(ss->class_total[k]);
                
                if (ss->class_word_pos[k][w] > 0) // word w appears positively for topic k
                {
                    model->log_prob_pos[k][w] = log(ss->class_word_pos[k][w]) - log(ss->class_word[k][w]);
                    if (model->log_prob_pos[k][w] > -0.000001) // percision consideration; effectively, log_prob_pos = 0
                    {
                        model->log_prob_pos[k][w] = -0.000001; // prob_pos = exp(-0.000001) = 0.999999, max prob_pos possible
                    }
                }
                else // word w appears for topic k, but never positively
                {
                    model->log_prob_pos[k][w] = -100; // still a small probability of being positive
                }
            }
            else // word w never appears for topic k
            {
                model->log_prob_w[k][w] = -100; // still a small probability of appearing
                
                model->log_prob_pos[k][w] = log(0.5); // 0.5 to be positive
            }
            
            // Sanity check on beta and rho
            // assert(exp(model->log_prob_w[k][w]) > 0 && exp(model->log_prob_w[k][w]) < 1);
            // assert(exp(model->log_prob_pos[k][w]) > 0 && exp(model->log_prob_pos[k][w]) < 1);
        }
    }

    if (estimate_alpha == 1)
    {
        model->alpha = opt_alpha(ss->alpha_suffstats, ss->num_docs, model->num_topics);

        printf("new alpha = %5.5f\n", model->alpha);
    }
}

/*
 * allocate sufficient statistics
 *
 */

polarlda_suffstats* new_polarlda_suffstats(polarlda_model* model)
{
    int num_topics = model->num_topics;
    int num_terms = model->num_terms;
    int i, j;

    polarlda_suffstats* ss = malloc(sizeof(polarlda_suffstats));
    ss->class_total = malloc(sizeof(double) * num_topics);
    ss->class_word = malloc(sizeof(double*) * num_topics);
    ss->class_word_pos = malloc(sizeof(double*) * num_topics);
    for (i = 0; i < num_topics; i++)
    {
        ss->class_total[i] = 0;
        ss->class_word[i] = malloc(sizeof(double) * num_terms);
        ss->class_word_pos[i] = malloc(sizeof(double) * num_terms);
        for (j = 0; j < num_terms; j++)
        {
	   ss->class_word[i][j] = 0;
           ss->class_word_pos[i][j] = 0;
        }
    }
    return(ss);
}


/*
 * various intializations for the sufficient statistics
 *
 */

void zero_initialize_ss(polarlda_suffstats* ss, polarlda_model* model)
{
    int k, w;
    for (k = 0; k < model->num_topics; k++)
    {
        ss->class_total[k] = 0;
        for (w = 0; w < model->num_terms; w++)
        {
            ss->class_word[k][w] = 0;
            ss->class_word_pos[k][w] = 0;
        }
    }
    ss->num_docs = 0;
    ss->alpha_suffstats = 0;
}


void random_initialize_ss(polarlda_suffstats* ss, polarlda_model* model)
{
    int num_topics = model->num_topics;
    int num_terms = model->num_terms;
    int k, n;
    double x, y;
    for (k = 0; k < num_topics; k++)
    {
        for (n = 0; n < num_terms; n++)
        {
            x = 0.5 * (1.0 / num_terms + myrand());
            y = 0.5 * (1.0 / num_terms + myrand());
            ss->class_word_pos[k][n] = x;
            ss->class_word[k][n] = x + y;
            ss->class_total[k] += ss->class_word[k][n];
        }
    }
}


void corpus_initialize_ss(polarlda_suffstats* ss, polarlda_model* model, corpus* c)
{
    int num_topics = model->num_topics;
    int i, k, d, n;
    document* doc;

    for (k = 0; k < num_topics; k++)
    {
        for (i = 0; i < NUM_INIT; i++)
        {
            d = floor(myrand() * c->num_docs);
            printf("initialized with document %d\n", d);
            doc = &(c->docs[d]);
            for (n = 0; n < doc->length; n++)
            {
                ss->class_word[k][doc->words[n]] += doc->counts[n];
                ss->class_word_pos[k][doc->words[n]] += doc->counts[n] * doc->polarities[n];
            }
        }
        for (n = 0; n < model->num_terms; n++)
        {
            ss->class_word[k][n] += 1.0;
            ss->class_total[k] = ss->class_total[k] + ss->class_word[k][n];
        }
    }
}

/*
 * allocate new polarlda model
 *
 */

polarlda_model* new_polarlda_model(int num_terms, int num_topics)
{
    int i,j;
    polarlda_model* model;

    model = malloc(sizeof(polarlda_model));
    model->num_topics = num_topics;
    model->num_terms = num_terms;
    model->alpha = 1.0;

    // initialize beta and rho
    model->log_prob_w = malloc(sizeof(double*) * num_topics);
    model->log_prob_pos = malloc(sizeof(double*) * num_topics);
    for (i = 0; i < num_topics; i++)
    {
        model->log_prob_w[i] = malloc(sizeof(double) * num_terms);
        model->log_prob_pos[i] = malloc(sizeof(double) * num_terms);
        for (j = 0; j < num_terms; j++)
        {
            model->log_prob_w[i][j] = 0;
            model->log_prob_pos[i][j] = 0;
        }
    }

    return(model);
}


/*
 * deallocate new polarlda model
 *
 */

void free_polarlda_model(polarlda_model* model)
{
    int i;

    for (i = 0; i < model->num_topics; i++)
    {
	free(model->log_prob_w[i]);
        free(model->log_prob_pos[i]);
    }
    free(model->log_prob_w);
    free(model->log_prob_pos);
}


/*
 * save a polarlda model
 *
 */

void save_polarlda_model(polarlda_model* model, char* model_root)
{
    char filename[200], filename1[200], filename2[200];
    FILE *fileptr, *fileptr1, *fileptr2;
    int i, j;

    sprintf(filename1, "%s.beta", model_root);
    sprintf(filename2, "%s.rho", model_root);
    fileptr1 = fopen(filename1, "w");
    fileptr2 = fopen(filename2, "w");
    for (i = 0; i < model->num_topics; i++)
    {
	   for (j = 0; j < model->num_terms; j++)
	   {
	   	fprintf(fileptr1, " %5.10f", model->log_prob_w[i][j]);
           	fprintf(fileptr2, " %5.10f", model->log_prob_pos[i][j]);
	   }
	   fprintf(fileptr1, "\n");
           fprintf(fileptr2, "\n");
    }
    fclose(fileptr1);
    fclose(fileptr2);

    sprintf(filename, "%s.other", model_root);
    fileptr = fopen(filename, "w");
    fprintf(fileptr, "num_topics %d\n", model->num_topics);
    fprintf(fileptr, "num_terms %d\n", model->num_terms);
    fprintf(fileptr, "alpha %5.10f\n", model->alpha);
    fclose(fileptr);
}


polarlda_model* load_polarlda_model(char* model_root)
{
    char filename[200], filename1[200], filename2[200];
    FILE *fileptr, *fileptr1, *fileptr2;
    int i, j, num_terms, num_topics;
    double x, alpha;

    sprintf(filename, "%s.other", model_root);
    printf("loading %s\n", filename);
    fileptr = fopen(filename, "r");
    fscanf(fileptr, "num_topics %d\n", &num_topics);
    fscanf(fileptr, "num_terms %d\n", &num_terms);
    fscanf(fileptr, "alpha %lf\n", &alpha);
    fclose(fileptr);

    polarlda_model* model = new_polarlda_model(num_terms, num_topics);
    model->alpha = alpha;

    sprintf(filename1, "%s.beta", model_root);
    sprintf(filename2, "%s.rho", model_root);
    printf("loading %s\n", filename);
    fileptr1 = fopen(filename1, "r");
    fileptr2 = fopen(filename2, "r");
    for (i = 0; i < num_topics; i++)
    {
        for (j = 0; j < num_terms; j++)
        {
            fscanf(fileptr1, "%lf", &x);
            model->log_prob_w[i][j] = x;
            fscanf(fileptr2, "%lf", &x);
            model->log_prob_pos[i][j] = x;
        }
    }
    fclose(fileptr1);
    fclose(fileptr2);

    return(model);
}
