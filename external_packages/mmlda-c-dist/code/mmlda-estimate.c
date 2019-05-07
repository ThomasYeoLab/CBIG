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

#include "mmlda-estimate.h"

/*
 * perform inference on a document and update sufficient statistics
 *
 */

double doc_e_step(document* doc1, document* doc2, double* gamma, double** phi1, double** phi2,
                  mmlda_model* model, mmlda_suffstats* ss)
{
    double likelihood;
    int n, k;

    // posterior inference

    likelihood = mmlda_inference(doc1, doc2, model, gamma, phi1, phi2);

    // update sufficient statistics

    double gamma_sum = 0;
    for (k = 0; k < model->num_topics; k++)
    {
        gamma_sum += gamma[k];
        ss->alpha_suffstats += digamma(gamma[k]);
    }
    ss->alpha_suffstats -= model->num_topics * digamma(gamma_sum);

    for (n = 0; n < doc1->length; n++)
    {
        for (k = 0; k < model->num_topics; k++)
        {
            ss->class_word1[k][doc1->words[n]] += doc1->counts[n]*phi1[n][k];
            ss->class_total1[k] += doc1->counts[n]*phi1[n][k];
        }
    }

    for (n = 0; n < doc2->length; n++)
    {
        for (k = 0; k < model->num_topics; k++)
        {
            ss->class_word2[k][doc2->words[n]] += doc2->counts[n]*phi2[n][k];
            ss->class_total2[k] += doc2->counts[n]*phi2[n][k];
        }
    }

    ss->num_docs = ss->num_docs + 1;

    return(likelihood);
}


/*
 * writes the word assignments line for a document to a file
 *
 */

void write_word_assignment(FILE* f, document* doc, double** phi, mmlda_model* model)
{
    int n;

    fprintf(f, "%03d", doc->length);
    for (n = 0; n < doc->length; n++)
    {
        fprintf(f, " %04d:%02d",
                doc->words[n], argmax(phi[n], model->num_topics));
    }
    fprintf(f, "\n");
    fflush(f);
}


/*
 * saves the gamma parameters of the current dataset
 *
 */

void save_gamma(char* filename, double** gamma, int num_docs, int num_topics)
{
    FILE* fileptr;
    int d, k;
    fileptr = fopen(filename, "w");

    for (d = 0; d < num_docs; d++)
    {
	fprintf(fileptr, "%5.10f", gamma[d][0]);
	for (k = 1; k < num_topics; k++)
	{
	    fprintf(fileptr, " %5.10f", gamma[d][k]);
	}
	fprintf(fileptr, "\n");
    }
    fclose(fileptr);
}


/*
 * run_em
 *
 */

void run_em(char* start, char* directory, corpus* corpus1, corpus* corpus2)
{

    int d, n;
    mmlda_model *model = NULL;
    double **var_gamma, **phi1, **phi2;

    // allocate variational parameters

    var_gamma = malloc(sizeof(double*)*(corpus1->num_docs));
    for (d = 0; d < corpus1->num_docs; d++)
	var_gamma[d] = malloc(sizeof(double) * NTOPICS);

    int max_length1 = max_corpus_length(corpus1);
    phi1 = malloc(sizeof(double*)*max_length1);
    for (n = 0; n < max_length1; n++)
	phi1[n] = malloc(sizeof(double) * NTOPICS);

    int max_length2 = max_corpus_length(corpus2);
    phi2 = malloc(sizeof(double*)*max_length2);
    for (n = 0; n < max_length2; n++)
    phi2[n] = malloc(sizeof(double) * NTOPICS);
    // initialize model

    char filename[200];
    char filename1[200];
    char filename2[200];

    mmlda_suffstats* ss = NULL;
    if (strcmp(start, "seeded")==0)
    {
        model = new_mmlda_model(corpus1->num_terms, corpus2->num_terms, NTOPICS);
        ss = new_mmlda_suffstats(model);
        corpus_initialize_ss(ss, model, corpus1, corpus2);
        mmlda_mle(model, ss, 0);
        model->alpha = INITIAL_ALPHA;
    }
    else if (strcmp(start, "random")==0)
    {
        model = new_mmlda_model(corpus1->num_terms, corpus2->num_terms, NTOPICS);
        ss = new_mmlda_suffstats(model);
        random_initialize_ss(ss, model);
        mmlda_mle(model, ss, 0);
        model->alpha = INITIAL_ALPHA;
    }
    else
    {
        model = load_mmlda_model(start);
        ss = new_mmlda_suffstats(model);
    }

    sprintf(filename,"%s/000",directory);
    save_mmlda_model(model, filename);

    // run expectation maximization

    int i = 0;
    double likelihood, likelihood_old = 0, converged = 1;
    sprintf(filename, "%s/likelihood.dat", directory);
    FILE* likelihood_file = fopen(filename, "w");

    while (((converged < 0) || (converged > EM_CONVERGED) || (i <= 2)) && (i <= EM_MAX_ITER))
    {
        i++; printf("**** em iteration %d ****\n", i);
        likelihood = 0;
        zero_initialize_ss(ss, model);

        // e-step

        for (d = 0; d < corpus1->num_docs; d++)
        {
            if ((d % 20) == 0) printf("document %d\n",d);
            likelihood += doc_e_step(&(corpus1->docs[d]),
                                     &(corpus2->docs[d]),
                                     var_gamma[d],
                                     phi1,
                                     phi2,
                                     model,
                                     ss);
        }

        // m-step

        mmlda_mle(model, ss, ESTIMATE_ALPHA);

        // check for convergence

        converged = (likelihood_old - likelihood) / (likelihood_old);
        if (converged < 0) VAR_MAX_ITER = VAR_MAX_ITER * 2;
        likelihood_old = likelihood;

        // output model and likelihood

        fprintf(likelihood_file, "%10.10f\t%5.5e\n", likelihood, converged);
        fflush(likelihood_file);
        if ((i % LAG) == 0)
        {
            sprintf(filename,"%s/%03d",directory, i);
            save_mmlda_model(model, filename);
            sprintf(filename,"%s/%03d.gamma",directory, i);
            save_gamma(filename, var_gamma, corpus1->num_docs, model->num_topics);
        }
    }

    // output the final model

    sprintf(filename,"%s/final",directory);
    save_mmlda_model(model, filename);
    sprintf(filename,"%s/final.gamma",directory);
    save_gamma(filename, var_gamma, corpus1->num_docs, model->num_topics);

    // output the word assignments (for visualization)

    sprintf(filename1, "%s/word-assignments1.dat", directory);
    FILE* w_asgn_file1 = fopen(filename1, "w");
    sprintf(filename2, "%s/word-assignments2.dat", directory);
    FILE* w_asgn_file2 = fopen(filename2, "w");
    for (d = 0; d < corpus1->num_docs; d++)
    {
        if ((d % 20) == 0) printf("final e step document %d\n",d);
        likelihood += mmlda_inference(&(corpus1->docs[d]), &(corpus2->docs[d]), model, var_gamma[d], phi1, phi2);
        write_word_assignment(w_asgn_file1, &(corpus1->docs[d]), phi1, model);
        write_word_assignment(w_asgn_file2, &(corpus2->docs[d]), phi2, model);
    }
    fclose(w_asgn_file1);
    fclose(w_asgn_file2);
    fclose(likelihood_file);
}


/*
 * read settings.
 *
 */

void read_settings(char* filename)
{
    FILE* fileptr;
    char alpha_action[100];
    fileptr = fopen(filename, "r");
    fscanf(fileptr, "var max iter %d\n", &VAR_MAX_ITER);
    fscanf(fileptr, "var convergence %f\n", &VAR_CONVERGED);
    fscanf(fileptr, "em max iter %d\n", &EM_MAX_ITER);
    fscanf(fileptr, "em convergence %f\n", &EM_CONVERGED);
    fscanf(fileptr, "alpha %s", alpha_action);
    if (strcmp(alpha_action, "fixed")==0)
    {
	ESTIMATE_ALPHA = 0;
    }
    else
    {
	ESTIMATE_ALPHA = 1;
    }
    fclose(fileptr);
}


/*
 * inference based on two modalities
 *
 */

void infer(char* model_root, char* save, corpus* corpus1, corpus* corpus2)
{
    FILE* fileptr;
    char filename[200];
    int i, d, n;
    mmlda_model *model;
    double **var_gamma, likelihood, **phi1, **phi2;
    document *doc1, *doc2;

    model = load_mmlda_model(model_root);
    var_gamma = malloc(sizeof(double*)*(corpus1->num_docs));
    for (i = 0; i < corpus1->num_docs; i++)
	var_gamma[i] = malloc(sizeof(double)*model->num_topics);
    sprintf(filename, "%s-lda-lhood.dat", save);
    fileptr = fopen(filename, "w");
    for (d = 0; d < corpus1->num_docs; d++)
    {
	if (((d % 1) == 0) && (d>0)) printf("document %d\n",d);

	doc1 = &(corpus1->docs[d]);
	phi1 = (double**) malloc(sizeof(double*) * doc1->length);
	for (n = 0; n < doc1->length; n++)
	    phi1[n] = (double*) malloc(sizeof(double) * model->num_topics);
	
    doc2 = &(corpus2->docs[d]);
    phi2 = (double**) malloc(sizeof(double*) * doc2->length);
    for (n = 0; n < doc2->length; n++)
        phi2[n] = (double*) malloc(sizeof(double) * model->num_topics);
    likelihood = mmlda_inference(doc1, doc2, model, var_gamma[d], phi1, phi2);

	fprintf(fileptr, "%5.5f\n", likelihood);
    }
    fclose(fileptr);
    sprintf(filename, "%s-gamma.dat", save);
    save_gamma(filename, var_gamma, corpus1->num_docs, model->num_topics);
}

/*
 * inference based on first modality only
 *
 */

void infer1(char* model_root, char* save, corpus* corpus1)
{
    FILE* fileptr;
    char filename[200];
    int i, d, n;
    mmlda_model *model;
    double **var_gamma, likelihood, **phi1;
    document *doc1;

    model = load_mmlda_model(model_root);
    var_gamma = malloc(sizeof(double*)*(corpus1->num_docs));
    for (i = 0; i < corpus1->num_docs; i++)
	var_gamma[i] = malloc(sizeof(double)*model->num_topics);
    sprintf(filename, "%s-lda-lhood.dat", save);
    fileptr = fopen(filename, "w");
    for (d = 0; d < corpus1->num_docs; d++)
    {
	if (((d % 1) == 0) && (d>0)) printf("document %d\n",d);

	doc1 = &(corpus1->docs[d]);
	phi1 = (double**) malloc(sizeof(double*) * doc1->length);
	for (n = 0; n < doc1->length; n++)
	    phi1[n] = (double*) malloc(sizeof(double) * model->num_topics);
	
    likelihood = mmlda_inference1(doc1, model, var_gamma[d], phi1);

	fprintf(fileptr, "%5.5f\n", likelihood);
    }
    fclose(fileptr);
    sprintf(filename, "%s-gamma.dat", save);
    save_gamma(filename, var_gamma, corpus1->num_docs, model->num_topics);
}
/*
 * inference based on second modality only
 *
 */

void infer2(char* model_root, char* save, corpus* corpus2)
{
    FILE* fileptr;
    char filename[200];
    int i, d, n;
    mmlda_model *model;
    double **var_gamma, likelihood, **phi2;
    document *doc2;

    model = load_mmlda_model(model_root);
    var_gamma = malloc(sizeof(double*)*(corpus2->num_docs));
    for (i = 0; i < corpus2->num_docs; i++)
	var_gamma[i] = malloc(sizeof(double)*model->num_topics);
    sprintf(filename, "%s-lda-lhood.dat", save);
    fileptr = fopen(filename, "w");
    for (d = 0; d < corpus2->num_docs; d++)
    {
	if (((d % 1) == 0) && (d>0)) printf("document %d\n",d);

    doc2 = &(corpus2->docs[d]);
    phi2 = (double**) malloc(sizeof(double*) * doc2->length);
    for (n = 0; n < doc2->length; n++)
        phi2[n] = (double*) malloc(sizeof(double) * model->num_topics);
    likelihood = mmlda_inference2(doc2, model, var_gamma[d], phi2);

	fprintf(fileptr, "%5.5f\n", likelihood);
    }
    fclose(fileptr);
    sprintf(filename, "%s-gamma.dat", save);
    save_gamma(filename, var_gamma, corpus2->num_docs, model->num_topics);
}


/*
 * main
 *
 */

int main(int argc, char* argv[])
{
    // (est / inf) alpha k settings data (random / seed/ model) (directory / out)

    corpus* corpus1;
    corpus* corpus2;

    long t1;

    if (argc == 10)
    {
        t1 = atol(argv[9]) * 2;
    }
    else
    {
        (void) time(&t1);
    }
    
    seedMT(t1);
    // seedMT(4357U);

    if (argc > 1)
    {
        if (strcmp(argv[1], "est")==0)
        {
            INITIAL_ALPHA = atof(argv[2]);
            NTOPICS = atoi(argv[3]);
            read_settings(argv[4]);
            corpus1 = read_data(argv[5]);
            corpus2 = read_data(argv[6]);
            assert(corpus1->num_docs==corpus2->num_docs);

            make_directory(argv[8]);
            run_em(argv[7], argv[8], corpus1, corpus2);
        }
        if (strcmp(argv[1], "inf")==0)
        {
            read_settings(argv[2]);
            corpus1 = read_data(argv[4]);
            corpus2 = read_data(argv[5]);
            infer(argv[3], argv[6], corpus1, corpus2);
        }
        if (strcmp(argv[1], "inf1")==0)
        {
        	read_settings(argv[2]);
        	corpus1 = read_data(argv[4]);
        	infer1(argv[3], argv[5], corpus1);
        }
        if (strcmp(argv[1], "inf2")==0)
        {
        	read_settings(argv[2]);
        	corpus2 = read_data(argv[4]);
        	infer2(argv[3], argv[5], corpus2);
        }
    }
    else
    {
        printf("usage : mmlda est [initial alpha] [k] [settings] [data1] [data2] [random/seeded/*] [directory] [seed]\n");
        printf("        mmlda inf [settings] [model] [data1] [data2] [name]\n");
    	printf("        mmlda inf1 [settings] [model] [data1] [name]\n");
    	printf("        mmlda inf2 [settings] [model] [data2] [name]\n");
    }
    return(0);
}
