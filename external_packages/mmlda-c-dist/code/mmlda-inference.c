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

#include "mmlda-inference.h"

/*
 * variational inference
 *
 */

double mmlda_inference(document* doc1, document* doc2, mmlda_model* model, double* var_gamma, double** phi1, double** phi2)
{
    double converged = 1;
    double phisum1 = 0, phisum2 = 0, likelihood = 0;
    double likelihood_old = 0, oldphi1[model->num_topics], oldphi2[model->num_topics];
    int k, n, var_iter;
    double digamma_gam[model->num_topics];

    // compute posterior dirichlet

    for (k = 0; k < model->num_topics; k++)
    {
        var_gamma[k] = model->alpha + ((doc1->total + doc2->total)/((double) model->num_topics));
        digamma_gam[k] = digamma(var_gamma[k]);
        for (n = 0; n < doc1->length; n++)
            phi1[n][k] = 1.0/model->num_topics;
        for (n = 0; n < doc2->length; n++)
            phi2[n][k] = 1.0/model->num_topics;

    }
    var_iter = 0;

    while ((converged > VAR_CONVERGED) &&
           ((var_iter < VAR_MAX_ITER) || (VAR_MAX_ITER == -1)))
    {
	var_iter++;
	   for (n = 0; n < doc1->length; n++)
	   {
            phisum1 = 0;
            for (k = 0; k < model->num_topics; k++)
            {
                oldphi1[k] = phi1[n][k];
                phi1[n][k] =
                    digamma_gam[k] +
                    model->log_prob_w1[k][doc1->words[n]];

                if (k > 0)
                    phisum1 = log_sum(phisum1, phi1[n][k]);
                else
                    phisum1 = phi1[n][k]; // note, phi is in log space
            }

            for (k = 0; k < model->num_topics; k++)
            {
                phi1[n][k] = exp(phi1[n][k] - phisum1);
                var_gamma[k] =
                    var_gamma[k] + doc1->counts[n]*(phi1[n][k] - oldphi1[k]);
                // !!! a lot of extra digamma's here because of how we're computing it
                // !!! but its more automatically updated too.
                //digamma_gam[k] = digamma(var_gamma[k]);
            }
        }

        for (n = 0; n < doc2->length; n++)
        {
            phisum2 = 0;
            for (k = 0; k < model->num_topics; k++)
            {
                oldphi2[k] = phi2[n][k];
                phi2[n][k] =
                    digamma_gam[k] +
                    model->log_prob_w2[k][doc2->words[n]];

                if (k > 0)
                    phisum2 = log_sum(phisum2, phi2[n][k]);
                else
                    phisum2 = phi2[n][k]; // note, phi is in log space
            }

            for (k = 0; k < model->num_topics; k++)
            {
                phi2[n][k] = exp(phi2[n][k] - phisum2);
                var_gamma[k] =
                    var_gamma[k] + doc2->counts[n]*(phi2[n][k] - oldphi2[k]);
                // !!! a lot of extra digamma's here because of how we're computing it
                // !!! but its more automatically updated too.
                //digamma_gam[k] = digamma(var_gamma[k]);
            }
        }
        for (k = 0; k < model->num_topics; k++)
        {
            digamma_gam[k] = digamma(var_gamma[k]);
        }

        likelihood = compute_likelihood(doc1, doc2, model, phi1, phi2, var_gamma);
        assert(!isnan(likelihood));
        converged = (likelihood_old - likelihood) / likelihood_old;
        likelihood_old = likelihood;

        // printf("[LDA INF] %8.5f %1.3e\n", likelihood, converged);
    }
    return(likelihood);
}


/*
 * compute likelihood bound
 *
 */

double
compute_likelihood(document* doc1, document* doc2, mmlda_model* model, double** phi1, double** phi2, double* var_gamma)
{
    double likelihood = 0, digsum = 0, var_gamma_sum = 0, dig[model->num_topics];
    int k, n;

    for (k = 0; k < model->num_topics; k++)
    {
	dig[k] = digamma(var_gamma[k]);
	var_gamma_sum += var_gamma[k];
    }
    digsum = digamma(var_gamma_sum);

    likelihood =
	lgamma(model->alpha * model -> num_topics)
	- model -> num_topics * lgamma(model->alpha)
	- (lgamma(var_gamma_sum));

    for (k = 0; k < model->num_topics; k++)
    {
	likelihood +=
	    (model->alpha - 1)*(dig[k] - digsum) + lgamma(var_gamma[k])
	    - (var_gamma[k] - 1)*(dig[k] - digsum);

	for (n = 0; n < doc1->length; n++)
	{
            if (phi1[n][k] > 0)
            {
                likelihood += doc1->counts[n]*
                    (phi1[n][k]*((dig[k] - digsum) - log(phi1[n][k])
                                + model->log_prob_w1[k][doc1->words[n]]));
            }
    }
    

    for (n = 0; n < doc2->length; n++)
    {
            if (phi2[n][k] > 0)
            {
                likelihood += doc2->counts[n]*
                    (phi2[n][k]*((dig[k] - digsum) - log(phi2[n][k])
                                + model->log_prob_w2[k][doc2->words[n]]));
            }
        
    }

    }   
    return(likelihood);
}

/*
 * variational inference2
 *
 */

double mmlda_inference2(document* doc2, mmlda_model* model, double* var_gamma, double** phi2)
{
    double converged = 1;
    double phisum2 = 0, likelihood = 0;
    double likelihood_old = 0, oldphi2[model->num_topics];
    int k, n, var_iter;
    double digamma_gam[model->num_topics];

    // compute posterior dirichlet

    for (k = 0; k < model->num_topics; k++)
    {
        var_gamma[k] = model->alpha + ((doc2->total)/((double) model->num_topics));
        digamma_gam[k] = digamma(var_gamma[k]);
        for (n = 0; n < doc2->length; n++)
            phi2[n][k] = 1.0/model->num_topics;        
    }
    var_iter = 0;

    while ((converged > VAR_CONVERGED) &&
           ((var_iter < VAR_MAX_ITER) || (VAR_MAX_ITER == -1)))
    {
    var_iter++;
       for (n = 0; n < doc2->length; n++)
       {
            phisum2 = 0;
            for (k = 0; k < model->num_topics; k++)
            {
                oldphi2[k] = phi2[n][k];
                phi2[n][k] =
                    digamma_gam[k] +
                    model->log_prob_w2[k][doc2->words[n]];

                if (k > 0)
                    phisum2 = log_sum(phisum2, phi2[n][k]);
                else
                    phisum2 = phi2[n][k]; // note, phi is in log space
            }

            for (k = 0; k < model->num_topics; k++)
            {
                phi2[n][k] = exp(phi2[n][k] - phisum2);
                var_gamma[k] =
                    var_gamma[k] + doc2->counts[n]*(phi2[n][k] - oldphi2[k]);
                // !!! a lot of extra digamma's here because of how we're computing it
                // !!! but its more automatically updated too.
                //digamma_gam[k] = digamma(var_gamma[k]);
            }
        }
         for (k = 0; k < model->num_topics; k++)
        {
            digamma_gam[k] = digamma(var_gamma[k]);
        }

        likelihood = compute_likelihood2(doc2, model, phi2, var_gamma);
        assert(!isnan(likelihood));
        converged = (likelihood_old - likelihood) / likelihood_old;
        likelihood_old = likelihood;

        // printf("[LDA INF] %8.5f %1.3e\n", likelihood, converged);
    }
    return(likelihood);
}

double
compute_likelihood2(document* doc2, mmlda_model* model, double** phi2, double* var_gamma)
{
    double likelihood = 0, digsum = 0, var_gamma_sum = 0, dig[model->num_topics];
    int k, n;

    for (k = 0; k < model->num_topics; k++)
    {
    dig[k] = digamma(var_gamma[k]);
    var_gamma_sum += var_gamma[k];
    }
    digsum = digamma(var_gamma_sum);

    likelihood =
    lgamma(model->alpha * model -> num_topics)
    - model -> num_topics * lgamma(model->alpha)
    - (lgamma(var_gamma_sum));

    for (k = 0; k < model->num_topics; k++)
    {
    likelihood +=
        (model->alpha - 1)*(dig[k] - digsum) + lgamma(var_gamma[k])
        - (var_gamma[k] - 1)*(dig[k] - digsum);

    for (n = 0; n < doc2->length; n++)
    {
            if (phi2[n][k] > 0)
            {
                likelihood += doc2->counts[n]*
                    (phi2[n][k]*((dig[k] - digsum) - log(phi2[n][k])
                                + model->log_prob_w2[k][doc2->words[n]]));
            }
    }
    
    }   
    return(likelihood);
}

/*
 * variational inference1
 *
 */

double mmlda_inference1(document* doc1, mmlda_model* model, double* var_gamma, double** phi1)
{
    double converged = 1;
    double phisum1 = 0, likelihood = 0;
    double likelihood_old = 0, oldphi1[model->num_topics];
    int k, n, var_iter;
    double digamma_gam[model->num_topics];

    // compute posterior dirichlet

    for (k = 0; k < model->num_topics; k++)
    {
        var_gamma[k] = model->alpha + ((doc1->total)/((double) model->num_topics));
        digamma_gam[k] = digamma(var_gamma[k]);
        for (n = 0; n < doc1->length; n++)
            phi1[n][k] = 1.0/model->num_topics;        
    }
    var_iter = 0;

    while ((converged > VAR_CONVERGED) &&
           ((var_iter < VAR_MAX_ITER) || (VAR_MAX_ITER == -1)))
    {
    var_iter++;
       for (n = 0; n < doc1->length; n++)
       {
            phisum1 = 0;
            for (k = 0; k < model->num_topics; k++)
            {
                oldphi1[k] = phi1[n][k];
                phi1[n][k] =
                    digamma_gam[k] +
                    model->log_prob_w1[k][doc1->words[n]];

                if (k > 0)
                    phisum1 = log_sum(phisum1, phi1[n][k]);
                else
                    phisum1 = phi1[n][k]; // note, phi is in log space
            }

            for (k = 0; k < model->num_topics; k++)
            {
                phi1[n][k] = exp(phi1[n][k] - phisum1);
                var_gamma[k] =
                    var_gamma[k] + doc1->counts[n]*(phi1[n][k] - oldphi1[k]);
                // !!! a lot of extra digamma's here because of how we're computing it
                // !!! but its more automatically updated too.
                //digamma_gam[k] = digamma(var_gamma[k]);
            }
        }

        
        for (k = 0; k < model->num_topics; k++)
        {
            digamma_gam[k] = digamma(var_gamma[k]);
        }

        likelihood = compute_likelihood1(doc1, model, phi1, var_gamma);
        assert(!isnan(likelihood));
        converged = (likelihood_old - likelihood) / likelihood_old;
        likelihood_old = likelihood;

        // printf("[LDA INF] %8.5f %1.3e\n", likelihood, converged);
    }
    return(likelihood);
}


/*
 * compute likelihood bound
 *
 */

double
compute_likelihood1(document* doc1, mmlda_model* model, double** phi1, double* var_gamma)
{
    double likelihood = 0, digsum = 0, var_gamma_sum = 0, dig[model->num_topics];
    int k, n;

    for (k = 0; k < model->num_topics; k++)
    {
    dig[k] = digamma(var_gamma[k]);
    var_gamma_sum += var_gamma[k];
    }
    digsum = digamma(var_gamma_sum);

    likelihood =
    lgamma(model->alpha * model -> num_topics)
    - model -> num_topics * lgamma(model->alpha)
    - (lgamma(var_gamma_sum));

    for (k = 0; k < model->num_topics; k++)
    {
    likelihood +=
        (model->alpha - 1)*(dig[k] - digsum) + lgamma(var_gamma[k])
        - (var_gamma[k] - 1)*(dig[k] - digsum);

    for (n = 0; n < doc1->length; n++)
    {
            if (phi1[n][k] > 0)
            {
                likelihood += doc1->counts[n]*
                    (phi1[n][k]*((dig[k] - digsum) - log(phi1[n][k])
                                + model->log_prob_w1[k][doc1->words[n]]));
            }
    }
    
    }   
    return(likelihood);
}