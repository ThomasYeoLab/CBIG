***************************
LATENT DIRICHLET ALLOCATION
***************************

David M. Blei
blei[at]cs.princeton.edu

(C) Copyright 2006, David M. Blei (blei [at] cs [dot] princeton [dot] edu)

This file is part of LDA-C.

LDA-C is free software; you can redistribute it and/or modify it under
the terms of the GNU General Public License as published by the Free
Software Foundation; either version 2 of the License, or (at your
option) any later version.

LDA-C is distributed in the hope that it will be useful, but WITHOUT
ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
for more details.

You should have received a copy of the GNU General Public License
along with this program; if not, write to the Free Software
Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
USA

------------------------------------------------------------------------

This is a C implementation of latent Dirichlet allocation (LDA), a
model of discrete data which is fully described in Blei et al. (2003)
(http://www.cs.berkeley.edu/~blei/papers/blei03a.pdf).

LDA is a hierarchical probabilistic model of documents.  Let \alpha be
a scalar and \beta_{1:K} be K distributions of words (called "topics").
As implemented here, a K topic LDA model assumes the following
generative process of an N word document:

          1. \theta | \alpha ~ Dirichlet(\alpha, ..., \alpha)

          2. for each word n = {1, ..., N}:

             a. Z_n | \theta ~ Mult(\theta)

             b. W_n | z_n, \beta ~ Mult(\beta_{z_n})

This code implements variational inference of \theta and z_{1:N} for a
document, and estimation of the topics \beta_{1:K} and Dirichlet
parameter \alpha.

------------------------------------------------------------------------


TABLE OF CONTENTS


A. COMPILING

B. TOPIC ESTIMATION

   1. SETTINGS FILE

   2. DATA FILE FORMAT

C. INFERENCE

D. PRINTING TOPICS

E. QUESTIONS, COMMENTS, PROBLEMS, UPDATE ANNOUNCEMENTS


------------------------------------------------------------------------

A. COMPILING

Type "make" in a shell.


------------------------------------------------------------------------

B. TOPIC ESTIMATION

Estimate the model by executing:

     lda est [alpha] [k] [settings] [data] [random/seeded/*] [directory]

The term [random/seeded/*] > describes how the topics will be
initialized.  "Random" initializes each topic randomly; "seeded"
initializes each topic to a distribution smoothed from a randomly
chosen document; or, you can specify a model name to load a
pre-existing model as the initial model (this is useful to continue EM
from where it left off).  To change the number of initial documents
used, edit lda-estimate.c.

The model (i.e., \alpha and \beta_{1:K}) and variational posterior
Dirichlet parameters will be saved in the specified directory every
ten iterations.  Additionally, there will be a log file for the
likelihood bound and convergence score at each iteration.  The
algorithm runs until that score is less than "em_convergence" (from
the settings file) or "em_max_iter" iterations are reached.  (To
change the lag between saved models, edit lda-estimate.c.)

The saved models are in two files:

     <iteration>.other contains alpha.

     <iteration>.beta contains the log of the topic distributions.
     Each line is a topic; in line k, each entry is log p(w | z=k)

The variational posterior Dirichlets are in:

     <iteration>.gamma

The settings file and data format are described below.


1. Settings file

See settings.txt for a sample.  See inf-settings.txt for an example of
a settings file for inference.  These are placeholder values; they
should be experimented with.

This is of the following form:

     var max iter [integer e.g., 10 or -1]
     var convergence [float e.g., 1e-8]
     em max iter [integer e.g., 100]
     em convergence [float e.g., 1e-5]
     alpha [fixed/estimate]

where the settings are

     [var max iter]

     The maximum number of iterations of coordinate ascent variational
     inference for a single document.  A value of -1 indicates "full"
     variational inference, until the variational convergence
     criterion is met.

     [var convergence]

     The convergence criteria for variational inference.  Stop if
     (score_old - score) / abs(score_old) is less than this value (or
     after the maximum number of iterations).  Note that the score is
     the lower bound on the likelihood for a particular document.

     [em max iter]

     The maximum number of iterations of variational EM.

     [em convergence]

     The convergence criteria for varitional EM.  Stop if (score_old -
     score) / abs(score_old) is less than this value (or after the
     maximum number of iterations).  Note that "score" is the lower
     bound on the likelihood for the whole corpus.

     [alpha]

     If set to [fixed] then alpha does not change from iteration to
     iteration.  If set to [estimate], then alpha is estimated along
     with the topic distributions.


2. Data format

Under LDA, the words of each document are assumed exchangeable.  Thus,
each document is succinctly represented as a sparse vector of word
counts. The data is a file where each line is of the form:

     [M] [term_1]:[count] [term_2]:[count] ...  [term_N]:[count]

where [M] is the number of unique terms in the document, and the
[count] associated with each term is how many times that term appeared
in the document.  Note that [term_1] is an integer which indexes the
term; it is not a string.


------------------------------------------------------------------------

C. INFERENCE

To perform inference on a different set of data (in the same format as
for estimation), execute:

     lda inf [settings] [model] [data] [name]

Variational inference is performed on the data using the model in
[model].* (see above).  Two files will be created : [name].gamma are
the variational Dirichlet parameters for each document;
[name].likelihood is the bound on the likelihood for each document.


------------------------------------------------------------------------

D. PRINTING TOPICS

The Python script topics.py lets you print out the top N
words from each topic in a .beta file.  Usage is:

     python topics.py <beta file> <vocab file> <n words>


------------------------------------------------------------------------

E. QUESTIONS, COMMENTS, PROBLEMS, AND UPDATE ANNOUNCEMENTS

Please join the topic-models mailing list,
topic-models@lists.cs.princeton.edu.

To join, go to http://lists.cs.princeton.edu and click on
"topic-models."
