## Polar Latent Dirichlet Allocation (polarLDA)
Polar latent dirichlet allocation (polarLDA) is an extension of latent dirichlet allocation (LDA; Blei et al., 2003). If you have not read the LDA paper, we suggest you read it before going through polarLDA. C implementation of polarLDA is based on Blei's C code for LDA, and is located in `code` folder. For details on polarLDA including the graphical model and detailed derivations, please refer to `doc/polarLDA.pdf`.

---

## Copyright
(C) Copyright 2018, Nanbo Sun (nanbosun@u.nus.edu), Xiuming Zhang

This file is part of POLARLDA-C.

POLARLDA-C is free software; you can redistribute it and/or modify it under
the terms of the GNU General Public License as published by the Free
Software Foundation; either version 2 of the License, or (at your
option) any later version.

POLARLDA-C is distributed in the hope that it will be useful, but WITHOUT
ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
for more details.

You should have received a copy of the GNU General Public License
along with this program; if not, write to the Free Software
Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
USA

---

## Usage
### Compiling
The C code of polarLDA is compiled by Makefile. The executable file `code/polarLDA`
was compiled in our lab server. If you want to use it in a different linux system,
we suggest you recompile it by the following commands in terminal:
```bash
# change to code directory
cd code
make
```
### Estimate latent factors/topics
To estimate latent factors/topics, run the following code in terminal:
```bash
polarLDA est [alpha] [k] [settings] [data] [random/seeded/*] [directory] [random_seed]
```

1. Input data format

    The input data format is the same as LDA C implementation. Thus, each document is succinctly represented as a sparse vector of word counts. `[data]` is a file where each line is of the form:
    ```
    [M] [term_1]:[count] [term_2]:[count] ... [term_N]:[count]
    ```
    where `[M]` is the number of unique terms in the document, and `[count]` associated with each term is how many times that term appears in the document. Note that `[term_1]` is an integer which indexes the term rather than a string.
    
    For example, let us assume that the document is associated with two languages, i.e., English and French. Let negative sign denote English, and positive sign denote French. Then `1:10` indicates that the 1st word is in French and it appears 10 times in the document; `8:-5` indicates that the 8th word is in English and it appears 5 times in the document.    

2. Input settings file

    See `code/inf-settings.txt` for an example. These are placeholder values; they should be experimented with.
    ```
    var max iter [integer; e.g., 10 or -1]
    var convergence [float; e.g., 1e-8]
    em max iter [integer; e.g., 100]
    em convergence [float; e.g., 1e-5]
    alpha [fixed/estimate]
    ```
    * `[var max iter]`
    
        The maximum number of iterations of coordinate ascent variational inference for a single document. A value of -1 indicates "full" variational inference, until the variational convergence criterion is met.
    * `[var convergence]`
    
        The convergence criteria for variational inference. Stop if (score_old - score) / abs(score_old) is less than this value (or after the maximum number of the iterations). Note that the score is the lower bound on the likelihood for a particular document.

    * `[em max iter]`
    
        The maximum number of iterations of variational EM.

    * `[em convergence]`

        The convergence criteria for varitional EM.  Stop if (score_old - score) / abs(score_old) is less than this value (or after the maximum number of iterations).  Note that the score is the lower bound on the likelihood for the whole corpus.
 
    * `[alpha]`
    
        If set to `[fixed]`, then alpha does not change from iteration to iteration. 
        If set to `[estimate]`, then alpha is estimated along with the topic distributions.

3. Initialization

    The term `[random/seeded/*]` describes how the topics will be initialized. 
    * `[random]` initializes each topic randomly and you can specify the random seed in `[random_seed]`, e.g. 1, 2, 3.
    * `[seeded]` initializes each topic to a distribution smoothed from a randomly chosen document. To change the number of initial documents used, edit `polarlda-estimate.c`. 
    * Alternatively, you can specify a model name to load a pre-existing model as the initial model (this is useful to continue EM from where it left off).

4. Output files

    The model (i.e., {alpha, beta, rho}) and variational posterior Dirichlet parameters will be saved in `[directory]` every 5 iterations. Additionally, there will be a log file for the likelihood bound and convergence score at each iteration. The algorithm runs until the convergence score is less than `[em convergence]` or until `[em max iter]` iterations (in the settings file) are reached. To change the lag between saved models, edit `polarlda-estimate.c`.

    * The saved model parameters are in 3 files:
        * `[iteration].other` contains alpha, where `final.other` is the final estimate;
        * `[iteration].beta` contains the log of the Multinomial distribution of topics. Each line is a topic; in line k, each entry is log p(w | z = k); and `final.beta` is the final estimate;
        * `[iteration].rho` contains the log of the Bernoulli distribution of polarity. Each line is a topic; in line k, each entry is log p(y | w, z = k); and `final.rho` is the final estimate;
    * The saved variational posterior Dirichlets are in:
        * `[iteration].gamma`, where `final.gamma` is the final estimate.
    * The saved variational posterior Categoricals are in:
        * `final.phi`. Phi is a M x K x V matrix, where M is number of documents, K is number of topics, V is number of unique words. To store it, we take the log and reshape it into (M x K) x V matrix. Each row is a 1 x V vector, and we only save non-zero word and its probability. Thus, the format is 
        ```
        [V] [term_1]:[prob] [term_2]:[prob] ... [term_N]:[prob]
        ```
        In `final.phi`, the first K rows are for M = 1, the second K rows are for K = 2 and so on. For example, if M = 10, K = 3, then you have 30 rows, the 1st row is for document 1 topic 1, the 2nd row is for document 1 topic 2, the 3nd row is for document 1 topic 3, the 4th row is for document 2 topic 1 and so on. 

### Infer factor/topic loadings
To infer loadings of the latent factors/topics, run the following code in terminal:
```bash
polarLDA inf [settings] [model] [data] [output_name]
```
Variational inference is performed on `[data]` using the estimated model parameters in `[model].*` (i.e., `[model].other`, `[model].beta`, `[model].rho`, see previous section "Output files"). We usually use "final" for `[model]` (i.e., use the final estimated model parameters for inference).

Two files will be created:
* `[output_name].gamma` contains the variational Dirichlet parameters for each document;
* `[output_name].likelihood` contains the bound on the likelihood for each document.

### Examples
An example of using the polarLDA C codes can be found in `$CBIG_CODE_DIR/stable_projects/disorder_subtypes/Tang2019_ASDFactors`, which includes bash wrapper scripts and unit tests.

