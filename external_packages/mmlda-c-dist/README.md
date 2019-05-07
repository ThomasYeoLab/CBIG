## Multi-modal Latent Dirichlet Allocation (MMLDA)
The multi-modal latent Dirichlet allocation (MMLDA) is an extension of latent 
Dirichlet allocation (LDA; Blei et al., 2003). If you have not read the LDA paper 
, I suggest you read it before you go through MMLDA. C implementation 
of MMLDA is based on Blei's C code for LDA, and is located in `code` folder. For now, 
it only supports 2 modalities but it is easy to extend it to more modalities. 
In addition, the graphical model and derivation of MMLDA is in `doc/MMLDA.pdf`.

---

## Copyright
(C) Copyright 2018, Nanbo Sun (nanbosun@u.nus.edu)

This file is part of MMLDA-C.

MMLDA-C is free software; you can redistribute it and/or modify it under
the terms of the GNU General Public License as published by the Free
Software Foundation; either version 2 of the License, or (at your
option) any later version.

MMLDA-C is distributed in the hope that it will be useful, but WITHOUT
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
The c code of MMLDA is compiled by Makefile. The executable file `code/MMLDA` 
has been compiled in CIRC server. If you want to use it in a different linux
system, I suggest you recompile it by following
```bash
# change to code directory
cd code
make
``` 

### Estimate latent factors/topics
To estimate latent factors/topics, Run the following code in terminal:
```bash
MMLDA est [alpha] [k] [settings] [data1] [data2] [random/seeded/*] [directory] [random_seed]
```

1. Input data format

    The input data format is the same as LDA C implementation.Thus,
    each document is succinctly represented as a sparse vector of word
    counts. The `[data1], [data2]` is a file where each line is of the form:
    ```
    [M] [term_1]:[count] [term_2]:[count] ...  [term_N]:[count]
    ```
    where `[M]` is the number of unique terms in the document, and the
    `[count]` associated with each term is how many times that term appeared
    in the document.  Note that `[term_1]` is an integer which indexes the
    term rather than a string.

2. Input setting file

    See `code/settings-100iter.txt` for a sample. These are placeholder values; 
    they should be experimented with.
    ```
    var max iter [integer e.g., 10 or -1]
    var convergence [float e.g., 1e-8]
    em max iter [integer e.g., 100]
    em convergence [float e.g., 1e-5]
    alpha [fixed/estimate]
    ```
    * `[var max iter]`

    The maximum number of iterations of coordinate ascent variational
    inference for a single document.  A value of -1 indicates "full"
    variational inference, until the variational convergence
    criterion is met.

    * `[var convergence]`

    The convergence criteria for variational inference.  Stop if
    (score_old - score) / abs(score_old) is less than this value (or
    after the maximum number of iterations).  Note that the score is
    the lower bound on the likelihood for a particular document.

    * `[em max iter]`

    The maximum number of iterations of variational EM.

    * `[em convergence]`

    The convergence criteria for varitional EM.  Stop if (score_old -
    score) / abs(score_old) is less than this value (or after the
    maximum number of iterations).  Note that "score" is the lower
    bound on the likelihood for the whole corpus.

    * `[alpha]`

    If set to `[fixed]` then alpha does not change from iteration to
    iteration.  
    If set to `[estimate]`, then alpha is estimated along
    with the topic distributions.

3. Initialization

    The term `[random/seeded/*]` describes how the topics will be
    initialized. `[random]` initializes each topic randomly and you can specify 
    `[random_seed]` e.g. 1, 2, 3; `[seeded]` initializes each topic to a distribution 
    smoothed from a randomly chosen document. Alternatively, you can specify a model
    name to load a pre-existing model as the initial model (this is useful to continue EM
    from where it left off). To change the number of initial documents
    used, edit `mmlda-estimate.c`.

4. Output files

    The model (i.e., alpha, beta1, beta2) and variational posterior
    Dirichlet parameters will be saved in `[directory]` every
    5 iterations. Additionally, there will be a log file for the
    likelihood bound and convergence score at each iteration.  The
    algorithm runs until that score is less than `[em_convergence]` (from
    the settings file) or `[em_max_iter]` iterations are reached. To
    change the lag between saved models, edit `mmlda-estimate.c`.

    * The saved models are in 3 files:
        * `[iteration].other` contains alpha.
        * `[iteration].beta1` contains the log of the topic distributions for the first modality. 
          Each line is a topic; in line k, each entry is log p(w | z=k)
        * `[iteration].beta2` contains the log of the topic distributions for the second modality. 
          Each line is a topic; in line k, each entry is log p(w | z=k)

    * The variational posterior Dirichlets are in:
        * `[iteration].gamma`

### Infer factor/topic loadings
For the inference, you can infer based on two modalities or only one modality.
```bash
# two modalities
MMLDA inf [settings] [model] [data1] [data2] [name]
# modality 1
MMLDA inf1 [settings] [model] [data1] [name]
# modality 2
MMLDA inf2 [settings] [model] [data2] [name]
```

Variational inference is performed on the data using the model in
`[model].*` (i.e., `[model].other`, `[model].beta1`, `[model].beta2`, see previous 
section "Output files"). For the `[model]` we usually use "final" (use the final 
estimated model parameters for inference). 

Two files will be created : 
* `[output_name].gamma` are the variational Dirichlet parameters for each document;
* `[output_name].likelihood` is the bound on the likelihood for each document.

### Examples
An example of using MMLDA C code can be found in 
`$CBIG_CODE_DIR/stable_projects/disorder_subtypes/Sun2019_ADJointFactors`  
project folder. You can find examples and unit tests there. 
