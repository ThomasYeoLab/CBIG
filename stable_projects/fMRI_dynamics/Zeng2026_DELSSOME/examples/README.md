# DELSSOME examples

Runnable demos for both pipelines. Every demo uses bundled **synthetic** data (see
[../README.md](../README.md), "Example data provenance"), runs on a single CPU core, and finishes
by comparing its output against a reference snapshot — so each one doubles as a regression test.

## Which demo do I want?

| Demo                                                                        | What it shows                                                                                            | Run time |
| --------------------------------------------------------------------------- | -------------------------------------------------------------------------------------------------------- | -------- |
| [indiv_level/CBIG_DELSSOME_indiv_EI_example.sh](indiv_level/)                | **Your own FC + FCD CDF → a per-region E/I ratio CSV.** The workflow most people are looking for. | ~3 min   |
| [indiv_level/CBIG_DELSSOME_indiv_gamlss_example.sh](indiv_level/)            | A cohort of E/I ratios → a normative lifespan trajectory, plus site harmonisation                       | ~15 s    |
| [group_level/scripts/CBIG_DELSSOME_CMAES_pFIC_example.sh](group_level/)      | Group-level pFIC inversion, Euler mode                                                                   | ~2 min   |
| [group_level/scripts/CBIG_DELSSOME_CMAES_pFIC_DL_example.sh](group_level/)   | Group-level pFIC inversion,**DELSSOME mode** (needs the pretrained predictor)                      | <1 min   |
| [group_level/scripts/CBIG_DELSSOME_CMAES_Hopf_example.sh](group_level/)      | Group-level Hopf inversion, Euler mode                                                                   | ~5 min   |
| [group_level/scripts/CBIG_DELSSOME_predictor_train_example.sh](group_level/) | Training a group-level DELSSOME cost predictor from scratch                                              | ~1 min   |

Group-level run times were measured on a single core of an **AMD EPYC 9555** (64-core, 3.20 GHz,
256 MB cache, 12-channel DDR5-6000); your hardware will differ.

## Prerequisites

The two pipelines use different conda environments:

```bash
# group-level demos
conda env create -f ../replication/config/CBIG_DELSSOME_python_env.yml
conda activate CBIG_DELSSOME

# individual-level demos (adds R + r-gamlss)
conda env create -f ../indiv_level/environment.yml
conda activate CBIG_DELSSOME_indiv
```

The reference outputs were generated on **CPU only**. Reproducing them bit-for-bit requires the
same package versions, which is what these environment files pin.

## Layout

```
examples/
├── README.md            <- this file
├── group_level/
│   ├── README.md
│   ├── scripts/         <- the four demo shell scripts
│   ├── config/          <- example .ini configs
│   ├── input/           <- group-level demo SC / FC / FCD-CDF
│   ├── reference_output/
│   ├── CBIG_DELSSOME_check_*.py
│   └── output/          <- created at run time (gitignored)
└── indiv_level/
    ├── README.md
    ├── input/           <- per-subject demo FC / FCD-CDF + GAMLSS tables
    ├── config/          <- demo-only speed overrides
    ├── reference_output/
    ├── CBIG_DELSSOME_indiv_EI_example.sh
    ├── CBIG_DELSSOME_indiv_gamlss_example.sh
    ├── CBIG_DELSSOME_indiv_check_*.py
    └── output/          <- created at run time (gitignored)
```

To re-run a demo from a clean slate, delete its `output/` directory. `input/` and
`reference_output/` are read-only release artefacts.
