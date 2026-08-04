# Replication

Reproducing the paper's full results requires the complete multi-dataset cohorts, which are
governed by their providers' Data Use Agreements and are not redistributed here. That material is
**deferred to a future release**; this folder currently holds only the project configuration:

- `config/CBIG_DELSSOME_python_env.yml` — conda environment for the group-level pipeline
- `config/CBIG_DELSSOME_tested_config.sh` — the environment this project was last tested with
- `config/CBIG_DELSSOME_generate_standalone.sh` — build a stand-alone copy of the project

Nothing here is needed to *use* DELSSOME:

- per-participant E/I ratios → [../indiv_level/README.md](../indiv_level/README.md)
- group-level model inversion → [../group_level/README.md](../group_level/README.md)
- runnable demos on bundled data → [../examples/](../examples/)
