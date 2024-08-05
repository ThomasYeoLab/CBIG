# Meta-matching to translate phenotypic predictive models from big to small data on structural
MRI


## References

+ Wulan, Naren, et al. "Translating phenotypic prediction models from big to small anatomical MRI data 
  using meta-matching." bioRxiv (2024): 2023-12.

----

## Background

Small sample size on structural MRI is inevitable in reality and significantly limits phenotypic 
prediction performance. Our goal was to improve prediction performance on small datasets for 
structural MRI brain imaging. We adapted the meta-matching framework from functional to structural MRI, 
and compared it with baseline methods (Elastic net and direct transfer learning). Our meta-matching-based 
approaches can greatly boost behavioral prediction performance for different small-scale structural MRI datasets.

![main_figures_from_paper](readme_figures/MMT1.png)

----

## Code Release

### Download whole repository
If you want to use the code from our lab's other stable projects (other than Naren2024_MMT1), 
you would need to download the whole CBIG repository.

- To download the version of the code that was last tested, you can either

    - visit this link:
    [https://github.com/ThomasYeoLab/CBIG/releases/tag/v0.32.0-Naren2024_MMT1](https://github.com/ThomasYeoLab/CBIG/releases/tag/v0.32.0-Naren2024_MMT1)

    or

    - run the following command, if you have Git installed
 
    ```
    git checkout -b Naren2024_MMT1 v0.32.0-Naren2024_MMT1
    ```
----

## Usage

### Environment setup
- Our code uses Python, here is the setup:
    1. Install [Miniconda](https://docs.conda.io/en/latest/miniconda.html) or 
       [Anaconda](https://www.anaconda.com/distribution/#download-section) with Python 3.x if you don't have conda
    2. Create conda environment from our `replication/config/CBIG_MMT1_python_env.yml` file by 
       `conda env create -f replication/config/CBIG_MMT1_python_env.yml`

### Example
- The example of our code is detailed in `examples/README.md`

### Replication
- If you have access to UK Biobank, HCP-YA, HCP-Aging datasets, you can replicate our result 
  using the instructions detailed in `replication/README.md`.


----

## Updates

- Release v0.32.0 (14/07/2024): Initial release of Naren2024_MMT1 project

----


## Bugs and Questions

Please contact Naren Wulan at wulannarenzhao@gmail.com and Thomas Yeo at yeoyeo02@gmail.com

