# Chord 
The chord diagram is used to visualize the 18x18 (Thomas' 17Networks + Subcortical) functional connectivity matrix.
The order of the 18 Networks are the same as 419x419 FC plot.
This code was originally written by Nanbo Sun and modified by Angela Tam.

# What does this folder contain?
* `bkr_colorscale.txt` : a file for a specific color scale used by `CBIG_chord_diagram_18x18.r`
* `CBIG_TRBPC_chord_diagram_18x18.r` : an R script to make the chord diagram
* `CBIG_TRBPC_chord.sh` : a bash script that calls the R script
* `CBIG_TRBPC_r_env.txt` : a specification file to build a conda environment that is identical to the environment that produced the original figures

# Installation
Because the code is written in R, you need to install R and the required packages first.

If you are using our `circuv`, I suggest the following way to install R.
1. Install `miniconda`
To install `miniconda`, you can check `setup/python_env_setup`
2. Install `r`

```bash
# Create and build a new conda environment called r_env
conda create --name r_env --file CBIG_TRBPC_r_env.txt
```

If the above procedure fails, please try the following:

```bash
# Create a new conda environment called r_env
conda create -n r_env

# Switch to r environment
source activate r_env

# Installs R
conda install -c r r
```

After the conda environment is set up, you need to install some R packages.

```bash
# Switch to the R environment
source activate r_env

# Call R
R

# Install the required packages in R in the following order
# Please note that you will be prompted to choose a mirror for the download

install.packages('circlize')
install.packages('igraph')
install.packages('rjson')
packageurl <- "https://cran.r-project.org/src/contrib/Archive/GetoptLong/GetoptLong_0.1.7.tar.gz"
install.packages(packageurl, repos=NULL, type="source")
source("https://bioconductor.org/biocLite.R")
biocLite("ComplexHeatmap")
install.packages('gridExtra')

# Please note that after the package 'ComplexHeatmap' is installed, you will encounter the following message
Old packages: 'boot', 'class', 'cluster', 'codetools', 'KernSmooth', 'lattice',
  'MASS', 'Matrix', 'mgcv', 'nlme', 'nnet', 'rpart', 'spatial', 'survival'
Update all/some/none? [a/s/n]:
# You can answer none (n)
```


# Usage
There are 4 arguments for the `CBIG_TRBPC_chord_diagram_18x18.r` script.
1. `mat18`: a csv file containing 18x18 matrix. It does not have a header. This is an output from ../matrix_plots/CBIG_TRBPC_chord_matrix_merge_pos_neg.m
2. `link_colorscale`: bkr or bwr. bkr means blue-black-red colorscale which matches our FC matrix colorscale.
                      bwr means blue-white-red colorscale.
3. `min_thre`: minimal absolute threshold for both positive and negative connectivity.
4. `max_thre`: maximal absolute threshold for both positive and negative connectivity.
5. `outname`: a eps file 
Please check `CBIG_TRBPC_chord.sh` for an example.
