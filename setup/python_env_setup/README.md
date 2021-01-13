# Python distribution used in CBIG
The Python distribution used in CBIG is based on [Anaconda](https://www.continuum.io/anaconda-overview). We use [conda](https://github.com/conda/conda) to manage Python packages, especially for packages that rely on external dependencies (e.g. Numpy). We use Python 3.5 in our lab, but packages that we use are available in both Python 2.7 and 3.5.

# Neuroimaging and Machine Learning Packages

Defined in `CBIG_pip_packages.txt`

- [nipy](https://github.com/nipy/nipy) version 0.4.0
- [nibabel](https://github.com/nipy/nibabel) version 2.0.2
- [nilearn](https://github.com/nipy/nibabel) version 0.2.5.1

Note that `scikit-learn` is included in Anaconda.

# Quick Installation (for Linux)

Before installing, please check if you have a `.condarc` file in your home directory. This is an optional configuration file of Conda that you might have (not a default configuration file installed by Conda). If so, please back up and rename this file before installing, e.g. `mv ~/.condarc ~/.condarc_bak`.

Open a terminal, switch to this directory and run the quick setup script:
```
cd $CBIG_CODE_DIR/setup/python_env_setup
bash CBIG_python_env_generic_setup.sh
```

We highly recommend installing Python 3.5 environment.

**[IMPORTANT]** To test if the setup is successful, **log out and log in again** before running: 
```bash
source activate CBIG_py3
python -m pytest $CBIG_CODE_DIR/setup/python_env_setup/tests/CBIG_python_env_setup_unit_test.py
```
# Install other packages for your own needs

The default packages installed with `CBIG_python_env_generic_setup.sh` are meant to maintain a common Python environment for users of the CBIG repo. If you wish to install other packages for your own test/ usage, you can install them using conda. See our guide below for [installing a package with conda](https://github.com/YeoPrivateLab/CBIG_private/tree/develop/setup/python_env_setup#install-a-certain-package).

If you deem the package to be useful for other people, please contact the admins so that we can have the package installed by default.

# Switch between conda environments
You can switch between CBIG's Python 3 environment and Python 2 environment by the command `source activate your_env_name`.
`your_env_name` is `CBIG_py3` for Python 3 environment and `CBIG_py2` for Python 2 enviroment.
```bash
username@yourpc:~> source activate CBIG_py2
(CBIG_py2) username@yourpc:~>
(CBIG_py2) username@yourpc:~> python my_python_2_script.py
(CBIG_py2) username@yourpc:~>
(CBIG_py2) username@yourpc:~> source activate CBIG_py3
(CBIG_py3) username@yourpc:~> python my_python_3_script.py
(CBIG_py3) username@yourpc:~>
(CBIG_py3) username@yourpc:~> source deactivate
username@yourpc:~>
```

# Conda crash course
Also see the official [conda](http://conda.pydata.org/docs/) documentation for details.
### List information of all available packages
```bash
username@userpc:~> conda list
```
### Get information of a certain package
```bash
username@userpc:~> conda list numpy
```

### Search for available versions of a certain package
```bash
username@userpc:~> conda search numpy
```

### Install a certain package
If no version number is given, conda will automatically install the latest version
```bash
username@yourpc:~> conda install numpy=1.9.3
```

### Update a certain package
```bash
username@yourpc:~> conda update numpy
```

### Update all packages
```bash
username@yourpc:~> conda update --all
```

### Create another conda environment
```bash
username@yourpc:~> conda create --name whatever_name what_software
```
#### Example
```bash
username@yourpc:~> conda create --name fancyname python=2.7
```

### Export conda environment
```bash
username@yourpc:~> conda env export --name your_env_name > sometextfile.txt
```

### Create an identical conda environment from that text file
```bash
username@yourpc:~> conda env create --name your_env_name --file sometextfile.txt
```

