# Python distribution used in CBIG

The Python distribution used in CBIG is based on [Miniconda](https://docs.conda.io/en/latest/miniconda.html). We use [conda](https://github.com/conda/conda) to manage Python packages, especially for packages that rely on external dependencies (e.g. Numpy). We use Python 3.6 in our lab.

# Neuroimaging and Machine Learning Packages

Defined in `CBIG_python3_env.yml`

-   [numpy](https://numpy.org/) version 1.19.5
-   [scipy](https://scipy.org/) version 1.5.4
-   [scikit-learn](https://scikit-learn.org/stable/) version 0.24.2
-   [tedana](https://tedana.readthedocs.io/en/stable/) version 0.0.11
-   [torch](https://pytorch.org/) version 1.9.0

# Quick Installation (for Linux)

Before installing, please check if you have a `.condarc` file in your home directory. This is an optional configuration file of Conda that you might have (not a default configuration file installed by Conda). If so, please back up and rename this file before installing, e.g. `mv ~/.condarc ~/.condarc_bak`.

To install the `CBIG_py3` environment, we provide a quick installation script under $CBIG_CODE_DIR/setup/python_env_setup. A common cause for errors when installing the CBIG_py3 environment is an outdated conda version. Please ensure that the conda version is equivalent to `4.10` or later. You can use `conda update conda` to update your conda to latest version.
If the `CBIG_py3` environment already exists, it will be renamed to `CBIG_py3_bak`as a backup.

To run the quick installation, open a terminal, switch to this directory and run following command:

```
cd $CBIG_CODE_DIR/setup/python_env_setup
bash CBIG_python_env_generic_setup.sh
```

**[IMPORTANT]** To test if CBIG_py3 is successful, **log out and log in again** before running:

```bash
source activate CBIG_py3
python -m pytest $CBIG_CODE_DIR/setup/python_env_setup/tests/CBIG_python_env_setup_unit_test.py
```

When `CBIG_python_env_setup_unit_test.py` is successful, you should see a message similar to :

```
================== 8 passed in 18.28s ==================
```

## Aws-cli package

Now we also provide the installation script in `CBIG_python_env_aws_setup.sh` for another environment called `CBIG_py3_aws`, which only contains package aws-cli for HCP data downloading. 

To use aws-cli, firstly you need follow the instructions [here](https://wiki.humanconnectome.org/display/PublicData/How+To+Connect+to+Connectome+Data+via+AWS) on how to create AWS credentials and connect to AWS data using these credentials. After that, you need to generate your configure [here](https://docs.aws.amazon.com/cli/latest/userguide/cli-configure-quickstart.html) and **save it (otherwise you need to regenerate it)**. 

After the above steps have been done, you can execute this script `CBIG_python_env_aws_setup.sh` and install aws-cli. The command for quick installation is as following:

```
cd $CBIG_CODE_DIR/setup/python_env_setup
sh CBIG_python_env_aws_setup.sh
```

**[IMPORTANT]** To test if CBIG_py3_aws is successful, **log out and log in again** before running:

```bash
source activate CBIG_py3_aws
aws configure
```

If successfully installed, it will ask you to input configuration basis. If you've obtained and saved your configuration basis, then you can input them. If you do not have configuration basis, please refer to the above paragraph for instructions.

# Install other packages for your own needs

The default packages installed with `CBIG_python_env_generic_setup.sh` are meant to maintain a common Python environment for users of the CBIG repo.

If you wish to install other packages for your own test/usage, you are **highly recommended** to create a new conda environment and install packages under your own conda environment.
See our guide for [installing a package with conda](https://github.com/YeoPrivateLab/CBIG_private/tree/develop/setup/python_env_setup#install-a-certain-package).

If you deem the package to be useful for whole lab, please contact the admins so that we can have the package installed by default.

# Switch between conda environments

You can switch between CBIG's Python 3 environment and your own environment by the command `source activate <Your_Environment>`.

```bash
username@yourpc:~> source activate <Your_Environment>
(<Your_Environment>) username@yourpc:~>
(<Your_Environment>) username@yourpc:~> python my_python3_script.py
(<Your_Environment>) username@yourpc:~>
(<Your_Environment>) username@yourpc:~> source activate CBIG_py3
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
username@yourpc:~> conda install numpy=1.21.2
```

If conda could not find the package (e.g. [tedana](https://tedana.readthedocs.io/en/stable/) installed in CBIG_py3), you could try [pip](https://pypi.org/) to install the package.

```bash
username@yourpc:~> pip install tedana=0.0.11
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
username@yourpc:~> conda create --name fancyname python=3.6
```

### Export conda environment

```bash
username@yourpc:~> conda env export --name your_env_name >  your_env_name_config.yml
```

### Create an identical conda environment from config file

```bash
username@yourpc:~> conda env create --file env_config.yml
```
